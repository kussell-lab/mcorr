package main

import (
	"fmt"
	"io"
	"log"
	"os"
	"runtime"

	"encoding/json"
	"math"

	"github.com/alecthomas/kingpin"
	"github.com/cheggaaa/pb"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

func main() {
	app := kingpin.New("mcorr", "Calculate mutation correlation from bacterial sequence data")
	app.Version("v0.2")

	alnFile := app.Arg("input", "Alignment file in XMFA format").Required().String()
	outFile := app.Arg("output", "Output file in CSV format").Required().String()

	jsonFile := app.Flag("json_file", "intermediate files in JSON format").Default("").String()
	mateFile := app.Flag("mate_file", "Experiment: calculate correlation between two clusters of strains").Default("").String()
	maxl := app.Flag("max_corr_len", "Maximum length of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("ncpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	numBoot := app.Flag("num_boot", "Number of bootstrapping on genes").Default("1000").Int()
	progress := app.Flag("progress", "Show progress").Default("false").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	setNumThreads(*ncpu)

	var numAlignment int
	var pbar *pb.ProgressBar
	if *progress {
		numAlignment = countAlignments(*alnFile)
		log.Printf("Total number of alignments: %d\n", numAlignment)
		pbar = pb.StartNew(numAlignment)
		defer pbar.Finish()
	}

	// prepare calculator.
	var calculator Calculator
	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3
	codonOffset := 0
	synonymous := true

	alnChan := readAlignments(*alnFile)
	var corrResChan chan CorrResults
	if *mateFile != "" {
		calculator = NewMateCalculator(codingTable, maxCodonLen, codonOffset, synonymous)
		mateAlnChan := readAlignments(*mateFile)
		corrResChan = calcTwoClade(alnChan, mateAlnChan, calculator)
	} else {
		calculator = NewCodingCalculator(codingTable, maxCodonLen, codonOffset, synonymous)
		corrResChan = calcSingleClade(alnChan, calculator)
	}

	if *jsonFile != "" {
		corrResChan = writeCorrResults(*jsonFile, corrResChan)
	}

	// prepare bootstrappers.
	bootstraps := []*Bootstrap{}
	notBootstrap := NewBootstrap("all", 1.0)
	notBootstrap.SetRandom(false)
	bootstraps = append(bootstraps, notBootstrap)
	for i := 0; i < *numBoot; i++ {
		id := fmt.Sprintf("boot_%d", i)
		sampleRatio := 1.0
		bootstraps = append(bootstraps, NewBootstrap(id, sampleRatio))
	}

	for corrResults := range corrResChan {
		if *progress {
			pbar.Increment()
		}
		for _, bs := range bootstraps {
			bs.Add(corrResults)
		}
	}

	w, err := os.Create(*outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("l,m,v,n,t,b\n")
	for _, bs := range bootstraps {
		results := bs.Results()
		for _, res := range results {
			w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, bs.ID))
		}
	}
}

// writeCorrResults
func writeCorrResults(file string, corrResChan chan CorrResults) (duplicateChan chan CorrResults) {
	duplicateChan = make(chan CorrResults)

	go func() {
		defer close(duplicateChan)
		f, err := os.Create(file)
		if err != nil {
			panic(err)
		}
		defer f.Close()

		encoder := json.NewEncoder(f)
		for corrResults := range corrResChan {
			var filteredCorrResults CorrResults
			filteredCorrResults.ID = corrResults.ID
			for _, corrRes := range corrResults.Results {
				if !math.IsNaN(corrRes.Mean) && !math.IsNaN(corrRes.Variance) {
					filteredCorrResults.Results = append(filteredCorrResults.Results, corrRes)
				}
			}
			if len(filteredCorrResults.Results) > 0 {
				if err := encoder.Encode(filteredCorrResults); err != nil {
					panic(err)
				}
				duplicateChan <- filteredCorrResults
			}
		}
	}()

	return
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

// calcSingleClade calculate correlation functions in a single cluster of sequence.
func calcSingleClade(alnChan chan Alignment, calculator Calculator) (corrResChan chan CorrResults) {
	corrResChan = make(chan CorrResults)
	done := make(chan bool)

	ncpu := runtime.GOMAXPROCS(0)
	for i := 0; i < ncpu; i++ {
		go func() {
			for aln := range alnChan {
				if len(aln.Sequences) > 1 {
					results := calculator.CalcP2(aln)
					corrResChan <- results
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(corrResChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()
	return
}

// calcTwoClade calculate correlation functions between two clades.
func calcTwoClade(alnChan, mateAlnChan chan Alignment, calculator Calculator) (corrResChan chan CorrResults) {
	type job struct {
		A, B Alignment
	}
	jobChan := make(chan job)
	go func() {
		defer close(jobChan)
		for aln := range alnChan {
			mateAln := <-mateAlnChan
			if len(aln.Sequences) > 1 && len(mateAln.Sequences) > 1 {
				j := job{A: aln, B: mateAln}
				jobChan <- j
			}
		}
	}()

	corrResChan = make(chan CorrResults)
	done := make(chan bool)

	ncpu := runtime.GOMAXPROCS(0)
	for i := 0; i < ncpu; i++ {
		go func() {
			for j := range jobChan {
				results := calculator.CalcP2(j.A, j.B)
				corrResChan <- results
			}
			done <- true
		}()
	}

	go func() {
		defer close(corrResChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()
	return
}

// setNumThreads sets number of CPU cores for using.
// if ncpu == 0, we will used all core avaible.
func setNumThreads(ncpu int) {
	if ncpu == 0 {
		ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(ncpu)
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(file string) (alnChan chan Alignment) {
	alnChan = make(chan Alignment)
	read := func() {
		defer close(alnChan)

		f := openFile(file)
		defer f.Close()
		xmfaReader := seq.NewXMFAReader(f)
		index := 0
		for {
			alignment, err := xmfaReader.Read()
			alnChan <- Alignment{ID: fmt.Sprintf("%d", index), Sequences: alignment}
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			alnChan <- Alignment{ID: fmt.Sprintf("%d", index), Sequences: alignment}
			index++
		}
	}
	go read()
	return
}
