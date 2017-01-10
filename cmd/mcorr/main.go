package main

import (
	"fmt"
	"io"
	"log"
	"os"
	"runtime"

	"github.com/alecthomas/kingpin"
	"github.com/cheggaaa/pb"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

func main() {
	app := kingpin.New("mcorr", "calculate correlated mutations for bacterial genomes")
	app.Version("v0.2")

	alnFile := app.Arg("align_file", "alignment file in XMFA format").Required().String()
	outFile := app.Arg("output_file", "output file in CSV format").Required().String()

	mateFile := app.Flag("mate_file", "mate file").Default("").String()
	maxl := app.Flag("max_len", "maximum length of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("threads", "number of threads (default: using all available cores)").Default("0").Int()
	numBoot := app.Flag("num_boot", "number of bootstrapping of genes").Default("1000").Int()
	progress := app.Flag("progress", "show progress").Default("false").Bool()

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
	var corrResChan chan []CorrResult
	if *mateFile != "" {
		calculator = NewMateCalculator(codingTable, maxCodonLen, codonOffset, synonymous)
		mateAlnChan := readAlignments(*mateFile)
		corrResChan = calcTwoClade(alnChan, mateAlnChan, calculator)
	} else {
		calculator = NewCodingCalculator(codingTable, maxCodonLen, codonOffset, synonymous)
		corrResChan = calcSingleClade(alnChan, calculator)
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

// Alignment is an array of mutliple sequences with same length.
type Alignment []seq.Sequence

// calcSingleClade calculate correlation functions in a single cluster of sequence.
func calcSingleClade(alnChan chan Alignment, calculator Calculator) (corrResChan chan []CorrResult) {
	corrResChan = make(chan []CorrResult)
	done := make(chan bool)

	ncpu := runtime.GOMAXPROCS(0)
	for i := 0; i < ncpu; i++ {
		go func() {
			for aln := range alnChan {
				results := calculator.CalcP2(aln)
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

// calcTwoClade calculate correlation functions between two clades.
func calcTwoClade(alnChan, mateAlnChan chan Alignment, calculator Calculator) (corrResChan chan []CorrResult) {
	type job struct {
		A, B []seq.Sequence
	}
	jobChan := make(chan job)
	go func() {
		defer close(jobChan)
		for aln := range alnChan {
			mateAln := <-mateAlnChan
			j := job{A: aln, B: mateAln}
			jobChan <- j
		}
	}()

	corrResChan = make(chan []CorrResult)
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
		for {
			alignment, err := xmfaReader.Read()
			if len(alignment) > 0 {
				alnChan <- alignment
			}

			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
		}
	}
	go read()
	return
}
