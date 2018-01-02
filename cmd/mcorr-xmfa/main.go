package main

import (
	"io"
	"log"
	"os"
	"runtime"
	"strings"

	"github.com/kussell-lab/biogo/seq"
	"github.com/kussell-lab/mcorr"
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"gopkg.in/alecthomas/kingpin.v2"
)

func main() {
	app := kingpin.New("mcorr-xmfa", "Calculate mutation correlation from bacterial sequence alignments in XMFA format.")
	app.Version("v20180102")

	alnFile := app.Arg("in", "Alignment file in XMFA format.").Required().String()
	outPrefix := app.Arg("out", "Output prefix.").Required().String()

	mateFile := app.Flag("second-alignment", "Second alignment file in XMFA format.").Default("").String()
	codonPos := app.Flag("codon-position", "Codon position (1: first codon position; 2: second codon position; 3: third codon position; 4: synoumous at third codon position.").Default("4").Int()
	maxl := app.Flag("max-corr-length", "Maximum distance of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	numBoot := app.Flag("num-boot", "Number of bootstrapping on genes").Default("1000").Int()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	if *ncpu <= 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	// prepare calculator.
	var calculator Calculator
	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3

	synonymous := false
	if *codonPos == 4 {
		synonymous = true
		*codonPos = 3
	}
	if *codonPos <= 0 || *codonPos > 4 {
		log.Fatalln("--codon-position should be in the range of 1 to 4.")
	}
	codonOffset := 0

	alnChan := readAlignments(*alnFile)
	var corrResChan chan mcorr.CorrResults
	if *mateFile != "" {
		calculator = NewMateCalculator(codingTable, maxCodonLen, codonOffset, *codonPos-1, synonymous)
		mateAlnChan := readAlignments(*mateFile)
		corrResChan = calcTwoClade(alnChan, mateAlnChan, calculator)
	} else {
		calculator = NewCodingCalculator(codingTable, maxCodonLen, codonOffset, *codonPos-1, synonymous)
		corrResChan = calcSingleClade(alnChan, calculator)
	}

	resChan := mcorr.PipeOutCorrResults(corrResChan, *outPrefix+".json")
	mcorr.CollectWrite(resChan, *outPrefix+".csv", *numBoot)
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

// calcSingleClade calculate correlation functions in a single cluster of sequence.
func calcSingleClade(alnChan chan Alignment, calculator Calculator) (corrResChan chan mcorr.CorrResults) {
	corrResChan = make(chan mcorr.CorrResults)
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
func calcTwoClade(alnChan, mateAlnChan chan Alignment, calculator Calculator) (corrResChan chan mcorr.CorrResults) {
	type job struct {
		A, B Alignment
	}
	jobChan := make(chan job)
	go func() {
		defer close(jobChan)
		for aln := range alnChan {
			mateAln := <-mateAlnChan
			if len(aln.Sequences) >= 1 && len(mateAln.Sequences) >= 1 {
				j := job{A: aln, B: mateAln}
				jobChan <- j
			}
		}
	}()

	corrResChan = make(chan mcorr.CorrResults)
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
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if len(alignment) > 0 {
				alnID := strings.Split(alignment[0].Id, " ")[0]
				alnChan <- Alignment{ID: alnID, Sequences: alignment}
			}
		}
	}
	go read()
	return
}
