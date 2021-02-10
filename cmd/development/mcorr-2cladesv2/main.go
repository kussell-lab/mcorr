package main

import (
	"io"
	"os"
	"runtime"
	"strings"

	"fmt"
	"github.com/apsteinberg/biogo/seq"
	"github.com/apsteinberg/mcorr"
	"github.com/apsteinberg/ncbiftp/taxonomy"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
)

// global variables.
func main() {
	fmt.Println("Couples dance late night") //assuming fmt is imported
	app := kingpin.New("mcorr-2clades-v2", "Calculate mutation correlation from bacterial sequence alignments from two clades in XMFA format.")
	app.Version("v20200807")

	alnFile := app.Arg("in-1", "Alignment file in XMFA format.").Required().String()
	//added by Asher
	mateAlnFile := app.Arg("in-2", "Alignment file in XMFA format.").Required().String()
	outPrefix := app.Arg("out", "Output prefix.").Required().String()

	maxl := app.Flag("max-corr-length", "Maximum distance of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	numBoot := app.Flag("num-boot", "Number of bootstrapping on genes").Default("1000").Int()
	showProgress := app.Flag("show-progress", "Show progress").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	if *ncpu <= 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	// show progress bar?
	var bar *pb.ProgressBar
	if *showProgress {
		max := getNumberOfAlignments(*alnFile)
		bar = pb.StartNew(max)
		defer bar.Finish()
	}

	// prepare calculator.
	var calculator PairCalculator
	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3

	synonymous := true
	codonPos := 3
	codonOffset := 0

	var alnChan chan Alignment
	if bar == nil {
		alnChan = readAlignments(*alnFile)
	} else {
		alnChan = make(chan Alignment)
		go func() {
			defer close(alnChan)
			count := 0
			c := readAlignments(*alnFile)
			for a := range c {
				alnChan <- a
				bar.Add(1)
				count++
			}
		}()
	}

	//added by Asher
	var mateAlnChan chan Alignment
	if bar == nil {
		mateAlnChan = readAlignments(*mateAlnFile)
	} else {
		mateAlnChan = make(chan Alignment)
		go func() {
			defer close(mateAlnChan)
			count := 0
			c := readAlignments(*mateAlnFile)
			for a := range c {
				mateAlnChan <- a
				bar.Add(1)
				count++
			}
		}()
	}
	// these are the lines that potentially need to change; mate calculator for mates
	// calctwoclades for clades
	//calculator = NewCodingCalculator(codingTable, maxCodonLen, codonOffset, codonPos-1, synonymous)
	//corrResChan := calcSingleClade(alnChan, calculator)

	calculator = NewMateCalculator(codingTable, maxCodonLen, codonOffset, codonPos-1, synonymous)
	corrResChan := calcTwoClade(alnChan, mateAlnChan, calculator)
	//what's in the json is actually Qs NOT P2!
	resChan := mcorr.PipeOutCorrResults(corrResChan, *outPrefix+".json")
	//division by d_sample or P2 is not until here!!!
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
func calcTwoClade(alnChan, mateAlnChan chan Alignment, calculator PairCalculator) (corrResChan chan mcorr.CorrResults) {
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
	go func() {
		defer close(alnChan)

		c := readXMFA(file)
		for alignment := range c {
			alnID := strings.Split(alignment[0].Id, " ")[0]
			alnChan <- Alignment{ID: alnID, Sequences: alignment}
		}
	}()

	return
}

// getNumberOfAlignments return total number of alignments in a xmfa file.
func getNumberOfAlignments(file string) (count int) {
	c := readXMFA(file)
	for a := range c {
		if len(a) >= 2 {
			count++
		}
	}
	return
}

// readXMFA reads a xmfa format file and returns a channel of []seq.Sequence.
func readXMFA(file string) chan []seq.Sequence {
	c := make(chan []seq.Sequence)
	go func() {
		defer close(c)

		f := mustOpen(file)
		defer f.Close()

		rd := seq.NewXMFAReader(f)
		for {
			a, err := rd.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if len(a) >= 2 {
				c <- a
			}
		}
	}()
	return c
}

// mustOpen is a helper function to open a file.
// and panic if error occurs.
func mustOpen(file string) (f *os.File) {
	var err error
	f, err = os.Open(file)
	if err != nil {
		panic(err)
	}
	return
}
