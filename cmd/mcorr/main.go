package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/alecthomas/kingpin"
	"github.com/cheggaaa/pb"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

func main() {
	kingpin.Version("v0.1")
	alnFile := kingpin.Arg("align_file", "input alignment file in FASTA format").Required().String()
	outFile := kingpin.Arg("out_file", "output prefix").Required().String()
	noncoding := kingpin.Flag("non_coding", "non_coding sequences?").Default("false").Bool()
	partial := kingpin.Flag("partial", "partial genes (like MLST data)").Default("false").Bool()
	maxl := kingpin.Flag("max_len", "maximum length of correlation to calculate").Default("300").Int()
	ncpu := kingpin.Flag("cpus", "number of threads (default 0, use all the cores)").Default("0").Int()
	numBoot := kingpin.Flag("num_boot", "number of bootstrapping").Default("100").Int()
	showProgress := kingpin.Flag("progress", "show progress?").Default("false").Bool()
	kingpin.Parse()

	var calculator Calculator
	if *noncoding {
		calculator = NewNoncodingCalculator(*maxl)
	} else {
		codingTable := taxonomy.GeneticCodes()["11"]
		maxCodonLen := *maxl / 3
		if *partial {
			calculator = NewPartialCalculator(codingTable, maxCodonLen)
		} else {
			codonOffset := 0
			synonymous := true
			calculator = NewCodingCalculator(codingTable, maxCodonLen, codonOffset, synonymous)
		}
	}

	setNumThreads(*ncpu)
	var pbar *pb.ProgressBar
	if *showProgress {
		count := countAlignments(*alnFile)
		pbar = pb.StartNew(count)
		defer pbar.Finish()
	}

	alnChan := readAlignments(*alnFile)
	corrResChan := calc(alnChan, calculator)

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
		if *showProgress {
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

func calc(alnChan chan []seq.Sequence, calculator Calculator) (corrResChan chan []CorrResult) {
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

// setNumThreads sets number of CPU cores for using.
// if ncpu == 0, we will used all core avaible.
func setNumThreads(ncpu int) {
	if ncpu == 0 {
		ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(ncpu)
}

// openFile is a helper function to open a file.
// and panic if error occurs.
func openFile(file string) (f *os.File) {
	var err error
	f, err = os.Open(file)
	if err != nil {
		panic(err)
	}
	return
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(file string) (alnChan chan []seq.Sequence) {
	alnChan = make(chan []seq.Sequence)
	read := func() {
		defer close(alnChan)

		f := openFile(file)
		defer f.Close()
		xmfaReader := NewXMFAReader(f)
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

// countAlignments return total number of alignments in a file.
func countAlignments(file string) (count int) {
	f := openFile(file)
	defer f.Close()
	rd := bufio.NewReader(f)
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		if line[0] == '=' {
			count++
		}
	}
	return
}
