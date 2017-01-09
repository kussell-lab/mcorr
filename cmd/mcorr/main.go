package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/alecthomas/kingpin"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

func main() {
	kingpin.Version("v0.1")
	alnFile := kingpin.Arg("align_file", "input alignment file in FASTA format").Required().String()
	outFile := kingpin.Arg("out_file", "output prefix").Required().String()
	mateFile := kingpin.Flag("mate_file", "mate file").Default("").String()
	noncoding := kingpin.Flag("non_coding", "non_coding sequences?").Default("false").Bool()
	partial := kingpin.Flag("partial", "partial genes (like MLST data)").Default("false").Bool()
	maxl := kingpin.Flag("max_len", "maximum length of correlation to calculate").Default("300").Int()
	ncpu := kingpin.Flag("cpus", "number of threads (default 0, use all the cores)").Default("0").Int()
	numBoot := kingpin.Flag("num_boot", "number of bootstrapping").Default("100").Int()
	splitGeneSize := kingpin.Flag("gene_size", "split genes").Default("0").Int()
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
			if *mateFile != "" {
				calculator = NewMateCalculator(codingTable, maxCodonLen, codonOffset, synonymous)
			} else {
				calculator = NewCodingCalculator(codingTable, maxCodonLen, codonOffset, synonymous)
			}
		}
	}

	setNumThreads(*ncpu)

	alnChan := readAlignments(*alnFile)
	var corrResChan chan []CorrResult
	if *mateFile != "" {
		mateAlnChan := readAlignments(*mateFile)
		corrResChan = calcWithMate(alnChan, mateAlnChan, calculator, *splitGeneSize)
	} else {
		corrResChan = calc(alnChan, calculator, *splitGeneSize)
	}

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

func calcWithMate(alnChan, mateAlnChan chan []seq.Sequence, calculator Calculator, geneSize int) (corrResChan chan []CorrResult) {
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
				aln := j.A
				mateAln := j.B
				subAln := [][]seq.Sequence{}
				mateSubAln := [][]seq.Sequence{}
				if geneSize <= 0 {
					subAln = append(subAln, aln)
					mateSubAln = append(mateSubAln, mateAln)
				} else {
					subAln = splitAlignment(aln, geneSize)
					mateSubAln = splitAlignment(mateAln, geneSize)
				}
				for i, aa := range subAln {
					results := calculator.CalcP2(aa, mateSubAln[i])
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

func calc(alnChan chan []seq.Sequence, calculator Calculator, geneSize int) (corrResChan chan []CorrResult) {
	corrResChan = make(chan []CorrResult)
	done := make(chan bool)

	ncpu := runtime.GOMAXPROCS(0)
	for i := 0; i < ncpu; i++ {
		go func() {
			for aln := range alnChan {
				subAln := [][]seq.Sequence{}
				if geneSize <= 0 {
					subAln = append(subAln, aln)
				} else {
					subAln = splitAlignment(aln, geneSize)
				}
				for _, aa := range subAln {
					results := calculator.CalcP2(aa)
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

// splitAlignment splits an alignment into mutiple parts.
func splitAlignment(sequences []seq.Sequence, size int) [][]seq.Sequence {
	genes := [][]seq.Sequence{}
	for _, gene := range sequences {
		s := gene.Seq
		parts, starts, ends := split(s, size)

		if len(genes) == 0 {
			genes = make([][]seq.Sequence, len(parts))
		}

		for i := 0; i < len(parts); i++ {
			id := fmt.Sprintf("%s_%d_%d", gene.Id, starts[i], ends[i])
			g := seq.Sequence{Id: id, Seq: parts[i]}
			genes[i] = append(genes[i], g)
		}
	}

	return genes
}

// split a sequence into multiple constant parts of constant size.
func split(s []byte, size int) (parts [][]byte, starts, ends []int) {
	numOfParts := len(s) / size
	mode := len(s) % size
	beginOfSeq := mode / 2

	for i := 0; i < numOfParts; i++ {
		part := []byte{}
		start := beginOfSeq + i*size
		end := beginOfSeq + (i+1)*size
		for j := start; j < end && j < len(s); j++ {
			part = append(part, s[j])
		}
		parts = append(parts, part)
		starts = append(starts, start)
		ends = append(ends, end)
	}

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
