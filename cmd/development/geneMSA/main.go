package main

import (
	"fmt"
	"github.com/apsteinberg/biogo/seq"
	"gopkg.in/alecthomas/kingpin.v2"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
)

func main() {
	app := kingpin.New("geneMSA", "split a multi-gene MSA file into MSAs for each gene")
	app.Version("v20201207")
	alnFile := app.Arg("master_MSA", "multi-sequence alignment file for all genes").Required().String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	numDigesters := app.Flag("threads", "Number of files read at a time (default: 20)").Default("20").Int()
	outdir := app.Flag("outdir", "output directory for gene MSAs").Default("genes").String()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	//timer

	start := time.Now()

	if *ncpu == 0 {
		*ncpu = runtime.NumCPU()
	}

	runtime.GOMAXPROCS(*ncpu)

	//make the folder and gene MSAs
	if _, err := os.Stat(*outdir); os.IsNotExist(err) {
		os.Mkdir(*outdir, 0755)
	}

	//reads all the genes in the master XMFA, and writes gene alignments for each gene

	done := make(chan struct{})
	//done := make(chan bool)
	//defer close(done)
	alignments, errc := readAlignments(done, *alnFile)

	//make alignment channel
	//var alnChan chan Alignment
	//alnChan = readAlignments(*alnFile)
	//done := make(chan bool)

	//start a fixed number of goroutines to read and digest files
	c := make(chan result)
	var wg sync.WaitGroup
	for i := 0; i < *numDigesters; i++ {
		wg.Add(1)
		go digester(done, alignments, c, i, &wg)
	}

	go func() {
		wg.Wait()
		close(c)
	}()
	//end of pipeline.
	for gene := range c {
		writeAln(gene.Alignment, *outdir)
	}
	//fmt.Printf("read %d alignments", c)
	//Check whether reading failed
	if err := <-errc; err != nil { // HLerrc
		panic(err)
	}

	duration := time.Since(start)
	fmt.Println("Time to write gene MSA files:", duration)

}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	num       int
	Sequences []seq.Sequence
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(done <-chan struct{}, file string) (<-chan Alignment, <-chan error) {
	alignments := make(chan Alignment)
	errc := make(chan error, 1)
	go func() {
		defer close(alignments)

		f, err := os.Open(file)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		xmfaReader := seq.NewXMFAReader(f)
		numAln := 0
		for {
			alignment, err := xmfaReader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if len(alignment) > 0 {
				numAln++
				alnID := strings.Split(alignment[0].Id, " ")[0]
				select {
				case alignments <- Alignment{alnID, numAln, alignment}:
				case <-done:
					fmt.Println("quit reading alignments")
				}
			}
		}
		errc <- err
	}()
	return alignments, errc
}

// digester reads path names from paths and sends digests of the corresponding
// files on alnChan until either paths or done is closed.
func digester(done <-chan struct{}, alignments <-chan Alignment, genes chan<- result, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	fmt.Printf("Worker %d starting\n", id)
	for aln := range alignments { // HLpaths
		//alnID := strings.Split(alignment[0].Id, " ")[0]
		gene := result{aln.num, aln}
		//writeAln(aln, outdir)
		select {
		//case c <- aln.num:
		case genes <- gene:
		//	writeAln(aln, outdir)
		case <-done:
			return
		}
	}
	fmt.Printf("Worker %d done\n", id)

}

// A result is a single gene alignment
type result struct {
	ID        int
	Alignment Alignment
}

func writeAln(aln Alignment, outdir string) {
	geneID := strconv.Itoa(aln.num)
	//geneName := strings.Split(aln.ID, "|")[1]
	geneName := "gene_" + geneID
	MSA := filepath.Join(outdir, geneName)
	f, err := os.Create(MSA)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	for _, s := range aln.Sequences {
		f.WriteString(">" + s.Id + "\n")
		f.Write(s.Seq)
		f.WriteString("\n")
	}
	f.WriteString("=\n")
	//f.Close()

}
