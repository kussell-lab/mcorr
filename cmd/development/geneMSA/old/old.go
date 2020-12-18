package old

import (
	"fmt"
	"github.com/apsteinberg/biogo/seq"
	"gopkg.in/alecthomas/kingpin.v2"
	"io"
	"os"
	"path/filepath"
	"runtime"
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
	defer close(done)
	alignments, errc := readAlignments(done, *alnFile)
	//start a fixed number of goroutines to read and digest files
	c := make(chan result)
	var wg sync.WaitGroup
	for i := 0; i < *numDigesters; i++ {
		wg.Add(1)
		go digester(*outdir, done, alignments, c, i, &wg)
	}

	////abort if done is closed
	//select {
	//case <- done:
	//}

	go func() {
		wg.Wait()
		close(c)
	}()
	//end of pipeline.
	//for gene := range c {
	//	writeAln(gene, *outdir)
	//}
	// Check whether reading failed
	if err := <-errc; err != nil { // HLerrc
		panic(err)
	}

	duration := time.Since(start)
	fmt.Println("Time to write gene MSA files:", duration)

}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
//func readAlignments(file string) (alnChan chan Alignment) {
//	alnChan = make(chan Alignment)
//	read := func() {
//		defer close(alnChan)
//
//		f, err := os.Open(file)
//		if err != nil {
//			panic(err)
//		}
//		defer f.Close()
//		xmfaReader := seq.NewXMFAReader(f)
//		numAln := 0
//		for {
//			alignment, err := xmfaReader.Read()
//			if err != nil {
//				if err != io.EOF {
//					panic(err)
//				}
//				break
//			}
//			if len(alignment) > 0 {
//				numAln++
//				alnID := strings.Split(alignment[0].Id, " ")[0]
//				alnChan <- Alignment{ID: alnID, Sequences: alignment}
//				fmt.Printf("\rRead %d alignments.", numAln)
//				fmt.Printf("\r alignment ID: %s", alnID)
//			}
//		}
//		fmt.Printf(" Total alignments %d\n", numAln)
//	}
//	go read()
//	return
//}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(done <-chan struct{}, file string) (<-chan []seq.Sequence, <-chan error) {
	alignments := make(chan []seq.Sequence)
	errc := make(chan error, 1)
	go func() {
		defer close(alignments)

		f, err := os.Open(file)
		if err != nil {
			panic(err)
		}
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
			select {
			case alignments <- alignment:
			case <-done:
				fmt.Println("quit reading alignments")
			}
		}
		errc <- err
	}()
	return alignments, errc
}

// digester reads path names from paths and sends digests of the corresponding
// files on alnChan until either paths or done is closed.
func digester(outdir string, done <-chan struct{}, alignments <-chan []seq.Sequence, genes chan<- result, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	fmt.Printf("Worker %d starting\n", id)

	for alignment := range alignments { // HLpaths
		alnID := strings.Split(alignment[0].Id, " ")[0]
		gene := result{alnID, alignment}
		writeAln(gene, outdir)
		//case <- done:
		//select {
		//case genes <- gene:
		//case <-done:
		//	return
		//}
	}
	fmt.Printf("Worker %d done\n", id)
}

// A result is a single gene alignment
type result struct {
	ID        string
	Sequences []seq.Sequence
}

//// processAlignments reads sequence alignment from a extended Multi-FASTA file,
//// and writes them to a separate gene MSA
//func processAlignments(done <- chan struct{}, file string, outdir string, id int, wg* sync.WaitGroup) (<- chan []seq.Sequence, <- chan error) {
//	alignments := make(chan []seq.Sequence)
//	errc := make(chan error, 1)
//	go func() {
//		defer close(alignments)
//
//		f, err := os.Open(file)
//		if err != nil {
//			panic(err)
//		}
//		defer f.Close()
//		xmfaReader := seq.NewXMFAReader(f)
//		for {
//			alignment, err := xmfaReader.Read()
//			if err != nil {
//				if err != io.EOF {
//					panic(err)
//				}
//				break
//			}
//			select {
//			case alignments <- alignment:
//				writeAln(alignment, outdir)
//			case <- done:
//				fmt.Println("quit writing alignment files")
//			}
//		}
//		errc <- err
//	}()
//	return alignments, errc
//}

func writeAln(aln result, outdir string) {

	geneName := strings.Split(aln.ID, "|")[1]
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
	f.Close()

}

func getNames(s string) (geneName, genomePos, genomeName string) {
	terms := strings.Split(s, " ")
	geneName = terms[0]
	genomePos = terms[1]
	genomeName = terms[2]
	return
}

// writeAll reads all the genes in the master XMFA, and writes gene alignments for each gene
func writeAll(alnFile string, outdir string, numDigesters int) error {
	done := make(chan struct{})
	defer close(done)
	alignments, errc := readAlignments(done, alnFile)
	//start a fixed number of goroutines to read and digest files
	c := make(chan result)
	var wg sync.WaitGroup
	for i := 0; i < numDigesters; i++ {
		wg.Add(1)
		go digester(outdir, done, alignments, c, i, &wg)
	}
	go func() {
		wg.Wait()
		close(c)
	}()
	//end of pipeline.
	// Check whether reading failed
	if err := <-errc; err != nil { // HLerrc
		return err
	}
	return nil
}
