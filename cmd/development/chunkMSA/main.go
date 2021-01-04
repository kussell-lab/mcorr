package main

import (
	"bufio"
	"fmt"
	"github.com/apsteinberg/biogo/seq"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"time"
)

func main() {

	app := kingpin.New("chunkMSA", "Splits an MSA file containing all sequences into smaller MSA file chunks'\n"+
		"for use with mcorr-pair-sync.")
	app.Version("v20201123")
	alnFile := app.Arg("master_MSA", "multi-sequence alignment file containing all sequences in the dataset").Required().String()
	strainList := app.Arg("strain_list", "list of all strains").Required().String()
	splits := app.Arg("splits", "number of MSA file chunks to split the Master MSA file into").Required().Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	chunkpath := app.Flag("chunk-folder", "folder name for chunked MSAs").Default("chunkedMSA").String()
	showProgress := app.Flag("show-progress", "Show progress").Default("false").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	//timer

	start := time.Now()

	if *ncpu == 0 {
		*ncpu = runtime.NumCPU()
	}

	runtime.GOMAXPROCS(*ncpu)

	// show progress bar
	var bar *pb.ProgressBar
	if *showProgress {
		max := countAlignments(*alnFile)
		bar = pb.StartNew(max)
		defer bar.Finish()
	}

	strains := readSamples(*strainList)
	fmt.Printf("\r %d total strains\n", len(strains))

	// chunk the strain slice into the desired number of strains
	// per split MSA
	numSplits := *splits
	strainsPer := len(strains) / numSplits
	strainChunks := chunkSlice(strains, strainsPer)

	//map which strain will go to which chunk
	//keys are strains, values are chunk number
	var chunkmap map[string]string
	chunkmap = make(map[string]string)
	chunkNum := 0
	// list of msa chunks
	var chunks []string
	//map strains to chunks
	for _, strainChunk := range strainChunks {
		msaChunk := strconv.Itoa(chunkNum)
		for _, strain := range strainChunk {
			chunkmap[strain] = msaChunk
		}
		chunks = append(chunks, msaChunk)
		chunkNum++
	}

	//make the folder and files for chunked MSAs
	if _, err := os.Stat(*chunkpath); os.IsNotExist(err) {
		os.Mkdir(*chunkpath, 0755)
	}
	//make MSAs
	for _, chunk := range chunks {
		MSA := filepath.Join(*chunkpath, "MSA_chunk"+chunk)
		f, err := os.Create(MSA)
		check(err)
		f.Close()
	}

	//define alignment channel
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

	numJob := *ncpu
	done := make(chan bool)
	chunkChan := make(chan cAlignment)
	for i := 0; i < numJob; i++ {
		go func() {
			for aln := range alnChan {
				alnMap := AssembleAlignments(aln, chunkmap)
				for ID, cAln := range alnMap {
					chunkChan <- cAlignment{chunkID: ID, Sequences: cAln}
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(chunkChan)
		for i := 0; i < numJob; i++ {
			<-done
		}
	}()

	WriteClusterMSA(chunkChan, *chunkpath)

	duration := time.Since(start)
	fmt.Println("Time to write chunked MSA files:", duration)
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

//channel structure to pump a gene alignment for a chunk into
type cAlignment struct {
	chunkID string //ID for chunk
	//	geneID    string         //gene name
	Sequences []seq.Sequence //the gene alignment for the cluster
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

//channel structure to pump clusters into
type geneSeq struct {
	seq       seq.Sequence
	clusterID string
}

func getNames(s string) (geneName, genomePos, genomeName string) {
	terms := strings.Split(s, " ")
	geneName = terms[0]
	genomePos = terms[1]
	genomeName = terms[2]
	return
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(file string) (alnChan chan Alignment) {
	alnChan = make(chan Alignment)
	read := func() {
		defer close(alnChan)

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
				alnChan <- Alignment{ID: alnID, Sequences: alignment}
				fmt.Printf("\rRead %d alignments.", numAln)
				fmt.Printf("\r alignment ID: %s", alnID)
			}
		}
		fmt.Printf(" Total alignments %d\n", numAln)
	}
	go read()
	return
}

// countAlignments return total number of alignments in a file.
func countAlignments(file string) (count int) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
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

//write the MSA file for the cluster
func WriteClusterMSA(chunkChan chan cAlignment, chunkpath string) {
	for c := range chunkChan {
		MSA := filepath.Join(chunkpath, "MSA_chunk"+c.chunkID)
		//f, err := os.OpenFile(filename, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0600)
		f, err := os.OpenFile(MSA, os.O_APPEND|os.O_WRONLY, 0600)
		if err != nil {
			panic(err)
		}
		for _, s := range c.Sequences {
			f.WriteString(">" + s.Id + "\n")
			f.Write(s.Seq)
			f.WriteString("\n")
		}
		f.WriteString("=\n")
		f.Close()
	}
}

func AssembleAlignments(aln Alignment, clustermap map[string]string) (alnMap map[string][]seq.Sequence) {
	Sequences := aln.Sequences
	alnMap = make(map[string][]seq.Sequence)
	//double check for duplicate sequences ....
	duplicates := make(map[string]bool)
	for _, s := range Sequences {
		_, _, strain := getNames(s.Id)
		_, found := duplicates[strain]
		if found {
			continue
		} else {
			duplicates[strain] = true
			cluster, found := clustermap[strain]
			if found {
				alnMap[cluster] = append(alnMap[cluster], s)
			}
		}
	}
	return
}

//channel structure to pump genes into
type clusterAln struct {
	clusterID string
	Sequences []seq.Sequence
}

//check for errors
func check(e error) {
	if e != nil {
		panic(e)
	}
}

func unique(intSlice []int) []string {
	keys := make(map[int]bool)
	list := []int{}
	strlist := []string{}
	for _, entry := range intSlice {
		if _, value := keys[entry]; !value {
			keys[entry] = true
			list = append(list, entry)
		}
	}
	sort.Ints(list)
	for _, num := range list {
		strnum := strconv.Itoa(num)
		strlist = append(strlist, strnum)
	}
	return strlist
}

// grabs the first sequence alignment from the XMFA and generates a list of strains
//as a string slice
func generateStrainSlice(file string) (strains []string) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	xmfaReader := seq.NewXMFAReader(f)
	numStrains := 0
	for {
		alignment, err := xmfaReader.Read()
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		if len(alignment) > 0 {
			for _, s := range alignment {
				_, _, strain := getNames(s.Id)
				strains = append(strains, strain)
				numStrains++
			}
			fmt.Printf("\r %d total strains", numStrains)
			break
		}
	}
	return
}

func chunkSlice(slice []string, chunkSize int) [][]string {
	var chunks [][]string
	for i := 0; i < len(slice); i += chunkSize {
		end := i + chunkSize

		// necessary check to avoid slicing beyond
		// slice capacity
		if end > len(slice) {
			end = len(slice)
		}

		chunks = append(chunks, slice[i:end])
	}

	return chunks
}

// readSamples return a list of samples from a sample file.
func readSamples(filename string) (samples []string) {
	f, err := os.Open(filename)
	if err != nil {
		log.Fatalf("Error when reading file %s:%v", filename, err)
	}
	defer f.Close()

	r := bufio.NewReader(f)
	for {
		line, err := r.ReadString('\n')

		if err != nil {
			if err != io.EOF {
				log.Fatalf("Error when reading file %s: %v", filename, err)
			}
			break
		}
		samples = append(samples, strings.TrimSpace(line))
	}
	return
}
