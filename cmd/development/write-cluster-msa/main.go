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

	app := kingpin.New("write-cluster-msa", "Write MSA for sequence clusters; option to split into core and flexible genomes")
	app.Version("v20201111")
	alnFile := app.Arg("master_MSA", "multi-sequence alignment file containing all sequences in the dataset").Required().String()
	clusterdict := app.Arg("cluster_dict", "hash table from makeCluster.py as a text file, relating cluster # to sequence name").Required().String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	CFsplit := app.Flag("core-flex-split", "If you want to split genomes into core and flexible genes").Default("true").Bool()
	threshold := app.Flag("core-cutoff", "Percentage above which to be considered a core gene").Default("90").Int()
	showProgress := app.Flag("show-progress", "Show progress; would not recommend; significantly slows down").Default("false").Bool()

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

	//load the cluster dict as a golang map
	f, err := os.Open(*clusterdict)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	rd := bufio.NewReader(f)
	var clustermap map[string]string
	clustermap = make(map[string]string)
	var intclusters []int
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		ln := strings.Split(line, ",")
		seq := ln[0]
		cluster := strings.Split(ln[1], "\n")
		clusternum, _ := strconv.Atoi(cluster[0])
		clustermap[seq] = cluster[0]
		intclusters = append(intclusters, clusternum)
		if err == io.EOF {
			break
		}
	}
	clusters := unique(intclusters)

	//make the cluster folders and MSAs
	for _, cluster := range clusters {
		clusterpath := "cluster" + cluster
		if _, err := os.Stat(clusterpath); os.IsNotExist(err) {
			os.Mkdir(clusterpath, 0755)
		}
		MSA := filepath.Join(clusterpath, "MSA_cluster"+cluster)
		f, err := os.Create(MSA)
		check(err)
		f.Close()
	}

	//determine if you're going to split the MSAs into core and flexible genomes
	var t float64
	var seqMap map[string][]string
	if *CFsplit {
		t, seqMap = splitPrep(*threshold, *clusterdict, clusters)
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

	//alnChan := readAlignments(*alnFile)

	numJob := *ncpu
	done := make(chan bool)
	//clusterChan := make(chan clusterAln)
	clusterChan := make(chan cAlignment)
	//percentChan := make(chan genePercents)
	var flex bool
	var CFgenes map[string]string
	var geneFrac map[string]float64

	for i := 0; i < numJob; i++ {
		go func() {
			for aln := range alnChan {
				alnMap := AssembleAlignments(aln, clustermap)
				if *CFsplit {
					flex, CFgenes, geneFrac = coreflexSplit(t, seqMap, alnMap)
					GetGenePercentages(aln.ID, geneFrac, clusters)
				}
				for ID, cAln := range alnMap {
					if *CFsplit {
						clusterChan <- cAlignment{clusterID: ID,
							geneID:    aln.ID,
							genetype:  CFgenes[ID],
							fraction:  geneFrac[ID],
							Sequences: cAln}
					} else {
						clusterChan <- cAlignment{clusterID: ID,
							geneID:    aln.ID,
							genetype:  "n/a",
							fraction:  float64(0),
							Sequences: cAln}
					}

					//clusterChan <- clusterAln{clusterID: ID, Sequences: cAln}
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(clusterChan)
		//defer close(percentChan)
		for i := 0; i < numJob; i++ {
			<-done
		}
	}()

	WriteClusterMSA(clusterChan, *CFsplit, flex)
	//WriteCFMSA(flex, CFgenes, clusterChan)
	//GetGenePercentages(percentChan, clusters)

	duration := time.Since(start)
	fmt.Println("Time to write cluster MSA files:", duration)
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

//channel structure to pump a gene alignment for a cluster into
type cAlignment struct {
	clusterID string         //ID for cluster
	geneID    string         //gene name
	genetype  string         //string saying if core or flex
	fraction  float64        //fraction of strains which have the gene
	Sequences []seq.Sequence //the gene alignment for the cluster
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

//channel structure to pump clusters into
type Cluster struct {
	geneseq   []byte
	header    string
	genename  string
	clusterID string
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
func WriteClusterMSA(clusterChan chan cAlignment, CFsplit bool, flex bool) {

	for c := range clusterChan {
		clusterpath := "cluster" + c.clusterID
		MSA := filepath.Join(clusterpath, "MSA_cluster"+c.clusterID)
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
		if CFsplit {
			WriteCFMSA(c, flex)
		}
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
