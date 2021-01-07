package old

import (
	"bufio"
	"fmt"
	"github.com/apsteinberg/biogo/seq"
	"io"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"
)

func main() {
	dir := "/Volumes/aps_timemachine/recombo/APS158_pangenomealignmentgenerator/write_msa_test"
	xmfa := "test_xmfa"
	xmfa = "1224_properheader"
	alnFile := filepath.Join(dir, xmfa)
	clusterdict := filepath.Join(dir, "clusterlist")
	CFsplit := true
	threshold := 90

	start := time.Now()
	//convert clusterdict to a golang map and a slice of strings
	clustermap, clusters := makeClusterMap(clusterdict)

	//make cluster folders and blank MSA files for writing
	initClusterMSAs(clusters)
	//prepare for splitting into core and flexible genomes
	var t float64
	var seqMap map[string][]string
	if CFsplit {
		t, seqMap = splitPrep(threshold, clusterdict, clusters)
	}
	//define alignment channel
	var alnChan chan Alignment
	alnChan = make(chan Alignment)
	go func() {
		defer close(alnChan)
		count := 0
		c := readAlignments(alnFile)
		for a := range c {
			alnChan <- a
			count++
		}
	}()

	//var flex bool
	var CFgenes map[string]string
	var geneFrac map[string]float64

	for aln := range alnChan {
		alnMap := AssembleAlignments(aln, clustermap)
		if CFsplit {
			_, CFgenes, geneFrac = coreflexSplit(t, seqMap, alnMap, clusters)
			GetGenePercentages(aln.ID, geneFrac, clusters)
			for ID, cAln := range alnMap {
				cluster := cAlignment{clusterID: ID,
					geneID:    aln.ID,
					genetype:  CFgenes[ID],
					fraction:  geneFrac[ID],
					Sequences: cAln}
				WriteClusterMSA(cluster, CFsplit)
			}
		} else {
			for ID, cAln := range alnMap {
				cluster := cAlignment{clusterID: ID,
					geneID:    aln.ID,
					genetype:  "n/a",
					fraction:  float64(0),
					Sequences: cAln}
				WriteClusterMSA(cluster, CFsplit)
			}
		}
	}
	duration := time.Since(start)
	fmt.Println("Time to write cluster MSA files:", duration)
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

//cAlignment is a channel structure to pump a gene alignment for a cluster into
type cAlignment struct {
	clusterID string         //ID for cluster
	geneID    string         //gene name
	genetype  string         //string saying if core or flex
	fraction  float64        //fraction of strains which have the gene
	Sequences []seq.Sequence //the gene alignment for the cluster
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

//makeClusterMap takes the clusterdict and returns both a map relating sequence to cluster ID
//and returns a vector of cluster IDs
func makeClusterMap(clusterdict string) (clustermap map[string]string, clusters []string) {
	f, err := os.Open(clusterdict)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	rd := bufio.NewReader(f)
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
	clusters = unique(intclusters)

	return
}

//unique returns a unique slice of strings from a slice of integers
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

//initClusterMSAs makes cluster folders and blank MSAs for writing
func initClusterMSAs(clusters []string) {
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
}

//check for errors
func check(e error) {
	if e != nil {
		panic(e)
	}
}

//AssembleAlignments maps strain sequences to their respective clusters
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

func getNames(s string) (geneName, genomePos, genomeName string) {
	terms := strings.Split(s, " ")
	geneName = terms[0]
	genomePos = terms[1]
	genomeName = terms[2]
	return
}

//WriteClusterMSA write the MSA file for the cluster
func WriteClusterMSA(c cAlignment, CFsplit bool) {
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
		WriteCFMSA(c)
	}
}
