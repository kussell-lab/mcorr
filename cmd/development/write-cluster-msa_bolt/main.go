package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"github.com/apsteinberg/biogo/seq"
	"github.com/boltdb/bolt"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
	"io"
	"io/ioutil"
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

	app := kingpin.New("writeclusters", "Write MSA for sequence clusters; option to split into core and flexible genomes")
	app.Version("v20201111")
	alnFile := app.Arg("master_MSA", "multi-sequence alignment file containing all sequences in the dataset").Required().String()
	clusterdict := app.Arg("cluster_dict", "hash table from makeCluster.py relating cluster # to sequence name").Required().String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	CFsplit := app.Flag("core-flex-split", "If you want to split genomes into core and flexible genes").Default("true").Bool()
	bufSize := app.Flag("buf_size", "gene load buf size (memory usage), default 1000 genes").Default("1000").Int()
	numGenes := app.Flag("num_genes", "number of genes per db").Default("1000").Int()
	threshold := app.Flag("core-cutoff", "Percentage above which to be considered a core gene").Default("90").Int()
	showProgress := app.Flag("show-progress", "Show progress; would not recommend; significantly slows down").Default("false").Bool()
	tmpDir := app.Flag("temp_dir", "temp dir").Default(".").String()

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

	//make cluster map, where keys are sequence names and values are cluster numbers
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
	//also make clusters slice, which is a string slice of the cluster names
	clusters := unique(intclusters)

	//set up databases for each cluster

	//make the cluster folders and MSA files
	//for _, cluster := range clusters {
	//	clusterpath := "cluster"+cluster
	//	if _, err := os.Stat(clusterpath); os.IsNotExist(err){
	//			os.Mkdir(clusterpath, os.ModeDir)
	//	}
	//	MSA := filepath.Join(clusterpath,"MSA_cluster"+cluster)
	//	f, err := os.Create(MSA)
	//	check(err)
	//	f.Close()
	//}

	//for if you're going to split the MSAs into core and flexible genomes
	var t float64
	var seqMap map[string][]string
	if *CFsplit {
		t, seqMap = splitPrep(*threshold, *clusterdict)
	}

	//define alignment channel

	alnChan := readAlignments(*alnFile)

	//define a channel to pump cluster alignments into
	clusterChan := make(chan cAlignment)
	//percentChan := make(chan genePercents)
	var flex bool
	var CFgenes map[string]string
	var geneFrac map[string]float64
	var genes []string

	//read through gene by gene and push each gene to a boltdb file for writing later
	//also make a key map of the genes to read through the boltdb files
	go func() {
		defer close(clusterChan)
		for aln := range alnChan {
			genes = append(genes, aln.ID)
			alnMap := AssembleAlignments(aln, clustermap)
			if *CFsplit {
				flex, CFgenes, geneFrac = coreflexSplit(t, seqMap, alnMap)
			}
			for ID, cAln := range alnMap {
				clusterChan <- cAlignment{clusterID: ID,
					geneID:    aln.ID,
					genetype:  CFgenes[ID],
					fraction:  geneFrac[ID],
					Sequences: cAln}
			}
		}
	}()

	jobChan := make(chan *LoadingJob)
	done := make(chan bool)
	for i := 0; i < *ncpu; i++ {
		go func() {
			for {
				dbFile, err := ioutil.TempFile(*tmpDir, fmt.Sprintf("boltdb_%d", i))
				if err != nil {
					panic(err)
				}
				job := &LoadingJob{
					dbFile:  dbFile.Name(),
					bufSize: *bufSize,
				}
				job.Load(clusterChan, *numGenes)
				jobChan <- job
			}
			done <- true
		}()
	}

	go func() {
		defer close(jobChan)
		for i := 0; i < *ncpu; i++ {
			<-done
		}
	}()

	//WriteClusterMSA(clusterChan)
	//WriteCFMSA(flex, CFgenes, clusterChan)
	//GetGenePercentages(percentChan, clusters)

	var jobs []*LoadingJob
	for job := range jobChan {
		jobs = append(jobs, job)
		//defer os.Remove(job.dbFile)
	}

	//make the cluster folders and MSA files
	for _, cluster := range clusters {
		clusterpath := "cluster" + cluster
		if _, err := os.Stat(clusterpath); os.IsNotExist(err) {
			os.Mkdir(clusterpath, os.ModeDir)
		}
		MSA := filepath.Join(clusterpath, "MSA_cluster"+cluster)
		f, err := os.Create(MSA)
		check(err)
		f.Close()
	}
	for cluster := range clusters {
		for gene := range genes {

		}
	}

	duration := time.Since(start)
	fmt.Println(duration)
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

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
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
func WriteClusterMSA(clusterChan chan clusterAln) {

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
	}
}

//assemble gene alignments for each cluster
func AssembleAlignments(aln Alignment, clustermap map[string]string) (alnMap map[string][]seq.Sequence) {
	Sequences := aln.Sequences
	//keys are the cluster names, the values are the
	alnMap = make(map[string][]seq.Sequence)
	for _, s := range Sequences {
		_, _, strain := getNames(s.Id)
		cluster, found := clustermap[strain]
		if found {
			alnMap[cluster] = append(alnMap[cluster], s)
		}
	}
	return
}

//channel structure to pump genes into
//type clusterAln struct {
//	clusterID	string
//	Sequences 	[]seq.Sequence
//}

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

func (job *LoadingJob) Load(clusterChan chan cAlignment, numGenes int) {
	//create boltdb
	job.db = createDB(job.dbFile)
	createBucket(job.db, "cluster_alignment")
	for cAln := range clusterChan {
		addClusterAln(job.db, cAln)
	}
}

func addClusterAln(db *bolt.DB, clusterAln cAlignment) error {
	confBytes, err := json.Marshal(clusterAln)
	if err != nil {
		return fmt.Errorf("could not marshal config json: %v", err)
	}
	cluster_gene := clusterAln.clusterID + "_" + clusterAln.geneID
	err = db.Update(func(tx *bolt.Tx) error {
		err = tx.Bucket([]byte("DB")).Bucket([]byte("cluster_alignment")).Put([]byte(cluster_gene), confBytes)
		if err != nil {
			return fmt.Errorf("could not set config: %v", err)
		}
		return nil
	})
	fmt.Println("Added ClusterAln")
	return err
}

// LoadingJob load genomes into a db.
type LoadingJob struct {
	db      *bolt.DB
	dbFile  string
	bufSize int
}

// createDB creates a bolt db.
func createDB(dbFile string) *bolt.DB {
	db, err := bolt.Open(dbFile, 0600, nil)
	if err != nil {
		log.Fatal(err)
	}
	return db
}

// createBucket creates a bucket.
func createBucket(db *bolt.DB, bucketName string) {
	fn := func(tx *bolt.Tx) error {
		_, err := tx.CreateBucketIfNotExists([]byte(bucketName))
		return err
	}

	err := db.Update(fn)
	if err != nil {
		log.Fatal(err)
	}
}
