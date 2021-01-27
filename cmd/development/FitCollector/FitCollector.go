package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"gonum.org/v1/gonum/stat"
	"gopkg.in/alecthomas/kingpin.v2"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"
)

func main() {

	app := kingpin.New("FitCollector", "Collect results from mcorr-fit or fitCorr for many sequence clusters")
	app.Version("v20210126")
	root := app.Arg("root", "root directory containing cluster out folders").Required().String()
	jsonSuffix := app.Flag("jsonSuffix", "suffix for json files").Default("_XMFA_OUT.json").String()
	lmfitSuffix := app.Flag("lmfitSuffix", "suffix for lmfit report output file").Default("_XMFA_OUT_lmfit_report.csv").String()
	ncpu := app.Flag("numCpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	numDigesters := app.Flag("numDigesters", "Number threads fetching results").Default("20").Int()
	out := app.Flag("outCsv", "prefix for output csv").Default("mcorrfit_results").String()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	if *ncpu <= 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	//timer
	start := time.Now()

	done := make(chan struct{})
	defer close(done)

	//root := "/Volumes/aps_timemachine/recombo/APS160.5_lmfit"
	//jsonSuffix := "_XMFA_OUT.json"
	//lmfitSuffix := "_XMFA_OUT_lmfit_report.csv"
	//numDigesters := 4
	//out := "blarrrrgh"
	clusterDirs := makeDirList(*root)
	//make a channel for cluster directories which closes when we're out of them
	clusters, fitFailed := clusterFileChan(done, *root, clusterDirs, *jsonSuffix, *lmfitSuffix)
	//start a fixed number of goroutines to send results on
	//make a results channel
	resChan := make(chan result)
	var wg sync.WaitGroup
	for i := 0; i < *numDigesters; i++ {
		wg.Add(1)
		go digester(done, clusters, resChan, i, &wg)
	}
	go func() {
		wg.Wait()
		close(resChan)
	}()
	//end of pipeline
	currentTime := time.Now()
	timeStamp := fmt.Sprintf(currentTime.Format("060102_1504"))
	outCsv := timeStamp + "_" + *out + ".csv"
	writeCSV(resChan, *root, outCsv)
	failOut := timeStamp + "_lmfitfailed.csv"
	writeFitFails(fitFailed, *root, failOut)

	duration := time.Since(start)
	fmt.Println("Time to collect results:", duration)

}

//clusterFiles is a struct for filepaths for all files associated with a cluster or pair of clusters
type clusterFiles struct {
	ID       string //cluster ID
	genome   string //type of genome
	json     string //json with gene corr profiles
	lmfitOut string //lmfit output csv
}

// result struct
type result struct {
	ID       string   //cluster ID
	genome   string   //type of genome
	numPairs float64  //average number of pairs
	StDev    float64  //St Dev of number of pairs
	numGenes int      //number of genes
	fitOut   []string //values from fitting
}

//Genes struct which contains
//
type Genes struct {
	Genes []CorrResults //'json:gene'
}

// CorrResults stores a list of CorrResult with an gene ID.
type CorrResults struct {
	ID      string
	Results []CorrResult
}

// CorrResult stores a correlation result.
type CorrResult struct {
	Lag      int
	Mean     float64
	Variance float64
	N        int
	Type     string
}

// makeDirList make a list of cluster directories
func makeDirList(root string) (DirList []string) {
	c, err := ioutil.ReadDir(root)
	if err != nil {
		panic(err)
	}
	for _, entry := range c {
		if entry.IsDir() && strings.HasPrefix(entry.Name(), "cluster") {
			DirList = append(DirList, entry.Name())
		}
	}
	return DirList
}

//clusterFileChan returns a channel of clusterFiles and a list of clusters where lmfit failed
func clusterFileChan(done <-chan struct{}, root string, DirList []string, jsonSuffix string, lmfitSuffix string) (<-chan clusterFiles, []clusterFiles) {
	clusterFileChan := make(chan clusterFiles)
	var lmfitFailed []clusterFiles
	go func() {
		defer close(clusterFileChan)
		for _, d := range DirList {
			//define the flex files
			core, flex := makeClusterFiles(root, d, jsonSuffix, lmfitSuffix)
			var cGenomes []clusterFiles
			if checkLmfit(root, core) {
				cGenomes = append(cGenomes, core)
			} else {
				lmfitFailed = append(lmfitFailed, core)
			}
			if checkLmfit(root, flex) {
				cGenomes = append(cGenomes, flex)
			} else {
				lmfitFailed = append(lmfitFailed, flex)
			}
			//cGenomes := []clusterFiles{core, flex}
			for _, c := range cGenomes {
				select {
				case clusterFileChan <- c:
				case <-done:
					return
				}
			}
		}
	}()
	return clusterFileChan, lmfitFailed
}

//makeClusterFiles returns ClusterFiles for core and flexible genomes
func makeClusterFiles(root string, d string, jsonSuffix string, lmfitSuffix string) (core clusterFiles, flex clusterFiles) {
	//define the flex files
	flexJson := d + "_FLEX" + jsonSuffix
	json := filepath.Join(root, d, flexJson)
	flexLmfit := d + "_FLEX" + lmfitSuffix
	lmfit := filepath.Join(root, d, flexLmfit)
	flex = clusterFiles{d, "FLEX", json, lmfit}
	//define the core files
	coreJson := d + "_CORE" + jsonSuffix
	json = filepath.Join(root, d, coreJson)
	coreLmfit := d + "_CORE" + lmfitSuffix
	lmfit = filepath.Join(root, d, coreLmfit)
	core = clusterFiles{d, "CORE", json, lmfit}
	return core, flex
}

//checkLmfit check to see if lmfit completed for the cluster
func checkLmfit(root string, c clusterFiles) bool {
	lmfit := filepath.Join(root, c.lmfitOut)
	_, err := os.Stat(lmfit)
	if os.IsNotExist(err) {
		return false
	} else {
		return true
	}
}

func getNumPairs(cluster clusterFiles) (avgPairs float64) {
	j, err := os.Open(cluster.json)
	if err != nil {
		fmt.Println(err)
	}
	defer j.Close()

	//initialize the Genes array
	var genes CorrResults
	//initialize the number of pairs
	N := 0
	//initialize the number of genes
	numGenes := 0
	r := bufio.NewReader(j)
	for {
		line, err := r.ReadString('\n')
		// read our opened xmlFile as a byte array.
		byteValue := []byte(line)
		//unmarshall into genes
		json.Unmarshal(byteValue, &genes)
		//get N for l = 0
		if len(genes.Results) > 0 {
			//fmt.Println("N: ", genes.Results[0].N)
			N = N + genes.Results[0].N
			numGenes++
			//if N != 0{
			//	numGenes++
			//}
		}
		if err != nil {
			if err != io.EOF {
				log.Fatalf("Error when reading file %s: %v", j, err)
			}
			break
		}
	}
	avgPairs = float64(N) / float64(numGenes)
	return avgPairs
}

func calcMeanStDev(cluster clusterFiles) (MeanVariance []float64, numGenes int) {
	j, err := os.Open(cluster.json)
	if err != nil {
		fmt.Println(err)
	}
	defer j.Close()

	//initialize the Genes array
	var genes CorrResults
	//initialize the num of pairs array
	var pairs []float64
	r := bufio.NewReader(j)
	for {
		line, err := r.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				log.Fatalf("Error when reading file %s: %v", j, err)
			}
			break
		}
		// read our opened xmlFile as a byte array.
		byteValue := []byte(line)
		//unmarshall into genes
		json.Unmarshal(byteValue, &genes)
		//get N for l = 0
		if len(genes.Results) > 0 {
			pairs = append(pairs, float64(genes.Results[0].N))
		}
	}
	m, v := stat.MeanStdDev(pairs, nil)
	MeanVariance = []float64{m, v}
	numGenes = len(pairs)
	return MeanVariance, numGenes
}

//digester aligns genomes and writes them to FASTA files until either SRAchan or done is closed
func digester(done <-chan struct{}, FileChan <-chan clusterFiles, c chan<- result, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	fmt.Printf("Digester %d starting\n", id)
	for cluster := range FileChan {
		MeanStDev, numGenes := calcMeanStDev(cluster)
		fitOut := getFitOut(cluster)
		select {
		case c <- result{cluster.ID, cluster.genome,
			MeanStDev[0], MeanStDev[1], numGenes, fitOut}:
		case <-done:
			return
		}
	}
	fmt.Printf("Digester %d done\n", id)
}

func getFitOut(cluster clusterFiles) (fitOut []string) {
	l, err := os.Open(cluster.lmfitOut)
	if err != nil {
		fmt.Println(err)
	}
	defer l.Close()

	reader := csv.NewReader(l)
	reader.FieldsPerRecord = -1
	i := 0
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		if i == 2 {
			fitOut = append(fitOut, record[1])
		}
		if i == 10 {
			fitOut = append(fitOut, record...)
			break
		}
		i++
	}
	return fitOut
}

func writeCSV(resChan chan result, root string, outName string) {
	path := filepath.Join(root, outName)
	recordFile, err := os.Create(path)
	if err != nil {
		fmt.Println("Error while creating the output csv ::", err)
		return
	}
	defer recordFile.Close()
	// Initialize the writer
	writer := csv.NewWriter(recordFile)
	defer writer.Flush()
	//write header
	header := []string{"ID", "genome", "avg_pairs", "stdev_pairs", "num_genes",
		"bp_analyzed", "ds", "thetaS", "f", "phiS", "thetaP",
		"phiP", "c", "dp", "dc", "chisq", "red-chisq"}
	err = writer.Write(header)
	for r := range resChan {
		// Write all the records
		var line []string
		mean := fmt.Sprintf("%f", r.numPairs)
		StDev := fmt.Sprintf("%f", r.StDev)
		numGenes := fmt.Sprintf("%d", r.numGenes)
		line = append(line, r.ID, r.genome, mean, StDev, numGenes)
		line = append(line, r.fitOut...)
		err = writer.Write(line)
		if err != nil {
			fmt.Println("Error while writing to the file ::", err)
			return
		}
	}
}

func writeFitFails(fitFailed []clusterFiles, root string, outName string) {
	path := filepath.Join(root, outName)
	recordFile, err := os.Create(path)
	if err != nil {
		fmt.Println("Error while creating the output csv ::", err)
		return
	}
	defer recordFile.Close()
	// Initialize the writer
	writer := csv.NewWriter(recordFile)
	defer writer.Flush()
	//write header
	header := []string{"ID", "genome"}
	err = writer.Write(header)
	for _, c := range fitFailed {
		// Write all the records
		var line []string
		line = append(line, c.ID, c.genome)
		err = writer.Write(line)
		if err != nil {
			fmt.Println("Error while writing to the file ::", err)
			return
		}
	}
}
