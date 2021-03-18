package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"gopkg.in/alecthomas/kingpin.v2"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"time"
)

func main() {
	app := kingpin.New("McorrJson2Csv", "converts gene json files to csv files which can be easily used with the python pandas package; need to divide output by d_sample to get corr profiles")
	app.Version("v20210218")
	currentDir, _ := os.Getwd()
	root := app.Flag("root", "folder containing single gene json file").Default(currentDir).String()
	geneJson := app.Arg("jsonName", "name of json file containing gene correlation profiles").Required().String()
	geneCsv := app.Arg("outputCsv", "name of output csv containing gene correlation profiles").Required().String()
	ncpu := app.Flag("numCpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	if *ncpu <= 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)
	//inputs for testing
	//define root directory
	//root := "/Volumes/GoogleDrive/My Drive/hpc/recombo/APS166_BF_Archive/E003188_34.0.periodic.trimmo.60.um"
	//geneJson := "E003188_34.0.periodic.trimmo.60.um.json"
	//geneCsv := "0218test.csv"

	//timer
	start := time.Now()

	jsonPath := filepath.Join(*root, *geneJson)
	j, err := os.Open(jsonPath)
	if err != nil {
		fmt.Println(err)
	}
	defer j.Close()
	//initialize the output csv
	initCsvOut(*root, *geneCsv)
	//initialize the Genes array
	var genes CorrResults
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
			write2Out(genes, *root, *geneCsv)
			numGenes++
		}
		if err != nil {
			if err != io.EOF {
				log.Fatalf("Error when reading file %s: %v", j, err)
			}
			break
		}
	}
	duration := time.Since(start)
	fmt.Println("Conversion time:", duration)
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

func initCsvOut(root string, outName string) {
	path := filepath.Join(root, outName)
	recordFile, err := os.Create(path)
	if err != nil {
		fmt.Println("Error while creating the fail output ::", err)
		return
	}
	defer recordFile.Close()
	// Initialize the writer
	writer := csv.NewWriter(recordFile)
	defer writer.Flush()
	//write header
	header := []string{"lag", "mean", "N", "variance", "type", "gene"}
	err = writer.Write(header)
	if err != nil {
		fmt.Println("Error while writing header for failure output ::", err)
		return
	}
}

func write2Out(genes CorrResults, root string, outName string) {
	path := filepath.Join(root, outName)
	recordFile, err := os.OpenFile(path, os.O_APPEND|os.O_WRONLY, 0644)
	if err != nil {
		fmt.Println("Error while writing to fail output ::", err)
		return
	}
	defer recordFile.Close()
	// Initialize the writer
	writer := csv.NewWriter(recordFile)
	defer writer.Flush()
	// Write all the records
	for i := 0; i < len(genes.Results); i++ {
		var line []string
		var Gene string
		l := fmt.Sprintf("%d", genes.Results[i].Lag)
		m := fmt.Sprintf("%e", genes.Results[i].Mean)
		N := fmt.Sprintf("%d", genes.Results[i].N)
		v := fmt.Sprintf("%f", genes.Results[i].Variance)
		Type := genes.Results[i].Type
		terms := strings.Split(genes.ID, "|")
		if len(terms) > 1 {
			Gene = terms[1]
		} else {
			Gene = terms[0]
		}

		line = []string{l, m, N, v, Type, Gene}
		err = writer.Write(line)
		if err != nil {
			fmt.Println("Error while writing to the file ::", err)
			return
		}
	}

}
