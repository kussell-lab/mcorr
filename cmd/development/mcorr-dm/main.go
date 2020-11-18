package main

import (
	"encoding/csv"
	"fmt"
	"github.com/sbinet/npyio"
	"gonum.org/v1/gonum/mat"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
	"io"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
)

func main() {
	app := kingpin.New("mcorr-dm", "convert mcorr-pair output to distance matrix")
	app.Version("v20201108")
	mcp := app.Arg("mcorr-pair_csv", "csv output from mcorr-pair as .csv file").Required().String()
	out := app.Arg("distmatrix", "prefix for distance matrix output as .npy file").Required().String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	showProgress := app.Flag("show-progress", "Show progress").Default("false").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	//timer

	start := time.Now()

	if *ncpu == 0 {
		*ncpu = runtime.NumCPU()
	}

	runtime.GOMAXPROCS(*ncpu)
	//outfile for distance matrix
	outFile := *out + ".npy"

	//make a map of strain names
	strainMap := make(map[string]int)
	var strainList []string
	f, err := os.Open(*mcp)
	if err != nil {
		log.Fatalln("Couldn't open csv file", err)
	}
	defer f.Close()
	r := csv.NewReader(f)
	//read the header
	if _, err := r.Read(); err != nil {
		panic(err)
	}
	//read a record and get a strain name
	j := 0
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		if record[5] == "all" {
			continue
		} else {
			record = strings.Split(record[5], "_vs_")
			if _, found := strainMap[record[0]]; found == false {
				strainMap[record[0]] = j
				j = j + 1
				strainList = append(strainList, record[0])
			}
			if _, found := strainMap[record[1]]; found == false {
				strainMap[record[1]] = j
				j = j + 1
				strainList = append(strainList, record[1])
			}
		}
	}

	//write strain list to file
	WriteStrainList(strainList)
	// show progress bar
	var bar *pb.ProgressBar
	if *showProgress {
		max := len(strainMap)
		bar = pb.StartNew(max)
		defer bar.Finish()
	}

	pairChan := readRecord(*mcp)

	numJob := *ncpu
	done := make(chan bool)
	dmChan := make(chan dmResults)
	for i := 0; i < numJob; i++ {
		go func() {
			count := 0
			for pair := range pairChan {
				pos_i, _ := strainMap[pair.strain1]
				pos_j, _ := strainMap[pair.strain2]
				dmChan <- dmResults{pos_i, pos_j, pair.dist}
				count = count + 1
				//ok := isIntegral(float64(count / 1000))
				//if ok {
				//	fmt.Printf("Added %c", count)
				//}
				//bar.Add(1)
			}
			done <- true
		}()
	}

	go func() {
		defer close(dmChan)
		for i := 0; i < numJob; i++ {
			<-done
		}
	}()

	WriteDM(dmChan, strainMap, outFile)
	duration := time.Since(start)
	fmt.Println("Time to convert from csv to distance matrix:", duration)

}

func readRecord(file string) (pairChan chan Pair) {
	pairChan = make(chan Pair)
	read := func() {
		defer close(pairChan)
		f, err := os.Open(file)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		r := csv.NewReader(f)
		//read the header
		if _, err := r.Read(); err != nil {
			panic(err)
		}
		for {
			record, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Fatal(err)
			}
			if record[5] == "all" {
				continue
			} else {
				strain1, strain2, dist := getRecord(record)
				pairChan <- Pair{strain1: strain1, strain2: strain2, dist: dist}
			}
		}
	}
	go read()
	return
}

type Pair struct {
	strain1 string
	strain2 string
	dist    float64
}

func getRecord(record []string) (strain1, strain2 string, dist float64) {
	strains := strings.Split(record[5], "_vs_")
	strain1, strain2 = strains[0], strains[1]
	dist, _ = strconv.ParseFloat(record[1], 64)
	return

}

//channel to pump distance matrix into
type dmResults struct {
	pos_i int
	pos_j int
	dist  float64
}

//write the distance matrix to a npy file
func WriteDM(dmChan chan dmResults, strainMap map[string]int, outFile string) {
	size := len(strainMap)
	distmatrix := mat.NewDense(size, size, nil)
	for dm := range dmChan {
		distmatrix.Set(dm.pos_i, dm.pos_j, dm.dist)
		distmatrix.Set(dm.pos_j, dm.pos_i, dm.dist)

	}
	f, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	err = npyio.Write(f, distmatrix)
}

//write the strain names to a npy array
func WriteStrainList(strainList []string) {
	w, err := os.Create("strains")
	if err != nil {
		panic(err)
	}
	defer w.Close()
	for _, name := range strainList {
		w.WriteString(name + "\n")
	}

}

func isIntegral(val float64) bool {
	return val == float64(int(val))
}
