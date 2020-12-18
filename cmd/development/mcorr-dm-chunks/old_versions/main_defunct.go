package old_versions

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"github.com/sbinet/npyio"
	"gonum.org/v1/gonum/mat"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
)

func main_defunct() {
	app := kingpin.New("mcorr-dm-chunks", "convert mcorr-pair-sync output to distance matrix")
	app.Version("v20201203")
	mpsoutdir := app.Arg("mcorr-pair-sync-dir", "directory containing mcorr-pair-sync output as .csv files").Required().String()
	strains := app.Arg("strain-list", "list of all strains").Required().String()
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

	//make a map of strain names from the strain list
	strainMap := make(map[string]int)
	var strainList []string
	f, err := os.Open(*strains)
	if err != nil {
		log.Fatalln("Couldn't open strain list", err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	j := 0
	for scanner.Scan() {
		strain := scanner.Text()
		strainMap[strain] = j
		strainList = append(strainList, strain)
		j++
	}
	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	// show progress bar
	var bar *pb.ProgressBar
	if *showProgress {
		max := len(strainMap)
		bar = pb.StartNew(max)
		defer bar.Finish()
	}

	var files []string
	//find dem csv files
	err = filepath.Walk(*mpsoutdir, func(path string, info os.FileInfo, err error) error {
		if filepath.Ext(path) == ".csv" {
			files = append(files, path)
		}
		return nil
	})
	if err != nil {
		panic(err)
	}

	//distmatrix := getDistances(strainMap, files...)
	//
	//res, err := os.Create(outFile)
	//if err != nil {
	//	panic(err)
	//}
	//defer res.Close()
	//err = npyio.Write(res, distmatrix)

	numJobs := len(files)
	jobs := make(chan int, numJobs)
	results := make(chan dmResults)

	//set up workers
	for w := 1; w <= *ncpu; w++ {
		go getDists(w, jobs, strainMap, files, results)
	}

	for j := 1; j <= numJobs; j++ {
		jobs <- j
	}
	close(jobs)
	for a := 1; a <= numJobs; a++ {
		<-results
	}

	res, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer res.Close()
	err = npyio.Write(res, results)

	//pairChan := readRecord(*mpsoutdir)
	//
	//numJob := *ncpu
	//done := make(chan bool)
	//dmChan := make(chan dmResults)
	//for i := 0; i < numJob; i++ {
	//	go func() {
	//		count := 0
	//		for pair := range pairChan {
	//			pos_i, _ := strainMap[pair.strain1]
	//			pos_j, _ := strainMap[pair.strain2]
	//			dmChan <- dmResults{pos_i, pos_j, pair.dist}
	//			count = count + 1
	//		}
	//		done <- true
	//	}()
	//}
	//
	//go func() {
	//	defer close(dmChan)
	//	for i := 0; i < numJob; i++ {
	//		<-done
	//	}
	//}()
	//
	//WriteDM(dmChan, strainMap, outFile)
	duration := time.Since(start)
	fmt.Println("Time to convert from csv to distance matrix:", duration)

}

func readRecord(mpsoutdir string) (pairChan chan Pair) {
	pairChan = make(chan Pair)
	read := func() {
		defer close(pairChan)
		var files []string
		//find dem csv files
		err := filepath.Walk(mpsoutdir, func(path string, info os.FileInfo, err error) error {
			if filepath.Ext(path) == ".csv" {
				files = append(files, path)
			}
			return nil
		})
		if err != nil {
			panic(err)
		}
		for _, file := range files {
			f, err := os.Open(file)
			if err != nil {
				panic(err)
			}
			//defer f.Close()
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
			//close the file
			f.Close()
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

func getDistances(strainMap map[string]int, wg *sync.WaitGroup, files ...string) *mat.Dense {
	size := len(strainMap)
	distmatrix := mat.NewDense(size, size, nil)
	defer wg.Done()
	//var wg sync.WaitGroup
	//var m sync.Mutex

	//filesLength := len(files)
	//contents := make(map[string][]byte, filesLength)
	//wg.Add(filesLength)

	for _, file := range files {
		go func(file string) {
			//content, err := ioutil.ReadFile(file)
			f, err := os.Open(file)
			if err != nil {
				log.Fatal(err)
			}
			r := csv.NewReader(f)

			//m.Lock()
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
					pos_i, _ := strainMap[strain1]
					pos_j, _ := strainMap[strain2]
					distmatrix.Set(pos_i, pos_j, dist)
					distmatrix.Set(pos_j, pos_i, dist)
				}

			}
			f.Close()
			//contents[file] = content
			//m.Unlock()
			//wg.Done()
		}(file)
	}

	//wg.Wait()

	return distmatrix
}

func getDists(id int, jobs <-chan int, strainMap map[string]int, files []string, results chan<- dmResults) {
	for j := range jobs {
		fmt.Println("worker", id, "started  job", j)
		file := files[j]
		//fmt.Println("worker", id, "started  job", j)
		f, err := os.Open(file)
		//defer f.Close()
		if err != nil {
			log.Fatal(err)
		}
		r := csv.NewReader(f)

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
				pos_i, _ := strainMap[strain1]
				pos_j, _ := strainMap[strain2]
				results <- dmResults{pos_i: pos_i, pos_j: pos_j, dist: dist}
			}

		}
		fmt.Println("worker", id, "finished job", j)
		f.Close()
	}

	return
}
