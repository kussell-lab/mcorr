package main

import (
	"bufio"
	"encoding/csv"
	"errors"
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

func main() {
	app := kingpin.New("mcorr-dm-chunks", "convert mcorr-pair-sync output to distance matrix")
	app.Version("v20201203")
	mpsoutdir := app.Arg("mcorr-pair-sync-dir", "directory containing mcorr-pair-sync output as .csv files").Required().String()
	strains := app.Arg("strain-list", "list of all strains").Required().String()
	out := app.Arg("distmatrix", "prefix for distance matrix output as .npy file").Required().String()
	ncpu := app.Flag("num-workers", "Number of files read at a time (default: 20)").Default("20").Int()
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

	dmMap, err := mapAll(*mpsoutdir, strainMap, *ncpu)

	res, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer res.Close()

	size := len(strainMap)
	distmatrix := mat.NewDense(size, size, nil)
	//dims := size*size
	//count := 0
	fmt.Printf("creating matrix ...")
	//loop over all rows in the matrix
	//set starting column for each row
	col_init := 1
	for i := 0; i < size; i++ {
		//loop over the columns
		for j := col_init; j < size; j++ {
			//get the dist
			//dmkey := strconv.Itoa(i)+strconv.Itoa(j)
			dist := dmMap[pos_key{i, j}]
			//dist := dmMap[dmkey]
			distmatrix.Set(i, j, dist)
			distmatrix.Set(j, i, dist)
		}
		//iterate the starting column
		col_init++
	}

	fmt.Printf("writing matrix ...")
	err = npyio.Write(res, distmatrix)
	duration := time.Since(start)
	fmt.Println("Time to convert from csv to distance matrix:", duration)

}

func getRecord(record []string) (strain1, strain2 string, dist float64) {
	strains := strings.Split(record[5], "_vs_")
	strain1, strain2 = strains[0], strains[1]
	dist, _ = strconv.ParseFloat(record[1], 64)
	return

}

//position key for distance matrix as a map
type pos_key struct {
	pos_i int
	pos_j int
}

// walkFiles starts a goroutine to walk the directory tree at root and send the
// path of each regular file on the string channel.  It sends the result of the
// walk on the error channel.  If done is closed, walkFiles abandons its work.
func walkFiles(done <-chan struct{}, root string) (<-chan string, <-chan error) {
	paths := make(chan string)
	errc := make(chan error, 1)
	go func() {
		// Close the paths channel after Walk returns.
		defer close(paths)
		// No select needed for this send, since errc is buffered.
		errc <- filepath.Walk(root, func(path string, info os.FileInfo, err error) error { // HL
			if filepath.Ext(path) == ".csv" {
				if err != nil {
					return err
				}
				if !info.Mode().IsRegular() {
					return nil
				}
				select {
				case paths <- path:
				case <-done: // HL
					return errors.New("walk canceled")
				}
			}
			return nil
		})
	}()
	return paths, errc
}

// A result is a map of the pairwise distances in the file
type result struct {
	path    string
	fileMap map[pos_key]float64
	err     error
}

// digester reads path names from paths and sends digests of the corresponding
// files on c until either paths or done is closed.
func digester(done <-chan struct{}, paths <-chan string, c chan<- result, strainMap map[string]int, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	fmt.Printf("Worker %d starting\n", id)
	for path := range paths { // HLpaths
		//data, err := ioutil.ReadFile(path)
		fileMap, err := mapFile(path, strainMap)
		select {
		case c <- result{path, fileMap, err}:
		case <-done:
			return
		}
	}
	fmt.Printf("Worker %d done\n", id)
}

//mapFile maps distances to matrix positions for the file
func mapFile(file string, strainMap map[string]int) (map[pos_key]float64, error) {

	f, err := os.Open(file)
	defer f.Close()
	if err != nil {
		log.Fatal(err)
	}
	r := csv.NewReader(f)
	//skip the header
	if _, err := r.Read(); err != nil {
		panic(err)
	}
	//initialize a map for this file
	fileMap := make(map[pos_key]float64)
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
			fileMap[pos_key{pos_i, pos_j}] = dist
			fileMap[pos_key{pos_j, pos_i}] = dist
		}

	}
	return fileMap, err
}

// mapAll reads all the files in the file tree rooted at root and returns a map
// of all the files.  If the directory walk
// fails or any read operation fails, mapAll returns an error.  In that case,
// mapAll does not wait for inflight read operations to complete.
func mapAll(root string, strainMap map[string]int, numDigesters int) (map[pos_key]float64, error) {
	// MD5All closes the done channel when it returns; it may do so before
	// receiving all the values from c and errc.
	done := make(chan struct{})
	defer close(done)

	paths, errc := walkFiles(done, root)

	// Start a fixed number of goroutines to read and digest files.
	c := make(chan result) // HLc
	var wg sync.WaitGroup

	for i := 0; i < numDigesters; i++ {
		wg.Add(1)
		go digester(done, paths, c, strainMap, i, &wg) // HLc
	}
	go func() {
		wg.Wait()
		close(c) // HLc
	}()
	// End of pipeline. it's business time

	dmMap := make(map[pos_key]float64)
	for r := range c {
		if r.err != nil {
			return nil, r.err
		}
		//m[r.path] = r.sum
		//add the values from the fileMap to the distance matrix map
		for k, dist := range r.fileMap {
			dmMap[k] = dist
			//dmkey := strconv.Itoa(k.pos_i) + strconv.Itoa(k.pos_j)
			//dmMap[dmkey] = dist
		}
	}
	// Check whether the Walk failed.
	if err := <-errc; err != nil { // HLerrc
		return nil, err
	}
	return dmMap, nil
}
