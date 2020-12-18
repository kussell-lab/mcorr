package old_versions

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"github.com/boltdb/bolt"
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

	var wg sync.WaitGroup
	var m sync.Mutex
	//create a map to store distance matrix values in
	//size := len(strainMap)
	//initialize boltdb databases
	//boltsize := int(size/4)
	//for s := 0; s <= 4; s++ {
	//	first := s*boltsize
	//	last := (s+1)*boltsize-1
	//	start := strconv.Itoa(first)
	//	end := strconv.Itoa(last)
	//	setupDB(start+end)
	//}
	//for s := 0; s < 4; s++ {
	//	boltname := strconv.Itoa(s)
	//	setupDB(boltname)
	//}
	//setup boltdb for distance matrix
	setupDB("distancematrix")
	dmMap := make(map[pos_key]float64)

	for i := 1; i < len(files); i++ {
		wg.Add(1)
		go csvtomap(i, &wg, &m, files[i], strainMap, dmMap)
	}

	wg.Wait()

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
			dist := dmMap[pos_key{i, j}]
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

//another equally valid way to do this
func getDistances(strainMap map[string]int, files ...string) map[pos_key]float64 {

	var wg sync.WaitGroup
	var m sync.Mutex

	filesLength := len(files)
	//contents := make(map[string][]byte, filesLength)
	dmMap := make(map[pos_key]float64)
	//basically the wait group will be released once we've gone through all the files
	wg.Add(filesLength)

	for _, file := range files {
		go func(file string) {
			f, err := os.Open(file)
			defer f.Close()
			if err != nil {
				log.Fatal(err)
			}
			r := csv.NewReader(f)

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
			//f.Close()
			//lock the map so only one thread at a time can add to the distance matrix map
			m.Lock()
			//add the values from the fileMap to the distance matrix map
			for k, dist := range fileMap {
				dmMap[k] = dist
			}
			m.Unlock()
			wg.Done()
		}(file)
	}

	wg.Wait()

	return dmMap
}

//position key for distance matrix as a map
type pos_key struct {
	Pos_i int
	Pos_j int
}

func csvtomap(id int, wg *sync.WaitGroup, m *sync.Mutex, file string, strainMap map[string]int, dmMap map[pos_key]float64) {
	defer wg.Done()
	fmt.Printf("Worker %d starting\n", id)
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
	//quartersize := int(len(strainMap)/4)
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
			//position := pos_key{pos_i, pos_j}
			fileMap[pos_key{pos_i, pos_j}] = dist
			//fileMap[pos_key{pos_j, pos_i}] = dist
		}

	}
	//lock the map so only one thread at a time can add to the distance matrix map
	m.Lock()
	//add the values from the fileMap to the distance matrix map
	db, err := bolt.Open("distancematrix.db", 0600, nil)
	if err != nil {
		log.Fatal(err)
	}
	defer db.Close()
	for k, dist := range fileMap {
		dmMap[k] = dist
		//quarter := checkquartile(k.Pos_i, quartersize)
		addEntry(db, k, dist)
	}
	m.Unlock()

	fmt.Printf("Worker %d done\n", id)
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

func setupDB(dbFile string) *bolt.DB {
	db, err := bolt.Open(dbFile+".db", 0600, nil)
	if err != nil {
		log.Fatal(err)
	}
	err = db.Update(func(tx *bolt.Tx) error {
		root, err := tx.CreateBucketIfNotExists([]byte("DB"))
		if err != nil {
			return fmt.Errorf("could not create root bucket: %v", err)
		}
		_, err = root.CreateBucketIfNotExists([]byte("ENTRIES"))
		if err != nil {
			return fmt.Errorf("could not create days bucket: %v", err)
		}
		return nil
	})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("DB Setup Done")
	return db
}

func addEntry(db *bolt.DB, position pos_key, distance float64) error {
	posBytes, err := json.Marshal(position)
	if err != nil {
		return fmt.Errorf("could not marshal position json: %v", err)
	}
	distBytes, err := json.Marshal(distance)
	if err != nil {
		return fmt.Errorf("could not marshal distance json: %v", err)
	}

	err = db.Update(func(tx *bolt.Tx) error {
		err := tx.Bucket([]byte("DB")).Bucket([]byte("ENTRIES")).Put(posBytes, distBytes)
		if err != nil {
			return fmt.Errorf("could not insert entry: %v", err)
		}
		return nil
	})
	return err
}

func checkquartile(pos int, quarter int) (quartile string) {
	//quarter := int(size/4)
	if pos < quarter {
		quartile = strconv.Itoa(0)
	} else if pos < 2*quarter {
		quartile = strconv.Itoa(1)
	} else if pos < 3*quarter {
		quartile = strconv.Itoa(2)
	} else {
		quartile = strconv.Itoa(3)
	}
	return
}
