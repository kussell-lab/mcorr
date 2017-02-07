package main

import (
	"encoding/json"
	"fmt"
	"github.com/alecthomas/kingpin"
	"io"
	"os"
)

func main() {
	app := kingpin.New("mcorr_collect", "Collect correlation results")
	app.Version("v0.1")

	jsonFile := app.Arg("json_file", "Input JSON file").Required().String()
	outFile := app.Arg("csv_file", "Output file in CSV file").Required().String()
	numBoot := app.Flag("num_boot", "Number of bootstrapping on genes").Default("1000").Int()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	corrResChan := readJSON(*jsonFile)
	// prepare bootstrappers.
	bootstraps := []*Bootstrap{}
	notBootstrap := NewBootstrap("all", 1.0)
	notBootstrap.SetRandom(false)
	bootstraps = append(bootstraps, notBootstrap)
	for i := 0; i < *numBoot; i++ {
		id := fmt.Sprintf("boot_%d", i)
		sampleRatio := 1.0
		bootstraps = append(bootstraps, NewBootstrap(id, sampleRatio))
	}

	for corrResults := range corrResChan {
		for _, bs := range bootstraps {
			bs.Add(corrResults)
		}
	}

	w, err := os.Create(*outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("l,m,v,n,t,b\n")
	for _, bs := range bootstraps {
		results := bs.Results()
		for _, res := range results {
			w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, bs.ID))
		}
	}
}

// CorrResult stores a correlation result.
type CorrResult struct {
	Lag      int
	Mean     float64
	Variance float64
	N        int
	Type     string
}

// CorrResults stores a list of CorrResult with an gene ID.
type CorrResults struct {
	ID      string
	Results []CorrResult
}

func readJSON(file string) chan CorrResults {
	c := make(chan CorrResults)
	go func() {
		defer close(c)
		f, err := os.Open(file)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		decoder := json.NewDecoder(f)
		for {
			var corrResults CorrResults
			if err := decoder.Decode(&corrResults); err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			c <- corrResults
		}
	}()
	return c
}
