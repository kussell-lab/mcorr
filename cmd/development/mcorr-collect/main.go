package main

import (
	"encoding/json"
	"fmt"
	"io"
	"os"
	"sort"

	"gopkg.in/alecthomas/kingpin.v2"
	"github.com/kussell-lab/mcorr"
)

func main() {
	app := kingpin.New("mcorr-collect", "Collect results.")
	app.Version("v0.1")

	alnFile := app.Arg("in", "json file").Required().String()
	outFile := app.Arg("out", "Output file in CSV format.").Required().String()
	numBoot := app.Flag("num-boot", "Number of bootstrapping on genes").Default("1000").Int()
	corrType := app.Flag("corr-type", "correlation type").Default("P4").String()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	resChan := readCorrRes(*alnFile)
	if *corrType == "P2" {
		mcorr.CollectWrite(resChan, *outFile, *numBoot)
	} else {
		bootstraps := mcorr.Collect(resChan, *numBoot)

		w, err := os.Create(*outFile)
		if err != nil {
			panic(err)
		}
		defer w.Close()

		w.WriteString("l,m,v,n,t,b\n")
		for _, bs := range bootstraps {
			results := bs.Results()
			qfactor := getQfactor(results)
			for _, res := range results {
				if res.Type == "Ks" || (res.Type == "P4" && res.Lag > 0) {
					if res.Type == "P4" {
						res.Mean *= qfactor
						res.Variance *= qfactor * qfactor
						res.Type = "P2"
					}
					w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n",
						res.Lag, res.Mean, res.Variance, res.N, res.Type, bs.ID))
				}
			}
		}
	}

}

// readCorrRes return a channel of CorrRes
func readCorrRes(filename string) chan mcorr.CorrResults {
	c := make(chan mcorr.CorrResults)
	go func() {
		defer close(c)
		f, err := os.Open(filename)
		if err != nil {
			panic(err)
		}
		defer f.Close()

		decoder := json.NewDecoder(f)
		for {
			var corrResults mcorr.CorrResults
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

// getQfactor return the q factor between p2 and p4.
func getQfactor(results []mcorr.CorrResult) float64 {
	p2values := make([]float64, 31)
	p4values := make([]float64, 31)
	for _, res := range results {
		if res.Lag <= 30 && res.Lag > 0 {
			if res.Type == "P2" {
				p2values[res.Lag] = res.Mean
			} else if res.Type == "P4" {
				p4values[res.Lag] = res.Mean
			}
		}
	}

	var factors []float64
	for i := range p2values {
		if p2values[i] > 0 && p4values[i] > 0 {
			factors = append(factors, p2values[i]/p4values[i])
		}
	}

	if len(factors) == 0 {
		return 0
	}

	sort.Float64s(factors)
	if len(factors)%2 == 0 {
		return (factors[len(factors)/2] + factors[len(factors)/2-1]) / 2
	}
	return (factors[len(factors)/2])
}
