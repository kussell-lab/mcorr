package main

import (
	"encoding/json"
	"io"
	"os"

	"github.com/alecthomas/kingpin"
	"github.com/mingzhi/mcorr"
)

func main() {
	app := kingpin.New("mcorr-collect", "Collect results.")
	app.Version("v0.1")

	alnFile := app.Arg("in", "json file").Required().String()
	outFile := app.Arg("out", "Output file in CSV format.").Required().String()
	numBoot := app.Flag("num-boot", "Number of bootstrapping on genes").Default("1000").Int()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	corrResChan := readCorrRes(*alnFile)

	mcorr.CollectWrite(corrResChan, *outFile, *numBoot)
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
