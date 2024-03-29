package mcorr

import (
	"encoding/json"
	"fmt"
	"os"
)

// PipeOutCorrResults pipe the the channel of CorrResults out to a file.
func PipeOutCorrResults(corrResChan chan CorrResults, outFile string) chan CorrResults {
	c := make(chan CorrResults)
	go func() {
		defer close(c)
		f, err := os.Create(outFile)
		if err != nil {
			panic(err)
		}
		defer f.Close()

		encoder := json.NewEncoder(f)
		for res := range corrResChan {
			if err := encoder.Encode(res); err != nil {
				panic(err)
			}
			c <- res
		}
	}()
	return c
}

// Collect feed correlation results into boostrappers and return them.
func Collect(corrResChan chan CorrResults, numBoot int) []*Bootstrap {
	// prepare bootstrappers.
	bootstraps := []*Bootstrap{}
	notBootstrap := NewBootstrap("all", 1.0)
	notBootstrap.SetRandom(false)
	bootstraps = append(bootstraps, notBootstrap)
	for i := 0; i < numBoot; i++ {
		id := fmt.Sprintf("boot_%d", i)
		sampleRatio := 1.0
		bootstraps = append(bootstraps, NewBootstrap(id, sampleRatio))
	}

	for corrResults := range corrResChan {
		for _, bs := range bootstraps {
			bs.Add(corrResults)
		}
	}
	return bootstraps
}

// CollectWrite collects and writes the correlation results.
func CollectWrite(corrResChan chan CorrResults, outFile string, numBoot int) {
	bootstraps := Collect(corrResChan, numBoot)

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("# l: the distance between two genomic positions\n")
	w.WriteString("# m: the mean value of correlatio profile\n")
	w.WriteString("# v: the variance of correlation profile\n")
	w.WriteString("# n: the total number of alignments used for calculation\n")
	w.WriteString("# t: the type of result: Ks is for d_sample, and P2 is for correlation profile\n")
	w.WriteString("# b: the bootstrap number (all means used all alignments).\n")

	w.WriteString("l,m,v,n,t,b\n")
	for _, bs := range bootstraps {
		//this is where division by P(l = 0) happens!!!!
		results := bs.Results()
		for _, res := range results {
			w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, bs.ID))
		}
	}
}
