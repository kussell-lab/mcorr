package mcorr

import (
	"fmt"
	"os"
)

// CollectWrite collects and writes the correlation results.
func CollectWrite(corrResChan chan CorrResults, outFile string, numBoot int) {
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

	w, err := os.Create(outFile)
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
