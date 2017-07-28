package mcorr

import "math/rand"
import "math"

// Bootstrap for one bootstrapping instance.
type Bootstrap struct {
	ID          string
	sampleRatio float64
	collector   *Collector
	isRandom    bool
}

// NewBootstrap creates a new Boot, given id and sample ratio.
// Sample ratio must be a float64 from 0 to 1.
// By default, bootstrap should do random sampling.
func NewBootstrap(id string, sampleRatio float64) *Bootstrap {
	b := Bootstrap{}
	b.ID = id
	if sampleRatio < 0 {
		sampleRatio = 0
	} else if sampleRatio > 1 {
		sampleRatio = 1
	}
	b.sampleRatio = sampleRatio
	b.collector = NewCollector()
	b.isRandom = true
	return &b
}

// SetRandom set random status
func (b *Bootstrap) SetRandom(r bool) {
	b.isRandom = r
}

// Add add one result into the Bootstrap.
func (b *Bootstrap) Add(results CorrResults) {
	if b.isRandom {
		k := poisson(b.sampleRatio)
		for i := 0; i < k; i++ {
			b.collector.Add(results)
		}
	} else {
		b.collector.Add(results)
	}

}

// Results return final results.
func (b *Bootstrap) Results() (results []CorrResult) {
	return b.collector.Results()
}

func poisson(lambda float64) int {
	L := math.Pow(math.E, -lambda)
	k := 0
	p := 1.0
	for p > L {
		k++
		p *= rand.Float64()
	}
	return k - 1
}
