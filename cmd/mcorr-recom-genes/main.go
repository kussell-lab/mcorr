package main

import (
	"encoding/json"
	"fmt"
	"io"
	"math"
	"os"
	"sort"

	"github.com/alecthomas/kingpin"
	"github.com/mingzhi/mcorr"
)

func main() {
	app := kingpin.New("mcorr-collect", "Collect results.")
	app.Version("v0.1")

	alnFile := app.Arg("in", "json file").Required().String()
	outFile := app.Arg("out", "Output file in CSV format.").Required().String()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	corrResChan := readCorrRes(*alnFile)

	var corrResList []mcorr.CorrResults
	for res := range corrResChan {
		corrResList = append(corrResList, res)
	}

	ksMean := Mean{}
	p2Means := &MeanList{}
	p4Means := &MeanList{}

	for _, corrResults := range corrResList {
		for _, res := range corrResults.Results {
			if res.Type == "P2" && res.Lag == 0 {
				ksMean.Value += res.Mean
				ksMean.Count++
			}
			if res.Type == "P2" && res.Lag > 0 {
				p2Means.Add(res.Lag/3, res.Mean)
			}
			if res.Type == "P4" && res.Lag > 0 {
				p4Means.Add(res.Lag/3, res.Mean)
			}
		}
	}
	rawP2values := estimateP2(p2Means, p4Means, ksMean)
	fmt.Println(rawP2values)
	w, err := os.Create(*outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()
	w.WriteString("id,r\n")
	for len(corrResList) > 0 {
		minResidual := math.Inf(1)
		minIdx := 0
		for i, corrResults := range corrResList {
			for _, res := range corrResults.Results {
				if res.Type == "P2" && res.Lag == 0 {
					ksMean.Value -= res.Mean
					ksMean.Count--
				}
				if res.Type == "P2" && res.Lag > 0 {
					p2Means.Remove(res.Lag/3, res.Mean)
				}
				if res.Type == "P4" && res.Lag > 0 {
					p4Means.Remove(res.Lag/3, res.Mean)
				}
			}
			p2values := estimateP2(p2Means, p4Means, ksMean)
			sum := 0.0
			for i := 1; i < 50; i++ {
				r := p2values[i] - rawP2values[i]
				sum += r * r
			}
			r := math.Sqrt(sum / 50.0)
			if r < minResidual {
				minIdx = i
				minResidual = r
			}
			for _, res := range corrResults.Results {
				if res.Type == "P2" && res.Lag == 0 {
					ksMean.Value += res.Mean
					ksMean.Count++
				}
				if res.Type == "P2" && res.Lag > 0 {
					p2Means.Add(res.Lag/3, res.Mean)
				}
				if res.Type == "P4" && res.Lag > 0 {
					p4Means.Add(res.Lag/3, res.Mean)
				}
			}
		}
		w.WriteString(fmt.Sprintf("%s,%g\n", corrResList[minIdx].ID, minResidual))
		temp := corrResList[:minIdx]
		temp = append(temp, corrResList[minIdx+1:]...)
		corrResList = temp
		// fmt.Println(len(corrResList))
	}

	// mcorr.CollectWrite(corrResChan, *outFile, *numBoot)
}

// Mean is a structure for store mean value.
type Mean struct {
	Value float64
	Count float64
}

// Mean return the mean.
func (m Mean) Mean() float64 {
	return m.Value / m.Count
}

// MeanList is a list of means.
type MeanList struct {
	Means []Mean
}

// Add adds a value at the idx.
func (m *MeanList) Add(idx int, val float64) {
	for len(m.Means) <= idx {
		m.Means = append(m.Means, Mean{})
	}
	m.Means[idx].Value += val
	m.Means[idx].Count++
}

// Remove removes a value at the idx.
func (m *MeanList) Remove(idx int, val float64) {
	m.Means[idx].Value -= val
	m.Means[idx].Count--
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
func getQfactor(p2values, p4values []float64) float64 {
	var factors []float64
	for i := range p2values {
		if p2values[i] > 0 && p4values[i] > 0 {
			factors = append(factors, p2values[i]/p4values[i])
		}
	}
	sort.Float64s(factors)
	if len(factors)%2 == 0 {
		return (factors[len(factors)/2] + factors[len(factors)/2-1]) / 2
	}
	return (factors[len(factors)/2])
}

func estimateP2(p2means, p4means *MeanList, ksmean Mean) []float64 {
	var p2values []float64
	for _, m := range p2means.Means {
		p2 := m.Value / m.Count
		p2values = append(p2values, p2)
	}

	var p4values []float64
	for _, m := range p4means.Means {
		p4 := m.Value / m.Count
		p4values = append(p4values, p4)
	}

	qfactor := getQfactor(p2values[1:10], p4values[1:10])

	for i := range p4values {
		p4values[i] *= qfactor / (ksmean.Value / ksmean.Count)
	}

	return p4values
}
