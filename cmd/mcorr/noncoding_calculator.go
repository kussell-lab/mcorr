package main

import "github.com/mingzhi/biogo/seq"

// NoncodingCalculator for calculating noncoding sequences.
type NoncodingCalculator struct {
	MaxLen int
}

// NewNoncodingCalculator return a NoncodingCalculator
func NewNoncodingCalculator(maxLen int) *NoncodingCalculator {
	return &NoncodingCalculator{
		MaxLen: maxLen,
	}
}

// CalcP2 calculate P2
func (cc *NoncodingCalculator) CalcP2(alignment []seq.Sequence) (results []CorrResult) {
	return calcP2Noncoding(alignment, cc.MaxLen)
}

func calcP2Noncoding(aln []seq.Sequence, maxLen int) (results []CorrResult) {
	for l := 0; l < maxLen; l++ {
		totalxy := 0.0
		totaln := 0
		for i := 0; i+l < len(aln[0].Seq); i++ {
			j := i + l
			basePairs := [][]byte{}
			for _, s := range aln {
				basePairs = append(basePairs, []byte{s.Seq[i], s.Seq[j]})
			}

			nc := doubleCounts(basePairs)
			xy, _, _, n := nc.Cov11()
			totalxy += xy
			totaln += n
		}
		res := CorrResult{
			Lag:  l,
			Mean: totalxy / float64(totaln),
			N:    totaln,
			Type: "P2"}
		results = append(results, res)
	}

	return
}

func doubleCounts(basePairs [][]byte) *NuclCov {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	c := NewNuclCov(alphabet)
	for _, basePair := range basePairs {
		a := basePair[0]
		b := basePair[1]
		if isATGC(a) && isATGC(b) {
			c.Add(a, b)
		}
	}
	return c
}

// ATGC is DNA alphabet.
const ATGC = "ATGC"
const atgc = "atgc"

func isATGC(b byte) bool {
	yes := false
	for i := 0; i < len(ATGC); i++ {
		if b == ATGC[i] {
			yes = true
			break
		}
	}
	return yes
}
