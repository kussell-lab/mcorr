package main

import (
	"math"

	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

// PartialCalculator for calculating partial genes.
type PartialCalculator struct {
	CodingTable *taxonomy.GeneticCode
	MaxCodonLen int
}

// NewPartialCalculator create a new PartialCalculator.
func NewPartialCalculator(codingTable *taxonomy.GeneticCode, maxCodonLen int) *PartialCalculator {
	return &PartialCalculator{CodingTable: codingTable, MaxCodonLen: maxCodonLen}
}

// CalcP2 calculates P2.
func (pc *PartialCalculator) CalcP2(alignment []seq.Sequence, others ...[]seq.Sequence) []CorrResult {
	offset := determineOffset(alignment, pc.CodingTable)
	results := calcP2Coding(alignment, offset, pc.MaxCodonLen, pc.CodingTable, true)
	ks := math.NaN()
	for _, res := range results {
		if res.Lag == 0 && res.Type == "P2" {
			ks = res.Mean
			break
		}
	}
	if !math.IsNaN(ks) && ks < 0.1 {
		return results
	}
	return results
}

func determineOffset(aln []seq.Sequence, codingTable *taxonomy.GeneticCode) int {
	max := -1.0
	maxOffset := -1
	for offset := 0; offset < 3; offset++ {
		results := calcP2Coding(aln, offset, 1, codingTable, true)
		m := results[0].Mean
		if m > max {
			max = m
			maxOffset = offset
		}
	}
	return maxOffset
}
