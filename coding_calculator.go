package main

import (
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

// CorrResult stores a correlation result.
type CorrResult struct {
	Lag      int
	Mean     float64
	Variance float64
	N        int
	Type     string
}

// Calculator define a interface for calculating correlations.
type Calculator interface {
	CalcP2(alignments []seq.Sequence) (corrResults []CorrResult)
}

// CodingCalculator for calculating coding sequences.
type CodingCalculator struct {
	CodingTable *taxonomy.GeneticCode
	MaxCodonLen int
	CodonOffset int
}

// NewCodingCalculator return a CodingCalculator
func NewCodingCalculator(codingTable *taxonomy.GeneticCode, maxCodonLen, codonOffset int) *CodingCalculator {
	return &CodingCalculator{
		CodingTable: codingTable,
		MaxCodonLen: maxCodonLen,
		CodonOffset: codonOffset,
	}
}

// CalcP2 calculate P2
func (cc *CodingCalculator) CalcP2(alignment []seq.Sequence) (results []CorrResult) {
	return calcP2Coding(alignment, cc.CodonOffset, cc.MaxCodonLen, cc.CodingTable)
}

func calcP2Coding(aln []seq.Sequence, codonOffset int, maxCodonLen int, codingTable *taxonomy.GeneticCode) (results []CorrResult) {
	codonSequences := [][]Codon{}
	for _, s := range aln {
		codons := extractCodons(s, codonOffset)
		codonSequences = append(codonSequences, codons)
	}

	for l := 0; l < maxCodonLen; l++ {
		totalxy := 0.0
		totaln := 0
		for i := 0; i+l < len(codonSequences[0]); i++ {
			codonPairs := [][]Codon{}
			j := i + l
			for _, cc := range codonSequences {
				codonPairs = append(codonPairs, []Codon{cc[i], cc[j]})
			}

			multiCodonPairs := synonymousSplit(codonPairs, codingTable)
			for _, codonPairs := range multiCodonPairs {
				if len(codonPairs) >= 2 {
					nc := doubleCodons(codonPairs)
					xy, _, _, n := nc.Cov11()
					totalxy += xy
					totaln += n
				}
			}
		}
		res := CorrResult{
			Lag:  l * 3,
			Mean: totalxy / float64(totaln),
			N:    totaln,
			Type: "P2"}
		results = append(results, res)
	}

	return
}

func doubleCodons(codonPairs [][]Codon) *NuclCov {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	c := NewNuclCov(alphabet)
	for _, codonPair := range codonPairs {
		a := codonPair[0][2]
		b := codonPair[1][2]
		c.Add(a, b)
	}
	return c
}

// Codon is a byte list of length 3
type Codon []byte

// extractCodons return a list of codons from a DNA sequence.
func extractCodons(s seq.Sequence, offset int) (codons []Codon) {
	for i := offset; i+3 <= len(s.Seq); i += 3 {
		c := s.Seq[i:(i + 3)]
		codons = append(codons, c)
	}
	return
}

// synonymousSplit split a list of codon pairs into multiple
// synonymous pairs.
func synonymousSplit(codonPairs [][]Codon, codingTable *taxonomy.GeneticCode) (multiCodonPairs [][][]Codon) {
	aaList := []string{}
	for _, codonPair := range codonPairs {
		// check gap.
		containsGap := false
		for _, codon := range codonPair {
			for i := 0; i < 3; i++ {
				if codon[i] == '-' || codon[i] == 'N' {
					containsGap = true
					break
				}
			}
		}
		if containsGap {
			continue
		}

		codonA := string(codonPair[0])
		codonB := string(codonPair[1])
		a := codingTable.Table[codonA]
		b := codingTable.Table[codonB]
		ab := string([]byte{a, b})
		index := -1
		for i := 0; i < len(aaList); i++ {
			if aaList[i] == ab {
				index = i
			}
		}
		if index == -1 {
			index = len(aaList)
			aaList = append(aaList, ab)
			multiCodonPairs = append(multiCodonPairs, [][]Codon{})
		}

		multiCodonPairs[index] = append(multiCodonPairs[index], codonPair)
	}

	return
}
