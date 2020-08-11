package main

import (
	"github.com/apsteinberg/mcorr"
	"github.com/apsteinberg/ncbiftp/taxonomy"
)

// Calculator define a interface for calculating correlations.
type PairCalculator interface {
	CalcP2(a Alignment, others Alignment) (corrResults mcorr.CorrResults)
}

// MateCalculator for calculating correlation for two clusters of sequences.
type MateCalculator struct {
	CodingTable   *taxonomy.GeneticCode
	MaxCodonLen   int
	CodonOffset   int
	CodonPosition int
	Synonymous    bool
}

// NewMateCalculator returns a MateCalculator
func NewMateCalculator(codingTable *taxonomy.GeneticCode, maxCodonLen, codonOffset, codonPos int, synonymous bool) *MateCalculator {
	return &MateCalculator{
		CodingTable:   codingTable,
		MaxCodonLen:   maxCodonLen,
		CodonOffset:   codonOffset,
		CodonPosition: codonPos,
		Synonymous:    synonymous,
	}
}

// CalcP2 calcualtes P2
// this was originally
// CalcP2(aln1 Alignment, mates ... Alignment)
func (cc *MateCalculator) CalcP2(aln1, mate Alignment) (corrResults mcorr.CorrResults) {
	//if len(mates) == 0 {
	//	return
	//}
	var results []mcorr.CorrResult

	cs1 := cc.extractCodonSequences(aln1)
	cs2 := cc.extractCodonSequences(mate)

	for l := 0; l < cc.MaxCodonLen; l++ {
		totalP2 := 0.0
		totaln := 0
		for pos := 0; pos+l < len(cs1[0]) && pos+l < len(cs2[0]); pos++ {
			//codon pairs are codons separated by pos+l basepairs on either sequence
			//cpList1 := cc.extractCodonPairs(cs1, pos, pos+l)
			//cpList2 := cc.extractCodonPairs(cs2, pos, pos+l)

			cpList1 := []CodonPair{}
			for _, cc := range cs1 {
				if pos+l < len(cc) {
					cpList1 = append(cpList1, CodonPair{A: cc[pos], B: cc[pos+l]})
				}
			}

			cpList2 := []CodonPair{}
			for _, cc := range cs2 {
				if pos+l < len(cc) {
					cpList2 = append(cpList2, CodonPair{A: cc[pos], B: cc[pos+l]})
				}
			}

			multiCP1 := [][]CodonPair{}
			multiCP2 := [][]CodonPair{}
			multiCP1 = synonymousSplit(cpList1, cc.CodingTable)
			multiCP2 = synonymousSplit(cpList2, cc.CodingTable)

			for _, cp1 := range multiCP1 {
				nc1 := doubleCodons(cp1, cc.CodonPosition)
				for _, cp2 := range multiCP2 {
					nc2 := doubleCodons(cp2, cc.CodonPosition)
					if cc.Synonymous {
						aa1 := cc.translateCodonPair(cp1[0])
						aa2 := cc.translateCodonPair(cp2[0])
						if aa1 == aa2 {
							xy, n := nc1.MateP11(nc2, 0)
							totalP2 += xy
							totaln += n
						}
					} else {
						xy, n := nc1.MateP11(nc2, 0)
						totalP2 += xy
						totaln += n
					}
				}
			}
		}
		if totaln > 0 {
			res1 := mcorr.CorrResult{
				Lag:  l * 3,
				Mean: totalP2 / float64(totaln),
				N:    totaln,
				Type: "P2",
			}
			results = append(results, res1)
		}
	}

	corrResults = mcorr.CorrResults{ID: aln1.ID, Results: results}

	return
}

func (cc *MateCalculator) translateCodonPair(cp CodonPair) string {
	a := cc.CodingTable.Table[string(cp.A)]
	b := cc.CodingTable.Table[string(cp.B)]
	return string([]byte{a, b})
}

func (cc *MateCalculator) extractCodonSequences(aln Alignment) (csList []CodonSequence) {
	for _, s := range aln.Sequences {
		csList = append(csList, extractCodons(s, cc.CodonOffset))
	}
	return
}

//maybe this is the issue?? doesn't get pairs that are separated by pos+l?
func (cc *MateCalculator) extractCodonPairs(codonSequences []CodonSequence, i, j int) [][]CodonPair {
	codonPairs := []CodonPair{}
	for _, cc := range codonSequences {
		if i < len(cc) && j < len(cc) {
			pair := CodonPair{A: cc[i], B: cc[j]}
			codonPairs = append(codonPairs, pair)
		}
	}

	if cc.Synonymous {
		return synonymousSplit(codonPairs, cc.CodingTable)
	}

	return [][]CodonPair{codonPairs}
}
