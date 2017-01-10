package main

import "github.com/mingzhi/ncbiftp/taxonomy"

// MateCalculator for calculating correlation for two clusters of sequences.
type MateCalculator struct {
	CodingTable *taxonomy.GeneticCode
	MaxCodonLen int
	CodonOffset int
	Synonymous  bool
}

// NewMateCalculator returns a MateCalculator
func NewMateCalculator(codingTable *taxonomy.GeneticCode, maxCodonLen, codonOffset int, synonymous bool) *MateCalculator {
	return &MateCalculator{
		CodingTable: codingTable,
		MaxCodonLen: maxCodonLen,
		CodonOffset: codonOffset,
		Synonymous:  synonymous,
	}
}

// CalcP2 calcualtes P2
func (cc *MateCalculator) CalcP2(aln1 Alignment, mates ...Alignment) (results []CorrResult) {
	if len(mates) == 0 {
		return
	}

	cs1 := cc.extractCodonSequences(aln1)
	cs2 := cc.extractCodonSequences(mates[0])

	for l := 0; l < cc.MaxCodonLen; l++ {
		totalxy := 0.0
		totaln := 0
		for pos := 0; pos+l < len(cs1[0]) && pos+l < len(cs2[0]); pos++ {
			cpList1 := cc.extractCodonPairs(cs1, pos, pos+l)
			cpList2 := cc.extractCodonPairs(cs2, pos, pos+l)
			for _, cp1 := range cpList1 {
				nc1 := doubleCodons(cp1)
				for _, cp2 := range cpList2 {
					nc2 := doubleCodons(cp2)
					if cc.Synonymous {
						aa1 := cc.translateCodonPair(cp1[0])
						aa2 := cc.translateCodonPair(cp2[0])
						if aa1 == aa2 {
							xy, _, _, n := nc1.CovMate11(nc2)
							totalxy += xy
							totaln += n
						}
					} else {
						xy, _, _, n := nc1.CovMate11(nc2)
						totalxy += xy
						totaln += n
					}
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

func (cc *MateCalculator) translateCodonPair(cp CodonPair) string {
	a := cc.CodingTable.Table[string(cp.A)]
	b := cc.CodingTable.Table[string(cp.B)]
	return string([]byte{a, b})
}

func (cc *MateCalculator) extractCodonSequences(aln Alignment) (csList []CodonSequence) {
	for _, s := range aln {
		csList = append(csList, extractCodons(s, cc.CodonOffset))
	}
	return
}

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
