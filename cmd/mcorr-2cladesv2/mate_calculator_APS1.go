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
func (cc *MateCalculator) CalcP2(aln1, mates Alignment) (corrResults mcorr.CorrResults) {
	//if len(mates) == 0 {
	//	return
	//}

	var results []mcorr.CorrResult
	var sigmaI_k float64
	var sigmaL_k float64
	cs1 := cc.extractCodonSequences(aln1)
	//original was:
	//cs2 := cc.extractCodonSequences(mates[0])
	cs2 := cc.extractCodonSequences(mates)
	//test
	for l := 0; l < cc.MaxCodonLen; l++ {
		//totalP2 := 0.0
		P2 := 0.0
		totaln := 0
		sigmaI := 0.0
		probAB := 0.0
		//sigmaL := 0.0
		//test
		for pos := 0; pos+l < len(cs1[0]) && pos+l < len(cs2[0]); pos++ {
			//codon pairs are codons separated by pos+l basepairs on either sequence
			cpList1 := cc.extractCodonPairs(cs1, pos, pos+l)
			cpList2 := cc.extractCodonPairs(cs2, pos, pos+l)
			for _, cp1 := range cpList1 {
				//nc1 := doubleCodons(cp1, cc.CodonPosition)
				for _, cp2 := range cpList2 {
					//nc2 := doubleCodons(cp2, cc.CodonPosition)
					if cc.Synonymous {
						aa1 := cc.translateCodonPair(cp1[0])
						aa2 := cc.translateCodonPair(cp2[0])
						if aa1 == aa2 {
							//store which codon is at i (A) and i + l (B)
							//for sequences from both lists
							seq1A := cp1[0].A[cc.CodonPosition]
							seq2A := cp2[0].A[cc.CodonPosition]
							seq1B := cp1[0].B[cc.CodonPosition]
							seq2B := cp2[0].B[cc.CodonPosition]
							if seq1A != seq2A {
								// get the pairwise state of the locus
								// sigma is 1 if there is a difference, zero otherwise
								sigmaI += 1.0
								//set sigma i for the kth pair for the dot product
								sigmaI_k = 1.0
							} else {
								sigmaI_k = 0.0
							}
							//set sigma i+l for the kth pair for the dot product
							if seq1B != seq2B {
								//actually sigma i + l, but labeling like this for convenience
								sigmaL_k = 1.0
							} else {
								sigmaL_k = 0.0
							}

							//xy, n := nc1.MateP11(nc2, 0)
							//totalP2 += xy
							//probability of both positions being different for the kth pair
							probAB += sigmaI_k * sigmaL_k
							//this is the total number of pairs
							totaln += 1
						}
					} //else {
					//	//xy, n := nc1.MateP11(nc2, 0)
					//	//totalP2 += xy
					//	//totaln += n
					//}
				}
			}
		}
		//sample diversity for the gene
		d_s := sigmaI / float64(totaln)
		//joint probability of substitutions for the gene
		q_s := probAB / float64(totaln)
		// not sure what to do if d_s is zero, ask Edo
		// i think assigning P2 to be zero is appropriate though
		if d_s != 0 {
			P2 = q_s / d_s
		} else {
			P2 = 0.0
		}

		if totaln > 0 {
			res1 := mcorr.CorrResult{
				Lag: l * 3,
				//Mean: totalP2 / float64(totaln),
				Mean: P2,
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
