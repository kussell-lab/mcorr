package main

import (
	"github.com/apsteinberg/ncbiftp/taxonomy"
)

// Codon stores a codon value, the position in a genome and the read id.
type Codon struct {
	Seq     string
	ReadID  string
	GenePos int
}

// ContainsGap return true if '-' in a sequence.
func (c Codon) ContainsGap() bool {
	for _, b := range c.Seq {
		if b == '-' {
			return true
		}
	}
	return false
}

// CodonPile stores a pile of Codon, which are at a particular genome position.
type CodonPile struct {
	genePos  int
	codonMap map[string]Codon
}

// NewCodonPile return a new CodonPile.
func NewCodonPile() *CodonPile {
	return &CodonPile{codonMap: make(map[string]Codon)}
}

// Add appends a new Codon.
func (cp *CodonPile) Add(c Codon) {
	cp.genePos = c.GenePos
	cp.codonMap[c.ReadID] = c
}

// LookUp search a codon by ReadName. If not found, it returns nil.
func (cp *CodonPile) LookUp(readID string) Codon {
	return cp.codonMap[readID]
}

// Len return the lenght of pileup Codons.
func (cp *CodonPile) Len() int {
	return len(cp.codonMap)
}

// GenePos return the gene position.
func (cp *CodonPile) GenePos() int {
	return cp.genePos
}

// CodonGene represents a gene with an array of CodonPile.
type CodonGene struct {
	CodonPiles []*CodonPile
}

// NewCodonGene return a new CodonGene.
func NewCodonGene() *CodonGene {
	return &CodonGene{}
}

// AddCodon add a codon.
func (cg *CodonGene) AddCodon(c Codon) {
	for len(cg.CodonPiles) <= c.GenePos {
		cg.CodonPiles = append(cg.CodonPiles, NewCodonPile())
	}
	cg.CodonPiles[c.GenePos].Add(c)
}

// DepthAt return the pile depth at position i.
func (cg *CodonGene) DepthAt(i int) int {
	if len(cg.CodonPiles) <= i {
		return 0
	}
	return cg.CodonPiles[i].Len()
}

// Len returns length of CodonPile array.
func (cg *CodonGene) Len() int {
	return len(cg.CodonPiles)
}

// CodonPair stores a pair of Codon
type CodonPair struct {
	A, B Codon
}

// PairCodonAt pairs codons at positions i and j.
func (cg *CodonGene) PairCodonAt(i, j int) (pairs []CodonPair) {
	if i >= len(cg.CodonPiles) || j >= len(cg.CodonPiles) {
		return
	}

	if i > j {
		j, i = i, j
	}

	pile1 := cg.CodonPiles[i]
	if i == j {
		for _, codon := range pile1.codonMap {
			pairs = append(pairs, CodonPair{A: codon, B: codon})
		}
	}
	pile2 := cg.CodonPiles[j]
	for readID, codon1 := range pile1.codonMap {
		codon2 := pile2.LookUp(readID)
		if codon2.ReadID != "" {
			pairs = append(pairs, CodonPair{A: codon1, B: codon2})
		}
	}
	return
}

// SynoumousSplitCodonPairs split codon pairs into synoumous pairs.
func SynoumousSplitCodonPairs(codonPairs []CodonPair, codeTable *taxonomy.GeneticCode) [][]CodonPair {
	var splittedPairs [][]CodonPair
	var aaArray []string
	for _, codonPair := range codonPairs {
		hasGap := false
		for _, codon := range []Codon{codonPair.A, codonPair.B} {
			for _, b := range codon.Seq {
				if !isATGC(byte(b)) {
					hasGap = true
					break
				}
			}
			if hasGap {
				break
			}
		}

		if hasGap {
			continue
		}

		a := codeTable.Table[codonPair.A.Seq]
		b := codeTable.Table[codonPair.B.Seq]
		ab := string([]byte{a, b})
		index := -1
		for i, aa := range aaArray {
			if aa == ab {
				index = i
			}
		}
		if index == -1 {
			index = len(aaArray)
			aaArray = append(aaArray, ab)
			splittedPairs = append(splittedPairs, []CodonPair{})
		}
		splittedPairs[index] = append(splittedPairs[index], codonPair)
	}
	return splittedPairs
}
