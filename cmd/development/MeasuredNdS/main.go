package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"github.com/apsteinberg/mcorr"
	"github.com/apsteinberg/ncbiftp/taxonomy"
	"github.com/kussell-lab/biogo/seq"
	"io"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"time"
)

// author: Asher Preska Steinberg

func main() {
	alnFile := "/Users/asherpreskasteinberg/Desktop/code/recombo/APS180_genebins/dnds_test/1224_properheader_GAPFILTERED"
	//strain list
	sampleFile := "/Users/asherpreskasteinberg/Desktop/code/recombo/APS180_genebins/dnds_test/strain_list"
	outdir := "/Users/asherpreskasteinberg/Desktop/code/recombo/APS180_genebins/dnds_test/two_strains"
	numSplitters := 4
	//timer
	start := time.Now()
	samples := readSamples(sampleFile)
	//get the total number of sequences
	totSeqs := len(samples)
	done := make(chan struct{})
	//read in alignments
	alignments, errc := readAlignments(done, alnFile)

	//start a fixed number of goroutines to read alignments and split into core/flex
	c := make(chan result)
	var wg sync.WaitGroup
	for i := 0; i < *numSplitters; i++ {
		wg.Add(1)
		go digester(done, alignments, c, totSeqs, minimum, maximum, i, &wg)
	}

	go func() {
		wg.Wait()
		close(c)
	}()
	//end of pipeline; write files
	for gene := range c {
		//if true, write to the out MSA
		if gene.bin {
			writeMSA(gene, *outdir, *min, *max)
		}
		getGenePercentage(gene, *outdir, *min, *max)
	}
	if err := <-errc; err != nil { // HLerrc
		panic(err)
	}
	//add the number of core and flex to the bottom of the spreadsheet

	duration := time.Since(start)
	fmt.Println("Time to make gene-binned MSA:", duration)

}

//makeGeneCSV initiates the gene percentage CSV
func makeGeneCSV(outdir string) {
	//prepare the gene percentage out csv
	name := "gene_percentages_dnds.csv"
	path := filepath.Join(outdir, name)
	w, err := os.Create(path)
	check(err)
	defer w.Close()
	csvwriter := csv.NewWriter(w)
	defer csvwriter.Flush()
	header := []string{"gene", "fraction of strains", "dN/dS"}
	err = csvwriter.Write(header)
	check(err)

	return
}

//check for errors
func check(e error) {
	if e != nil {
		panic(e)
	}
}

// readSamples return a list of samples from a sample file.
func readSamples(filename string) (samples []string) {
	f, err := os.Open(filename)
	if err != nil {
		log.Fatalf("Error when reading file %s:%v", filename, err)
	}
	defer f.Close()

	r := bufio.NewReader(f)
	for {
		line, err := r.ReadString('\n')

		if err != nil {
			if err != io.EOF {
				log.Fatalf("Error when reading file %s: %v", filename, err)
			}
			break
		}
		samples = append(samples, strings.TrimSpace(line))
	}
	return
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(done <-chan struct{}, file string) (<-chan Alignment, <-chan error) {
	alignments := make(chan Alignment)
	errc := make(chan error, 1)
	go func() {
		defer close(alignments)

		f, err := os.Open(file)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		xmfaReader := seq.NewXMFAReader(f)
		numAln := 0
		for {
			alignment, err := xmfaReader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if len(alignment) > 0 {
				numAln++
				alnID := strings.Split(alignment[0].Id, " ")[0]
				select {
				case alignments <- Alignment{alnID, numAln, alignment}:
					fmt.Printf("\rRead %d alignments.\n", numAln)
					fmt.Printf("\r alignment ID: %s\n", alnID)
				case <-done:
					fmt.Printf(" Total alignments %d\n", numAln)
				}
			}
		}
		errc <- err
	}()
	return alignments, errc
}

// Alignment is an array of multiple sequences with same length.
type Alignment struct {
	ID        string
	num       int
	Sequences []seq.Sequence
}

// A result is a single gene alignment which either belongs to the bin or not
type result struct {
	Alignment Alignment
	dNdS      float64 // average dN/dS ratio for the gene
	frac      float64 //fraction of strains that have the gene
}

// digester figures out what fraction of strains have the gene and the dN/dS ratio for the gene
// then sends these processed results on alnChan until either the master MSA or done channel is closed.
func digester(done <-chan struct{}, alignments <-chan Alignment, genes chan<- result, totSeqs int, min float64,
	max float64, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	//fmt.Printf("Worker %d starting\n", id)
	for aln := range alignments { // HLpaths
		//get the fraction of sequences which have the gene
		var frac float64
		//boolean for if it's in the bin or not
		var bin bool
		//count number of strains with the gene; which for gap-filtered MSAs filtered with FilterGaps
		//is just the number of strains which have the sequence
		var count int
		count = len(aln.Sequences)
		//get fraction of strains which have the gene
		frac = float64(count) / float64(totSeqs)
		//is it part of the bin or not
		if frac > min && frac <= max {
			bin = true
		} else {
			bin = false
		}
		gene := result{aln, bin, frac}
		select {
		case genes <- gene:
		case <-done:
			return
		}
	}
	//fmt.Printf("Worker %d done\n", id)

}

//base this off of calcP2 in mcorr-pair ...
func calcdNdS(aln Alignment, codingTable *taxonomy.GeneticCode) (dN float64, dS float64) {
	codonSequences := [][]Codon{}
	sequences := []seq.Sequence{}
	sequences = append(sequences, aln.Sequences...)
	for _, s := range sequences {
		codons := extractCodons(s, 0)
		codonSequences = append(codonSequences, codons)
		for i, seq1 := range codonSequences {
			for j := i + 1; j < len(codonSequences); j++ {
				seq2 := codonSequences[j]

			}
		}
	}
}

//func calcdNdS (aln Alignment, codingTable *taxonomy.GeneticCode)(dN float64, dS float64){
//	codonSequences := [][]Codon{}
//	for _, s := range aln.Sequences {
//		codons := extractCodons(s, 0)
//		codonSequences = append(codonSequences, codons)
//	}
//
//	totalS := 0.0
//	totalnS := 0
//
//	for i := 0; i < len(codonSequences[0]); i++ {
//		codonPairs := []CodonPair{}
//		for _, cc := range codonSequences {
//			if i < len(cc) {
//				codonPairs = append(codonPairs, CodonPair{A: cc[i], B: cc[i]})
//			}
//		}
//
//		multiCodonPairs := [][]CodonPair{}
//		//get average number of synonymous mutations for the gene
//		multiCodonPairs = synonymousSplit(codonPairs, codingTable)
//		for _, codonPairs := range multiCodonPairs {
//			if len(codonPairs) >= 2 {
//				for i := 1; i < 4; i++ {
//					codonPosition := i
//					nc := doubleCodons(codonPairs, codonPosition)
//					xy, n := nc.P11(0)
//					totalS += xy
//					totalnS += n
//				}
//
//			}
//		}
//		// get the average number of nonsynonymous mutations for the gene
//
//	}
//}

// Codon is a byte list of length 3
type Codon []byte

// CodonSequence is a sequence of codons.
type CodonSequence []Codon

// CodonPair is a pair of Codons.
type CodonPair struct {
	A, B Codon
}

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
func synonymousSplit(codonPairs []CodonPair, codingTable *taxonomy.GeneticCode) (multiCodonPairs [][]CodonPair) {
	aaList := []string{}
	for _, codonPair := range codonPairs {
		// check gap.
		containsGap := false
		for _, codon := range []Codon{codonPair.A, codonPair.B} {
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

		codonA := string(codonPair.A)
		codonB := string(codonPair.B)
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
			multiCodonPairs = append(multiCodonPairs, []CodonPair{})
		}

		multiCodonPairs[index] = append(multiCodonPairs[index], codonPair)
	}

	return
}

func doubleCodons(codonPairs []CodonPair, codonPosition int) *mcorr.NuclCov {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	c := mcorr.NewNuclCov(alphabet)
	for _, codonPair := range codonPairs {
		a := codonPair.A[codonPosition]
		b := codonPair.B[codonPosition]
		c.Add(a, b)
	}
	return c
}
