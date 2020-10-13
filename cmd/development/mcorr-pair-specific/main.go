package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"strings"
	"time"

	"github.com/apsteinberg/biogo/seq"
	"github.com/apsteinberg/mcorr"
	"github.com/apsteinberg/ncbiftp/taxonomy"
	"github.com/tobgu/qframe"
	"gopkg.in/alecthomas/kingpin.v2"
)

func main() {
	app := kingpin.New("mcorr-pair-specific", "Calculate mutation correlation for each pair of isolates from a specific list. THIS IS A TEST")
	app.Version("v20200917")

	alnFile := app.Arg("in", "Alignment file in XMFA format").Required().String()
	outFile := app.Arg("out", "Output file in CSV format").Required().String()

	mateFile := app.Flag("second-alignment", "Second alignment file in XMFA format").Default("").String()
	maxl := app.Flag("max-corr-length", "Maximum length of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	codonPos := app.Flag("codon-position", "Codon position (1: first codon position; 2: second codon position; 3: third codon position; 4: synonymous at third codon position.").Default("4").Int()
	pairList := app.Flag("pair-list", "list of isolate pairs you want to run").Default("").String()
	kingpin.MustParse(app.Parse(os.Args[1:]))
	start := time.Now()

	if *ncpu == 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	synonymous := true
	if *codonPos == 4 {
		synonymous = true
		*codonPos = 3
	}
	if *codonPos <= 0 || *codonPos > 4 {
		log.Fatalln("--codon-position should be in the range of 1 to 4.")
	}

	var mateMap map[string]*seq.Sequence
	if *mateFile != "" {
		f, err := os.Open(*mateFile)
		if err != nil {
			panic(err)
		}
		rd := seq.NewFastaReader(f)
		sequences, err := rd.ReadAll()
		if err != nil {
			panic(err)
		}

		mateMap = make(map[string]*seq.Sequence)
		for _, s := range sequences {
			geneid := strings.Split(s.Id, " ")[0]
			mateMap[geneid] = s
		}

		f.Close()
	}

	//define the list of isolate pairs
	var isoPairs qframe.QFrame
	csvfile, err := os.Open(*pairList)
	if err != nil {
		log.Fatal(err)
	}
	isoPairs = qframe.ReadCSV(csvfile)

	alnChan := readAlignments(*alnFile)

	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3
	codonOffset := 0

	numJob := *ncpu
	done := make(chan bool)
	resChan := make(chan mcorr.CorrResults)
	for i := 0; i < numJob; i++ {
		go func() {
			for aln := range alnChan {
				var mateSequence *seq.Sequence
				if mateMap != nil {
					geneid, _ := getNames(aln.Sequences[0].Id)
					s, found := mateMap[geneid]
					if found {
						mateSequence = s
					}
				}

				corrRes := calcP2Coding(aln, isoPairs, codonOffset, maxCodonLen, codingTable, synonymous, *codonPos-1, mateSequence)
				for _, res := range corrRes {
					resChan <- res
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(resChan)
		for i := 0; i < numJob; i++ {
			<-done
		}
	}()

	CollectWrite(resChan, *outFile)
	//time it
	duration := time.Since(start)
	fmt.Println(duration)
}

// Alignment is an array of multiple sequences with same length
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(file string) (alnChan chan Alignment) {
	alnChan = make(chan Alignment)
	read := func() {
		defer close(alnChan)

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
				alnChan <- Alignment{ID: alnID, Sequences: alignment}
				fmt.Printf("\rRead %d alignments.", numAln)
				fmt.Printf("\r alignment ID: %s", alnID)
			}
		}
		fmt.Printf(" Total alignments %d\n", numAln)
	}
	go read()
	return
}

// isoPairs was pairList before
func calcP2Coding(aln Alignment, isoPairs qframe.QFrame, codonOffset int, maxCodonLen int, codingTable *taxonomy.GeneticCode, synonymous bool, codonPos int, mateSequence *seq.Sequence) (results []mcorr.CorrResults) {
	codonSequences := [][]Codon{}
	sequences := []seq.Sequence{}
	if mateSequence != nil {
		sequences = append(sequences, *mateSequence)
	}
	sequences = append(sequences, aln.Sequences...)
	for _, s := range sequences {
		codons := extractCodons(s, codonOffset)
		codonSequences = append(codonSequences, codons)
	}

	//define the list of isolate pairs
	//csvfile, err := os.Open(*pairList)
	//if err != nil {
	//	log.Fatal(err)
	//}
	//isoPairs := qframe.ReadCSV(csvfile)
	////fmt.Println(isoPairs)

	count := 0
	//breaks := 0
	for i, seq1 := range codonSequences {
		totPairs := isoPairs.Len()
		////if neither of the columns has one of the isolates, skip the whole
		//inner loop
		_, genomeName1 := getNames(aln.Sequences[i].Id)
		check := isoPairs.Filter(qframe.Or(
			qframe.Filter{Column: "pair_0", Comparator: "=", Arg: genomeName1},
			qframe.Filter{Column: "pair_1", Comparator: "=", Arg: genomeName1}))
		if check.Len() == 0 {
			continue
		}

		for j := i + 1; j < len(codonSequences); j++ {
			_, genomeName1 := getNames(aln.Sequences[i].Id)
			_, genomeName2 := getNames(aln.Sequences[j].Id)

			// check if this pair is part of our list, in either order
			check1 := isoPairs.Filter(qframe.And(
				qframe.Filter{Column: "pair_0", Comparator: "=", Arg: genomeName1},
				qframe.Filter{Column: "pair_1", Comparator: "=", Arg: genomeName2}))
			check2 := isoPairs.Filter(qframe.And(
				qframe.Filter{Column: "pair_0", Comparator: "=", Arg: genomeName2},
				qframe.Filter{Column: "pair_1", Comparator: "=", Arg: genomeName1}))

			if check1.Len()+check2.Len() == 0 {
				//fmt.Printf("I'm skipping again")
				continue
			}
			id := genomeName1 + "_vs_" + genomeName2
			seq2 := codonSequences[j]
			crRes := mcorr.CorrResults{ID: id}
			for l := 0; l < maxCodonLen; l++ {
				d := 0.0
				t := 0
				for k := 0; k < len(seq1)-l; k++ {
					c1 := seq1[k]
					c2 := seq2[k]
					a1, found1 := codingTable.Table[string(c1)]
					a2, found2 := codingTable.Table[string(c2)]
					if found1 && found2 && a1 == a2 {
						b1 := seq1[k+l]
						b2 := seq2[k+l]

						good := true
						if synonymous {
							d1, found1 := codingTable.Table[string(c1)]
							d2, found2 := codingTable.Table[string(c2)]
							if found1 && found2 && d1 == d2 {
								good = true
							} else {
								good = false
							}
						}
						if good {
							var codonPositions []int
							if codonPos < 0 || codonPos > 2 {
								codonPositions = []int{0, 1, 2}
							} else {
								codonPositions = append(codonPositions, codonPos)
							}
							for _, codonP := range codonPositions {
								if c1[codonP] != c2[codonP] {
									if b1[codonP] != b2[codonP] {
										d++
									}
								}
								t++
							}
						}
					}
				}
				cr := mcorr.CorrResult{}
				cr.Lag = l * 3
				cr.Mean = d / float64(t)
				cr.N = t
				cr.Type = "P2"
				crRes.Results = append(crRes.Results, cr)
			}
			results = append(results, crRes)
			// break the loop if we've gotten all isolate pairs
			count = count + 1
			//fmt.Printf("\rCurrent count is %d pairs.", count)
			if count == totPairs {
				//breaks = breaks + 1
				//fmt.Printf("\ronly broke %d times.", breaks)
				break
			}
		}
		if mateSequence != nil {
			break
		}
		if count == totPairs {
			break
		}
	}
	fmt.Printf("\rFinal count is %d pairs.", count)
	return

}

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

// countAlignments return total number of alignments in a file.
func countAlignments(file string) (count int) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	rd := bufio.NewReader(f)
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		if line[0] == '=' {
			count++
		}
	}
	return
}

func getNames(s string) (geneName, genomeName string) {
	//double check this for the sra files
	terms := strings.Split(s, " ")
	geneName = terms[0]
	//genomeName = terms[1]
	//for the output of referencealignmentgenerator
	genomeName = terms[2]

	// for h pylori BIGSdb file
	//	terms := strings.Split(s, " ")
	//	geneName = terms[2]
	//	temp := terms[0]
	//	strainTerms := strings.Split(temp, ":")
	//	genomeName = strainTerms[0]

	return
}

// CollectWrite collects and writes the correlation results.
func CollectWrite(corrResChan chan mcorr.CorrResults, outFile string) {
	// prepare bootstrappers.
	bootstraps := make(map[string]*mcorr.Bootstrap)
	notBootstrap := mcorr.NewBootstrap("all", 1.0)
	notBootstrap.SetRandom(false)
	bootstraps["all"] = notBootstrap

	for corrResults := range corrResChan {
		id := corrResults.ID
		if _, found := bootstraps[id]; !found {
			bootstraps[id] = mcorr.NewBootstrap(id, 1.0)
			bootstraps[id].SetRandom(false)
		}
		bootstraps[id].Add(corrResults)
		bootstraps["all"].Add(corrResults)
	}

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("l,m,v,n,t,b\n")
	for _, bs := range bootstraps {
		results := bs.Results()
		for _, res := range results {
			w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, bs.ID))
		}
	}
}