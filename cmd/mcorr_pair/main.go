package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"runtime"

	"github.com/alecthomas/kingpin"
	"github.com/cheggaaa/pb"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

func main() {
	app := kingpin.New("mcorr", "Calculate mutation correlation for each pair of isolates")
	app.Version("v0.2")

	alnFile := app.Arg("input", "Alignment file in XMFA format").Required().String()
	outFile := app.Arg("output", "Output file in CSV format").Required().String()

	mateFile := app.Flag("mate_file", "mate file").Default("").String()
	maxl := app.Flag("max_corr_len", "Maximum length of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("ncpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	progress := app.Flag("progress", "show progress").Default("false").Bool()
	codonPos := app.Flag("codon_pos", "codon position").Default("3").Int()
	synonymous := app.Flag("synonymous", "synonymous").Default("false").Bool()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	if *ncpu == 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	var pbar *pb.ProgressBar
	if *progress {
		numAln := countAlignments(*alnFile)
		pbar = pb.StartNew(numAln)
		defer pbar.Finish()
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

	alnChan := readAlignments(*alnFile)

	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3
	codonOffset := 0

	if *codonPos == 4 {
		*synonymous = true
		*codonPos = 3
	}

	numJob := *ncpu
	done := make(chan bool)
	resChan := make(chan CorrResults)
	for i := 0; i < numJob; i++ {
		go func() {
			for aln := range alnChan {
				var mateSequence *seq.Sequence
				if mateMap != nil {
					geneid := strings.Split(aln.Sequences[0].Id, " ")[0]
					s, found := mateMap[geneid]
					if found {
						mateSequence = s
					}
				}
				corrRes := calcP2Coding(aln, codonOffset, maxCodonLen, codingTable, *synonymous, *codonPos-1, mateSequence)
				if pbar != nil {
					pbar.Increment()
				}
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

	collectors := collect(resChan)

	w, err := os.Create(*outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("l,m,v,n,t,b\n")
	for _, cc := range collectors {
		results := cc.Results()
		for _, res := range results {
			w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, cc.ID))
		}
	}
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
		for {
			alignment, err := xmfaReader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if len(alignment) > 0 {
				alnID := strings.Split(alignment[0].Id, " ")[0]
				alnChan <- Alignment{ID: alnID, Sequences: alignment}
			}
		}
	}
	go read()
	return
}

// CorrResult stores a correlation result.
type CorrResult struct {
	Lag      int
	Mean     float64
	Variance float64
	N        int
	Type     string
}

// CorrResults stores a list of CorrResult with an gene ID.
type CorrResults struct {
	ID      string
	Results []CorrResult
}

func calcP2Coding(aln Alignment, codonOffset int, maxCodonLen int, codingTable *taxonomy.GeneticCode, synonymous bool, codonPos int, mateSequence *seq.Sequence) (results []CorrResults) {
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

	for i, seq1 := range codonSequences {
		for j := i + 1; j < len(codonSequences); j++ {
			id := aln.Sequences[i].Id + "_vs_" + aln.Sequences[j].Id
			seq2 := codonSequences[j]
			crRes := CorrResults{ID: id}
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
				cr := CorrResult{}
				cr.Lag = l * 3
				cr.Mean = d / float64(t)
				cr.N = t
				cr.Type = "P2"
				crRes.Results = append(crRes.Results, cr)
			}
			results = append(results, crRes)
		}
		if mateSequence != nil {
			break
		}
	}

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

func collect(resChan chan CorrResults) (collectors []*Collector) {
	cmap := make(map[string]*Collector)
	for res := range resChan {
		terms := strings.Split(res.ID, "_vs_")
		isolate1 := getGenome(terms[0])
		isolate2 := getGenome(terms[1])
		if isolate1 > isolate2 {
			isolate1, isolate2 = isolate2, isolate1
		}
		id := strings.Join([]string{isolate1, isolate2}, "_vs_")
		_, found := cmap[id]
		if !found {
			cmap[id] = NewCollector()
		}
		cmap[id].Add(res)
	}

	for id, c := range cmap {
		c.ID = id
		collectors = append(collectors, c)
	}

	return
}

func getGenome(s string) string {
	return strings.Split(strings.Split(s, "|")[0], "_")[0]
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
