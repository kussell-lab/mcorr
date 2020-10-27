package main

import (
	"bufio"
	"fmt"
	//"github.com/kussell-lab/biogo/seq"
	"io"
	"log"
	"os"
	"runtime"
	"strings"
	"time"

	"github.com/apsteinberg/biogo/seq"
	"github.com/apsteinberg/mcorr"
	"github.com/apsteinberg/ncbiftp/taxonomy"
	//	"github.com/tobgu/qframe"
	"gopkg.in/alecthomas/kingpin.v2"
)

func main() {
	app := kingpin.New("mcorr-parallel", "Calculate mutation correlation between a list of sequences and all other sequences in an MSA")
	app.Version("v20201025")

	alnFile := app.Arg("in", "Alignment file in XMFA format").Required().String()
	outFile := app.Arg("out", "Output file in CSV format").Required().String()

	//mateFile := app.Flag("second-alignment", "Second alignment file in XMFA format").Default("").String()
	maxl := app.Flag("max-corr-length", "Maximum length of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	codonPos := app.Flag("codon-position", "Codon position (1: first codon position; 2: second codon position; 3: third codon position; 4: synonymous at third codon position.").Default("4").Int()
	seqList := app.Flag("sequence-list", "list of sequences for which you'll compute correlations with all sequences in MSA").Default("").String()
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

	////define the list of isolate pairs
	//var mateSeqs qframe.QFrame
	//csvfile, err := os.Open(*seqList)
	//if err != nil {
	//	log.Fatal(err)
	//}
	//mateSeqs = qframe.ReadCSV(csvfile)

	mateSeqs := []string{"SRR3343540", "1", "B", "2", "C", "3"}

	var mateMap map[string]*seq.Sequence
	if *seqList != "" {
		f, err := os.Open(*alnFile)
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
			genomeID := strings.Split(s.Name, " ")[2]

			//check := mateSeqs.Filter(
			//	qframe.Filter{Column: "Sequences", Comparator: "=", Arg: genomeID})
			//if check.Len() == 1 {
			//	geneid := strings.Split(s.Id, " ")[0]
			//	mateMap[geneid] = s
			//}
			_, found := Find(mateSeqs, genomeID)
			if found {
				geneid := strings.Split(s.Id, " ")[0]
				mateMap[geneid] = s
				fmt.Printf("found %s\n", genomeID)
			}

		}

		f.Close()
	}

	alnChan := readAlignments(*alnFile)

	//mateAlnChan :=

	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3
	codonOffset := 0

	numJob := *ncpu
	done := make(chan bool)
	resChan := make(chan mcorr.CorrResults)
	for i := 0; i < numJob; i++ {
		go func() {
			for mate := range mateSeqs {
				for aln := range alnChan {
					_, genomeID := getNames(aln.Sequences[0].Id)
					var mateSequence *seq.Sequence
					_, found := Find(mateSeqs, genomeID)
					if mateMap != nil && found {
						geneid, _ := getNames(aln.Sequences[0].Id)
						s, found := mateMap[geneid]
						if found {
							mateSequence = s
						}
					}
					//mateSeqs := []string{"SRR3343540", "1", "B", "2", "C", "3"}
					//find the same gene from our sequence list
					mateAlnChan := findMateAln(*alnFile, mateSeqs[mate], aln.ID)
					mateAln := <-mateAlnChan
					corrRes := calcP2Coding(aln, mateAln, codonOffset, maxCodonLen, codingTable, synonymous, *codonPos-1, mateSequence)
					for _, res := range corrRes {
						resChan <- res
					}
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

func findMateAln(file string, mate string, alnID string) (mateAln chan singleAln) {
	mateAln = make(chan singleAln)
	go func() {
		defer close(mateAln)

		//c := readXMFA(file)
		//for aln := range c {
		//	mateAlnID, genomeID := getNames(aln[0].Id)
		//	//mateAlnID, _, genomeID := strings.Split(alignment[0].Id, " ")
		//	//_, found := Find(mateSeqs, genomeID)
		//	if mateAlnID == alnID && mate == genomeID {
		//		mateAln <- Alignment{ID: mateAlnID, Sequences: aln}
		//	}
		//}
		f, err := os.Open(file)
		if err != nil {
			panic(err)
		}
		rd := seq.NewFastaReader(f)
		sequences, err := rd.ReadAll()
		if err != nil {
			panic(err)
		}
		for _, s := range sequences {
			genomeId := strings.Split(s.Name, " ")[2]
			geneId := strings.Split(s.Id, " ")[0]
			if genomeId == mate && geneId == alnID {
				mateAln <- singleAln{ID: alnID, Sequences: s}
			}
		}

	}()

	return
}

// Alignment is an array of multiple sequences with same length
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

// one sequence
type singleAln struct {
	ID        string
	Sequences *seq.Sequence
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
func calcP2Coding(aln Alignment, mateAln singleAln, codonOffset int, maxCodonLen int, codingTable *taxonomy.GeneticCode, synonymous bool, codonPos int, mateSequence *seq.Sequence) (results []mcorr.CorrResults) {
	codonSequences := [][]Codon{}
	//mateCodonSequences := [][]Codon{}
	sequences := []seq.Sequence{}
	//mateSequences := seq.Sequence{}
	mateCodonSequence := []Codon{}
	//if mateSequence != nil {
	//	mateSequences = append(sequences, *mateSequence)
	//	for _, s := range mateSequences {
	//		codons := extractCodons(s, codonOffset)
	//		mateCodonSequences = append(codonSequences, codons)
	//	}
	//}

	//mateSequences = append(mateSequences, mateAln.Sequences ...)
	//for _, s := range mateSequences {
	//	codons := extractCodons(s, codonOffset)
	//	mateCodonSequences = append(mateCodonSequences, codons)
	//}

	mateCodons := extractCodonsSingle(mateAln.Sequences, codonOffset)
	mateCodonSequence = mateCodons

	sequences = append(sequences, aln.Sequences...)
	for _, s := range sequences {
		codons := extractCodons(s, codonOffset)
		codonSequences = append(codonSequences, codons)
	}

	count := 0
	//breaks := 0
	seq1 := mateCodonSequence
	for i := 0; i < len(mateCodonSequence); i++ {
		//for i, seq1 := range mateCodonSequence {
		//seq1 := mateCodonSequence
		genomeName1 := strings.Split(mateAln.ID, " ")[2]

		for j := i + 1; j < len(codonSequences); j++ {
			//_, genomeName1 := getNames(mateAln.Sequences[i].Id)

			_, genomeName2 := getNames(aln.Sequences[j].Id)

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
		}
		if mateSequence != nil {
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

// extractCodons return a list of codons from a single DNA sequence.
func extractCodonsSingle(s *seq.Sequence, offset int) (codons []Codon) {
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

//function for finding the same gene in the list of sequences for genome 1

//func findMateAln(file string, alnID string) (mateAln chan Alignment) {
//	mateAln = make(chan Alignment)
//	go func() {
//		defer close(mateAln)
//
//		c := readAlignments(file)
//		for alignment := range c {
//			mateAlnID := strings.Split(alignment[0].Id, " ")[0]
//			if mateAlnID == alnID {
//				mateAln <- Alignment{ID: mateAlnID, Sequences: alignment}
//			}
//		}
//	}()
//
//	return
//}

//func (r *FastaReader) ReadSeqList() (seqs []*seq.Sequence, err error) {
//	for {
//		seq, err := r.Read()
//		if err == io.EOF {
//			seqs = append(seqs, seq)
//			return seqs, nil
//		}
//		if err != nil {
//			return nil, err
//		}
//		seqs = append(seqs, seq)
//	}
//	panic("unreachable")
//}
//
//// A reader for reading sequences in FASTA format.
//type FastaReader struct {
//	DeflineParser    func(string) string            // parse the definition line to get the seq id
//	AnnotationParser func(string) map[string]string // parse the defline to get the annotations.
//	r                *bufio.Reader
//}

//// Read reads one record from r. The record is a SeqRecord.
//func (r *FastaReader) Read() (seq *seq.Sequence, err error) {
//	seq, err = r.parseRecord()
//	seq.Seq = bytes.Replace(seq.Seq, []byte(" "), []byte(""), -1)
//	return
//}

// Find takes a slice and looks for an element in it. If found it will
// return it's key, otherwise it will return -1 and a bool of false.
func Find(slice []string, val string) (int, bool) {
	for i, item := range slice {
		if item == val {
			return i, true
		}
	}
	return -1, false
}

// readXMFA reads a xmfa format file and returns a channel of []seq.Sequence.
func readXMFA(file string) chan []seq.Sequence {
	c := make(chan []seq.Sequence)
	go func() {
		defer close(c)

		f := mustOpen(file)
		defer f.Close()

		rd := seq.NewXMFAReader(f)
		for {
			a, err := rd.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			//temporary change from 2 to 1
			if len(a) >= 1 {
				c <- a
			}
		}
	}()
	return c
}

// mustOpen is a helper function to open a file.
// and panic if error occurs.
func mustOpen(file string) (f *os.File) {
	var err error
	f, err = os.Open(file)
	if err != nil {
		panic(err)
	}
	return
}
