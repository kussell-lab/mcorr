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
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
)

func main() {
	app := kingpin.New("mcorr-pair-sync", "Calculate mutation correlation between individual pairs of isolates from two MSA files.")
	app.Version("v20201026")

	alnFile := app.Arg("first-alignment", "Alignment file in XMFA format").Required().String()
	mateFile := app.Arg("second-alignment", "Second alignment file in XMFA format").Required().String()
	outFile := app.Arg("out", "Output file in CSV format").Required().String()

	//mateFile := app.Flag("second-alignment", "Second alignment file in XMFA format").Default("").String()
	maxl := app.Flag("max-corr-length", "Maximum length of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	codonPos := app.Flag("codon-position", "Codon position (1: first codon position; 2: second codon position; 3: third codon position; 4: synonymous at third codon position.").Default("4").Int()
	showProgress := app.Flag("show-progress", "Show progress").Default("true").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))
	//timer

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

	//var mateMap map[string]Alignment
	var mateMap map[Key]Alignment
	f, err := os.Open(*mateFile)
	if err != nil {
		panic(err)
	}
	defer f.Close() //MIGHT NEED TO REMOVE!
	xmfaReader := seq.NewXMFAReader(f)
	//mateMap = make(map[string]Alignment)
	mateMap = make(map[Key]Alignment)
	NumSeq := 0
	for {
		s, err := xmfaReader.Read()
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		if len(s) > 0 {
			NumSeq++
			//grab the header for the gene alignment
			header := strings.Split(s[0].Id, " ")
			//here, we need both the gene name AND the genome position in
			//case there are duplicate genes on the genome
			//geneid := header[0]+"_"+header[1]
			geneid := header[0]
			genomePos := header[1]
			mateMap[Key{geneid, genomePos}] = Alignment{ID: geneid, Pos: genomePos, Sequences: s}

		}
	}
	//fmt.Print("length of mate map is %s", NumSeq)

	// show progress bar
	var bar *pb.ProgressBar
	if *showProgress {
		max := NumSeq
		bar = pb.StartNew(max)
		defer bar.Finish()
	}

	//special error log
	//If the file doesn't exist, create it or append to the file
	//file, err := os.OpenFile("201030-1523-mps_logs.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
	//if err != nil {
	//	log.Fatal(err)
	//}
	//log.SetOutput(file)

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
				var mateSequence Alignment
				if mateMap != nil {
					geneid, pos, _ := getNames(aln.Sequences[0].Id)
					s, found := mateMap[Key{geneid, pos}]
					if found {
						mateSequence = s
						corrRes := calcP2Coding(aln, codonOffset, maxCodonLen, codingTable, synonymous, *codonPos-1, mateSequence)
						for _, res := range corrRes {
							resChan <- res
						}
					}
					//if !found {
					//	fmt.Print("gotcha gotcha")
					//}
				}
				bar.Add(1)
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
	//was originally string
	ID        string // gene ID
	Pos       string //this is new; the position on the genome, so we don't compare alleles
	Sequences []seq.Sequence
}

//this is for the key map
type Key struct {
	ID, Pos string
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
				//alnID := strings.Split(alignment[0].Id, " ")[0]
				header := strings.Split(alignment[0].Id, " ")
				//here, we need both the gene name AND the genome position in
				//case there are duplicate genes on the genome
				//alnID := header[0]+"_"+header[1]
				geneID := header[0]
				genomePos := header[1]
				alnChan <- Alignment{ID: geneID, Pos: genomePos, Sequences: alignment}
				//bar.Add(1)
				//fmt.Printf("\rRead %d alignments.", numAln)
				//fmt.Printf("\r alignment ID: %s", alnID)
			}
		}
		//fmt.Printf(" Total alignments %d\n", numAln)
	}
	go read()
	return
}

func calcP2Coding(aln Alignment, codonOffset int, maxCodonLen int, codingTable *taxonomy.GeneticCode, synonymous bool, codonPos int, mateSequence Alignment) (results []mcorr.CorrResults) {
	codonSequences := [][]Codon{}
	sequences := []seq.Sequence{}
	//for the mate file
	mateCodonSequences := [][]Codon{}
	mateSequences := []seq.Sequence{}

	//if mateSequence != nil {
	//	sequences = append(sequences, *mateSequence)
	//}

	//define the sequences of the mate file
	mateSequences = append(mateSequences, mateSequence.Sequences...)
	//mateSequences = mateSequence
	for _, s := range mateSequences {
		mateCodons := extractCodons(s, codonOffset)
		mateCodonSequences = append(mateCodonSequences, mateCodons)
	}
	//define those of the alignment file
	sequences = append(sequences, aln.Sequences...)
	for _, s := range sequences {
		codons := extractCodons(s, codonOffset)
		codonSequences = append(codonSequences, codons)
	}

	//make sure that we don't run a pair of genomes where the genome
	//from the aln.Sequences is also part of
	//mateSequence list; this will result in a duplicate
	var mates []string
	for i, _ := range mateCodonSequences {
		var mate string
		_, _, mate = getNames(mateSequence.Sequences[i].Id)
		mates = append(mates, mate)
	}

	for i := 0; i < len(mateCodonSequences); i++ {
		//for i, seq1 := range mateCodonSequences {

		//mateName := strings.Split(mateSequences[i].Name, " ")
		//genomeName1 := mateName[2]

		for j := 0; j < len(codonSequences); j++ {
			//for j := i + 1; j < len(codonSequences); j++ {
			//_, genomeName1 := getNames(aln.Sequences[i].Id)
			//_, genomeName1 := getNames(mateSequence.Sequences[i].Id)
			seq1 := mateCodonSequences[i]
			_, _, genomeName1 := getNames(mateSequence.Sequences[i].Id)
			_, _, genomeName2 := getNames(aln.Sequences[j].Id)

			//don't run self vs self
			if genomeName1 == genomeName2 {
				continue
			}
			//make sure we don't run a genome pair twice
			if genomeName1 > genomeName2 {
				_, found := Find(mates, genomeName2)
				if found {
					continue
				} else {
					genomeName1, genomeName2 = genomeName2, genomeName1
				}
			}

			id := genomeName1 + "_vs_" + genomeName2
			seq2 := codonSequences[j]
			crRes := mcorr.CorrResults{ID: id}
			for l := 0; l < maxCodonLen; l++ {
				d := 0.0
				t := 0
				//took out constraint of k+l < len(seq2)
				for k := 0; k+l < len(seq1) && k+l < len(seq2); k++ {
					c1 := seq1[k]
					c2 := seq2[k]
					a1, found1 := codingTable.Table[string(c1)]
					a2, found2 := codingTable.Table[string(c2)]
					if found1 && found2 && a1 == a2 {
						b1 := seq1[k+l]
						b2 := seq2[k+l]

						good := true
						if synonymous {
							d1, found1 := codingTable.Table[string(b1)]
							d2, found2 := codingTable.Table[string(b2)]
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
				//log.Printf("%s %d", aln.ID, t)
				cr.Type = "P2"
				crRes.Results = append(crRes.Results, cr)
			}
			results = append(results, crRes)
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

func getNames(s string) (geneName, genomePos, genomeName string) {
	terms := strings.Split(s, " ")
	//this is for the helicobacter test files
	//geneName = terms[0]
	//genomeName = terms[1]
	//this is the genomeName for the MSA files assembled from ReferenceAlignmentGenerator
	//geneName = terms[0]+"_"+terms[1]
	geneName = terms[0]
	genomePos = terms[1]
	genomeName = terms[2]
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

//check error log
func check(e error) {
	if e != nil {
		panic(e)
	}
}
