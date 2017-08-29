// Calculate correlation functions (P2 and P4) from read mapping results.
package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"runtime"
	"sort"
	"strings"

	"github.com/biogo/hts/sam"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/mcorr"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"gopkg.in/alecthomas/kingpin.v2"
)

// SubProfile Substitution/mutation profile.
type SubProfile struct {
	Pos     int
	Profile []float64
}

// MinBaseQuality min base quality
var MinBaseQuality int

// MinMapQuality min map quality
var MinMapQuality int

// MinReadLength minimal read length
var MinReadLength int

// MinAlleleNumber minimal allele (pairs) number
var MinAlleleNumber int

func main() {
	// Command variables.
	var bamFile string      // bam or sam file
	var outFile string      // output file
	var maxl int            // max length of correlation
	var ncpu int            // number of CPUs
	var minDepth int        // min depth
	var minCoverage float64 // min coveage
	var gffFile string      // gff file
	var corrChanFile string // CorrResults channel results.

	// Parse command arguments.
	app := kingpin.New("mcorr-bam", "Calculate mutation correlation from bacterial metagenomic sequence in BAM read files.")
	app.Version("v20170728")
	gffFileArg := app.Arg("gff", "Gff3 file").Required().String()
	bamFileArg := app.Arg("in", "input file.").Required().String()
	outFileArg := app.Arg("out", "output file.").Required().String()

	maxlFlag := app.Flag("max-corr-length", "Max len of correlations (base pairs).").Default("300").Int()
	ncpuFlag := app.Flag("num-cpu", "Number of CPUs (default: using all available cores).").Default("0").Int()
	minDepthFlag := app.Flag("min-depth", "Minimal depth at each position.").Default("2").Int()
	minCoverageFlag := app.Flag("min-coverage", "Minimal coverage of a gene.").Default("0.5").Float64()
	minBaseQFlag := app.Flag("min-base-qual", "Minimal base quality").Default("30").Int()
	minMapQFlag := app.Flag("min-map-qual", "Minimal mapping quality").Default("30").Int()
	minReadLenFlag := app.Flag("min-read-length", "Minimal read length").Default("60").Int()
	codonPosition := app.Flag("codon-position", "Codon position (1: first codon position; 2: second codon position; 3: third codon position; 4: synoumous at third codon position.").Default("4").Int()
	numBoot := app.Flag("num-boot", "Number of bootstrapping on genes").Default("1000").Int()
	corrChanFileFlag := app.Flag("temp", "Temp results").Default("").String()
	minAlleleNumber := app.Flag("min-allele-number", "Minimal number of alleles").Default("0").Int()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	bamFile = *bamFileArg
	outFile = *outFileArg
	maxl = *maxlFlag / 3
	if *ncpuFlag <= 0 {
		ncpu = runtime.NumCPU()
	} else {
		ncpu = *ncpuFlag
	}
	runtime.GOMAXPROCS(ncpu)
	minDepth = *minDepthFlag
	minCoverage = *minCoverageFlag
	gffFile = *gffFileArg
	MinBaseQuality = *minBaseQFlag
	MinMapQuality = *minMapQFlag
	MinReadLength = *minReadLenFlag
	MinAlleleNumber = *minAlleleNumber
	corrChanFile = *corrChanFileFlag

	synoumous := false
	if *codonPosition == 4 {
		synoumous = true
		*codonPosition = 3
	}
	if *codonPosition <= 0 || *codonPosition > 4 {
		log.Fatalln("--codon-position should be in the range of 1 to 4.")
	}

	// Read sequence reads.
	var recordsChan chan GeneSamRecords
	gffRecMap := readGffs(gffFile)
	_, recordsChan = readStrainBamFile(bamFile, gffRecMap)

	codeTable := taxonomy.GeneticCodes()["11"]

	done := make(chan bool)
	p2Chan := make(chan mcorr.CorrResults)
	for i := 0; i < ncpu; i++ {
		go func() {
			for geneRecords := range recordsChan {
				geneLen := geneRecords.End - geneRecords.Start
				gene := pileupCodons(geneRecords)
				ok := checkCoverage(gene, geneLen, minDepth, minCoverage)
				if ok {
					p2 := calcP2(gene, 10, minDepth, codeTable, *codonPosition-1, synoumous)
					p4 := calcP4(gene, maxl, minDepth, codeTable, *codonPosition-1, synoumous)
					p2 = append(p2, p4...)
					p2Chan <- mcorr.CorrResults{Results: p2, ID: geneRecords.ID}
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(p2Chan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	var resChan chan mcorr.CorrResults
	if corrChanFile != "" {
		resChan = mcorr.PipeOutCorrResults(p2Chan, corrChanFile)
	} else {
		resChan = p2Chan
	}

	bootstraps := mcorr.Collect(resChan, *numBoot)

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("l,m,v,n,t,b\n")
	for _, bs := range bootstraps {
		results := bs.Results()
		qfactor := getQfactor(results)
		for _, res := range results {
			if res.Type == "Ks" || (res.Type == "P4" && res.Lag > 0) {
				if res.Type == "P4" {
					res.Mean *= qfactor
					res.Type = "P2"
				}
				w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n",
					res.Lag, res.Mean, res.Variance, res.N, res.Type, bs.ID))
			}
		}
	}
}

// getQfactor return the q factor between p2 and p4.
func getQfactor(results []mcorr.CorrResult) float64 {
	p2values := make([]float64, 31)
	p4values := make([]float64, 31)
	for _, res := range results {
		if res.Lag <= 30 && res.Lag > 0 {
			if res.Type == "P2" {
				p2values[res.Lag] = res.Mean
			} else if res.Type == "P4" {
				p4values[res.Lag] = res.Mean
			}
		}
	}

	var factors []float64
	for i := range p2values {
		if p2values[i] > 0 && p4values[i] > 0 {
			factors = append(factors, p2values[i]/p4values[i])
		}
	}

	if len(factors) == 0 {
		return 0
	}

	sort.Float64s(factors)
	if len(factors)%2 == 0 {
		return (factors[len(factors)/2] + factors[len(factors)/2-1]) / 2
	}
	return (factors[len(factors)/2])
}

// pileupCodons pileup codons of a list of reads at a gene.
func pileupCodons(geneRecords GeneSamRecords) (codonGene *CodonGene) {
	codonGene = NewCodonGene()
	for _, read := range geneRecords.Records {
		codonArray := getCodons(read, geneRecords.Start, geneRecords.Strand)
		for _, codon := range codonArray {
			if !codon.ContainsGap() {
				codonGene.AddCodon(codon)
			}
		}
	}

	return
}

// getCodons split a read into a list of Codon.
func getCodons(read *sam.Record, offset, strand int) (codonArray []Codon) {
	// get the mapped sequence of the read onto the reference.
	mappedSeq, _ := Map2Ref(read)
	for i := 2; i < len(mappedSeq); {
		if (read.Pos+i-offset+1)%3 == 0 {
			codonSeq := mappedSeq[i-2 : i+1]
			genePos := (read.Pos+i-offset+1)/3 - 1
			if genePos >= 0 {
				if strand == -1 {
					codonSeq = seq.Reverse(seq.Complement(codonSeq))
				}
				codon := Codon{ReadID: read.Name, Seq: string(codonSeq), GenePos: genePos}
				codonArray = append(codonArray, codon)
			}
			i += 3
		} else {
			i++
		}
	}

	return
}

func isATGC(b byte) bool {
	if b == 'A' {
		return true
	} else if b == 'T' {
		return true
	} else if b == 'C' {
		return true
	} else if b == 'G' {
		return true
	}

	return false
}

// P2 stores p2 calculation results.
type P2 struct {
	Total float64
	Count int
}

// doubleCount count codon pairs.
func doubleCount(nc *mcorr.NuclCov, codonPairArray []CodonPair, position int) {
	for _, cp := range codonPairArray {
		a := cp.A.Seq[position]
		b := cp.B.Seq[position]
		nc.Add(a, b)
	}
}

func calcP2(gene *CodonGene, maxl, minDepth int, codeTable *taxonomy.GeneticCode, codonPosition int, synoumous bool) (p2Res []mcorr.CorrResult) {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	for i := 0; i < gene.Len(); i++ {
		for j := i; j < gene.Len(); j++ {
			codonPairRaw := gene.PairCodonAt(i, j)
			if len(codonPairRaw) < 2 {
				continue
			}
			lag := codonPairRaw[0].B.GenePos - codonPairRaw[0].A.GenePos
			if lag < 0 {
				lag = -lag
			}
			if lag >= maxl {
				break
			}

			var splittedCodonPairs [][]CodonPair
			if synoumous {
				splittedCodonPairs = SynoumousSplitCodonPairs(codonPairRaw, codeTable)
			} else {
				splittedCodonPairs = [][]CodonPair{codonPairRaw}
			}

			for _, synPairs := range splittedCodonPairs {
				if len(synPairs) > minDepth {
					nc := mcorr.NewNuclCov(alphabet)
					doubleCount(nc, synPairs, codonPosition)

					for len(p2Res) <= lag {
						p2Res = append(p2Res, mcorr.CorrResult{Type: "P2", Lag: len(p2Res)})
					}
					xy, n := nc.P11(MinAlleleNumber)
					p2Res[lag].N += n
					p2Res[lag].Mean += xy
				}
			}
		}
	}

	for i := 0; i < len(p2Res); {
		if p2Res[i].N == 0 {
			p2Res = append(p2Res[:i], p2Res[i+1:]...)
		} else {
			p2Res[i].Mean /= float64(p2Res[i].N)
			p2Res[i].Lag *= 3
			i++
		}
	}

	return
}

func calcP4(gene *CodonGene, maxl, minDepth int, codeTable *taxonomy.GeneticCode, codonPosition int, synoumous bool) (p4Res []mcorr.CorrResult) {
	var valueArray []float64
	var countArray []int
	var posArray []int
	for i := 0; i < gene.Len(); i++ {
		value, count := autoCov(gene, i, minDepth, codeTable, codonPosition, synoumous)
		if count > 0 {
			pos := gene.CodonPiles[i].GenePos()
			valueArray = append(valueArray, value)
			countArray = append(countArray, count)
			posArray = append(posArray, pos)
		}
	}
	for i := 0; i < len(valueArray); i++ {
		value1 := valueArray[i]
		count1 := countArray[i]
		xbar := value1 / float64(count1)
		for j := i; j < len(valueArray); j++ {
			value2 := valueArray[j]
			count2 := countArray[j]
			ybar := value2 / float64(count2)
			lag := posArray[j] - posArray[i]
			if lag < 0 {
				lag = -lag
			}
			if lag >= maxl {
				break
			}
			for len(p4Res) <= lag {
				p4Res = append(p4Res, mcorr.CorrResult{Type: "P4", Lag: len(p4Res)})
			}
			p4Res[lag].Mean += xbar * ybar
			p4Res[lag].N++
		}
	}

	for i := 0; i < len(p4Res); {
		if p4Res[i].N == 0 {
			p4Res = append(p4Res[:i], p4Res[i+1:]...)
		} else {
			p4Res[i].Mean /= float64(p4Res[i].N)
			p4Res[i].Lag *= 3
			i++
		}
	}

	return
}

func autoCov(gene *CodonGene, i, minDepth int, codeTable *taxonomy.GeneticCode, codonPosition int, synoumous bool) (value float64, count int) {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	codonPairRaw := gene.PairCodonAt(i, i)
	if len(codonPairRaw) < 2 {
		return
	}
	lag := codonPairRaw[0].B.GenePos - codonPairRaw[0].A.GenePos
	if lag < 0 {
		lag = -lag
	}

	var splittedCodonPairs [][]CodonPair
	if synoumous {
		splittedCodonPairs = SynoumousSplitCodonPairs(codonPairRaw, codeTable)
	} else {
		splittedCodonPairs = [][]CodonPair{codonPairRaw}
	}

	for _, synPairs := range splittedCodonPairs {
		if len(synPairs) > minDepth {
			nc := mcorr.NewNuclCov(alphabet)
			doubleCount(nc, synPairs, codonPosition)

			xy, n := nc.P11(MinAlleleNumber)
			value += xy
			count += n
		}
	}
	return
}

// Map2Ref Obtains a read mapping to the reference genome.
func Map2Ref(r *sam.Record) (s []byte, q []byte) {
	p := 0                 // position in the read sequence.
	read := r.Seq.Expand() // read sequence.
	qual := r.Qual
	length := 0
	for _, c := range r.Cigar {
		switch c.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual, sam.CigarSoftClipped:
			length += c.Len()
		}
	}
	if length != len(read) || len(read) != len(qual) {
		return
	}

	for _, c := range r.Cigar {
		switch c.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual:
			s = append(s, read[p:p+c.Len()]...)
			q = append(q, qual[p:p+c.Len()]...)
			p += c.Len()
		case sam.CigarInsertion, sam.CigarSoftClipped:
			p += c.Len()
		case sam.CigarDeletion, sam.CigarSkipped:
			for i := 0; i < c.Len(); i++ {
				s = append(s, '-')
				q = append(q, 0)
			}
		}
	}

	s = bytes.ToUpper(s)

	for i, a := range q {
		if int(a) < MinBaseQuality {
			s[i] = '-'
		}
	}

	return
}

func checkCoverage(gene *CodonGene, geneLen, minDepth int, minCoverage float64) (ok bool) {
	num := 0
	for _, pile := range gene.CodonPiles {
		if pile.Len() > minDepth {
			num++
		}
	}
	coverage := float64(num) / float64(geneLen) * 3.0 // codon pile is in unit of codons (3)
	ok = coverage > minCoverage
	return
}

// readLines return all trimmed lines.
func readLines(filename string) []string {
	f, err := os.Open(filename)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	rd := bufio.NewReader(f)
	var lines []string
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			}
			break
		}
		lines = append(lines, strings.TrimSpace(line))
	}
	return lines
}

// subsample
func subsample(geneRecords GeneSamRecords, maxDepth float64) GeneSamRecords {
	length := float64(geneRecords.End - geneRecords.Start)
	readNum := len(geneRecords.Records)
	readLen := float64(geneRecords.Records[0].Len())
	maxReadNum := int(length * maxDepth / readLen)
	if readNum <= maxReadNum {
		return geneRecords
	}

	oldRecords := geneRecords.Records
	geneRecords.Records = []*sam.Record{}
	ratio := float64(maxReadNum) / float64(readNum)
	for _, read := range oldRecords {
		if rand.Float64() < ratio {
			geneRecords.Records = append(geneRecords.Records, read)
		}
	}

	return geneRecords
}
