package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/kussell-lab/mcorr"

	"gopkg.in/alecthomas/kingpin.v2"
)

func main() {
	app := kingpin.New("mcorr-vcf", "Calculate mutational correlation from VCF files.")
	app.Version("v20171020")
	vcfFileArg := app.Arg("vcf-file", "VCF input file.").Required().String()
	outFileArg := app.Arg("out-prefix", "output prefix.").Required().String()
	maxlFlag := app.Flag("max-corr-length", "max length of correlations (bp).").Default("300").Int()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	vcfChan := readVCF(*vcfFileArg)

	p2arr := make([]float64, *maxlFlag)
	p2counts := make([]int64, *maxlFlag)
	var buffer []VCFRecord
	for rec := range vcfChan {
		if len(buffer) == 0 || rec.Pos-buffer[0].Pos < *maxlFlag {
			buffer = append(buffer, rec)
		} else {
			for i := 0; i < len(buffer); i++ {
				nc := mcorr.NewNuclCov([]byte{'0', '1'})
				for k := 0; k < len(buffer[0].GTs); k++ {
					nc.Add(buffer[0].GTs[k], buffer[i].GTs[k])
				}
				lag := buffer[i].Pos - buffer[0].Pos
				xy, n := nc.P11(0)
				p2arr[lag] += xy / float64(n)
				p2counts[lag]++
			}
			buffer = buffer[1:]
		}
	}

	w, err := os.Create(*outFileArg)
	if err != nil {
		panic(err)
	}
	defer w.Close()
	w.WriteString("l,m,n,v,t,b\n")
	var ks float64
	for k := 0; k < len(p2arr); k++ {
		m := p2arr[k] / float64(p2counts[k])
		n := p2counts[k]
		t := "P2"
		b := "all"
		v := 0.0
		l := k
		if l == 0 {
			ks = m
			t = "Ks"
		} else {
			m /= ks
		}
		if n > 0 {
			w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", l, m, v, n, t, b))
		}
	}
}

// readVCF return a channel of VCF record.
func readVCF(filename string) (c chan VCFRecord) {
	go func() {
		defer close(c)
		f, err := os.Open(filename)
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
			if line[0] == '#' {
				continue
			}

			line = strings.TrimSpace(line)
			terms := strings.Split(line, "\t")
			var rec VCFRecord
			rec.Chrom = terms[0]
			rec.Pos = atoi(terms[1])
			rec.Ref = terms[3]
			rec.Alt = terms[4]
			if len(rec.Alt) == 1 && len(rec.Ref) == 1 {
				inGT := false
				for _, t := range terms {
					if t == "GT" {
						inGT = true
					} else if inGT {
						for _, gt := range t {
							if gt != '|' {
								rec.GTs = append(rec.GTs, byte(gt))
							}
						}
					}
				}
				c <- rec
			}
		}
	}()
	return
}

func atoi(s string) int {
	v, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return v
}
