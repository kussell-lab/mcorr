package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"github.com/kussell-lab/biogo/seq"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
)

func splitPrep(cutoff int, clusterdict string, clusters []string) (threshold float64, seqMap map[string][]string) {
	//define the threshold as a fraction
	threshold = float64(cutoff) / 100
	seqMap = make(map[string][]string)
	//load the cluster dict as a golang map
	f, err := os.Open(clusterdict)
	if err != nil {
		log.Fatal(err)
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
		ln := strings.Split(line, ",")
		seq := ln[0]
		cluster := strings.Split(ln[1], "\n")
		seqMap[cluster[0]] = append(seqMap[cluster[0]], seq)
		//make MSA files for core and flexible genes ...
		clusterpath := "cluster" + cluster[0]
		coreMSA := filepath.Join(clusterpath, "MSA_CORE_cluster"+cluster[0])
		core, err := os.Create(coreMSA)
		check(err)
		core.Close()
		flexMSA := filepath.Join(clusterpath, "MSA_FLEX_cluster"+cluster[0])
		flex, err := os.Create(flexMSA)
		check(err)
		flex.Close()
		if err == io.EOF {
			break
		}
	}

	//make a gene percentage csv file
	//get the cluster names for the headers
	var clusternames []string
	for _, cluster := range clusters {
		clusternames = append(clusternames, "cluster"+cluster)
	}

	w, err := os.Create("gene_percentages.csv")
	if err != nil {
		panic(err)
	}
	defer w.Close()
	csvwriter := csv.NewWriter(w)
	defer csvwriter.Flush()
	header := []string{"gene"}
	header = append(header, clusternames...)
	csvwriter.Write(header)

	return
}

func coreflexSplit(threshold float64,
	seqMap map[string][]string,
	alnMap map[string][]seq.Sequence,
	clusters []string) (flex bool, CFgenes map[string]string, geneFrac map[string]float64) {
	//map where keys are the genes, and values are "core" and "flex"
	//map where keys are cluster IDs, and return values are the percentage of strains in cluster which have the gene
	geneFrac = make(map[string]float64)
	// map where keys are cluster IDs, and return values are whether this gene is core or flex
	CFgenes = make(map[string]string)
	//boolean which will tell us if the gene is not in one of clusters,
	//meaning its a flexible gene for all clusters
	flex = false
	//map out percentages for genes, and determine if its core or flexible
	for _, cluster := range clusters {
		var frac float64
		cAln, found := alnMap[cluster]
		// if there's no sequence there, then no strains have the gene
		if !found {
			flex = true
			//set gene frac to zero
			frac = 0
			geneFrac[cluster] = frac
		} else {
			strains, _ := seqMap[cluster]
			//# of strains in the cluster
			totstrains := float64(len(strains))
			//# of strains in cluster with the gene
			//count of strains in cluster with gene
			var count int
			for _, s := range cAln {
				//check if the strain has at least one full codon
				NumFullCodons := extractFullCodons(s)
				if NumFullCodons > 0 {
					count++
				}
			}
			frac = float64(count) / totstrains
			geneFrac[cluster] = frac
			//if the gene isn't present in one cluster, it's considered a flexible gene
			if frac == 0 {
				flex = true
			}
		}
	}
	//map out if its a core or flexible gene
	for cluster, frac := range geneFrac {
		if flex {
			CFgenes[cluster] = "FLEX"
		} else {
			if frac > threshold {
				CFgenes[cluster] = "CORE"
			} else {
				CFgenes[cluster] = "FLEX"
			}
		}
	}

	return
}

//func coreflexSplit(threshold float64,
//	seqMap map[string][]string,
//	alnMap map[string][]seq.Sequence) (flex bool, CFgenes map[string]string, geneFrac map[string]float64) {
//	//map where keys are the genes, and values are "core" and "flex"
//	//map where keys are cluster IDs, and return values are the percentage of strains in cluster which have the gene
//	geneFrac = make(map[string]float64)
//	// map where keys are cluster IDs, and return values are whether this gene is core or flex
//	CFgenes = make(map[string]string)
//	//boolean which will tell us if the gene is not in one of clusters,
//	//meaning its a flexible gene for all clusters
//	flex = false
//	//map out percentages for genes, and determine if its core or flexible
//	for cluster, cAln := range alnMap {
//		strains, _ := seqMap[cluster]
//		//# of strains in the cluster
//		totstrains := float64(len(strains))
//		//# of strains in cluster with the gene
//		//count of strains in cluster with gene
//		var count int
//		for _, s := range cAln {
//			//check if the strain has at least one full codon
//			NumFullCodons := extractFullCodons(s)
//			if NumFullCodons > 0 {
//				count++
//			}
//		}
//		frac := float64(count) / totstrains
//		//if frac != 1.0 && frac > 0{
//		//	fmt.Printf("not one!")
//		//}
//		geneFrac[cluster] = frac
//		if frac == 0 {
//			flex = true
//		}
//		//map out if its a core or flexible gene
//		if frac > threshold {
//			CFgenes[cluster] = "CORE"
//		} else {
//			CFgenes[cluster] = "FLEX"
//		}
//	}
//
//	return
//}

//write the core or flexible MSA files for the cluster
func WriteCFMSA(cAln cAlignment) {
	var cfMSA string

	clusterpath := "cluster" + cAln.clusterID
	//if flex {
	//	cfMSA = filepath.Join(clusterpath, "MSA_FLEX_cluster"+cAln.clusterID)
	//} else {
	//	cfMSA = filepath.Join(clusterpath, "MSA_"+cAln.genetype+"_cluster"+cAln.clusterID)
	//}
	cfMSA = filepath.Join(clusterpath, "MSA_"+cAln.genetype+"_cluster"+cAln.clusterID)
	//f, err := os.OpenFile(filename, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0600)
	f, err := os.OpenFile(cfMSA, os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
		panic(err)
	}
	for _, s := range cAln.Sequences {
		f.WriteString(">" + s.Id + "\n")
		f.Write(s.Seq)
		f.WriteString("\n")
	}
	f.WriteString("=\n")
	f.Close()
}

//func GetGenePercentages(PercentChan chan genePercents, clusters []string){
//	//get the cluster names for the headers
//	var clusternames []string
//	for _, cluster := range clusters{
//		clusternames = append(clusternames,"cluster"+cluster)
//	}
//
//	w, err := os.Create("gene_percentages.csv")
//	if err != nil {
//		panic(err)
//	}
//	defer w.Close()
//	csvwriter := csv.NewWriter(w)
//	header := []string{"gene"}
//	header = append(header, clusternames ...)
//	csvwriter.Write(header)
//	for gene := range PercentChan{
//		clusterPercent := []string{}
//		for _, cluster := range clusters{
//			p := fmt.Sprintf("%f", gene.Percent[cluster])
//			clusterPercent = append(clusterPercent, p)
//		}
//		genePercent := []string{gene.geneID}
//		genePercent = append(genePercent, clusterPercent ...)
//		csvwriter.Write(genePercent)
//
//	}
//}

func GetGenePercentages(geneID string, Percent map[string]float64, clusters []string) {

	w, err := os.OpenFile("gene_percentages.csv", os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
		panic(err)
	}
	defer w.Close()
	csvwriter := csv.NewWriter(w)
	defer csvwriter.Flush()
	//get percentages and write a line
	clusterPercent := []string{}
	for _, cluster := range clusters {
		p := fmt.Sprintf("%f", Percent[cluster])
		clusterPercent = append(clusterPercent, p)
	}
	genePercent := []string{geneID}
	genePercent = append(genePercent, clusterPercent...)
	csvwriter.Write(genePercent)

}

type genePercents struct {
	geneID  string
	Percent map[string]float64
}

//channel structure to pump genes into
type coreflex struct {
	clusterID string
	Sequences []seq.Sequence
}

// extractFullCodons returns the number of full codons
//there needs to be at least 1 full codon for us to say the strain "has the gene"
func extractFullCodons(s seq.Sequence) (NumFullCodons int) {
	var codons []Codon
	for i := 0; i+3 <= len(s.Seq); i += 3 {
		c := s.Seq[i:(i + 3)]
		//check for gaps
		containsGap := false
		for k := 0; k < 3; k++ {
			if c[k] == '-' || c[k] == 'N' {
				containsGap = true
				break
			}
		}
		if containsGap {
			continue
		} else {
			codons = append(codons, c)
		}

	}
	NumFullCodons = len(codons)
	return
}

// Codon is a byte list of length 3
type Codon []byte
