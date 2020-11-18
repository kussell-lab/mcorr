package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"github.com/apsteinberg/biogo/seq"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
)

func splitPrep(cutoff int, clusterdict string) (threshold float64, seqMap map[string][]string) {
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
		//clusterpath := "cluster"+cluster[0]
		//coreMSA := filepath.Join(clusterpath,"MSA_CORE_cluster"+cluster[0])
		//core, err := os.Create(coreMSA)
		//check(err)
		//core.Close()
		//flexMSA := filepath.Join(clusterpath,"MSA_FLEX_cluster"+cluster[0])
		//flex, err := os.Create(flexMSA)
		//check(err)
		//flex.Close()
		if err == io.EOF {
			break
		}
	}

	return
}

func coreflexSplit(threshold float64,
	seqMap map[string][]string,
	alnMap map[string][]seq.Sequence) (flex bool, CFgenes map[string]string, geneFrac map[string]float64) {
	//map where keys are the genes, and values are "core" and "flex"
	//map where keys are cluster IDs, and return values are the percentage of strains in cluster which have the gene
	geneFrac = make(map[string]float64)
	// map where keys are cluster IDs, and return values are whether this gene is core or flex
	CFgenes = make(map[string]string)
	//boolean which will tell us if the gene is not in one of clusters,
	//meaning its a flexible gene for all clusters
	flex = false
	//map out percentages for genes, and determine if its core or flexible
	for cluster, cAln := range alnMap {
		strains, _ := seqMap[cluster]
		//# of strains in the cluster
		totstrains := len(strains)
		//# of strains in cluster with the gene
		strainswgene := len(cAln)
		frac := float64(strainswgene / totstrains)
		geneFrac[cluster] = frac
		if frac == 0 {
			flex = true
		}
		//map out if its a core or flexible gene
		if frac > threshold {
			CFgenes[cluster] = "CORE"
		} else {
			CFgenes[cluster] = "FLEX"
		}
	}

	return
}

//write the core and flexible MSA files for the cluster
func WriteCFMSA(flex bool, CFgenes map[string]string, clusterChan chan clusterAln) {
	var cfMSA string
	for c := range clusterChan {
		clusterpath := "cluster" + c.clusterID
		if flex {
			cfMSA = filepath.Join(clusterpath, "MSA_FLEX_cluster"+c.clusterID)
		} else {
			cf := CFgenes[c.clusterID]
			cfMSA = filepath.Join(clusterpath, "MSA_"+cf+"_cluster"+c.clusterID)
		}
		//f, err := os.OpenFile(filename, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0600)
		f, err := os.OpenFile(cfMSA, os.O_APPEND|os.O_WRONLY, 0600)
		if err != nil {
			panic(err)
		}
		for _, s := range c.Sequences {
			f.WriteString(">" + s.Id + "\n")
			f.Write(s.Seq)
			f.WriteString("\n")
		}
		f.WriteString("=\n")
		f.Close()
	}
}

func GetGenePercentages(PercentChan chan genePercents, clusters []string) {
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
	header := []string{"gene"}
	header = append(header, clusternames...)
	csvwriter.Write(header)
	for gene := range PercentChan {
		clusterPercent := []string{}
		for _, cluster := range clusters {
			p := fmt.Sprintf("%f", gene.Percent[cluster])
			clusterPercent = append(clusterPercent, p)
		}
		genePercent := []string{gene.geneID}
		genePercent = append(genePercent, clusterPercent...)
		csvwriter.Write(genePercent)

	}
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

//channel structure to pump a gene alignment for a cluster into
type cAlignment struct {
	clusterID string         //ID for cluster
	geneID    string         //gene name
	genetype  string         //string saying if core or flex
	fraction  float64        //fraction of strains which have the gene
	Sequences []seq.Sequence //the gene alignment for the cluster
}
