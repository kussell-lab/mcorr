package main

// VCFRecord store information for a VCF.
type VCFRecord struct {
	Chrom    string // chromosome name
	Pos      int    // position in the chromosome
	Ref, Alt string // reference and alternative allels
	GTs      []byte
}
