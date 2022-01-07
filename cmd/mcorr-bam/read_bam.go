package main

import (
	"io"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/kussell-lab/biogo/feat/gff"
)

// SamReader is an interface for sam or bam reader.
type SamReader interface {
	Header() *sam.Header
	Read() (*sam.Record, error)
}

func readSamRecords(fileName string) (headerChan chan *sam.Header, samRecChan chan *sam.Record) {
	headerChan = make(chan *sam.Header)
	samRecChan = make(chan *sam.Record)
	go func() {
		defer close(headerChan)
		defer close(samRecChan)

		// Open file stream, and close it when finished.
		f, err := os.Open(fileName)
		if err != nil {
			panic(err)
		}
		defer f.Close()

		// Decide if it is a .sam or .bam file.
		var reader SamReader
		if fileName[len(fileName)-3:] == "bam" {
			bamReader, err := bam.NewReader(f, 0)
			if err != nil {
				panic(err)
			}
			defer bamReader.Close()
			reader = bamReader
		} else {
			reader, err = sam.NewReader(f)
			if err != nil {
				panic(err)
			}
		}

		header := reader.Header()
		headerChan <- header

		// Read sam records and send them to the channel,
		// until it hit an error, which raises a panic
		// if it is not a IO EOF.
		for {
			rec, err := reader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			samRecChan <- rec
		}
	}()
	return
}

// GeneSamRecords stores Sam Records.
type GeneSamRecords struct {
	ID      string
	Start   int
	End     int
	Strand  int
	Records []*sam.Record
}

// readPanGenomeBamFile reads bam file, and return the header and a channel of sam records.
func readPanGenomeBamFile(fileName string) (header *sam.Header, recordsChan chan GeneSamRecords) {
	headerChan, samRecChan := readSamRecords(fileName)
	header = <-headerChan
	recordsChan = make(chan GeneSamRecords)
	go func() {
		defer close(recordsChan)
		currentRefID := ""
		var records []*sam.Record
		for rec := range samRecChan {
			if currentRefID == "" {
				currentRefID = rec.Ref.Name()
			}
			if rec.Ref.Name() != currentRefID {
				if len(records) > 0 {
					recordsChan <- GeneSamRecords{Start: 0, Records: records, End: records[0].Ref.Len(), ID: currentRefID}
					records = []*sam.Record{}
				}
				currentRefID = rec.Ref.Name()
			}
			records = append(records, rec)
		}
		if len(records) > 0 {
			recordsChan <- GeneSamRecords{Start: 0, Records: records, End: records[0].Ref.Len()}
		}
	}()

	return
}

//readStrainBamFile read []sam.Record from a bam file of mapping reads to a strain genome file.
func readStrainBamFile(fileName string, gffMap map[string][]*gff.Record) (header *sam.Header, resultChan chan GeneSamRecords) {
	headerChan, samRecChan := readSamRecords(fileName)
	header = <-headerChan
	resultChan = make(chan GeneSamRecords)
	go func() {
		defer close(resultChan)
		var currentRecord GeneSamRecords
		var currentIdx int
		var currentRef string
		for read := range samRecChan {
			// skip if read quality is low.
			if !checkReadQuality(read) {
				continue
			}

			// skip if read is not in the reference genomes.
			if _, found := gffMap[read.Ref.Name()]; !found {
				continue
			}

			// update current reference and current idx
			if currentRef != read.Ref.Name() {
				currentRef = read.Ref.Name()
				currentIdx = 0
			}

			// update current record.
			if currentRecord.End <= read.Pos {
				if currentRecord.Start < currentRecord.End && len(currentRecord.Records) > 0 {
					resultChan <- currentRecord
				}
				records := gffMap[currentRef]
				for currentIdx < len(records) && records[currentIdx].End <= read.Pos {
					currentIdx++
				}
				if currentIdx < len(records) {
					rec := records[currentIdx]
					currentRecord = GeneSamRecords{Start: rec.Start - 1, End: rec.End, ID: rec.ID(), Strand: rec.Strand}
				}
			}

			if read.Pos >= currentRecord.Start && read.Pos < currentRecord.End {
				currentRecord.Records = append(currentRecord.Records, read)
			}
		}
		if currentRecord.End > currentRecord.Start && len(currentRecord.Records) > 0 {
			resultChan <- currentRecord
		}

	}()
	return
}

func readGffs(fileName string) map[string][]*gff.Record {
	m := make(map[string][]*gff.Record)
	f, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	gffReader := gff.NewReader(f)

	records, err := gffReader.ReadAll()
	if err != nil {
		panic(err)
	}

	for _, rec := range records {
		if rec.Feature == "CDS" {
			m[rec.SeqName] = append(m[rec.SeqName], rec)
		}
	}
	return m
}

// checkReadQuality return false if the read fails quality check.
func checkReadQuality(read *sam.Record) bool {
	if int(read.MapQ) < MinMapQuality || read.Len() < MinReadLength {
		return false
	}

	// contains only match or mismatch
	for _, cigar := range read.Cigar {
		if cigar.Type() != sam.CigarMatch && cigar.Type() != sam.CigarSoftClipped {
			return false
		}
	}

	return true
}
