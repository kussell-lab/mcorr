package main

import (
	"bufio"
	"io"

	"bytes"

	"github.com/mingzhi/biogo/seq"
)

// XMFAReader for extended Multi-FASTA format
type XMFAReader struct {
	r *bufio.Reader
}

// NewXMFAReader return a XMFAReader.
func NewXMFAReader(rd io.Reader) *XMFAReader {
	xmfa := XMFAReader{}
	xmfa.r = bufio.NewReader(rd)
	return &xmfa
}

// Read returns a alignment (block).
func (r XMFAReader) Read() (alignment []seq.Sequence, err error) {
	var seqID, seqBytes []byte
	for {
		var line []byte
		line, err = r.r.ReadBytes('\n')
		if err != nil {
			return
		}

		line = bytes.TrimSpace(line)
		if line[0] == '=' {
			break
		} else if line[0] == '#' {
			continue
		} else if line[0] == '>' {
			if len(seqBytes) > 0 {
				ss := seq.Sequence{}
				ss.Id = string(seqID)
				ss.Seq = seqBytes
				alignment = append(alignment, ss)
			}
			seqID = line[1:]
			seqBytes = []byte{}
		} else {
			seqBytes = append(seqBytes, line...)
		}

	}

	if len(seqBytes) > 0 {
		ss := seq.Sequence{}
		ss.Id = string(seqID)
		ss.Seq = seqBytes
		alignment = append(alignment, ss)
	}

	return
}
