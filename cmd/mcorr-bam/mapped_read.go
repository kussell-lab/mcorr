package main

// MappedRead contains the section of a read mapped to a reference genome.
type MappedRead struct {
	Pos  int
	Seq  []byte
	Qual []byte
}

// Len return the lenght of a sequence.
func (m MappedRead) Len() int {
	return len(m.Seq)
}
