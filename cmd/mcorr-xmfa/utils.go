package main

import (
	"bufio"
	"io"
	"os"
)

// openFile is a helper function to open a file.
// and panic if error occurs.
func openFile(file string) (f *os.File) {
	var err error
	f, err = os.Open(file)
	if err != nil {
		panic(err)
	}
	return
}

// countAlignments return total number of alignments in a file.
func countAlignments(file string) (count int) {
	f := openFile(file)
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
