package mcorr

import (
	"bytes"
	"fmt"
)

// NuclCov contains covariance of nucleotide acid in a DNA sequence.
type NuclCov struct {
	Doublets []int
	Alphabet []byte
}

// NewNuclCov return a NuclCov given the alphabet.
func NewNuclCov(alphabet []byte) *NuclCov {
	sizeOfAlphabet := len(alphabet)
	nc := NuclCov{Alphabet: alphabet}
	nc.Doublets = make([]int, sizeOfAlphabet*sizeOfAlphabet)
	return &nc
}

// Add insert a pair of nucliotide acids.
// It returns error when the nucliotide acid is not in the alphabet.
func (nc *NuclCov) Add(a, b byte) error {
	indexA := bytes.IndexByte(nc.Alphabet, a)
	indexB := bytes.IndexByte(nc.Alphabet, b)
	sizeOfAlphabet := len(nc.Alphabet)
	if indexA >= 0 && indexB >= 0 {
		nc.Doublets[indexA*sizeOfAlphabet+indexB]++
		return nil
	}

	var err error
	if indexA < 0 && indexB < 0 {
		err = fmt.Errorf("%c and %c are not in Alphabet: %s", a, b, string(nc.Alphabet))
	} else if indexA < 0 {
		err = fmt.Errorf("%c is not in Alphabet: %s", a, string(nc.Alphabet))
	} else {
		err = fmt.Errorf("%c is not in Alphabet: %s", b, string(nc.Alphabet))
	}

	return err
}

// Count returns the total number of pairs.
func (nc *NuclCov) Count() int {
	n := 0
	for _, a := range nc.Doublets {
		n += a
	}
	return n
}

// XBar return the count of 1 at the first position
func (nc *NuclCov) XBar() (x float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > 0 {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > 0 {
					c := float64(nc.Doublets[i] * nc.Doublets[j])

					if i/sizeOfAlphabet != j/sizeOfAlphabet {
						x += c
					}
					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
		}
	}
	return
}

// YBar return the count of 1 at the first position
func (nc *NuclCov) YBar() (x float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > 0 {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > 0 {
					c := float64(nc.Doublets[i] * nc.Doublets[j])
					if i%sizeOfAlphabet != j%sizeOfAlphabet {
						x += c
					}
					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
		}
	}
	return
}

// P00 returns the covariance.
func (nc *NuclCov) P00() (xy float64, n int) {
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > 0 {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > 0 {
					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
			xy += float64(nc.Doublets[i] * (nc.Doublets[i] - 1) / 2)
		}
	}
	return
}

// P01 returns the covariance.
func (nc *NuclCov) P01() (xy float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > 0 {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > 0 {
					c := float64(nc.Doublets[i] * nc.Doublets[j])
					if i%sizeOfAlphabet != j%sizeOfAlphabet {
						xy += c
					}
					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
		}
	}
	return
}

// P10 returns the covariance.
func (nc *NuclCov) P10() (xy float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > 0 {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > 0 {
					c := float64(nc.Doublets[i] * nc.Doublets[j])
					if i/sizeOfAlphabet != j/sizeOfAlphabet {
						xy += c
					}
					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
		}
	}
	return
}

// P11 returns the covariance.
func (nc *NuclCov) P11() (xy float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > 0 {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > 0 {
					c := float64(nc.Doublets[i] * nc.Doublets[j])
					if i%sizeOfAlphabet != j%sizeOfAlphabet && i/sizeOfAlphabet != j/sizeOfAlphabet {
						xy += c
					}

					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
		}
	}
	return
}

// MateP11 calculate covariance between two clusters.
func (nc *NuclCov) MateP11(nc2 *NuclCov) (xy float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > 0 {
			for j := 0; j < len(nc2.Doublets); j++ {
				if i != j && nc2.Doublets[j] > 0 {
					c := float64(nc.Doublets[i] * nc2.Doublets[j])
					if i%sizeOfAlphabet != j%sizeOfAlphabet && i/sizeOfAlphabet != j/sizeOfAlphabet {
						xy += c
					}
				}
			}
		}
	}
	n1 := 0
	n2 := 0
	for i := 0; i < len(nc.Doublets); i++ {
		n1 += nc.Doublets[i]
		n2 += nc2.Doublets[i]
	}
	n = n1 * n2
	return
}

// MateP00 calculate covariance between two clusters.
func (nc *NuclCov) MateP00(nc2 *NuclCov) (xy float64, n int) {
	n1, n2 := 0, 0
	for i := 0; i < len(nc.Doublets); i++ {
		xy += float64(nc.Doublets[i] * nc2.Doublets[i])
		n1 += nc.Doublets[i]
		n2 += nc2.Doublets[i]
	}
	n = n1 * n2
	return
}

// Append another NuclCov.
func (nc *NuclCov) Append(nc2 *NuclCov) error {
	// Check alphabet
	diffAlphabetError := fmt.Errorf("Different alphbet %s, %s", string(nc.Alphabet), string(nc2.Alphabet))
	if len(nc.Alphabet) != len(nc2.Alphabet) {
		return diffAlphabetError
	}
	for i, a := range nc.Alphabet {
		b := nc2.Alphabet[i]
		if a != b {
			return diffAlphabetError
		}
	}

	for i := 0; i < len(nc.Doublets); i++ {
		nc.Doublets[i] += nc2.Doublets[i]
	}

	return nil
}
