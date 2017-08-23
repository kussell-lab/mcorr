package mcorr

import (
	"testing"
)

func TestXBar(t *testing.T) {
	doublets := [][]byte{[]byte{'A', 'T'}, []byte{'A', 'C'}, []byte{'C', 'C'}}
	nc := NewNuclCov([]byte{'A', 'T', 'G', 'C'})
	for _, db := range doublets {
		nc.Add(db[0], db[1])
	}
	expectedX := 2.0
	expectedN := 3
	gotX, gotN := nc.XBar()
	if expectedX != gotX {
		t.Errorf("Expect %g, got %g", expectedX, gotX)
	}
	if expectedN != gotN {
		t.Errorf("Expect %d, got %d", expectedN, gotN)
	}
}

func TestYBar(t *testing.T) {
	doublets := [][]byte{[]byte{'A', 'T'}, []byte{'A', 'C'}, []byte{'C', 'C'}, []byte{'C', 'G'}}
	nc := NewNuclCov([]byte{'A', 'T', 'G', 'C'})
	for _, db := range doublets {
		nc.Add(db[0], db[1])
	}
	expectedX := 5.0
	expectedN := 4 * 3 / 2
	gotX, gotN := nc.YBar()
	if expectedX != gotX {
		t.Errorf("Expect %g, got %g", expectedX, gotX)
	}
	if expectedN != gotN {
		t.Errorf("Expect %d, got %d", expectedN, gotN)
	}
}

func TestP11(t *testing.T) {
	doublets := [][]byte{[]byte{'A', 'T'}, []byte{'A', 'C'}, []byte{'C', 'C'}}
	nc := NewNuclCov([]byte{'A', 'T', 'G', 'C'})
	for _, db := range doublets {
		nc.Add(db[0], db[1])
	}
	expectedX := 1.0
	expectedN := 3
	gotX, gotN := nc.P11()
	if expectedX != gotX {
		t.Errorf("Expect %g, got %g", expectedX, gotX)
	}
	if expectedN != gotN {
		t.Errorf("Expect %d, got %d", expectedN, gotN)
	}
}

func TextP00(t *testing.T) {
	doublets := [][]byte{[]byte{'A', 'T'}, []byte{'A', 'C'}, []byte{'C', 'C'}}
	nc := NewNuclCov([]byte{'A', 'T', 'G', 'C'})
	for _, db := range doublets {
		nc.Add(db[0], db[1])
	}
	expectedX := 0.0
	expectedN := 3
	gotX, gotN := nc.P11()
	if expectedX != gotX {
		t.Errorf("Expect %g, got %g", expectedX, gotX)
	}
	if expectedN != gotN {
		t.Errorf("Expect %d, got %d", expectedN, gotN)
	}
}

func TextP10(t *testing.T) {
	doublets := [][]byte{[]byte{'A', 'T'}, []byte{'A', 'C'}, []byte{'C', 'C'}}
	nc := NewNuclCov([]byte{'A', 'T', 'G', 'C'})
	for _, db := range doublets {
		nc.Add(db[0], db[1])
	}
	expectedX := 1.0
	expectedN := 3
	gotX, gotN := nc.P10()
	if expectedX != gotX {
		t.Errorf("Expect %g, got %g", expectedX, gotX)
	}
	if expectedN != gotN {
		t.Errorf("Expect %d, got %d", expectedN, gotN)
	}
}

func TextP01(t *testing.T) {
	doublets := [][]byte{[]byte{'A', 'T'}, []byte{'A', 'C'}, []byte{'C', 'C'}}
	nc := NewNuclCov([]byte{'A', 'T', 'G', 'C'})
	for _, db := range doublets {
		nc.Add(db[0], db[1])
	}
	expectedX := 1.0
	expectedN := 3
	gotX, gotN := nc.P01()
	if expectedX != gotX {
		t.Errorf("Expect %g, got %g", expectedX, gotX)
	}
	if expectedN != gotN {
		t.Errorf("Expect %d, got %d", expectedN, gotN)
	}
}
