package main

import (
	"math"
)

// MeanVar is for calculate mean and variance in the increment way.
type MeanVar struct {
	N             int     // number of values.
	M1            float64 // first moment.
	Dev           float64
	NDev          float64
	M2            float64 // second moment.
	BiasCorrected bool
}

// NewMeanVar return a new MeanVar.
func NewMeanVar() *MeanVar {
	return &MeanVar{}
}

// Add adds a value.
func (m *MeanVar) Add(v float64) {
	if m.N < 1 {
		m.M1 = 0
		m.M2 = 0
	}

	m.N++
	n0 := m.N
	m.Dev = v - m.M1
	m.NDev = m.Dev / float64(n0)
	m.M1 += m.NDev
	m.M2 += float64(m.N-1) * m.Dev * m.NDev
}

// Mean returns the mean result.
func (m *MeanVar) Mean() float64 {
	return m.M1
}

// Variance returns the variance.
func (m *MeanVar) Variance() float64 {
	if m.N < 2 {
		return math.NaN()
	}

	if m.BiasCorrected {
		return m.M2 / float64(m.N-1)
	}

	return m.M2 / float64(m.N)
}

// Append add another result.
func (m *MeanVar) Append(m2 *MeanVar) {
	if m.N == 0 {
		m.N = m2.N
		m.M1 = m2.M1
		m.Dev = m2.Dev
		m.NDev = m2.NDev
		m.M2 = m2.M2
	} else {
		if m2.N > 0 {
			total1 := m.M1 * float64(m.N)
			total2 := m2.M1 * float64(m2.N)
			newMean := (total1 + total2) / float64(m.N+m2.N)
			delta1 := m.Mean() - newMean
			delta2 := m2.Mean() - newMean
			sm := (m.M2 + m2.M2) + float64(m.N)*delta1*delta1 + float64(m2.N)*delta2*delta2
			m.M1 = newMean
			m.M2 = sm
			m.N = m.N + m2.N
		}
	}
}
