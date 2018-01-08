package mcorr

import (
	"math"
)

// MeanVar is for calculate mean and variance in the increment way.
type MeanVar struct {
	n             int     // number of values.
	m1            float64 // first moment.
	dev           float64
	nDev          float64
	m2            float64 // second moment.
	biasCorrected bool
}

// NewMeanVar return a new MeanVar.
func NewMeanVar() *MeanVar {
	return &MeanVar{}
}

// Add adds a value.
func (m *MeanVar) Add(v float64) {
	m.n++
	m.dev = v - m.m1
	m.nDev = m.dev / float64(m.n)
	m.m1 += m.nDev
	m.m2 += float64(m.n-1) * m.dev * m.nDev
}

// Mean returns the mean result.
func (m *MeanVar) Mean() float64 {
	return m.m1
}

// Variance returns the variance.
func (m *MeanVar) Variance() float64 {
	if m.n < 2 {
		return math.NaN()
	}

	if m.biasCorrected {
		return m.m2 / float64(m.n-1)
	}

	return m.m2 / float64(m.n)
}

// N returns the number of values.
func (m *MeanVar) N() int {
    return m.n
}

// IsBiasCorrected return true if the variance will be bias corrected.
func (m *MeanVar) IsBiasCorrected() bool {
    return m.biasCorrected
}

// SetBiasCorrected sets if bias corrected.
func (m *MeanVar) SetBiasCorrected(biasCorrected bool) {
    m.biasCorrected = biasCorrected
}

// Append add another result.
func (m *MeanVar) Append(m2 *MeanVar) {
	if m.n == 0 {
		m.n = m2.n
		m.m1 = m2.m1
		m.dev = m2.dev
		m.nDev = m2.nDev
		m.m2 = m2.m2
	} else {
		if m2.n > 0 {
			total1 := m.m1 * float64(m.n)
			total2 := m2.m1 * float64(m2.n)
			newMean := (total1 + total2) / float64(m.n+m2.n)
			delta1 := m.Mean() - newMean
			delta2 := m2.Mean() - newMean
			sm := (m.m2 + m2.m2) + float64(m.n)*delta1*delta1 + float64(m2.n)*delta2*delta2
			m.m1 = newMean
			m.m2 = sm
			m.n = m.n + m2.n
		}
	}
}
