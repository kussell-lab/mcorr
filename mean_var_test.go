package mcorr

import (
    "testing"
    "math"
)

func TestMeanAndVariance(t *testing.T) {
    mv := NewMeanVar();
    if mv.Mean() != 0 {
        t.Error("Empty MeanVar should return zero for mean\n")
    }

    if !math.IsNaN(mv.Variance()) {
        t.Error("Empty MeanVar should return NaN for variance\n")
    }

    mv.Add(1.0)
    if mv.Mean() != 1.0 {
        t.Errorf("Expected 1.0, but got %g\n", mv.Mean())
    }

    if !math.IsNaN(mv.Variance()) {
        t.Errorf("Expected NaN, but got %g\n", mv.Variance())
    }

    resValues := []float64{1.0, 2.0, 4.0, 7.0}
    sum := 1.0
    for _, val := range resValues {
        sum += val
    }
    expectedMean := sum / float64(len(resValues) + 1)
    expectedVariance := (1.0 - expectedMean) * (1.0 - expectedMean)
    for _, val := range resValues {
        expectedVariance += (val - expectedMean) * (val - expectedMean)
    }
    expectedVariance /= float64(len(resValues) + 1)
    for _, val := range resValues {
        mv.Add(val)
    }
    if mv.Mean() != expectedMean {
        t.Errorf("Expected %g, but got %g\n", expectedMean, mv.Mean())
    }
    if mv.Variance() != expectedVariance {
        t.Errorf("Expected %g, but got %g\n", expectedVariance, mv.Variance())
    }

}

func TestN(t *testing.T) {
    mv := NewMeanVar()
    if mv.N() != 0 {
        t.Error("Empty MeanVariance should return zero for N()\n")
    }
    values := []float64{1, 2, 3, 4}
    for _, v := range values {
        mv.Add(v)
    }
    if mv.N() != len(values) {
        t.Errorf("Expected %d, but got %d\n", len(values), mv.N())
    }
}
