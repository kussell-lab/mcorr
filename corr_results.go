package mcorr

// CorrResult stores a correlation result.
type CorrResult struct {
	Lag      int
	Mean     float64
	Variance float64
	N        int
	Type     string
}

// CorrResults stores a list of CorrResult with an gene ID.
type CorrResults struct {
	ID      string
	Results []CorrResult
}