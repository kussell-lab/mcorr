package mcorr

// Collector collect correlation results.
type Collector struct {
	m    map[string][]*MeanVar
	minN int
}

// NewCollector return a new Collector.
func NewCollector() *Collector {
	c := Collector{}
	c.m = make(map[string][]*MeanVar)
	return &c
}

// Add add an array of CorrResult.
func (c *Collector) Add(results CorrResults) {
	for _, res := range results.Results {
		for len(c.m[res.Type]) <= res.Lag {
			c.m[res.Type] = append(c.m[res.Type], NewMeanVar())
		}
		if res.N > c.minN {
			c.m[res.Type][res.Lag].Add(res.Mean)
		}
	}
}

// Means return means of a particular type.
func (c *Collector) Means(corrType string) (values []float64) {
	for _, mv := range c.MeanVars(corrType) {
		values = append(values, mv.Mean())
	}
	return
}

// Vars return variances of a particular type.
func (c *Collector) Vars(corrType string) (values []float64) {
	for _, mv := range c.MeanVars(corrType) {
		values = append(values, mv.Variance())
	}
	return
}

// Ns return variances of a particular type.
func (c *Collector) Ns(corrType string) (nums []int) {
	for _, mv := range c.MeanVars(corrType) {
		nums = append(nums, mv.N())
	}
	return
}

// MeanVars return a list of meanvar.MeanVar.
func (c *Collector) MeanVars(corrType string) (values []*MeanVar) {
	return c.m[corrType]
}

// CorrTypes return all corr types.
func (c *Collector) CorrTypes() (corrTypes []string) {
	for key := range c.m {
		corrTypes = append(corrTypes, key)
	}
	return
}

// Results get results
func (c *Collector) Results() (results []CorrResult) {
	// Failed fitting.
	if len(c.Means("P2")) == 0 {
		return nil
	}

	// calculate ks first
	ks := c.Means("P2")[0]
	results = append(results,
		CorrResult{
			Lag:      0,
			N:        c.Ns("P2")[0],
			Type:     "Ks",
			Mean:     c.Means("P2")[0],
			Variance: c.Vars("P2")[0],
		})
	if ks == 0 {
		return
	}

	for _, ctype := range c.CorrTypes() {
		means := c.Means(ctype)
		vars := c.Vars(ctype)
		ns := c.Ns(ctype)
		for i := 0; i < len(means); i++ {
			if !(ctype == "P2" && i == 0) && ns[i] > 0 {
				res := CorrResult{}
				res.Lag = i
				res.N = ns[i]
				res.Type = ctype
				res.Mean = means[i] / ks
				res.Variance = vars[i] / (ks * ks)
				results = append(results, res)
			}
		}
	}

	return
}
