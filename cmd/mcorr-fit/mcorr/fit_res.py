class FitRes(object):
    """Fitting results"""
    def __init__(self, group, fit_res, d_sample):
        self.group = group
        self.d_sample = d_sample
        self.residual = fit_res.residual
        params = fit_res.params.valuesdict()
        if "thetaP" in params:
            self.theta_pool = params['thetaP']
        if 'phiP' in params:
            self.phi_pool = params['phiP']
        if 'f' in params:
            self.fbar = params['f']
        if 'phiP' in params:
            self.ratio = self.phi_pool / self.theta_pool
            if 'f' in params:
                self.rho = self.phi_pool * self.fbar
        if 'c' in params:
            self.c = params['c']
        if 'dc' in params:
            self.d_clonal = params['dc']
        if 'dp' in params:
            self.d_pool = params['dp']
        if 'phiS' in params:
            self.phi_clonal = params['phiS']
        if 'thetaS' in params:
            self.theta_clonal = params['thetaS']

    def get_values(self, attributes):
        """Get attribute values"""
        values = []
        for name in attributes:
            if hasattr(self, name):
                values.append(getattr(self, name))
            else:
                values.append("NA")
        return values


