class FitRes(object):
    """Fitting results"""
    def __init__(self, group, fit_res, d_sample):
        self.group = group
        self.d_sample = d_sample
        self.residual = fit_res.residual
        params = fit_res.params.valuesdict()
        if "theta" in params:
            self.theta_pool = params['theta']
        if 'phi' in params:
            self.phi_pool = params['phi']
        if 'fbar' in params:
            self.fbar = params['fbar']
        if 'phi' in params:
            self.ratio = self.phi_pool / self.theta_pool
            if 'fbar' in params:
                self.rho = self.phi_pool * self.fbar
        if 'c' in params:
            self.c = params['c']
        if 'dclonal' in params:
            self.d_clonal = params['dclonal']
        if 'dpool' in params:
            self.d_pool = params['dpool']
        if 'phi_clonal' in params:
            self.phi_clonal = params['phi_clonal']
        if 'theta_clonal' in params:
            self.theta_clonal = params['theta_clonal']

    def get_values(self, attributes):
        """Get attribute values"""
        values = []
        for name in attributes:
            if hasattr(self, name):
                values.append(getattr(self, name))
            else:
                values.append("NA")
        return values


