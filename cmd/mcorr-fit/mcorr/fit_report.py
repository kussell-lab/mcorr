import numpy
class FitReport(object):
    """statistics report of fitting results"""
    def __init__(self, fit_results, param_name, label_name=None):
        """generate FitReport from fit_results of the param"""
        self.param_name = param_name
        if label_name is None:
            self.label_name = param_name
        else:
            self.label_name = label_name

        self.boot_data = []
        self.raw_value = None
        for res in fit_results:
            if hasattr(res, param_name):
                value = getattr(res, param_name)
                group = res.group
                if group == "all":
                    self.raw_value = value
                else:
                    self.boot_data.append(value)
    def get_param_name(self):
        return self.param_name

    def get_label_name(self):
        return self.label_name

    def get_raw_value(self):
        return self.raw_value

    def get_boot_size(self):
        """return the size of the bootstrapping data"""
        return len(self.boot_data)

    def get_boot_mean(self):
        """return mean of the bootstrapping data"""
        return numpy.mean(self.boot_data)

    def get_boot_std(self):
        """return standard deviation of the bootstrapping data"""
        return numpy.std(self.boot_data)

    def get_boot_median(self):
        """return median of the bootstrapping data"""
        return numpy.median(self.boot_data)

    def get_boot_lower_bound(self):
        """return bootstrapping lower bound"""
        return numpy.percentile(self.boot_data, 5)

    def get_boot_upper_bound(self):
        """return bootstrapping upper bound"""
        return numpy.percentile(self.boot_data, 95)

    def report(self):
        value = ""
        value += "[%s]\n" % self.get_label_name()
        if self.get_raw_value():
            value += "value = %g\n" % self.get_raw_value()
        if len(self.boot_data) >= 10:
            value += "bootstrapping mean = %g\n" % self.get_boot_mean()
            value += "bootstrapping standard deviation = %g\n" % self.get_boot_std()
            value += "bootstrapping size = %d\n" % self.get_boot_size()
            value += "bootstrapping median = %g\n" % self.get_boot_median()
            value += "bootstrapping lower bound (5%%) = %g\n" % \
                self.get_boot_lower_bound()
            value += "bootstrapping upper bound (95%%) = %g\n" % \
                self.get_boot_upper_bound()
        return value
