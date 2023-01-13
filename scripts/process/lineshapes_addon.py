from nmrglue.analysis.lineshapes1d import location_2params, sim_lorentz_fwhm, center_fwhm
import numpy as np

class split_lorentz_fwhm(location_2params):
    """
    Lorentzian lineshape class representing a split peak, with height at center of
    each multiplet, full-width half-maximum, and J-constant parameters.
    """
    name = "split-lorentz"

    def sim(self, M, p):
        x = np.arange(M)
        x0, fwhm, j = p
        return sim_lorentz_fwhm(x, x0 - int(j/2), fwhm) + sim_lorentz_fwhm(x, x0 + int(j/2), fwhm)

    def guessp(self, sig):
        """Probably won't work."""
        c, fwhm = center_fwhm(sig)
        return (c, fwhm, 250)

    def pnames(self, M):
        return ("x0", "fwhm", "J-constant")
