import numpy as np
import bilby
from .model import FaradayRotation

class RMLikelihood(bilby.likelihood.Likelihood):
    def __init__(self,
        stokes_q,
        stokes_u,
        freq,
        freq_cen
    ):
        """
        The Gaussian likelihood from Bannister et al. (2019) - used for
        measuring pulsar/fast radio burst rotation measures.

        Parameters
        ----------
        stokes_q, stokes_u: array_like
            The polarisation data to analyse
        freq: array_like
            Corresponding frequencies the data covers (Hz)
        freq_cen: float
            Centre frequency of the archive (Hz)
        """

        super().__init__()
        self.stokes_q = stokes_u
        self.stokes_u = stokes_u
        self.freq = freq
        self.freq_cen = freq_cen
        self.parameters = dict(pa_zero=None, rm=None, sigma=None)

    def log_likelihood(self):
        self.sigma = self.parameters["sigma"]

        self.RM_model = FaradayRotation(self.stokes_q, self.stokes_u, self.freq,
            self.freq_cen, self.parameters["pa_zero"], self.parameters["rm"]
            )

        self.residual = self.RM_model.fit_QU()

        ln_l = np.sum(-(self.residual/(2*(self.sigma**2))) -
            np.log(2*np.pi*(self.sigma**2)) / 2)
        return ln_l


class GFRLikelihood(bilby.likelihood.Likelihood):
    def __init__(self,
        stokes_q,
        stokes_u,
        stokes_v,
        stokes_freq,
        freq_cen
    ):
        """
        Modified Gaussian likelihood for measuring the generalised Faraday
        effect in pulsars and fast radio bursts.

        Parameters
        ----------
        stokes_q, stokes_u, stokes_v: array_like
            The polarisation data to analyse.
        freq: array_like
            Corresponding frequencies the data covers.
        freq_cen: float
            Centre frequency of the archive (Hz)
        """
        super().__init__()
        self.stokes_q = stokes_q
        self.stokes_u = stokes_u
        self.stokes_v = stokes_v
        self.freq = freq
        self.freq_cen = freq_cen
        self.parameters = dict(psi_zero=None, RM=None, alpha=None, chi=None,
            phi=None, theta=None, sigma=None)

    def log_likelihood(self):
        self.sigma = self.parameters["sigma"]

        # Total polarisation
        self.p = np.sqrt(self.stokes_q**2 + self.stokes_u**2 + self.stokes_v**2)

        stokes_q_m, stokes_u_m, stokes_v_m = self.model = stokes_model(
            self.freq,
            freq_cen,
            DEG2RAD*self.parameters["psi_zero"],
            self.parameters["RM"],
            self.parameters["alpha"],
            DEG2RAD*self.parameters["chi"],
            DEG2RAD*self.parameters["phi"],
            DEG2RAD*self.parameters["theta"]
        )

        self.residual = (
            -((self.stokes_q/self.p) - stokes_q_m)**2
            - ((self.stokes_u/self.p) - stokes_u_m)**2
            - ((self.stokes_v/self.p) - stokes_v_m)**2
            )

        ln_l = np.sum((self.residual/(2*(self.sigma**2)))
            - np.log(2*np.pi*(self.sigma**2)) / 2)
        return ln_l
