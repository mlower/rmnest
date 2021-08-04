import numpy as np
import bilby
from rmnest.model import GeneralisedFaradayRotation


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
        self.stokes_q = stokes_q
        self.stokes_u = stokes_u
        self.freq = freq
        self.freq_cen = freq_cen

        # Linear polarisation
        self.l = np.sqrt(self.stokes_q**2 + self.stokes_u**2)

        self.parameters = dict.fromkeys(["psi_zero", "rm", "sigma"], None)

    def log_likelihood(self):
        self.sigma = self.parameters["sigma"]
        fr_model = GeneralisedFaradayRotation(
            self.freq,
            self.freq_cen,
            self.parameters["psi_zero"],
            2,
            self.parameters["rm"],
            0,
            0,
            0
        )

        self.residual = (
            ((self.stokes_q/self.l) - fr_model.m_q)**2
            + ((self.stokes_u/self.l) - fr_model.m_u)**2
            )

        ln_l = np.sum(-(self.residual/(2*(self.sigma**2)))
            - np.log(2*np.pi*(self.sigma**2)) / 2)
        return ln_l



class GFRLikelihood(bilby.likelihood.Likelihood):
    def __init__(self,
        stokes_q,
        stokes_u,
        stokes_v,
        freq,
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

        # Total polarisation
        self.p = np.sqrt(self.stokes_q**2 + self.stokes_u**2 + self.stokes_v**2)

        self.parameters = dict.fromkeys(
            ["psi_zero", "grm", "alpha", "chi", "phi", "theta", "sigma"], None
        )

    def log_likelihood(self):
        self.sigma = self.parameters["sigma"]
        gfr_model = GeneralisedFaradayRotation(
            self.freq,
            self.freq_cen,
            self.parameters["psi_zero"],
            self.parameters["alpha"],
            self.parameters["grm"],
            self.parameters["chi"],
            self.parameters["phi"],
            self.parameters["theta"]
        )

        self.residual = (
            ((self.stokes_q/self.p) - gfr_model.m_q)**2
            + ((self.stokes_u/self.p) - gfr_model.m_u)**2
            + ((self.stokes_v/self.p) - gfr_model.m_v)**2
            )

        ln_l = np.sum(-(self.residual/(2*(self.sigma**2)))
            - np.log(2*np.pi*(self.sigma**2)) / 2)
        return ln_l
