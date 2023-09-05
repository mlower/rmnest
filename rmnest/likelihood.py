from __future__ import annotations
import numpy as np
import bilby
from uncertainties import unumpy
from rmnest.model import FaradayRotation, GeneralisedFaradayRotation


class FRLikelihood(bilby.likelihood.Likelihood):
    """
    Faraday rotation likelihood to measure pulsar/fast radio burst rotation measures.

    Parameters
    ----------
    freq: array_like
        Corresponding frequencies the data covers (Hz)
    freq_cen: float
        Centre frequency of the archive (Hz)
    s_q: array_like
        Stokes Q flux or intensity measurements.
    s_u: array_like
        Stokes U flux or intensity measurements.
    rms_q: array_like
        RMS of the Stokes Q flux or intensity measurements.
    rms_u: array_like
        RMS of the Stokes U flux or intensity measurements.

    Notes
    -----
    The Gaussian likelihood is defined as:
    See supplementary materials of Bannister et al. (2019) for details
    """

    def __init__(
        self,
        freq: np.ndarray,
        freq_cen: float,
        s_q: np.ndarray,
        s_u: np.ndarray,
    ) -> None:
        super().__init__()
        self.freq = freq
        self.freq_cen = freq_cen

        self.s_q = s_q
        self.s_u = s_u

        self.parameters = dict.fromkeys(["psi_zero", "rm", "sigma"], 0.0)

    def log_likelihood(self) -> float:
        fr_model = FaradayRotation(
            self.freq,
            self.freq_cen,
            self.parameters["psi_zero"],
            self.parameters["rm"],
        )

        self.residuals = (
            (self.s_q**2) + (self.s_u**2) 
             - ((self.s_q * fr_model.m_q) + (self.s_u * fr_model.m_u))**2
        )

        ln_l = np.sum(- (self.residuals / (self.parameters["sigma"]**2)) / 2 -
                       np.log(2 * np.pi * self.parameters["sigma"]**2) / 2)

        return -0.5 * ln_l


class GFRLikelihood(bilby.likelihood.Likelihood):
    def __init__(
        self,
        freq: np.ndarray,
        freq_cen: float,
        s_q: np.ndarray,
        s_u: np.ndarray,
        s_v: np.ndarray,
        rms_q: np.ndarray | None = None,
        rms_u: np.ndarray | None = None,
        rms_v: np.ndarray | None = None,
    ) -> None:
        """Modified Gaussian likelihood for measuring the generalised Faraday effect in pulsars and fast radio bursts.

        Parameters
        ----------
        freq : np.ndarray
            Corresponding frequencies the data covers (MHz)
        freq_cen : float
            Centre frequency of the archive (MHz)
        s_q : np.ndarray
            Stokes Q flux or intensity measurements.
        s_u : np.ndarray
            Stokes U flux or intensity measurements.
        s_v : np.ndarray
            Stokes V flux or intensity measurements.
        rms_q : np.ndarray | None, optional
            RMS of the Stokes Q flux or intensity measurements, by default None.
        rms_u : np.ndarray | None, optional
            RMS of the Stokes U flux or intensity measurements, by default None.
        rms_v : np.ndarray | None, optional
            RMS of the Stokes V flux or intensity measurements, by default None.
        """
        super().__init__()
        self.freq = freq
        self.freq_cen = freq_cen

        if rms_q is None:
            rms_q = np.zeros_like(s_q)

        if rms_u is None:
            rms_u = np.zeros_like(s_u)

        if rms_v is None:
            rms_v = np.zeros_like(s_v)

        s_q = unumpy.uarray(s_q, rms_q)
        s_u = unumpy.uarray(s_u, rms_u)
        s_v = unumpy.uarray(s_v, rms_v)
        s_p = unumpy.sqrt(s_q**2 + s_u**2 + s_v**2)

        self._norm_s_q = s_q / s_p
        self._norm_s_u = s_u / s_p
        self._norm_s_v = s_v / s_p

        self.parameters = dict.fromkeys(
            ["psi_zero", "grm", "alpha", "chi", "phi", "theta", "sigma"], 0.0
        )

    @property
    def norm_s_q(self) -> np.ndarray:
        return unumpy.nominal_values(self._norm_s_q)

    @property
    def norm_s_q_rms(self) -> np.ndarray:
        return unumpy.std_devs(self._norm_s_q)

    @property
    def norm_s_q_sigma(self) -> np.ndarray:
        return np.sqrt(self.norm_s_q_rms**2 + self.parameters["sigma"] ** 2)

    @property
    def norm_s_u(self) -> np.ndarray:
        return unumpy.nominal_values(self._norm_s_u)

    @property
    def norm_s_u_rms(self) -> np.ndarray:
        return unumpy.std_devs(self._norm_s_u)

    @property
    def norm_s_u_sigma(self) -> np.ndarray:
        return np.sqrt(self.norm_s_u_rms**2 + self.parameters["sigma"] ** 2)

    @property
    def norm_s_v(self) -> np.ndarray:
        return unumpy.nominal_values(self._norm_s_v)

    @property
    def norm_s_v_rms(self) -> np.ndarray:
        return unumpy.std_devs(self._norm_s_v)

    @property
    def norm_s_v_sigma(self) -> np.ndarray:
        return np.sqrt(self.norm_s_v_rms**2 + self.parameters["sigma"] ** 2)

    def log_likelihood(self):
        gfr_model = GeneralisedFaradayRotation(
            self.freq,
            self.freq_cen,
            self.parameters["psi_zero"],
            self.parameters["grm"],
            alpha=self.parameters["alpha"],
            chi=self.parameters["chi"],
            phi=self.parameters["phi"],
            theta=self.parameters["theta"],
        )

        residual_q = (self.norm_s_q - gfr_model.m_q) / self.norm_s_q_sigma
        residual_u = (self.norm_s_u - gfr_model.m_u) / self.norm_s_u_sigma
        residual_v = (self.norm_s_v - gfr_model.m_v) / self.norm_s_v_sigma
        ln_l_q = np.sum(residual_q**2 + np.log(2 * np.pi * self.norm_s_q_sigma**2))
        ln_l_u = np.sum(residual_u**2 + np.log(2 * np.pi * self.norm_s_u_sigma**2))
        ln_l_v = np.sum(residual_v**2 + np.log(2 * np.pi * self.norm_s_v_sigma**2))
        return -0.5 * (ln_l_q + ln_l_u + ln_l_v)
