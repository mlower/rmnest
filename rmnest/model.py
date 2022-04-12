import numpy as np

from scipy import constants
from scipy.spatial.transform import Rotation


class FaradayRotation(object):
    """
    Fits a Faraday rotation model directly to the input Stokes Q and U spectra.

    See supplementary materials of Bannister et al. (2019) for details
    (arXiv:1906.11476)
    """

    def __init__(self, stokes_q, stokes_u, freq, freq_cen, pa_0, rm):
        """
        Parameters
        ----------
        stokes_q: array_like
            A list of Stokes Q flux or intensity measurements.
        stokes_u: array_like
            A list of Stokes U flux or intensity measurements.
        freq: array_like
            List of observing frequencies corresponding to each element in the
            input Stokes Q/U lists. (Hz)
        freq_cen: float
            Centre frequency of the observing band. (Hz)
        pa_0: float
            Linear polarisation position position angle corrsponding to
            the centre frequency. (deg)
        rm: float
            Rotation measure that encodes the strength of the Faraday effect.
            (rad m^-2)
        """
        self.stokes_q = stokes_q
        self.stokes_u = stokes_u
        self.freq = freq
        self.freq_cen = freq_cen
        self.pa_0 = pa_0
        self.rm = rm

    def compute_position_angle(self):
        """Polarisation position angle as a function of freq."""
        pa_freq = self.pa_0 + self.rm * (
            ((c / self.freq) ** 2) - ((c / (self.freq_cen)) ** 2)
        )
        return pa_freq

    def fit_QU(self):
        """Fits the linearly polarised emission"""
        pa = self.compute_position_angle()

        residuals = (
            (self.stokes_q**2)
            + (self.stokes_u**2)
            - (self.stokes_q * np.cos(2 * pa) + self.stokes_u * np.sin(2 * pa)) ** 2
        )

        return residuals


class GeneralisedFaradayRotation(object):
    """A phenomenological generalised Faraday rotation model.

    Parameters
    ----------
    freq : np.ndarray
        List of observing frequencies. (MHz)
    freq_cen : float
        Centre frequency of the observing band. (MHz)
    psi_0 : float
        Linear polarisation position angle corrsponding to the centre frequency (deg).
    grm : float
        Generalised rotation measure that encodes the strength of the Faraday effect (rad m^-alpha)
    alpha : float, optional
        Frequency scaling index, by default 2
    chi : float, optional
        Offset in the ellipticity angle, i.e the latitude of the poalrisation vector
        on the Poincaré sphere. (deg), by default 0
    phi : float, optional
        Rotation of the polarisation vector about the Stokes V axis. (deg) by default 0
    theta : float, optional
        Rotation of the poalrisation vector about the Stokes U axis. (deg) by default 0

    Notes
    -----
    Fits directly to the input Stokes Q, U and V spectra.
    See Lower (2020) for details (arXiv:2108.09429)
    """

    def __init__(
        self,
        freq: np.ndarray,
        freq_cen: float,
        psi_0: float,
        grm: float,
        alpha: float = 2,
        chi: float = 0,
        phi: float = 0,
        theta: float = 0,
    ) -> None:
        self.freq = freq
        self.freq_cen = freq_cen
        self.psi_0 = psi_0
        self.grm = grm
        self.alpha = alpha
        self.chi = chi
        self.phi = phi
        self.theta = theta

        # Model linear position angle
        psi = np.deg2rad(psi_0) + grm * (
            ((constants.c / (freq * constants.mega)) ** alpha)
            - ((constants.c / (freq_cen * constants.mega)) ** alpha)
        )

        # Model Stokes components
        stokes_q = np.cos(2 * psi) * np.cos(2 * np.deg2rad(chi))
        stokes_u = np.sin(2 * psi) * np.cos(2 * np.deg2rad(chi))
        stokes_v = np.sin(2 * np.deg2rad(chi))

        stokes_params = np.array(
            [stokes_q, stokes_u, np.repeat(stokes_v, len(self.freq))], dtype=float
        )
        # Rotation about the V and U axis
        rot = Rotation.from_euler("zy", [phi, theta], degrees=True)
        self._rotated_stokes = rot.apply(stokes_params.T, inverse=True).T

        self._m_psi = psi

    @property
    def rotated_stokes(self) -> np.ndarray:
        return self._rotated_stokes

    @property
    def m_q(self) -> np.ndarray:
        return self.rotated_stokes[0]

    @property
    def m_u(self) -> np.ndarray:
        return self.rotated_stokes[1]

    @property
    def m_v(self) -> np.ndarray:
        return self.rotated_stokes[2]

    @property
    def m_psi(self) -> np.ndarray:
        return self._m_psi
