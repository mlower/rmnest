import numpy as np

from scipy import constants
from scipy.spatial.transform import Rotation


class FaradayRotation(object):
    """
    Fits a Faraday rotation model directly to the input Stokes Q and U spectra.

    See supplementary materials of Bannister et al. (2019) for details
    (arXiv:)
    """
    def __init__(self,
        stokes_q,
        stokes_u,
        freq,
        freq_cen,
        pa_0,
        rm
    ):
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
        """ Polarisation position angle as a function of freq."""
        pa_freq = (self.pa_0 + self.rm*(((c/self.freq)**2)
            - ((c/(self.freq_cen))**2)))
        return pa_freq


    def fit_QU(self):
        """ Fits the linearly polarised emission """
        pa = self.compute_position_angle()

        residuals = ((self.stokes_q**2) + (self.stokes_u**2)
            - (self.stokes_q*np.cos(2*pa) + self.stokes_u*np.sin(2*pa))**2)

        return residuals


class GeneralisedFaradayRotation(object):
    """
    Fits a Faraday rotation model directly to the input Stokes Q and U spectra.

    See supplementary materials of Bannister et al. (2019) for details
    (arXiv:)
    """
    def __init__(self,
        freq,
        freq_cen,
        psi_0,
        alpha,
        grm,
        chi,
        phi,
        theta
    ):
        """
        Parameters
        ----------
        freq: array_like
            List of observing frequencies. (Hz)
        freq_cen: float
            Centre frequency of the observing band. (Hz)
        psi_0: float
            Linear polarisation position position angle corrsponding to
            the centre frequency. (deg)
        alpha: float
            Frequency scaling index.
        grm: float
            Generalised rotation measure that encodes the strength of the
            generalised Faraday effect. (rad m^-alpha)
        chi: float
            Offset in the ellipticity angle , i.e the latitude of the
            poalrisation vector on the Poincaré sphere. (deg)
        phi: float
            Rotation of the polarisation vector about the Stokes V axis. (deg)
        theta: float
            Rotation of the poalrisation vector about the Stokes U axis. (deg)
        """
        self.freq = freq
        self.freq_cen = freq_cen
        self.psi_0 = psi_0
        self.alpha = alpha
        self.grm = grm
        self.chi = chi
        self.phi = phi
        self.theta = theta

        # Model linear position angle
        psi = np.deg2rad(psi_0) + grm * (
            ((constants.c / (freq * constants.mega)) ** alpha)
            - ((constants.c / (freq_cen * constants.mega)) ** alpha)
        )
        self._psi = psi

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

    @property
    def psi(self):
        return self._psi

    @property
    def rotated_stokes(self):
        return self._rotated_stokes

    @property
    def m_q(self):
        return self.rotated_stokes[0]

    @property
    def m_u(self):
        return self.rotated_stokes[1]

    @property
    def m_v(self):
        return self.rotated_stokes[2]
