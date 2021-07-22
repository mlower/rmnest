import numpy as np

c = 2.998e8 # vacuum speed of light in m/s


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


    def get_Psi(self):
        """ Model linear position angle """
        pa_freq = (self.psi_0 + self.grm*(((c/self.freq)**self.alpha)
            - ((c/(self.freq_cen))**self.alpha)))
        return pa_freq

    
    def Q_model(self):
        """ Stokes Q model """
        return np.cos(2*self.psi) * np.cos(2*self.chi)


    def U_model(self):
        """ Stokes U model """
        return np.sin(2*self.psi) * np.cos(2*self.chi)

    def V_model(self):
        """ Stokes V model """
        return np.sin(2*self.chi)


    def rotate_V(self):
        """ Rotation about the V axis """
        return np.array([
            [np.cos(self.phi), -np.sin(self.phi), 0],
            [np.sin(self.phi), np.cos(self.phi), 0],
            [0, 0, 1]])


    def rotate_U(self):
        """ Rotation about the U axis """
        return np.array([
            [np.cos(self.theta), 0, np.sin(self.theta)],
            [0, 1, 0],
            [-np.sin(self.theta), 0, np.cos(self.theta)]])


    def generate_gfr_model(self):
        """ Generate the generalised Faraday rotation model """

        self.psi = self.get_Psi()
        stokes_params = np.transpose(
            np.array([self.Q_model(), self.U_model(), self.V_model()])
            )

        rotated_stokes = np.transpose(
            stokes_params.dot(self.rotate_U()).dot(self.rotate_V())
        )

        return rotated_stokes[0], rotated_stokes[1], rotated_stokes[2]
