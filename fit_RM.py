import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=False)
import argparse
import os
import sys
import bilby
import psrchive

def get_input_arguments(parser):
    parser.add_argument("-a", dest="archive", type=str, default=None, 
        help="Input data, must be a PSRCHIVE format archive.")
    parser.add_argument("-o", dest="outdir", type=str, default="./", 
        help="Output destination.")
    parser.add_argument("-f", dest="fscrunch", type=int, default=None,
        help="Frequency scrunch data to this many channels.")
    parser.add_argument("-d", dest="dedisperse", type=str, default="False",
        help="Tell psrchive to dedisperse the data, default = False")
    parser.add_argument("-l", dest="label", type=str, default="RM_fit",
        help="Label added to output files.")
    parser.add_argument("--window", type=str, default="0.0:1.0", 
        help="Window to place around the pulse, default = 0.0:1.0")
    return parser.parse_args()


def get_median_and_bounds(posterior, nbins=80):
    pdf, vals = np.histogram(posterior, nbins)
    pdf = list(np.float_(pdf))

    pdf_normalised = pdf/np.sum(pdf)
    cdf = np.cumsum(pdf_normalised)

    median = vals[np.argmin(np.abs(cdf - 0.5))]
    low_bound = vals[np.argmin(np.abs(cdf - 0.16))]
    upp_bound = vals[np.argmin(np.abs(cdf - 0.84))]

    return median, low_bound, upp_bound


def get_rms(data, nbin):
    """ Returns the off pulse RMS """
    off_pulse = data[:int(0.2*nbin)]
    return np.sqrt((1/(0.2*nbin)) * (np.sum(off_pulse**2)))


def find_good_bins(data, nbin, thresh):
    """ Finds which bins contain signal > 3 * RMS """
    rms = get_rms(data, nbin)
    return np.argwhere(data > thresh*rms)


def position_angle(freq, pa_zero, rm, freq_cen):
    """ Polarisation position angle as a function of freq."""
    c = 2.998e8 # vacuum speed of light in m/s
    return pa_zero + rm*(((c/freq)**2) - ((c/(freq_cen*1e6))**2))


def fit_QU(freq, q, u, pa_zero, rm, freq_cen):
    """ Fits the linearly polarised emission """
    pa = position_angle(freq, pa_zero, rm, freq_cen)
    return (q**2) + (u**2) - (q*np.cos(2*pa) + u*np.sin(2*pa))**2


class RMLikelihood(bilby.likelihood.Likelihood):
    def __init__(self, q, u, freq, freq_cen, func):
        """
        The Gaussian likelihood from Bannister et al. (2019) - used for
        measuring pulsar/FRB rotation measures.

        Parameters
        ----------
        q, u: array_like
            The polarisation data to analyse
        freq: array_like
            Corresponding frequencies the data covers (Hz)
        freq_cen: float
            Centre frequency of the archive (Hz)
        func:
            The python function to fit to the data. Note, this must take the
            dependent variable as its first argument. The other arguments
            will require a prior and will be sampled over (unless a fixed
            value is given)
        sigma: None, float, array_like
            If None, the standard deviation of the noise is unknown and will be
            estimated (note: this requires a prior to be given for sigma). If
            not None, this defines the standard-deviation of the data points.
            This can either be a single float, or an array with length equal
            to that for `q` and `u`
        """

        super().__init__()
        self.q = q
        self.u = u
        self.freq = freq
        self.freq_cen = freq_cen
        self._func = func
        self.parameters = dict(pa_zero=None, rm=None, sigma=None)

    def log_likelihood(self):
        self.sigma = self.parameters["sigma"]

        self.residual = self._func(self.freq, self.q, self.u, 
            self.parameters["pa_zero"], self.parameters["rm"],
            self.freq_cen)

        ln_l = np.sum(-(self.residual/(2*(self.sigma**2))) -
            np.log(2*np.pi*(self.sigma**2)) / 2)
        return ln_l


def fit_rotation_measure(archive, outdir, label, nchan, nbin, window):
    """ Runs the rotation measure fitting routine. """

    # Extract the on-pulse data from the archive 
    data = archive.get_data()[0,:,:,:]
    on_pulse = np.mean(data[:,:,window[0]:window[1]], axis=2)

    # Extract Stokes I and find bad frequency channels
    Stokes_I = on_pulse[0,:]
    zeroed_chans = np.argwhere(Stokes_I <= 0.0)

    # Extract Stokes Q & U
    Stokes_Q = np.delete(on_pulse[1,:], zeroed_chans)
    Stokes_U = np.delete(on_pulse[2,:], zeroed_chans)

    # Get channel frequencies and centre frequency
    freqs = np.delete(archive.get_frequencies(), zeroed_chans)
    freq_cen = archive.get_centre_frequency()

    # Set bilby priors
    priors = dict()
    priors["rm"] = bilby.core.prior.Uniform(-2000, 2000, "RM")
    priors["pa_zero"] = bilby.core.prior.Uniform(-np.pi/2, np.pi, r"$\Psi_{0}$")
    priors["sigma"] = bilby.core.prior.Uniform(0, 1e4, r"$\sigma$")

    likelihood = RMLikelihood(Stokes_Q, Stokes_U, freqs*1e6, freq_cen*1e6, fit_QU)

    result = bilby.run_sampler(
        likelihood=likelihood, priors=priors, sampler="dynesty",
        nlive=1024, outdir=outdir, plot=False, label=label)

    result.plot_corner(dpi=100)
    plt.close()

    posterior = result.posterior.rm
    median, low_bound, upp_bound = get_median_and_bounds(posterior)
    rm = median
    rm_upp = upp_bound - median
    rm_low = median - low_bound

    print("RM = {0} +{1}/-{2} (68% CI)".format(rm, rm_upp, rm_low))


# If run directly
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    args = get_input_arguments(parser)

    if args.archive == None:
        raise ValueError('No archive specified.')
    else:
        archive = psrchive.Archive_load(args.archive)

    archive.remove_baseline()
    
    if args.dedisperse == None:
        pass
    else:
        archive.dedisperse()
    
    if args.fscrunch == None:
        pass
    else:
        archive.fscrunch_to_nchan(args.fscrunch)

    nchan = archive.get_nchan()
    nbin = archive.get_nbin()

    window_start = int(float(args.window.split(":")[0]) * nbin)
    window_end = int(float(args.window.split(":")[1]) * nbin)

    window = [window_start, window_end]

    fit_rotation_measure(archive, args.outdir, args.label, nchan, 
        nbin, window)

    print("Done!")
