import argparse
import bilby
import psrchive
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=False)

from utils import *
from likelihood import *


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
    parser.add_argument("--gfr", type=str2bool, default=False,
        help="Fit for generalised Faraday rotation (GFR).")
    parser.add_argument("--alpha", dest="free_alpha", type=str2bool,
        default=False, help="Use a free spectral dependence for GFR fitting.")
    return parser.parse_args()

class RMNest(object):
    def __init__(self, stokes_q, stokes_u, stokes_v, freqs, freq_cen):
        self.stokes_q = stokes_q
        self.stokes_u = stokes_u
        self.stokes_v = stokes_v
        self.freqs = freqs
        self.freq_cen = freq_cen


    def fit_rm(self, label="RM_Nest", outdir="./", **kwargs):
        """ Runs the rotation measure fitting routine. """
        self.priors = self._get_rm_priors()
        self.likelihood = RMLikelihood(
            self.stokes_q,
            self.stokes_u,
            self.freqs*1e6,
            self.freq_cen*1e6,
        )
        bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)
        result = bilby.run_sampler(
            likelihood=self.likelihood,
            priors=self.priors,
            sampler="dynesty",
            nlive=1024,
            outdir=outdir,
            plot=False,
            label=label,
            **kwargs
        )

        self.result = result
        self.post_json_file = bilby.result.result_file_name(outdir, result.label)

        posterior = result.posterior.rm
        median, low_bound, upp_bound = get_median_and_bounds(posterior)
        rm = median
        rm_upp = upp_bound - median
        rm_low = median - low_bound

        print("RM = {0} +{1}/-{2} rad/m^2 (68% CI)".format(rm, rm_upp, rm_low))


    def fit_gfr(self, label="RM_Nest", outdir="./", free_alpha=False, **kwargs):
        """ Runs the rotation measure fitting routine. """
        self.priors = self._get_gfr_priors(free_alpha)
        self.likelihood = GFRLikelihood(
            self.stokes_q,
            self.stokes_u,
            self.stokes_v,
            self.freq*1e6,
            self.freq_cen*1e6,
        )
        bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)
        result = bilby.run_sampler(
            likelihood=self.likelihood,
            priors=self.priors,
            sampler="dynesty",
            nlive=1024,
            outdir=outdir,
            plot=False,
            label=label,
            **kwargs
        )

        self.result = result
        self.post_json_file = bilby.result.result_file_name(outdir, result.label)

        posterior = result.posterior.rm
        median, low_bound, upp_bound = get_median_and_bounds(posterior)
        rm = median
        rm_upp = upp_bound - median
        rm_low = median - low_bound

        print("RM = {0} +{1}/-{2} rad/m^2 (68% CI)".format(rm, rm_upp, rm_low))


    def plot_corner(self):
        self.result.plot_corner(dpi=100)
        plt.close()


    @classmethod
    def from_psrchive(cls, ar_file, window, dedisperse=None, fscrunch=None):
        archive = psrchive.Archive_load(ar_file)
        archive.remove_baseline()
        if dedisperse is not None:
            archive.dedisperse()
        if fscrunch is not None:
            archive.fscrunch_to_nchan(fscrunch)

        nchan = archive.get_nchan()
        nbin = archive.get_nbin()

        window_start = int(float(window.split(":")[0]) * nbin)
        window_end = int(float(window.split(":")[1]) * nbin)

        window = [window_start, window_end]

        # Get weights and extract the on-pulse data from the archive
        data = apply_weights(archive.get_data()[0,:,:,:], archive.get_weights())
        on_pulse = np.mean(data[:,:,window[0]:window[1]], axis=2)

        # Extract Stokes I and find bad frequency channels
        Stokes_I = on_pulse[0,:]
        zeroed_chans = np.argwhere(Stokes_I <= 0.0)

        # Extract Stokes Q & U
        stokes_q = np.delete(on_pulse[1,:], zeroed_chans)
        stokes_u = np.delete(on_pulse[2,:], zeroed_chans)
        stokes_v = np.delete(on_pulse[3,:], zeroed_chans)

        # Get channel frequencies and centre frequency
        freqs = np.delete(archive.get_frequencies(), zeroed_chans)
        freq_cen = archive.get_centre_frequency()

        return cls(stokes_q, stokes_u, stokes_v, freqs, freq_cen)


    def _get_rm_priors(self):
        # Set bilby priors
        priors = bilby.prior.PriorDict()
        priors["rm"] = bilby.core.prior.Uniform(-2000, 2000,
            r"RM (rad m$^{-2}$)")
        priors["pa_zero"] = bilby.core.prior.Uniform(-np.pi/2, np.pi,
            r"$\Psi_{0}$ (deg)")
        priors["sigma"] = bilby.core.prior.Uniform(0, 1e4, r"$\sigma$")
        return priors


    def _get_gfr_priors(self, free_alpha):
        priors = bilby.prior.PriorDict()
        priors["psi_zero"] = bilby.core.prior.Uniform(-45, 45,
            r"$\Psi_{0} (deg)$")
        priors["chi"] = bilby.core.prior.Uniform(-90, 90, r"$\chi_{0} (deg)$")
        priors["phi"] = bilby.core.prior.Uniform(-180, 0, r"$\varphi (deg)$")
        priors["theta"] = bilby.core.prior.Uniform(0, 360, r"$\theta (deg)$")
        priors["sigma"] = bilby.core.prior.Uniform(0, 100, r"$\sigma$")

        # Set spectral dependency to be free, or fixed at freq^-3
        if free_alpha == "True":
            priors["alpha"] = bilby.core.prior.Uniform(0, 10, r"$\alpha$")
            priors["grm"] = bilby.core.prior.Uniform(0, 200,
                r"GRM (rad m$^{-\alpha}$)")
        else:
            priors["alpha"] = bilby.core.prior.DeltaFunction(3, r"$\alpha$")
            priors["grm"] = bilby.core.prior.Uniform(0, 200,
                r"GRM (rad m$^{-3}$)")

        return priors


def main():
    parser = argparse.ArgumentParser()
    args = get_input_arguments(parser)

    if args.archive is None:
        raise ValueError('No archive specified.')

    rmnest = RMNest.from_psrchive(
        args.archive,
        args.window,
        dedisperse=args.dedisperse,
        fscrunch=args.fscrunch
    )

    if args.gfr == True:
        if args.free_alpha == True:
            rmnest.fit_gfr(label=args.label, outdir=args.outdir,
                free_alpha=True)
        else:
            rmnest.fit_gfr(label=args.label, outdir=args.outdir)
    else:
        rmnest.fit_rm(label=args.label, outdir=args.outdir)

    rmnest.plot_corner()

    print("Done!")


# If run directly
if __name__ == "__main__":
    main()
