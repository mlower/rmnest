import argparse
import bilby

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=False)

from rmnest.utils import *
from rmnest.likelihood import RMLikelihood, GFRLikelihood


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

    def fit(self, gfr=False, free_alpha=False, label="RM_Nest", outdir="./", sampler="dynesty", **kwargs):
        """ Runs the rotation measure fitting routine. """
        if gfr:
            self.priors = self._get_gfr_priors(free_alpha)
            self.likelihood = GFRLikelihood(
                self.stokes_q,
                self.stokes_u,
                self.stokes_v,
                self.freqs,
                self.freq_cen,
            )
        else:
            self.priors = self._get_rm_priors()
            self.likelihood = RMLikelihood(self.stokes_q, self.stokes_u, self.freqs, self.freq_cen)

        bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)
        result = bilby.run_sampler(
            likelihood=self.likelihood,
            priors=self.priors,
            sampler=sampler,
            nlive=512,
            outdir=outdir,
            plot=False,
            label=label,
            **kwargs
        )

        self.result = result
        self.post_json_file = bilby.result.result_file_name(outdir, result.label)

    def print_summary(self):
        for iparam, param in enumerate(self.result.search_parameter_keys):
            posterior = self.result.posterior[param]
            median, low_bound, upp_bound = get_median_and_bounds(posterior)
            param_label = self.result.parameter_labels[iparam]
            print(f"{param_label} = {0} +{upp_bound - median}/-{median - low_bound} (68% CI)")

    def print_bilby_summary(self, quantiles=(0.16, 0.84)):
        for iparam, param in enumerate(self.result.search_parameter_keys):
            posterior = self.result.posterior[param]
            param_summary = self.result.get_one_dimensional_median_and_error_bar(param, fmt=".2f", quantiles=quantiles)
            param_label = self.result.parameter_labels[iparam]
            conf_interval = quantiles[1] - quantiles[0]
            print(f"{param_label} = {param_summary.median} +{param_summary.plus}/-{param_summary.minus} ({int(conf_interval*100):d}% CI)")

    def plot_corner(self):
        self.result.plot_corner(dpi=100)
        plt.close()


    @classmethod
    def from_psrchive(cls, ar_file, window, dedisperse=None, fscrunch=None):
        import psrchive
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

    @classmethod
    def from_stokesfile(cls, stokes_file):
        spec_data = np.loadtxt(stokes_file)
        freqs = spec_data[:,0]
        freq_cen = np.median(freqs)
        print(f"Using freq_cen = {freq_cen}")
        return cls(spec_data[:,1], spec_data[:,2], spec_data[:,3], spec_data[:,4], freqs, freq_cen)

    def _get_rm_priors(self):
        # Set bilby priors
        priors = bilby.prior.PriorDict()
        priors["rm"] = bilby.core.prior.Uniform(-2000, 2000,
            r"RM (rad m$^{-2}$)")
        priors["psi_zero"] = bilby.core.prior.Uniform(-90, 90,
            r"$\Psi_{0}$ (deg)")
        priors["sigma"] = bilby.core.prior.Uniform(0, 1e4, r"$\sigma$")
        return priors


    def _get_gfr_priors(self, free_alpha=False):
        priors = bilby.prior.PriorDict()

        # Set spectral dependency to be free, or fixed at freq^-3
        if free_alpha == True:
            priors["grm"] = bilby.core.prior.Uniform(0, 200,
                r"GRM (rad m$^{-\alpha}$)")
            priors["alpha"] = bilby.core.prior.Uniform(0, 10, r"$\alpha$")
        else:
            priors["grm"] = bilby.core.prior.Uniform(0, 200,
                r"GRM (rad m$^{-3}$)")
            priors["alpha"] = bilby.core.prior.DeltaFunction(3, r"$\alpha$")

        priors["psi_zero"] = bilby.core.prior.Uniform(-90, 90,
            r"$\Psi_{0} (deg)$")
        priors["chi"] = bilby.core.prior.Uniform(-45, 45, r"$\chi_{0} (deg)$")
        priors["phi"] = bilby.core.prior.Uniform(-180, 180, r"$\varphi (deg)$")
        priors["theta"] = bilby.core.prior.Uniform(0, 180, r"$\theta (deg)$")
        priors["sigma"] = bilby.core.prior.Uniform(0, 100, r"$\sigma$")

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
    rmnest.fit(gfr=args.gfr, free_alpha=args.free_alpha, label=args.label, outdir=args.outdir)
    rmnest.print_summary()
    rmnest.plot_corner()

    print("Done!")


# If run directly
if __name__ == "__main__":
    main()
