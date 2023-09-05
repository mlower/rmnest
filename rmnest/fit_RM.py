import bilby
import numpy as np

from rmnest import utils
from rmnest.likelihood import FRLikelihood, GFRLikelihood


class RMNest(object):
    def __init__(
        self, freqs, freq_cen, s_q, s_u, s_v, rms_q=None, rms_u=None, rms_v=None
    ):
        self.freqs = freqs
        self.freq_cen = freq_cen
        self.s_q = s_q
        self.s_u = s_u
        self.s_v = s_v
        self.rms_q = rms_q
        self.rms_u = rms_u
        self.rms_v = rms_v

    def fit(
        self,
        gfr=False,
        free_alpha=False,
        label="RM_Nest",
        outdir="./",
        sampler="dynesty",
        **kwargs,
    ):
        """Runs the rotation measure fitting routine."""
        if gfr:
            self.priors = self._get_gfr_priors(free_alpha)
            self.likelihood = GFRLikelihood(
                self.freqs,
                self.freq_cen,
                self.s_q,
                self.s_u,
                self.s_v,
                self.rms_q,
                self.rms_u,
                self.rms_v,
            )
        else:
            self.priors = self._get_fr_priors()
            self.likelihood = FRLikelihood(
                self.freqs, 
                self.freq_cen, 
                self.s_q, 
                self.s_u
            )

        bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)
        result = bilby.run_sampler(
            likelihood=self.likelihood,
            priors=self.priors,
            sampler=sampler,
            nlive=512,
            outdir=outdir,
            plot=False,
            label=label,
            **kwargs,
        )

        self.result = result
        self.post_json_file = bilby.result.result_file_name(outdir, result.label)

    def print_summary(self):
        for iparam, param in enumerate(self.result.search_parameter_keys):
            posterior = self.result.posterior[param]
            median, low_bound, upp_bound = utils.get_median_and_bounds(posterior)
            param_label = self.result.parameter_labels[iparam]
            print(
                f"{param_label} = {median} +{upp_bound - median}/-{median - low_bound} (68% CI)"
            )

    def print_bilby_summary(self, quantiles=(0.16, 0.84)):
        for iparam, param in enumerate(self.result.search_parameter_keys):
            param_summary = self.result.get_one_dimensional_median_and_error_bar(
                param, fmt=".2f", quantiles=quantiles
            )
            param_label = self.result.parameter_labels[iparam]
            conf_interval = quantiles[1] - quantiles[0]
            print_str = (
                f"{param_label} = {param_summary.median} +{param_summary.plus}/-{param_summary.minus}"
                f" ({int(conf_interval*100):d}% CI)"
            )
            print(print_str)

    def plot_corner(self):
        self.result.plot_corner(dpi=100)

    @classmethod
    def from_psrchive(cls, ar_file, window, dedisperse=False, fscrunch=None):
        import psrchive

        archive = psrchive.Archive_load(ar_file)
        archive.remove_baseline()
        if dedisperse:
            archive.dedisperse()
        if fscrunch is not None:
            archive.fscrunch_to_nchan(fscrunch)

        nbin = archive.get_nbin()

        window_start = int(float(window.split(":")[0]) * nbin)
        window_end = int(float(window.split(":")[1]) * nbin)

        window = [window_start, window_end]

        # Get weights and extract the on-pulse data from the archive
        data = utils.apply_weights(archive.get_data()[0, :, :, :], archive.get_weights())
        on_pulse = np.mean(data[:, :, window[0] : window[1]], axis=2)

        # Extract Stokes I and find bad frequency channels
        stokes_i = on_pulse[0, :]
        zeroed_chans = np.argwhere(stokes_i <= 0.0)

        # Extract Stokes Q & U
        stokes_q = np.delete(on_pulse[1, :], zeroed_chans)
        stokes_u = np.delete(on_pulse[2, :], zeroed_chans)
        stokes_v = np.delete(on_pulse[3, :], zeroed_chans)

        # Get channel frequencies and centre frequency
        freqs = np.delete(archive.get_frequencies(), zeroed_chans)
        freq_cen = archive.get_centre_frequency()

        return cls(freqs, freq_cen, stokes_q, stokes_u, stokes_v)

    @classmethod
    def from_stokesfile(cls, filename):
        spec = np.loadtxt(filename, unpack=True)
        if len(spec) == 9:
            freqs, s_i, rms_i, s_q, rms_q, s_u, rms_u, s_v, rms_v = spec
        elif len(spec) == 5:
            freqs, s_i, s_q, s_u, s_v = spec
            rms_q = rms_u = rms_v = None
        else:
            raise ValueError("Invalid number of columns in Stokes file.")
        freq_cen = np.median(freqs)
        print(f"Using freq_cen = {freq_cen}")
        return cls(freqs, freq_cen, s_q, s_u, s_v, rms_q, rms_u, rms_v)

    def _get_fr_priors(self):
        # Set bilby priors
        priors = bilby.prior.PriorDict()
        priors["rm"] = bilby.core.prior.Uniform(-2000, 2000, r"RM (rad m$^{-2}$)")
        priors["psi_zero"] = bilby.core.prior.Uniform(-90, 90, r"$\Psi_{0}$ (deg)")
        priors["sigma"] = bilby.core.prior.Uniform(0, 1e4, r"$\sigma$")
        return priors

    def _get_gfr_priors(self, free_alpha=False):
        priors = bilby.prior.PriorDict()

        # Set spectral dependency to be free, or fixed at freq^-3
        if free_alpha:
            priors["grm"] = bilby.core.prior.Uniform(0, 200, r"GRM (rad m$^{-\alpha}$)")
            priors["alpha"] = bilby.core.prior.Uniform(0, 10, r"$\alpha$")
        else:
            priors["grm"] = bilby.core.prior.Uniform(0, 200, "GRM (rad m$^{-3}$)")
            priors["alpha"] = bilby.core.prior.DeltaFunction(3, r"$\alpha$")

        priors["psi_zero"] = bilby.core.prior.Uniform(-90, 90, r"$\Psi_{0} (deg)$")
        priors["chi"] = bilby.core.prior.Uniform(-45, 45, r"$\chi (deg)$")
        priors["phi"] = bilby.core.prior.Uniform(-180, 180, r"$\varphi (deg)$")
        priors["theta"] = bilby.core.prior.Uniform(0, 180, r"$\vartheta (deg)$")
        priors["sigma"] = bilby.core.prior.Uniform(0, 100, r"$\sigma$")

        return priors
