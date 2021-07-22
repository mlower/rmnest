import numpy as np


def get_rc_params():
    """
    Get the rcParams that will be used in all the plots.
    Run: plt.rcParams.update(get_rc_params())
    """

    rc_params = {
        "text.usetex": True,
        "font.family": "serif",

        "figure.dpi": 125,
        "legend.fontsize": 12,
        "legend.frameon": False,
        "legend.markerscale": 1.0,

        "axes.labelsize": 14,

        "xtick.direction": 'in',
        "xtick.labelsize": 14,
        "xtick.minor.visible": True,
        "xtick.top": True,
        "xtick.major.width": 1,

        "ytick.direction": 'in',
        "ytick.labelsize": 14,
        "ytick.minor.visible": True,
        "ytick.right": True,
        "ytick.major.width": 1,

        "savefig.transparent": False,
    }

    return rc_params


def apply_weights(data, weights, pol=True):
    """ Apply weights to zero RFI affected channels """

    bins = int(np.shape(data)[-1])
    mask = weights.T * np.ones((1,bins))

    if pol == True:
        for i in range(0,4):
            data[i,:,:] = np.multiply(mask, data[i,:,:])
        return data
    else:
        return np.multiply(mask, data)


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


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
