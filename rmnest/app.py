import click

from rmnest.fit_RM import RMNest


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "show_default": True}
)
@click.argument(
    "archive",
    type=click.Path(exists=True),
    help="Input data, must be a PSRCHIVE format archive.",
)
@click.option(
    "-o",
    "--outdir",
    type=click.Path(exists=True),
    default="./",
    help="Output destination.",
)
@click.option(
    "-f", "--fscrunch", type=int, help="Frequency scrunch data to this many channels."
)
@click.option(
    "--window",
    type=str,
    default="0.0:1.0",
    help="Window to place around the pulse, default = 0.0:1.0.",
)
@click.option("--label", type=str, default="RM_fit", help="Label added to output files.")
@click.option("--dedisperse", is_flag=True, help="Tell psrchive to dedisperse the data.")
@click.option("--gfr", is_flag=True, help="Fit for generalised Faraday rotation (GFR).")
@click.option(
    "--free_alpha", is_flag=True, help="Use a free spectral dependence for GFR fitting."
)
def main(archive, outdir, fscrunch, window, label, dedisperse, gfr, free_alpha):
    if archive is None:
        raise ValueError("No archive specified.")

    rmnest = RMNest.from_psrchive(
        archive, window, dedisperse=dedisperse, fscrunch=fscrunch
    )
    rmnest.fit(gfr=gfr, free_alpha=free_alpha, label=label, outdir=outdir)
    rmnest.print_summary()
    rmnest.plot_corner()

    print("Done!")


if __name__ == "__main__":
    main()
