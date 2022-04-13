import click

from rmnest.fit_RM import RMNest


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "show_default": True}
)
def main() -> None:
    pass


@main.command()
@click.argument("ar_file", type=click.Path(exists=True))
@click.option(
    "-o", "--outdir", type=click.Path(), default="./", help="Output destination."
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
@click.option("--label", type=str, default="RM_Nest", help="Label added to output files.")
@click.option("--dedisperse", is_flag=True, help="Tell psrchive to dedisperse the data.")
@click.option("--gfr", is_flag=True, help="Fit for generalised Faraday rotation (GFR).")
@click.option(
    "--free_alpha", is_flag=True, help="Use a free spectral dependence for GFR fitting."
)
def archive(ar_file, outdir, fscrunch, window, label, dedisperse, gfr, free_alpha):
    """Fit RMNest to an archive file."""
    rmnest = RMNest.from_psrchive(
        ar_file, window, dedisperse=dedisperse, fscrunch=fscrunch
    )
    rmnest.fit(gfr=gfr, free_alpha=free_alpha, label=label, outdir=outdir)
    rmnest.print_summary()
    rmnest.plot_corner()

    print("Done!")


@main.command()
@click.argument("stokes_file", type=click.Path(exists=True))
@click.option(
    "-o", "--outdir", type=click.Path(), default="./", help="Output destination."
)
@click.option("--label", type=str, default="RM_Nest", help="Label added to output files.")
@click.option("--gfr", is_flag=True, help="Fit for generalised Faraday rotation (GFR).")
@click.option(
    "--free_alpha", is_flag=True, help="Use a free spectral dependence for GFR fitting."
)
def txtfile(stokes_file, outdir, label, gfr, free_alpha):
    """Fit RMNest to a text file with columns: freq, I, Q, U, V."""
    rmnest = RMNest.from_stokesfile(stokes_file)
    rmnest.fit(gfr=gfr, free_alpha=free_alpha, label=label, outdir=outdir)
    rmnest.print_summary()
    rmnest.plot_corner()

    print("Done!")


if __name__ == "__main__":
    main()
