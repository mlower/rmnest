<a href="https://ascl.net/2204.008"><img src="https://img.shields.io/badge/ascl-2204.008-blue.svg?colorB=262255" alt="ascl:2204.008" /></a>
[![PyPI](https://img.shields.io/pypi/v/rmnest.svg?label=PyPI)](https://pypi.python.org/pypi/rmnest)
[![Python](https://img.shields.io/pypi/pyversions/rmnest.svg?label=Python)](https://pypi.python.org/pypi/rmnest)
[![License](https://img.shields.io/pypi/l/rmnest.svg?colorB=purple&label=License)](https://github.com/mlower/rmnest/blob/main/LICENSE)

# RMNest

*RMNest* is an open source python package for estimating both standard and generalised
rotation measures via direct fits to Stokes *Q*, *U* and *V* spectra.

## Installation

The latest release of *RMNest* can be installed from [PyPi](https://pypi.python.org/pypi/rmnest) by running
the following

```bash
pip install rmnest
```

Note that while a working installation of the PSRCHIVE Python-3 bindings is
not necessary for using *RMNest*, it is strongly recommended.

## Requirements

The following packages are required to running *RMNest*.

- numpy: Array manipulation

- matplotlib: Modules for plotting

- bilby: Inference calculations framework

- dynesty: Modules for nested sampling

## Usage

*RMNest* can be run directly from the command line within using `rmnest`.
As an example, the below command would run a standard rotation-measure fit on the provided test data after frequency-averaging to 128 channels
within a [pulse] phase window between phase = 0.45 to 0.55

```bash
rmnest archive test/2020-03-16-18\:12\:00.calib.ST -o test/output/ -l testrun --window 0.45:0.55 -f 128
```

Alternatively, fitting for the generalised form of Faraday rotation, sometimes referred to as Faraday conversion
(see e.g. [Kennett & Melrose 1998](https://ui.adsabs.harvard.edu/abs/1998PASA...15..211K/abstract)), can be performed
by adding the ``--gfr`` and ``--free_alpha`` flags as

```bash
rmnest <archive>.ar -o <outdir> -l testrun --window 0.45:0.55 --gfr --free_alpha
```

Omitting the `--free_alpha` flag will result in the spectral exponent being fixed to 3. Details of the underlying phenomenological model can be
found in a technical document by [Lower (2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv210809429L).

The likelihood and Faraday rotation models, as well as the general `RMFit` class in `fit_RM.py`, can also be imported like any other API.

## Issues and Contributing

If you encounter any issues with *RMNest*, or have in mind a feature that
currently does not exist, then you can contribute by openning a
[Github Issue](https://github.com/mlower/rmnest/issues) and outlining the feature.

## Referencing RMNest

If you make use of *RMNest* in your research, we would greatly appreciate it if you
cite both the ASCL entry ([Lower et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ascl.soft04008L))
and the papers behind its development.

```
@software{2022ascl.soft04008L,
       author = {{Lower}, Marcus E. and {Kumar}, Pravir and {Shannon}, Ryan M.},
        title = "{RMNest: Bayesian approach to measuring Faraday rotation and conversion in radio signals}",
     keywords = {Software},
 howpublished = {Astrophysics Source Code Library, record ascl:2204.008},
         year = 2022,
        month = apr,
          eid = {ascl:2204.008},
        pages = {ascl:2204.008},
archivePrefix = {ascl},
       eprint = {2204.008},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022ascl.soft04008L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

For standard rotation measure fitting, then
please cite [Bannister et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B).

```
@ARTICLE{2019Sci...365..565B,
       author = {{Bannister}, K.~W. and {Deller}, A.~T. and {Phillips}, C. and {Macquart}, J. -P. and {Prochaska}, J.~X. and {Tejos}, N. and {Ryder}, S.~D. and {Sadler}, E.~M. and {Shannon}, R.~M. and {Simha}, S. and {Day}, C.~K. and {McQuinn}, M. and {North-Hickey}, F.~O. and {Bhandari}, S. and {Arcus}, W.~R. and {Bennert}, V.~N. and {Burchett}, J. and {Bouwhuis}, M. and {Dodson}, R. and {Ekers}, R.~D. and {Farah}, W. and {Flynn}, C. and {James}, C.~W. and {Kerr}, M. and {Lenc}, E. and {Mahony}, E.~K. and {O'Meara}, J. and {Os{\l}owski}, S. and {Qiu}, H. and {Treu}, T. and {U}, V. and {Bateman}, T.~J. and {Bock}, D.~C. -J. and {Bolton}, R.~J. and {Brown}, A. and {Bunton}, J.~D. and {Chippendale}, A.~P. and {Cooray}, F.~R. and {Cornwell}, T. and {Gupta}, N. and {Hayman}, D.~B. and {Kesteven}, M. and {Koribalski}, B.~S. and {MacLeod}, A. and {McClure-Griffiths}, N.~M. and {Neuhold}, S. and {Norris}, R.~P. and {Pilawa}, M.~A. and {Qiao}, R. -Y. and {Reynolds}, J. and {Roxby}, D.~N. and {Shimwell}, T.~W. and {Voronkov}, M.~A. and {Wilson}, C.~D.},
        title = "{A single fast radio burst localized to a massive galaxy at cosmological distance}",
      journal = {Science},
     keywords = {ASTRONOMY, Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2019,
        month = aug,
       volume = {365},
       number = {6453},
        pages = {565-570},
          doi = {10.1126/science.aaw5903},
archivePrefix = {arXiv},
       eprint = {1906.11476},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

If you used *RMNest* for generalised Faraday rotation measure fitting, please include
a citation to [Lower (2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv210809429L).

```
@ARTICLE{2021arXiv210809429L,
       author = {{Lower}, Marcus E.},
        title = "{A phenomenological model for measuring generalised Faraday rotation}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - High Energy Astrophysical Phenomena},
         year = 2021,
        month = aug,
          eid = {arXiv:2108.09429},
        pages = {arXiv:2108.09429},
archivePrefix = {arXiv},
       eprint = {2108.09429},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv210809429L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
