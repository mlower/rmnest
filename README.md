# RMNest
A simple Bayesian inference package for estimating both standard and generalised rotation measures via direct fits to Stokes Q, U &amp; V spectra.

Originally based on the methodology outlined in the online supplementary materials of [Bannister et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B/abstract).

Installation:
```
pip install rmnest
```

Dependencies: bilby (PSRCHIVE Python-3 bindings strongly recommended).

Example usage on the test data:
```
python ./rmnest/fit_RM.py -a test/2020-03-16-18\:12\:00.calib.ST -o test/output/ -l testrun --window 0.45:0.55 -f 128
```

Fitting for generalised Faraday rotation (see e.g. [Kennett &amp; Melrose (1998)](https://ui.adsabs.harvard.edu/abs/1998PASA...15..211K/abstract) for a primer on this phenomenon) can be performed by adding the `--gfr` and `--alpha` flags as:
```
python ./rmnest/fit_RM.py -a <archive>.ar -o <outdir> -l testrun --window 0.45:0.55 --alpha True
```
Omitting the `--alpha` flag will result in the spectral exponent being fixed to 3. Details of the underlying phenomenological model will be presented in a forthcoming work.

The likelihood and Faraday rotation models, as well as the general `RMFit` class in `fit_RM.py`, can also be imported like any other API. 

If the standard RM-fitting is used, please cite
```
@ARTICLE{2020ApJ...896L..37L,
       author = {{Lower}, Marcus E. and {Shannon}, Ryan M. and {Johnston}, Simon and {Bailes}, Matthew},
        title = "{Spectropolarimetric Properties of Swift J1818.0-1607: A 1.4 s Radio Magnetar}",
      journal = {\apjl},
     keywords = {992, 1108, 1306, 1353, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2020,
        month = jun,
       volume = {896},
       number = {2},
          eid = {L37},
        pages = {L37},
          doi = {10.3847/2041-8213/ab9898},
archivePrefix = {arXiv},
       eprint = {2004.11522},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020ApJ...896L..37L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
and 
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
and include a [link](https://github.com/mlower/rmnest) to this repository.


