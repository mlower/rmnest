# rmnest
Bayesian method for estimating rotation measures via direct fits to Stokes Q &amp; U

Based on the methodology outlined in the online supplementary materials of [Bannister et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B/abstract).

Dependencies: bilby &amp; PSRCHIVE Python-3 bindings.

Example usage on the test data:
```
python fit_RM.py -a test/2020-03-16-18\:12\:00.calib.ST -o test/output/ -l testrun --window 0.45:0.55 -f 128
```

If used, please cite
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
and add a link to this repository.
