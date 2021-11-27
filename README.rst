|PyPI| |Python| |License|


*RMNest* is an open source python package for estimating both standard and generalised 
rotation measures via direct fits to Stokes *Q*, *U* and *V* spectra.

Installation
------------
The latest release of *RMNest* can be installed from `PyPi`_ by running 
the following::

    pip install rmnest

Note that while a working installation of the PSRCHIVE Python-3 bindings is
not necessary for using *RMNest*, it is strongly recommended.

.. _PyPI: https://pypi.python.org/pypi/rmnest

Requirements
------------
The following packages are required to running *RMNest*.

 - numpy: Array manipulation

 - matplotlib: Modules for plotting

 - bilby: Inference calculations framework
 
 - dynesty: Modules for nested sampling
 
 
Usage
-----
*RMNest* can be run directly from the command line within using the `fit_RM.py` script. 
As an example, the below command would run a standard rotation-measure fit on the provided test data after frequency-averaging to 128 channels 
within a [pulse] phase window between phase = 0.45 to 0.55:
::

  python ./rmnest/fit_RM.py -a test/2020-03-16-18\:12\:00.calib.ST -o test/output/ -l testrun --window 0.45:0.55 -f 128

Alternatively, fitting for the generalised form of Faraday rotation, sometimes referred to as Faraday conversion 
(see e.g. `Kennett & Melrose (1998)`_), can be performed by adding the ``--gfr`` and ``--alpha`` flags as:
::

  python ./rmnest/fit_RM.py -a <archive>.ar -o <outdir> -l testrun --window 0.45:0.55 --alpha True

Omitting the `--alpha` flag will result in the spectral exponent being fixed to 3. Details of the underlying phenomenological model can be found in a technical document by [Lower (2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv210809429L).

The likelihood and Faraday rotation models, as well as the general `RMFit` class in `fit_RM.py`, can also be imported like any other API.
    

.. _Kennett & Melrose 1998: https://ui.adsabs.harvard.edu/abs/1998PASA...15..211K/abstract


Issues and Contributing
-----------------------
If you encounter any issues with *RMNest*, or have in mind a feature that 
currently does not exist, then you can contribute by openning a `Github Issue`_ and 
outlining the feature. 

.. _Github Issue: https://github.com/mlower/rmnest/issues

Referencing RMNest
------------------
If you make use of *RMNest* in your research, we would greatly appreciate it if you 
cite the papers behind its development.

For instance, if you only make use of the the standard rotation measure fitting, then
please cite both `Bannister et al. (2019)`_ and `Lower et al. (2020)`_, and include
a `link to this repository`_.

.. _Bannister et al. (2019): https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B
.. _Lower et al. (2020): https://ui.adsabs.harvard.edu/abs/2020ApJ...896L..37L
.. _link to this repository: https://github.com/mlower/rmnest

::

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

and
::

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


Alternatively, if the generalised Faraday rotation fitting is used, please include 
a citation to `Lower (2021)`_.

.. _Lower (2021): https://ui.adsabs.harvard.edu/abs/2021arXiv210809429L

::

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


.. |Logo| image:: logo/fruitbat_logo.svg
    :alt: Fruitbat Logo

.. |PyPI| image:: https://img.shields.io/pypi/v/rmnest.svg?label=PyPI
    :target: https://pypi.python.org/pypi/rmnest
    :alt: PyPI - Latest Release

.. |Python| image:: https://img.shields.io/pypi/pyversions/rmnest.svg?label=Python
    :target: https://pypi.python.org/pypi/rmnest
    :alt: PyPI - Python Versions

.. |License| image:: https://img.shields.io/pypi/l/rmnest.svg?colorB=purple&label=License
    :target: https://github.com/abatten/rmnest/raw/master/LICENSE
    :alt: PyPI - License
