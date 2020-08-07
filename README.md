# expression-profiles
Gene expression profiles and analytics.
Developed by UNM and others of Team Helium in the NIH Data Commons Pilot Phase Consortium as part of Key Capability 8, Scientific Use Cases.
GTEx RNAseq pre-processing and analysis workflows initially developed by Oleg Ursu and Giovanni Bocci.
Jupyter notebooks and Python/Pandas port by Jeremy Yang.
Dockerization and AWS deployments by Mike Garcia.
See R Shiny app at http://unmtid-shinyapps.net/exfiles/.

## Dependencies

* R 3.6+
* R packages: readr, data.table, labdsv, wCorr
* [BioClients](https://github.com/jeremyjyang/BioClients) for `BioClients.lincs.Client_cmap`.
