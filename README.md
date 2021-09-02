# `Exfiles` (expression-profiles)

Gene expression profiles and analytics.
Developed by UNM for the NIH Data Commons Pilot Phase Consortium as part of Key Capability 8, Scientific Use Cases.
GTEx RNAseq pre-processing and analysis workflows initially developed by Oleg Ursu and Giovanni Bocci.
Jupyter notebooks and Python/Pandas port by Jeremy Yang.
Dockerization and Cloud/AWS tooling by Mike Garcia.
See R Shiny app at https://unmtid-shinyapps.net/exfiles/.

## Dependencies

* Python 3.6+; Pandas, [BioClients](https://github.com/jeremyjyang/BioClients)
* R 3.6+; readr, data.table, labdsv, wCorr

## Workflow

See [WORKFLOW.md](doc/WORKFLOW.md)

## Glossary

| Term | Description |
|---:|:---|
| TPM | Transcripts per million |
| Log2FoldChange | LOG2(TPM\_F/TPM\_M) |
| Tau | Tissue specificity index (Yanai et al., 2004) |

## See also:

* [GTEx Documentation](https://www.gtexportal.org/home/documentationPage)
* [GTEx Datasets](https://www.gtexportal.org/home/datasets) (with data dictionaries).
