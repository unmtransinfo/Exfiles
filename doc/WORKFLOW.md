# Exfiles Workflow

Steps for updating the Exfiles dataset from sources.

## Dependencies

* R 3.6+
* readr, wCorr, data.table

## Steps

1. Download RNA-Seq files from [GTEx portal](https://www.gtexportal.org/)
  1. RNA-Seq: [Gene TPMs](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz)
  1. Annotations: [Sample Attributes](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)
  1. Annotations: [Subject Phenotypes](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt)
1. Generate gene xrefs, from Ensembl/Biomart and HUGO.
    * [Go_gtex_GeneXref.sh](sh/Go_gtex_GeneXref.sh)
    * Download [protein-coding genes](ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding _gene.txt) from EBI.
    * Interactively download BIOMART ENSG2NCBI (human) mapping from Ensembl.org/biomart, with NCBI and HUGO IDs and HUGO symbols.  Select only those with Ensembl Protein Family IDs, for protein-encoding genes.
    * [gtex_gene_xref.R](R/gtex_gene_xref.R)
    * Get TCRD targets from IDG.
1. Process RNA-Seq data to expression profiles (Exfiles). ___WARNING: Requires big memory computer (120+GB).___
    * [Go_gtex_prep.sh](sh/Go_gtex_prep.sh)
    * [gtex_prep_app.py](python/gtex_prep_app.py)
        1. READ: GTEx Subjects data, 1-row/subject.
        1. READ: GTEx Samples data, 1-row/sample.
        1. READ: GTEx RNAseq expression TPM data, 1-row/gene, 1-col/sample.
        1. READ: gene IDs file, from GTEx/Ensembl/HGNC, via gtex_gene_xref.R. 
        1. REMOVE: samples with Hardy score >2 (prefer healthier).
        1. REMOVE: samples with high degree of autolysis (self-digestion).
        1. MERGE: Samples and subjects, to 1-row/sample.
        1. RESHAPE: RNAseq data from 1-col/sample, to 3 cols: gene, sample, TPM.
        1. REMOVE: genes in pseudoautosomal regions (PAR) of chromosome Y.
        1. AGGREGATE: samples, computing median TPM by gene+tissue.
        1. AGGREGATE: samples, computing median TPM by gene+tissue+sex.
        1. OUTPUT: median TPMs, 1-row/gene+tissue+sex: 
        1. OUTPUT: expression profiles, 1-row/gene+sex: exfiles_eps.tsv.gz
    * Requires ~3hr.
1. Compute pairwise similarity coefficients between all profiles.
    * [Go_exfiles_sim.sh](sh/Go_exfiles_sim.sh)
    * [exfiles_similarity_ruzicka.R](R/exfiles_similarity_ruzicka.R)
    * Requires ~20GB RAM, ~10 min.
1. Compute pairwise correlation coefficients between all profiles.
    * [Go_exfiles_cor.sh](sh/Go_exfiles_cor.sh)
    * [exfiles_similarity_wcorr.R](R/exfiles_similarity_wcorr.R)
    * Requires ~6hr.
1. Combine correlations and similarity into one file.
    * [Go_exfiles_post.sh](sh/Go_exfiles_post.sh)
    * [exfiles_similarity_post.py](python/exfiles_similarity_post.py)
