#!/usr/bin/env python3
"""gtex_prep_app.py

GTEx RNAseq Preprocessing, for Sex As Biological Variable (SABV) analyses.

Command-line version; see also Jupyter notebook gtex_prep.ipynb

 - Author: Jeremy Yang
 - Based on R code by Oleg Ursu.
 - Required: Python 3.6+, Pandas

Workflow (prep):
 - READ: GTEx Subjects data, 1-row/subject.
 - READ: GTEx Samples data, 1-row/sample.
 - READ: GTEx RNAseq expression TPM data, 1-row/gene, 1-col/sample.
 - READ: gene IDs file, from GTEx/Ensembl/HGNC, via gtex_gene_xref.R. 
 - REMOVE: samples with Hardy score >2 (prefer healthier).
 - REMOVE: samples with high degree of autolysis (self-digestion).
 - MERGE: Samples and subjects, to 1-row/sample.
 - RESHAPE: RNAseq data from 1-col/sample, to 3 cols: gene, sample, TPM.
 - REMOVE: genes in pseudoautosomal regions (PAR) of chromosome Y.
 - AGGREGATE: samples, computing median TPM by gene+tissue.
 - AGGREGATE: samples, computing median TPM by gene+tissue+sex.

 - *TODO*: AGGREGATE: samples, computing median TPM by gene+tissue+sex+age.

 - OUTPUT: median TPMs, 1-row/gene+tissue+sex: 
 - OUTPUT: expression profiles, 1-row/gene+sex: exfiles_eps.tsv.gz

"""
#############################################################################
import sys,os,io,re,time,argparse,logging
import numpy,scipy,scipy.stats
import pandas as pd

#############################################################################
### (GTEx_v7_Annotations_SubjectPhenotypesDS.txt)
#############################################################################
def ReadSubjects(ifile):
  fin = open(ifile)
  logging.info(f"=== GTEx Subjects datafile: {fin.name}")
  subjects = pd.read_csv(fin, sep="\t")
  logging.info(f"Subjects dataset nrows: {subjects.shape[0]} ; ncols: {subjects.shape[1]}:")
  return subjects

#############################################################################
### Format: one line per tissue name, in preferred order.
#############################################################################
def ReadTissues(ifile):
  fin = open(ifile)
  logging.info(f"=== GTEx Tissues datafile: {fin.name}")
  tissues = pd.read_csv(fin, sep=";", index_col=False, header=None, names=["name"])
  tissues = tissues.name.str.strip()
  logging.info(f"n_tissues: {tissues.size}:")
  logging.info(f"tissues:\n{str(tissues)}")
  return tissues

#############################################################################
### Keep only healthier subjects: 
### (DTHHRDY = 4-point Hardy Scale Death Classification.)
### Keep 0, 1, 2 and reject 3, 4 and NA.
#############################################################################
def CleanSubjects(subjects):
  logging.info(f"=== Removing subjects with Hardy score > 2 or NA: {subjects[~(subjects.DTHHRDY<=2)].shape[0]}")
  subjects = subjects[subjects.DTHHRDY<=2]
  logging.info(f"Subjects dataset nrows: {subjects.shape[0]} ; ncols: {subjects.shape[1]}:")
  DescribeDf(subjects)
  return subjects

#############################################################################
def DescribeSubjects(subjects):
  logging.info("=== DescribeSubjects:")
  for name,val in subjects.AGE.value_counts().sort_index().iteritems():
    logging.info(f"\tAGE {name}: {val:4d}")
  for name,val in subjects.DTHHRDY.value_counts(sort=True, dropna=False).sort_index().iteritems():
    logging.info(f"\tDTHHRDY {name}: {val:4d}")

#############################################################################
def DescribeDf(df):
  buff = io.StringIO()
  df.info(buf=buff, verbose=True, null_counts=True)
  logging.info(re.sub(re.compile("^", re.M), "\t", buff.getvalue()))

#############################################################################
### (GTEx_v7_Annotations_SampleAttributesDS.txt)
#############################################################################
def ReadSamples(ifile):
  logging.info("=== ReadSamples:")
  fin = open(ifile)
  logging.info(f"GTEx Samples datafile: {fin.name}")
  samples = pd.read_csv(fin, sep="\t")
  samples = samples[["SAMPID", "SMATSSCR", "SMTS", "SMTSD"]]
  logging.info(f"Samples dataset nrows: {samples.shape[0]} ; ncols: {samples.shape[1]}:")
  ### SUBJID is first two hyphen-delimted fields of SAMPID.
  samples["SUBJID"] = samples.SAMPID.str.extract("^([^-]+-[^-]+)-", expand=False)
  DescribeDf(samples)
  return samples

#############################################################################
### Clean & tidy cols.
#############################################################################
def CleanSamples(samples):
  logging.info("=== CleanSamples:")
  samples_pre = samples.SMTSD.unique()
  logging.info(f"\tSamples: {samples.shape[0]}; tissues: {samples.SMTSD.nunique()}")
  logging.debug("\tsamples.SEX.apply()...")
  samples.SEX = samples.SEX.apply(lambda x: "F" if x==2 else "M" if x==1 else None)
  logging.debug("\tsamples.dropna(subset=['SEX'])...")
  if (samples.SEX.isna().sum()>0):
    samples.dropna(subset=["SEX"], inplace=True)
  logging.info(f"\tSamples: {samples.shape[0]}; tissues: {samples.SMTSD.nunique()}")
  ### Remove samples with severe degree of autolysis (self-digestion).
  ### NOTE that we keep SMATSSCR NAs.
  logging.debug("\t(samples.SMATSSCR!=3)&(samples.SMATSSCR!=2)...")
  samples = samples.loc[(samples.SMATSSCR!=3)&(samples.SMATSSCR!=2)]
  logging.info(f"\tSamples: {samples.shape[0]}; tissues: {samples.SMTSD.nunique()}")
  samples.loc[(samples.SMTS.str.strip()=="") & samples.SMTSD.str.startswith("Skin -"), "SMTS"] = "Skin"
  logging.info(f"\tSamples: {samples.shape[0]}; tissues: {samples.SMTSD.nunique()}")
  samples_post = samples.SMTSD.unique()
  samples_removed = set(samples_pre) - set(samples_post)
  if len(samples_removed)>0:
    logging.info("\tSamples removed: \n\t"+("\n\t".join(sorted(list(samples_removed)))))
  else:
    logging.info("\tSamples removed: (none)")
  return samples

#############################################################################
def DescribeSamples(samples):
  logging.info("=== DescribeSamples:")
  logging.info(f"Samples dataset nrows: {samples.shape[0]} ; ncols: {samples.shape[1]}:")
  for name,val in samples.SEX.value_counts().sort_index().iteritems():
    logging.info(f"\tSEX {name}: {val:4d}")
  logging.info(f"\tSamples: {samples.shape[0]}; tissues: {samples.SMTSD.nunique()}")
  #i=0
  #for name,val in samples.SMTSD.value_counts().sort_index().iteritems():
  #  i+=1
  #  logging.debug(f"\t{i}. '{name}': {val:4d}")

#############################################################################
def ReadRnaseq(ifile):
  """
READ GENE TPMs (full or demo subset)
 Top 2 rows, format:
	#1.2
	nrow	ncol
 Full file is ~56k rows, 2.6GB uncompressed.  Demo ~1k rows.
 *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
 *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz

 Truncate ENSGV version, use un-versioned ENSG for mapping. Ok?
  """
  logging.info("=== ReadRnaseq:")
  fin = open(ifile, "rb")
  logging.info(f"GTEx RNAseq TPM datafile: {fin.name}")
  rnaseq = pd.read_table(fin, compression="gzip", sep="\t", skiprows=2)
  logging.info(f"RNAseq dataset nrows: {rnaseq.shape[0]} ; ncols: {rnaseq.shape[1]}:")
  rnaseq = rnaseq.drop(columns=["Description"])
  rnaseq = rnaseq.rename(columns={"Name":"ENSG"})
  rnaseq.ENSG = rnaseq.ENSG.str.extract("^([^\.]+)\..*$", expand=False)
  samples = rnaseq.columns[1:]
  logging.info(f"RNAseq samples count: {samples.size}:")
  logging.info(f"RNAseq unique samples count: {samples.nunique()}:")
  logging.info(f"RNAseq genes (ENSG) count: {rnaseq.ENSG.size}:")
  logging.info(f"RNAseq unique genes (ENSG) count: {rnaseq.ENSG.nunique()}:")
  return rnaseq

#############################################################################
def ReadGenes(ifile):
  """Read gene IDs, etc.: ENSG,NCBI,HGNCID,symbol,name"""
  logging.info("=== ReadGenes:")
  fin = open(ifile)
  logging.info(f"GTEx/Ensembl/HGNC genes datafile: {fin.name}")
  genes = pd.read_csv(fin, sep="\t", na_values=[""], dtype={2:str})
  logging.info(f"Genes dataset nrows: {genes.shape[0]} ; ncols: {genes.shape[1]}:")
  #genes.columns = ["ENSG","NCBI","HGNC"]
  #genes.dropna(inplace=True)
  return genes

#############################################################################
def CleanRnaseq(rnaseq, keep_all_tissues):
  """
Memory intensive. Divide task to manage memory use.
For each tissue, group and concatenate results.
  """
  logging.info(f"CleanRnaseq IN: nrows = {rnaseq.shape[0]}, cols: {str(rnaseq.columns.tolist())}")
  logging.info("=== CleanRnaseq:")
  if not keep_all_tissues:
    logging.info("Remove sex-specific tissues...")
    tissues=list(rnaseq.SMTSD.sort_values().unique())
    sex_specific_tissues=[]
    for smtsd in tissues:
      if rnaseq[rnaseq.SMTSD==smtsd].SEX.nunique()<2:
        sex_specific_tissues.append(smtsd)
        logging.info(f"\tRemoving sex-specific tissue: \"{smtsd}\"")
        rnaseq = rnaseq[rnaseq.SMTSD!=smtsd]
    smtsd_breast = "Breast - Mammary Tissue"
    logging.info(f"Remove manually, not 100%% sex-specific: \"{smtsd_breast}\"...")
    rnaseq = rnaseq[rnaseq.SMTSD!=smtsd_breast] 

  logging.info("For each tissue, remove genes with TPMs all zero...")
  tissues=list(rnaseq.SMTSD.sort_values().unique())
  for i,smtsd in enumerate(tissues):
    rnaseq_this = rnaseq[rnaseq.SMTSD==smtsd]
    tpm_all0  = (rnaseq_this[["ENSG","TPM"]].groupby(by=["ENSG"], as_index=True).max()==0).rename(columns={"TPM":"tpm_all0"})
    n_all0 = (tpm_all0.tpm_all0.value_counts()[True] if True in tpm_all0.tpm_all0.value_counts() else 0)
    if n_all0>0:
      logging.info(f"\t{smtsd}: removing TPMs-all-zero genes: {n_all0}")
      rnaseq_this = pd.merge(rnaseq_this, tpm_all0, left_on=["ENSG"], right_index=True)
      rnaseq_this = rnaseq_this[~rnaseq_this["tpm_all0"]]
      rnaseq_this.drop(columns=["tpm_all0"], inplace=True)
    if i==0:
      rnaseq_out=rnaseq_this
    else:
      rnaseq_out=pd.concat([rnaseq_out,rnaseq_this])
  rnaseq = rnaseq_out

  rnaseq = rnaseq[["ENSG","SMTSD","SAMPID","SMATSSCR","SEX","AGE","DTHHRDY","TPM"]]
  rnaseq = rnaseq.sort_values(by=["ENSG","SMTSD","SAMPID"])
  rnaseq = rnaseq.reset_index(drop=True)
  logging.info(f"RNAseq unique samples count: {rnaseq.SAMPID.nunique()}:")
  logging.info(f"RNAseq unique tissues count: {rnaseq.SMTSD.nunique()}:")
  logging.info(f"RNAseq unique gene count: {rnaseq.ENSG.nunique()}")
  logging.info(f"CleanRnaseq OUT: nrows = {rnaseq.shape[0]}, cols: {str(rnaseq.columns.tolist())}")
  return rnaseq

#############################################################################
def Aggregate_median_SABV(rnaseq):
  """Compute median TPM by gene+tissue+sex."""
  logging.info("=== Aggregate_median_SABV:")
  logging.info(f"Aggregate_median_SABV IN: nrows = {rnaseq.shape[0]}, cols: {str(rnaseq.columns.tolist())}")
  rnaseq = rnaseq[["ENSG", "SMTSD", "SEX", "TPM"]].groupby(by=["ENSG", "SMTSD", "SEX"], as_index=False).median()
  logging.info(f"Aggregate_median_SABV OUT: nrows = {rnaseq.shape[0]}, cols: {str(rnaseq.columns.tolist())}")
  return rnaseq

#############################################################################
def Aggregate_median_SABV_AGE(rnaseq):
  """Compute median TPM by gene+tissue+sex+age."""
  logging.info("=== Aggregate_median_SABV_AGE:")
  logging.info(f"Aggregate_median_SABV_AGE IN: nrows = {rnaseq.shape[0]}, cols: {str(rnaseq.columns.tolist())}")
  rnaseq = rnaseq[["ENSG", "SMTSD", "SEX", "AGE", "TPM"]].groupby(by=["ENSG", "SMTSD", "SEX", "AGE"], as_index=False).median()
  logging.info(f"Aggregate_median_SABV_AGE OUT: nrows = {rnaseq.shape[0]}, cols: {str(rnaseq.columns.tolist())}")
  return rnaseq

#############################################################################
def PivotToProfiles(rnaseq):
  """
	Reshape to one-row-per-gene format.
	From:   ENSG,SMTSD,SEX,TPM,LOG_TPM
	To:	    ENSG,SEX,TPM_1,TPM_2,...TPM_N (N tissues)
	Preserve tissue order.
	Some missing TPM values, written as "" to TSV output.
  """
  logging.info(f"PivotToProfiles IN: nrows = {rnaseq.shape[0]}, cols: {str(rnaseq.columns.tolist())}")
  logging.info(f"PivotToProfiles tissue count: {rnaseq.SMTSD.nunique()}")
  tissues = pd.Series(pd.unique(rnaseq.SMTSD.sort_values()))

  # Assure only 1-row per unique (ensg,smtsd) tuple (or pivot will fail).
  #rnaseq = rnaseq.drop_duplicates(subset=["ENSG","SMTSD"], keep="first")

  rnaseq_f = rnaseq[rnaseq.SEX=="F"].drop(columns=["SEX"])
  rnaseq_m = rnaseq[rnaseq.SEX=="M"].drop(columns=["SEX"])
  rnaseq_f = rnaseq_f[["ENSG","SMTSD","TPM"]]
  rnaseq_m = rnaseq_m[["ENSG","SMTSD","TPM"]]
  exfiles_f = rnaseq_f.pivot(index="ENSG", columns="SMTSD")
  exfiles_f.columns = exfiles_f.columns.get_level_values(1)
  exfiles_f = exfiles_f.reset_index(drop=False)
  exfiles_f["SEX"] = "F"
  exfiles_m = rnaseq_m.pivot(index="ENSG", columns="SMTSD")
  exfiles_m.columns = exfiles_m.columns.get_level_values(1)
  exfiles_m = exfiles_m.reset_index(drop=False)
  exfiles_m["SEX"] = "M"
  exfiles = pd.concat([exfiles_f,exfiles_m])
  cols = ["ENSG","SEX"]+tissues.tolist()
  exfiles = exfiles[cols]
  DescribeDf(exfiles)
  logging.info(f"PivotToProfiles OUT: nrows = {exfiles.shape[0]}, cols: {str(exfiles.columns.tolist())}")
  return exfiles

#############################################################################
if __name__=="__main__":
  parser = argparse.ArgumentParser(description="GTEx Exfiles/SABV preprocessor")
  parser.add_argument("--i_subject", dest="ifile_subject", help="input subjects file")
  parser.add_argument("--i_sample", dest="ifile_sample", help="input samples file")
  parser.add_argument("--i_rnaseq", dest="ifile_rnaseq", help="input rnaseq file")
  parser.add_argument("--i_gene", dest="ifile_gene", help="input gene file")
  parser.add_argument("--o_median", dest="ofile_median", help="output median TPM, 1-row/gene+tissue+sex (TSV)")
  parser.add_argument("--o_median_sexage", dest="ofile_median_sexage", help="output median TPM, 1-row/gene+tissue+sex+age (TSV)")
  parser.add_argument("--o_sample", dest="ofile_sample", help="output sample TPM, 1-row/gene+sample (TSV)")
  parser.add_argument("--o_profiles", dest="ofile_profiles", help="output profiles, 1-row/gene+sex (TSV)")
  parser.add_argument("--o_tissue", dest="ofile_tissue", help="output tissues (TSV)")
  parser.add_argument("--decimals", type=int, default=3, help="output decimal places")
  parser.add_argument("--keep_all_tissues", action="store_true", help="normally remove reproductive+breast")
  parser.add_argument("-v", "--verbose", default=0, action="count")
  args = parser.parse_args()

  logging.basicConfig(format="%(levelname)s:%(message)s", level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  t0 = time.time()

  logging.info(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()))

  if args.verbose:
    logging.info(f"Python: {sys.version.split()[0]}; Pandas: {pd.__version__}; Scipy: {scipy.__version__} ; Numpy: {numpy.__version__}")

  if not args.ifile_subject:
    parser.error("Input subject file required.")

  subjects = ReadSubjects(args.ifile_subject)

  if args.verbose:
    DescribeSubjects(subjects)

  subjects = CleanSubjects(subjects)

  if not args.ifile_sample:
    parser.error("Input sample file required.")
  samples = ReadSamples(args.ifile_sample)

  logging.info("=== MERGE samples and subjects:")
  samples = pd.merge(samples, subjects, how="inner", on="SUBJID")

  if args.verbose:
    DescribeSamples(samples)

  if args.ofile_tissue:
    sample_tissues = samples[["SMTS", "SMTSD"]].reset_index(drop=True)
    sample_tissues = sample_tissues.drop_duplicates().sort_values(["SMTS", "SMTSD"])
    logging.info(f"=== Output tissues file: {args.ofile_tissue}")
    sample_tissues.round(args.decimals).to_csv(args.ofile_tissue, sep="\t", index=False)

  samples = CleanSamples(samples)

  if not args.ifile_gene:
    parser.error("Input gene file required.")
  genes = ReadGenes(args.ifile_gene)

  if not args.ifile_rnaseq:
    parser.error("Input RNAseq file required.")
  t1 = time.time()
  rnaseq = ReadRnaseq(args.ifile_rnaseq)
  logging.info(f"ReadRnaseq elapsed: {time.time()-t1}s")

  # Merge/inner with gene IDs file, to retain only protein-coding genes.
  rnaseq = pd.merge(rnaseq, genes[["ENSG"]], on="ENSG", how="inner")
  logging.info(f"RNAseq unique gene count (inner join with protein-coding gene ENSGs): {rnaseq.ENSG.nunique()}")

  logging.info("=== Remove genes in pseudoautosomal regions (PAR) of chromosome Y (ENSGR):")
  n_ensgr = rnaseq.ENSG.str.startswith("ENSGR").sum()
  logging.info(f"ENSGR gene TPMs: {n_ensgr} ({100*n_ensgr/rnaseq.shape[0]:.2f}%)")
  rnaseq = rnaseq[~rnaseq.ENSG.str.startswith("ENSGR")]
  logging.info(f"RNAseq unique gene count (after PAR removal): {rnaseq.ENSG.nunique()}")

  logging.info("=== MELT: One row per ENSG+SAMPID+TPM triplet:")
  ### Easier to handle but ~3x storage.
  rnaseq = rnaseq.melt(id_vars = "ENSG", var_name = "SAMPID", value_name = "TPM")
  DescribeDf(rnaseq)
  logging.info(f"RNAseq unique gene count (after melt): {rnaseq.ENSG.nunique()}")

  # Merge/inner with gene IDs file. This time to add IDs, names.
  rnaseq = pd.merge(rnaseq, genes, on="ENSG", how="left")
  logging.info(f"RNAseq unique gene count (after merge with gene IDs): {rnaseq.ENSG.nunique()}")

  logging.info("=== Merge with samples:")
  rnaseq = pd.merge(rnaseq, samples, how="inner", on="SAMPID")
  logging.info(f"RNAseq unique gene count (after merge with samples): {rnaseq.ENSG.nunique()}")
  logging.info(f"RNAseq unique tissue count (after merge with samples): {rnaseq.SMTSD.nunique()}")

  rnaseq = CleanRnaseq(rnaseq, args.keep_all_tissues)

  if args.ofile_sample:
    logging.info(f"=== Output sample TPM file: {args.ofile_sample}")
    rnaseq.round(args.decimals).to_csv(args.ofile_sample, sep="\t", index=False)

  if args.ofile_median_sexage:
    logging.info("=== Compute median TPM by gene+tissue+sex+age:")
    rnaseq_sexage = Aggregate_median_SABV_AGE(rnaseq)
    logging.info(f"SABV TPM median unique counts: genes: {rnaseq_sexage.ENSG.nunique()}")
    logging.info(f"=== Output median (by gene+tissue+sex+age) TPM file: {args.ofile_median_sexage}")
    rnaseq_sexage.round(args.decimals).to_csv(args.ofile_median_sexage, sep="\t", index=False)

  logging.info("=== Compute median TPM by gene+tissue+sex:")
  rnaseq = Aggregate_median_SABV(rnaseq)

  logging.info(f"SABV TPM median unique counts: genes: {rnaseq.ENSG.nunique()}")

  if args.ofile_median:
    logging.info(f"=== Output median (by gene+tissue+sex) TPM file: {args.ofile_median}")
    rnaseq.round(args.decimals).to_csv(args.ofile_median, sep="\t", index=False)

  logging.info("=== Pivot to one-row-per-gene format (profiles).")
  rnaseq_profiles = PivotToProfiles(rnaseq)
  if args.ofile_profiles:
    logging.info(f"=== Output profiles file: {args.ofile_profiles}")
    rnaseq_profiles.round(args.decimals).to_csv(args.ofile_profiles, sep="\t", index=False)

  logging.info(f"Elapsed: {time.time()-t0}s")
  logging.info(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()))
