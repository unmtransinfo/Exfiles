#!/usr/bin/env python3
"""gtex_rnaseq_prep_app.py

GTEx RNAseq Preprocessing

Command-line version; see also Jupyter notebook gtex_rnaseq_prep.ipynb

 - Author: Jeremy Yang
 - Based on R code by Oleg Ursu.
 - Required: Python3, Pandas 0.22+

Workflow:
 - READ: GTEx Subjects data
 - READ: GTEx Samples data
 - READ: GTEx RNAseq expression TPM data
 - READ: Ensembl gene IDs file ENSG2NCBI, from Ensembl.
 - REMOVE: samples with Hardy score >2 (prefer healthier).
 - MERGE: Samples and subjects, to one row per sample.
 - RESHAPE: RNAseq data from 1 col/sample, to 3 cols: gene, sample, TPM.
 - REMOVE: genes in pseudoautosomal regions (PAR) of chromosome Y.
 - COMPUTE: median TPM by gene+tissue.
 - COMPUTE: LOG10(TPM+1)
 - COMPUTE: TAU, tissue specificity index (Yanai et al., 2004).
 - STRATIFY: RNAseq data by sex.
 - RESHAPE: into one row per gene+tissue, TPM_F, TPM_M.
 - COMPUTE: median TPM by gene+tissue+sex.
 - COMPUTE: Log fold-change, log of ratio (F/M).

 - SAVE: files for downstream SABV expression profile analytics.

"""
#############################################################################
import sys,os,io,re,time,argparse
import numpy,scipy,scipy.stats
import pandas

#############################################################################
### (GTEx_v7_Annotations_SubjectPhenotypesDS.txt)
#############################################################################
def ReadSubjects(ifile, verbose):
  fin = open(ifile)
  print('=== GTEx Subjects datafile: %s'%fin.name, file=sys.stdout)
  subjects = pandas.read_csv(fin, sep='\t', index_col='SUBJID')
  print("Subjects dataset nrows: %d ; ncols: %d:"%(subjects.shape[0],subjects.shape[1]), file=sys.stdout)
  return subjects

#############################################################################
### Keep only healthier subjects: 
### (DTHHRDY = 4-point Hardy Scale Death Classification.)
#############################################################################
def CleanSubjects(subjects, verbose):
  print("=== Subjects with Hardy score > 2 or NA: %d (removing)"%(subjects.query('DTHHRDY > 2').shape[0]), file=sys.stdout)
  subjects = subjects.query('DTHHRDY <= 2')
  print("Subjects dataset nrows: %d ; ncols: %d:"%(subjects.shape[0],subjects.shape[1]), file=sys.stdout)
  DescribeDf(subjects, verbose)
  return subjects

#############################################################################
def DescribeSubjects(subjects):
  print("=== DescribeSubjects:", file=sys.stdout)
  for name,val in subjects.AGE.value_counts().sort_index().iteritems():
    print('\tAGE %s: %4d'%(name,val), file=sys.stdout)
  for name,val in subjects.DTHHRDY.value_counts(sort=True, dropna=False).sort_index().iteritems():
    print('\tDTHHRDY %s: %4d'%(name,val), file=sys.stdout)

#############################################################################
def DescribeDf(df, verbose):
  buff = io.StringIO()
  df.info(buf=buff,verbose=bool(verbose>0),null_counts=bool(verbose>0))
  print(re.sub(re.compile('^', re.M), '\t', buff.getvalue()), file=sys.stdout)

#############################################################################
### (GTEx_v7_Annotations_SampleAttributesDS.txt)
#############################################################################
def ReadSamples(ifile, verbose):
  print("=== ReadSamples:", file=sys.stdout)
  fin = open(ifile)
  print('GTEx Samples datafile: %s'%fin.name, file=sys.stdout)
  samples = pandas.read_csv(fin, sep='\t', index_col='SAMPID')
  samples = samples[['SMATSSCR', 'SMTS', 'SMTSD']]
  print("Samples dataset nrows: %d ; ncols: %d:"%(samples.shape[0],samples.shape[1]), file=sys.stdout)
  DescribeDf(samples, verbose)
  ### SUBJID is first two hyphen-delimted fields of SAMPID.
  samples['SUBJID'] = samples.index
  samples['SUBJID'] = samples.SUBJID.str.extract('^([^-]+-[^-]+)-', expand=True)
  return samples

#############################################################################
### Clean & tidy cols. Remove samples with high degree of autolysis (self-digestion)."""
#############################################################################
def CleanSamples(samples, verbose):
  print("=== CleanSamples:", file=sys.stdout)
  samples.dropna(how='any', inplace=True)
  samples.SEX = samples.SEX.apply(lambda x: 'female' if x==2 else 'male' if x==1 else None)
  samples = samples[samples.SMATSSCR < 2]
  samples.loc[(samples.SMTS.str.strip()=='') & samples.SMTSD.str.startswith("Skin -"), 'SMTS'] = 'Skin'
  print("Samples dataset nrows: %d ; ncols: %d:"%(samples.shape[0],samples.shape[1]), file=sys.stdout)
  return samples

#############################################################################
def DescribeSamples(samples):
  print("=== DescribeSamples:", file=sys.stdout)
  for name,val in samples.SEX.value_counts().sort_index().iteritems():
    print('\tSEX %s: %4d'%(name,val), file=sys.stdout)
  for name,val in (samples.SMTS+" : "+samples.SMTSD).value_counts().sort_index().iteritems():
    print('\tSMTS:SMTSD %s: %4d'%(name,val), file=sys.stdout)

#############################################################################
### READ GENE TPMs (full or demo subset)
### Full file is ~56k rows, 2.6GB uncompressed.  Demo ~1k rows.
### *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
### *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz
#############################################################################
def ReadRnaseq(ifile, verbose):
  print("=== ReadRnaseq:", file=sys.stdout)
  fin = open(ifile, "rb")
  print('GTEx RNAseq TPM datafile: %s'%fin.name, file=sys.stdout)
  t0 = time.time()
  rnaseq = pandas.read_table(fin, compression='gzip', sep='\t', skiprows=2)
  print("RNAseq dataset nrows: %d ; ncols: %d:"%(rnaseq.shape[0],rnaseq.shape[1]), file=sys.stdout)
  print("ReadRnaseq elapsed: %ds"%(time.time()-t0), file=sys.stdout)
  rnaseq = rnaseq.drop(columns=['Description'])
  rnaseq = rnaseq.rename(columns={'Name':'ENSG'})
  return rnaseq

#############################################################################
### Read gene symbols.
#############################################################################
def ReadGenes(ifile, verbose):
  print("=== ReadGenes:", file=sys.stdout)
  fin = open(ifile)
  print('Biomart ENSG2NCBI genes datafile: %s'%fin.name, file=sys.stdout)
  genes = pandas.read_csv(fin, sep='\t', usecols=[1,2,3], na_values=[''], dtype={2:str})
  print("Genes dataset nrows: %d ; ncols: %d:"%(genes.shape[0],genes.shape[1]), file=sys.stdout)
  genes.columns = ['ENSG','NCBI','HGNC']
  genes.dropna(inplace=True)
  return genes

#############################################################################
### Remove data for gene-tissue pairs with all zero expression.
### Remove data for gene-tissue pairs not present in both sexes.
#############################################################################
def CleanRnaseq(rnaseq, verbose):
  print("=== CleanRnaseq:", file=sys.stdout)
  max_tpm_0 = (rnaseq[['ENSG', 'SMTSD', 'TPM']].groupby(by=['ENSG','SMTSD'], as_index=True).max() == 0).rename(columns={'TPM':'max_tpm_0'})
  #print(max_tpm_0.max_tpm_0.value_counts(), file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, max_tpm_0, left_on=['ENSG', 'SMTSD'], right_index=True)
  rnaseq = rnaseq[~rnaseq['max_tpm_0']]
  rnaseq.drop(columns=['max_tpm_0'], inplace=True)

  sex_count = (rnaseq[['ENSG', 'SMTSD', 'SEX']].groupby(by=['ENSG','SMTSD'], as_index=True).nunique()).rename(columns={'SEX':'sex_count'})
  #print(sex_count.sex_count.value_counts(), file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, sex_count, left_on=['ENSG', 'SMTSD'], right_index=True)
  rnaseq = rnaseq[rnaseq['sex_count'] == 2]
  rnaseq.drop(columns=['sex_count'], inplace=True)
  return rnaseq

#############################################################################
def SABV_stratify(rnaseq_level, verbose):
  print("=== SABV_stratify:", file=sys.stdout)
  print("=== Assign levels (Low-Med-High) to ranks, cutoff quantiles .25, .75:", file=sys.stdout)
  rnaseq_level['LEVEL'] = rnaseq_level.TPM_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')

  rnaseq_level['AGE'] = 'ALL'
  rnaseq_level['SEX'] = 'ALL'

  ### Compute median TPM by gene+tissue+sex:
  rnaseq_sex = rnaseq[['ENSG', 'SMTSD', 'SEX', 'TPM']].groupby(by=['ENSG','SMTSD','SEX'], as_index=False).median()
  #print(rnaseq_sex.shape)

  ### Combine rows into one row per gene+tissue, cols for M and F TPM.
  rnaseq_sex_f = rnaseq_sex.loc[rnaseq_sex['SEX'] == 'female']
  rnaseq_sex_f = rnaseq_sex_f[['ENSG', 'SMTSD', 'TPM']].rename(columns={'TPM':'TPM_F'})
  rnaseq_sex_m = rnaseq_sex.loc[rnaseq_sex['SEX'] == 'male']
  rnaseq_sex_m = rnaseq_sex_m[['ENSG', 'SMTSD', 'TPM']].rename(columns={'TPM':'TPM_M'})
  rnaseq_sex = pandas.merge(rnaseq_sex_f, rnaseq_sex_m, how='inner', on=['ENSG','SMTSD'])
  return rnaseq_sex

#############################################################################
def SABV_analyze(rnaseq_sex, verbose):
  print("=== SABV_analyze:", file=sys.stdout)
  print("=== Assign gene-tissue rank (quantile) among tissues (F):", file=sys.stdout)
  t0 = time.time()
  rnaseq_level_f = GTRanks(rnaseq_sex[['ENSG','SMTSD','TPM_F']].copy(), 'TPM_F')
  print("GTRanks (F) elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  rnaseq_level_f['LEVEL_F'] = rnaseq_level_f.TPM_F_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')
  print("=== Assign gene-tissue rank (quantile) among tissues (M):", file=sys.stdout)
  t0 = time.time()
  rnaseq_level_m = GTRanks(rnaseq_sex[['ENSG','SMTSD','TPM_M']].copy(), 'TPM_M')
  print("GTRanks (M) Elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  rnaseq_level_m['LEVEL_M'] = rnaseq_level_m.TPM_M_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')

  rnaseq_level_sex = pandas.merge(rnaseq_level_f, rnaseq_level_m, on = ['ENSG','SMTSD'], how = 'inner')

  print("=== Compute Log fold-change, log of ratio (F/M):", file=sys.stdout)
  rnaseq_level_sex['log2foldchange'] = ((rnaseq_level_sex.TPM_F+1) / (rnaseq_level_sex.TPM_M+1)).apply(lambda x: numpy.log2(max(x, 1/x)))

  #WilcoxonSignedRank(rnaseq_level_sex, rnaseq_sex, verbose)

  return rnaseq_level_sex

#############################################################################
def WilcoxonSignedRank(rnaseq_level, rnaseq_sex, verbose):
  print("=== WilcoxonSignedRank:", file=sys.stdout)
  ### For each gene, compute sex difference via Wilcoxon signed-rank test,
  ### with Wilcox treatment, discarding all zero-differences.

  wilcox = pandas.DataFrame({'ENSG':rnaseq_sex.ENSG.drop_duplicates().sort_values(), 'stat':None, 'pval':None}).reset_index(drop=True)

  for i in range(wilcox.shape[0]):
    tpm_f_this = rnaseq_sex.TPM_F[rnaseq_level.ENSG == wilcox.ENSG[i]]
    tpm_m_this = rnaseq_sex.TPM_M[rnaseq_level.ENSG == wilcox.ENSG[i]]
    stat, pval = scipy.stats.wilcoxon(x=tpm_f_this, y=tpm_m_this, zero_method='wilcox')
    wilcox.stat.iloc[i] = stat
    wilcox.pval.iloc[i] = pval 

  # Output?

#############################################################################
### Compute tissue specificity index (Yanai et al., 2004).
### $ \tau = \frac{\sum_{i=0}^N (1 - x_i)}{N - 1} $
### * N = number of tissues
### * x = expression profile component normalized by the maximal component value
###
### Validate with example vector from paper.  Should be 0.95.
### print('%.2f'%TAU([0,8,0,0,0,2,0,2,0,0,0,0]), file=sys.stdout)
#############################################################################
def TAU(X):
  N = len(X)
  xmax = max(X)
  if xmax==0: return(0.0)
  tau = 0.0
  for x in X:
    tau += (1 - x/xmax)
  tau /= (N - 1)
  return(tau)

#############################################################################
### Assign gene-tissue rank (quantile) among tissues.
### Ranks, for given gene, tissue expression from gene-tissue TPMs.
#############################################################################
def GTRanks(rnaseq, tpm_col):
  tpm_rank = pandas.Series(dtype="float", index=range(rnaseq.shape[0]))
  for i in rnaseq.index:
    ensg = rnaseq.ENSG[i]
    val_this = rnaseq[tpm_col][i]
    vals_ensg = rnaseq[tpm_col][rnaseq.ENSG==ensg]
    vals_ensg = vals_ensg.sort_values().reset_index(drop=True)
    j = vals_ensg[vals_ensg == val_this].index[0]
    tpm_rank.iloc[i] = j/vals_ensg.size 

  rnaseq[tpm_col+'_RANK'] = tpm_rank
  return(rnaseq)

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='GTEx RNAseq Exfiles/SABV preprocessor')
  parser.add_argument("--i_subject",dest="ifile_subject",help="input subjects file")
  parser.add_argument("--i_sample",dest="ifile_sample",help="input samples file")
  parser.add_argument("--i_rnaseq",dest="ifile_rnaseq",help="input rnaseq file")
  parser.add_argument("--i_gene",dest="ifile_gene",help="input gene file")
  parser.add_argument("--o_median",dest="ofile_median",help="output median TPM (TSV)")
  parser.add_argument("--o_tau",dest="ofile_tau",help="output TAU (TSV)")
  parser.add_argument("--o_level",dest="ofile_level",help="output level TPM (TSV)")
  parser.add_argument("--o_median_sex",dest="ofile_median_sex",help="output SABV median TPM (TSV)")
  parser.add_argument("--o_level_sex",dest="ofile_level_sex",help="output SABV level TPM (TSV)")
  parser.add_argument("--decimals",type=int,default=3,help="output decimal places")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  if not args.ifile_subject:
    parser.error('Input subject file required.')

  if args.verbose:
    print('Python: %s\nPandas: %s'%(sys.version,pandas.__version__), file=sys.stdout)

  subjects = ReadSubjects(args.ifile_subject, args.verbose)

  if args.verbose:
    DescribeSubjects(subjects)

  subjects = CleanSubjects(subjects, args.verbose)

  if not args.ifile_sample:
    parser.error('Input sample file required.')
  samples = ReadSamples(args.ifile_sample, args.verbose)

  print('=== MERGE samples and subjects:', file=sys.stdout)
  samples = pandas.merge(samples, subjects, how='inner', left_on='SUBJID', right_index=True)

  if args.verbose:
    DescribeSamples(samples)

  samples = CleanSamples(samples, args.verbose)

  if not args.ifile_gene:
    parser.error('Input gene file required.')
  genes = ReadGenes(args.ifile_gene, args.verbose)

  if not args.ifile_rnaseq:
    parser.error('Input RNAseq file required.')
  rnaseq = ReadRnaseq(args.ifile_rnaseq, args.verbose)

  print('=== MELT: One row per ENSG+SAMPID+TPM triplet:', file=sys.stdout)
  ### Easier to handle but ~3x storage.
  rnaseq = rnaseq.melt(id_vars = "ENSG", var_name = "SAMPID", value_name = "TPM")
  DescribeDf(rnaseq,args.verbose)

  rnaseq = pandas.merge(genes, rnaseq, on='ENSG', how='inner')

  print('=== Remove genes in pseudoautosomal regions (PAR) of chromosome Y ("ENSGR"):', file=sys.stdout)
  n_ensgr = rnaseq.ENSG.str.startswith('ENSGR').sum()
  print('ENSGR gene TPMs: %d (%.2f%%)'%(n_ensgr,100*n_ensgr/rnaseq.shape[0]), file=sys.stdout)
  rnaseq = rnaseq[~rnaseq.ENSG.str.startswith('ENSGR')]

  print('=== Merge with samples:', file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, samples, how="inner", left_on="SAMPID", right_index=True)

  rnaseq = CleanRnaseq(rnaseq, args.verbose)

  print('=== Compute median TPM by gene+tissue:', file=sys.stdout)
  rnaseq_med = rnaseq[['ENSG', 'SMTSD', 'TPM']].groupby(by=['ENSG','SMTSD'], as_index=False).median()
  print("RNAseq unique counts: genes: %d ; tissues: %d ; gene-tissue pairs: %d"%(rnaseq_med.ENSG.nunique(), rnaseq_med.SMTSD.nunique(), rnaseq_med.shape[0]), file=sys.stdout)

  if args.ofile_median:
    print("=== Output TPM median file: %s"%args.ofile_median, file=sys.stdout)
    rnaseq_med.round(parser.decimals).to_csv(args.ofile_median, sep='\t', index=False)

  print('=== Compute LOG10(TPM+1):', file=sys.stdout)
  rnaseq_med['LOG_TPM'] = rnaseq_med.TPM.apply(lambda x: numpy.log10(x+1))

  print('=== Compute TAU (tissue specificity, Yanai et al., 2004):', file=sys.stdout)
  rnaseq_tau = rnaseq_med.groupby(['ENSG']).TPM.agg(TAU)
  rnaseq_tau = pandas.DataFrame(rnaseq_tau).rename(columns={'TPM':'TAU'})
  if args.ofile_tau:
    print("=== Output TAU file: %s"%args.ofile_tau, file=sys.stdout)
    rnaseq_tau.round(parser.decimals).to_csv(args.ofile_tau, sep='\t') #Index is ENSG

  print('=== Compute gene-tissue ranks (GTRanks):', file=sys.stdout)
  t0 = time.time()
  rnaseq_level = GTRanks(rnaseq_med.copy(), 'TPM')
  print("GTRanks elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  print("TPM level unique counts: genes: %d"%(rnaseq_level.ENSG.nunique()), file=sys.stdout)

  if args.ofile_level:
    print("=== Output TPM level file: %s"%args.ofile_level, file=sys.stdout)
    rnaseq_level.round(parser.decimals).to_csv(args.ofile_level, sep='\t', index=False)

  print('=== Compute median TPM by gene+tissue+sex:', file=sys.stdout)
  rnaseq_sex = SABV_stratify(rnaseq_level, args.verbose)
  print("SABV TPM median unique counts: genes: %d"%(rnaseq_sex.ENSG.nunique()), file=sys.stdout)

  if args.ofile_median_sex:
    print("=== Output SABV TPM median file: %s"%args.ofile_median_sex, file=sys.stdout)
    rnaseq_sex.round(parser.decimals).to_csv(args.ofile_median_sex, sep='\t', index=False)

  rnaseq_level_sex = SABV_analyze(rnaseq_sex, args.verbose)
  print("SABV TPM level unique counts: genes: %d"%(rnaseq_level_sex.ENSG.nunique()), file=sys.stdout)

  if args.ofile_level_sex:
    print("=== Output SABV TPM level file: %s"%args.ofile_level_sex, file=sys.stdout)
    rnaseq_level_sex.round(parser.decimals).to_csv(args.ofile_level_sex, sep='\t', index=False)
