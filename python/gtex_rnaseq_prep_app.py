#!/usr/bin/env python3
"""gtex_rnaseq_prep_app.py

GTEx RNAseq Preprocessing

Command-line version; see also Jupyter notebook gtex_rnaseq_prep.ipynb

 - Author: Jeremy Yang
 - Based on R code by Oleg Ursu.
 - Required: Python3, Pandas 0.22+
 - Clean, tidy, reshape RNAseq expression data.
 - Compute tissue specificity index (Yanai et al., 2004).
 - Save files for downstream SABV expression profile analytics.

"""
#############################################################################
import sys,os,re,time,io,argparse
#import urllib.request
import numpy,scipy,scipy.stats
import pandas

#############################################################################
### (GTEx_v7_Annotations_SubjectPhenotypesDS.txt)
#############################################################################
def ReadSubjects(ifile, verbose):
  fin = open(ifile)
  print('GTEx Subjects datafile: %s'%fin.name, file=sys.stderr)
  subjects = pandas.read_csv(fin, sep='\t', index_col='SUBJID')
  print("Subjects dataset nrows: %d ; ncols: %d:"%(subjects.shape[0],subjects.shape[1]), file=sys.stderr)
  return subjects

#############################################################################
### Keep only healthier subjects: 
### (DTHHRDY = 4-point Hardy Scale Death Classification.)
#############################################################################
def CleanSubjects(subjects, verbose):
  print("Subjects with Hardy score > 2 or NA: %d (removing)"%(subjects.query('DTHHRDY > 2').shape[0]), file=sys.stderr)
  subjects = subjects.query('DTHHRDY <= 2')
  print("Subjects dataset nrows: %d ; ncols: %d:"%(subjects.shape[0],subjects.shape[1]), file=sys.stderr)
  subjects.info(buf=sys.stderr,verbose=False)
  return subjects

#############################################################################
def DescribeSubjects(subjects):
    for name,val in subjects.AGE.value_counts().sort_index().iteritems():
      print('\tAGE %s: %4d'%(name,val), file=sys.stderr)
    for name,val in subjects.DTHHRDY.value_counts(sort=True, dropna=False).sort_index().iteritems():
      print('\tDTHHRDY %s: %4d'%(name,val), file=sys.stderr)

#############################################################################
### (GTEx_v7_Annotations_SampleAttributesDS.txt)
#############################################################################
def ReadSamples(ifile, verbose):
  print("DEBUG: in ReadSamples...", file=sys.stderr)
  fin = open(ifile)
  print('GTEx Samples datafile: %s'%fin.name, file=sys.stderr)
  samples = pandas.read_csv(fin, sep='\t', index_col='SAMPID')
  samples = samples[['SMATSSCR', 'SMTS', 'SMTSD']]
  print("Samples dataset nrows: %d ; ncols: %d:"%(samples.shape[0],samples.shape[1]), file=sys.stderr)
  samples.info(buf=sys.stderr,verbose=False)
  ### SUBJID is first two hyphen-delimted fields of SAMPID.
  samples['SUBJID'] = samples.index
  samples['SUBJID'] = samples.SUBJID.str.extract('^([^-]+-[^-]+)-', expand=True)
  return samples

#############################################################################
### Clean & tidy cols. Remove samples with high degree of autolysis (self-digestion)."""
#############################################################################
def CleanSamples(samples, verbose):
  print("DEBUG: in CleanSamples...", file=sys.stderr)
  samples.dropna(how='any', inplace=True)
  samples.SEX = samples.SEX.apply(lambda x: 'female' if x==2 else 'male' if x==1 else None)
  samples = samples[samples.SMATSSCR < 2]

  print("DEBUG: CleanSamples, fixing Skin - ...", file=sys.stderr)
  samples.loc[(samples.SMTS.str.strip()=='') & samples.SMTSD.str.startswith("Skin -"), 'SMTS'] = 'Skin'


  print("Samples dataset nrows: %d ; ncols: %d:"%(samples.shape[0],samples.shape[1]), file=sys.stderr)
  return samples

#############################################################################
def DescribeSamples(samples):
  for name,val in samples.SEX.value_counts().sort_index().iteritems():
    print('\tSEX %s: %4d'%(name,val), file=sys.stderr)
  for name,val in (samples.SMTS+" : "+samples.SMTSD).value_counts().sort_index().iteritems():
    print('\tSMTS:SMTSD %s: %4d'%(name,val), file=sys.stderr)

#############################################################################
### READ GENE TPMs (full or demo subset)
### Full file is ~56k rows, 2.6GB uncompressed.  Demo ~1k rows.
### *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
### *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz
#############################################################################
def ReadRnaseq(ifile, verbose):
  fin = open(ifile, "rb")
  print('GTEx RNAseq TPM datafile: %s'%fin.name, file=sys.stderr)
  t0 = time.time()
  rnaseq = pandas.read_table(fin, compression='gzip', sep='\t', skiprows=2)
  print("RNAseq dataset nrows: %d ; ncols: %d:"%(rnaseq.shape[0],rnaseq.shape[1]), file=sys.stderr)
  print("ReadRnaseq elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  rnaseq.info(buf=sys.stderr,verbose=False)
  rnaseq = rnaseq.drop(columns=['Description'])
  rnaseq = rnaseq.rename(columns={'Name':'ENSG'})
  return rnaseq

#############################################################################
### Read gene symbols.
#############################################################################
def ReadGenes(ifile, verbose):
  fin = open(ifile)
  print('Biomart ENSG2NCBI genes datafile: %s'%fin.name, file=sys.stderr)
  genes = pandas.read_csv(fin, sep='\t', usecols=[1,2,3], na_values=[''], dtype={2:str})
  print("Genes dataset nrows: %d ; ncols: %d:"%(genes.shape[0],genes.shape[1]), file=sys.stderr)
  genes.columns = ['ENSG','NCBI','HGNC']
  genes.dropna(inplace=True)
  return genes

#############################################################################
### Remove data for gene-tissue pairs with all zero expression.
### Remove data for gene-tissue pairs not present in both sexes.
#############################################################################
def CleanRnaseq(rnaseq, verbose):
  max_tpm_0 = (rnaseq[['ENSG', 'SMTSD', 'TPM']].groupby(by=['ENSG','SMTSD'], as_index=True).max() == 0).rename(columns={'TPM':'max_tpm_0'})
  print(max_tpm_0.max_tpm_0.value_counts(), file=sys.stderr)
  rnaseq = pandas.merge(rnaseq, max_tpm_0, left_on=['ENSG', 'SMTSD'], right_index=True)
  rnaseq = rnaseq[~rnaseq['max_tpm_0']]
  rnaseq.drop(columns=['max_tpm_0'], inplace=True)

  sex_count = (rnaseq[['ENSG', 'SMTSD', 'SEX']].groupby(by=['ENSG','SMTSD'], as_index=True).nunique()).rename(columns={'SEX':'sex_count'})
  print(sex_count.sex_count.value_counts(), file=sys.stderr)
  rnaseq = pandas.merge(rnaseq, sex_count, left_on=['ENSG', 'SMTSD'], right_index=True)
  rnaseq = rnaseq[rnaseq['sex_count'] == 2]
  rnaseq.drop(columns=['sex_count'], inplace=True)
  return rnaseq

#############################################################################
def SABV(rnaseq_level, verbose):
  rnaseq_level['LEVEL'] = rnaseq_level.TPM_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')
  rnaseq_level['AGE'] = 'ALL'
  rnaseq_level['SEX'] = 'ALL'

  ### Compute median TPM by gene+tissue+SEX:
  rnaseq_med_sex = rnaseq[['ENSG', 'SMTSD', 'SEX', 'TPM']].groupby(by=['ENSG','SMTSD','SEX'], as_index=False).median()
  print(rnaseq_med_sex.shape)

  ### Combine rows into one row per gene+tissue, cols for M and F TPM.
  rnaseq_med_sex_f = rnaseq_med_sex.loc[rnaseq_med_sex['SEX'] == 'female']
  rnaseq_med_sex_f = rnaseq_med_sex_f[['ENSG', 'SMTSD', 'TPM']].rename(columns={'TPM':'TPM_F'})
  rnaseq_med_sex_m = rnaseq_med_sex.loc[rnaseq_med_sex['SEX'] == 'male']
  rnaseq_med_sex_m = rnaseq_med_sex_m[['ENSG', 'SMTSD', 'TPM']].rename(columns={'TPM':'TPM_M'})
  rnaseq_med_sex = pandas.merge(rnaseq_med_sex_f, rnaseq_med_sex_m, how='inner', on=['ENSG','SMTSD'])

  t0 = time.time()
  rnaseq_level_f = GTRanks(rnaseq_med_sex[['ENSG','SMTSD','TPM_F']].copy(), 'TPM_F')
  print("GTRanks elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  rnaseq_level_f['LEVEL_F'] = rnaseq_level_f.TPM_F_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')

  t0 = time.time()
  rnaseq_level_m = GTRanks(rnaseq_med_sex[['ENSG','SMTSD','TPM_M']].copy(), 'TPM_M')
  print("Elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  rnaseq_level_m['LEVEL_M'] = rnaseq_level_m.TPM_M_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')

  rnaseq_level = pandas.merge(rnaseq_level_f, rnaseq_level_m, on = ['ENSG','SMTSD'], how = 'inner')

  ### For each gene, compute sex difference via Wilcox test:
  ### (Wilcoxon signed-rank test, with Wilcox treatment, discarding all zero-differences.)


  wilcox = pandas.DataFrame({'ENSG':rnaseq_med_sex.ENSG.drop_duplicates().sort_values(), 'stat':None, 'pval':None}).reset_index(drop=True)

  for i in range(wilcox.shape[0]):
    tpm_f_this = rnaseq_med_sex.TPM_F[rnaseq_level.ENSG == wilcox.ENSG[i]]
    tpm_m_this = rnaseq_med_sex.TPM_M[rnaseq_level.ENSG == wilcox.ENSG[i]]
    stat, pval = scipy.stats.wilcoxon(x=tpm_f_this, y=tpm_m_this, zero_method='wilcox')
    wilcox.stat.iloc[i] = stat
    wilcox.pval.iloc[i] = pval 

  ### Log fold-change is log of ratio.

  rnaseq_level['log2foldchange'] = ((rnaseq_level.TPM_F+1) / (rnaseq_level.TPM_M+1)).apply(lambda x: numpy.log2(max(x, 1/x)))


#############################################################################
### Compute tissue specificity index (Yanai et al., 2004).
### $ \tau = \frac{\sum_{i=0}^N (1 - x_i)}{N - 1} $
### * N = number of tissues
### * x = expression profile component normalized by the maximal component value
###
### Validate with example vector from paper.  Should be 0.95.
### print('%.2f'%TAU([0,8,0,0,0,2,0,2,0,0,0,0]), file=sys.stderr)
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
### Assign gene-tissue rank (quantile of median) among tissues.
### Low-Med-High cutoff quantiles: .25 and .75.  These ranks measure, for a
### given gene, relative tissue expression from the gene-tissue TPMs.
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

  parser = argparse.ArgumentParser(
        description='GTEx RNAseq Exfiles/SABV preprocessor')
  #ops = ['cid2Activity', 'tid2Targetcomponents']
  #parser.add_argument("op",choices=ops,help='operation')
  parser.add_argument("--i_subject",dest="ifile_subject",help="input subjects file")
  parser.add_argument("--i_sample",dest="ifile_sample",help="input samples file")
  parser.add_argument("--i_rnaseq",dest="ifile_rnaseq",help="input rnaseq file")
  parser.add_argument("--i_gene",dest="ifile_gene",help="input gene file")
  parser.add_argument("--o_median",dest="ofile_median",help="output median TPM (TSV)")
  parser.add_argument("--o_level",dest="ofile_level",help="output level TPM (TSV)")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  if not args.ifile_subject:
    parser.error('Input subject file required.')

  print('Python: %s\nPandas: %s'%(sys.version,pandas.__version__), file=sys.stderr)

  subjects = ReadSubjects(args.ifile_subject, args.verbose)

  if args.verbose:
    DescribeSubjects(subjects)

  subjects = CleanSubjects(subjects, args.verbose)

  if not args.ifile_sample:
    parser.error('Input sample file required.')
  samples = ReadSamples(args.ifile_sample, args.verbose)

  ### MERGE samples and subjects:"""
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

  ### MELT: One row per ENSG+SAMPID+TPM triplet:
  ### Easier to handle but ~3x storage.
  rnaseq = rnaseq.melt(id_vars = "ENSG", var_name = "SAMPID", value_name = "TPM")
  rnaseq.info(buf=sys.stderr,verbose=False)

  rnaseq = pandas.merge(genes, rnaseq, on='ENSG', how='inner')

  ### Remove genes in pseudoautosomal regions (PAR) of chromosome Y ("ENSGR").
  n_ensgr = rnaseq.ENSG.str.startswith('ENSGR').sum()
  print('ENSGR gene TPMs: %d (%.2f%%)'%(n_ensgr,100*n_ensgr/rnaseq.shape[0]), file=sys.stderr)

  rnaseq = rnaseq[~rnaseq.ENSG.str.startswith('ENSGR')]

  ### Merge with samples:
  rnaseq = pandas.merge(rnaseq, samples, how="inner", left_on="SAMPID", right_index=True)

  rnaseq = CleanRnaseq(rnaseq, args.verbose)


  ### Compute median TPM by gene+tissue:

  rnaseq_med = rnaseq[['ENSG', 'SMTSD', 'TPM']].groupby(by=['ENSG','SMTSD'], as_index=False).median()
  print("Unique counts: genes: %d ; tissues: %d ; gene-tissue pairs: %d"%(rnaseq_med.ENSG.nunique(), rnaseq_med.SMTSD.nunique(), rnaseq_med.shape[0]), file=sys.stderr)

  if args.ofile_median:
    print("Output TPM median file: %s"%args.ofile_median, file=sys.stderr)
    rnaseq_med.round(3).to_csv(args.ofile_median, sep='\t', index=False)

  ### LOG10(TPM+1) useful transformation.

  rnaseq_med['LOG_TPM'] = rnaseq_med.TPM.apply(lambda x: numpy.log10(x+1))

  rnaseq_tau = rnaseq_med.groupby(['ENSG']).TPM.agg(TAU)
  rnaseq_tau = pandas.DataFrame(rnaseq_tau).rename(columns={'TPM':'TAU'})


  t0 = time.time()
  rnaseq_level = GTRanks(rnaseq_med.copy(), 'TPM')
  print("GTRanks elapsed: %ds"%(time.time()-t0), file=sys.stderr)


  if args.ofile_level:
    print("Output TPM level file: %s"%args.ofile_level, file=sys.stderr)
    rnaseq_level.round(3).to_csv(args.ofile_level, sep='\t', index=False)


  SABV(rnaseq_level, args.verbose)

