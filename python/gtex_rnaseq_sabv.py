#!/usr/bin/env python3
"""gtex_rnaseq_sabv.py

GTEx RNAseq Sex As Biological Variable (SABV) Analysis

Command-line version; see also Jupyter notebook gtex_rnaseq_sabv.ipynb

 - Author: Jeremy Yang
 - Based on R code by Oleg Ursu.
 - Required: Python3, Pandas 0.22+

Workflow (analysis):
 - READ: median TPMs, 1-row/gene+tissue+sex (from gtex_rnaseq_prep_app.py).
 - READ: sample TPMs, 1-row/gene+sample (from gtex_rnaseq_prep_app.py).
 - READ: gene IDs and names mapping file.

 - COMPUTE: TAU, TAU_F, TAU_M, tissue specificity index (Yanai et al., 2004)., by gene and gene+sex.
 - COMPUTE: TPM rank (quantile) among tissues, including by sex.
 - COMPUTE: Log fold-change, log of ratio (TPM_F/TPM_M).

 - COMPUTE: Wilcoxon rank sum test, F vs M, each gene+tissue.

 - OUTPUT: SABV metrics, TAU Log fold-change, 1-row/gene+tissue+sex.
 - OUTPUT: Wilcoxon rank Sum stats and pvals, 1-row/gene+tissue.

"""
#############################################################################
import sys,os,io,re,time,argparse,logging
import numpy,scipy,scipy.stats
import pandas

#############################################################################
def ReadMedianTPMs(ifile):
  fin = open(ifile)
  logging.info('=== GTEx Median-TPMs datafile: %s'%fin.name)
  tpms = pandas.read_csv(fin, sep='\t')
  logging.info("Median-TPMs dataset nrows: %d ; ncols: %d:"%(tpms.shape[0],tpms.shape[1]))
  return tpms

#############################################################################
def ReadSampleTPMs(ifile):
  fin = open(ifile)
  logging.info('=== GTEx Sample-TPMs datafile: %s'%fin.name)
  tpms = pandas.read_csv(fin, sep='\t')
  logging.info("Sample-TPMs dataset nrows: %d ; ncols: %d:"%(tpms.shape[0],tpms.shape[1]))
  return tpms

#############################################################################
### Format: one line per tissue name, in preferred order.
#############################################################################
def ReadTissues(ifile):
  fin = open(ifile)
  logging.info('=== GTEx Tissues datafile: %s'%fin.name)
  tissues = pandas.read_csv(fin, sep=';', index_col=False, header=None, names=['name'])
  tissues = tissues.name.str.strip()
  logging.info("n_tissues: %d:"%(tissues.size))
  logging.debug("tissues:\n%s"%(str(tissues)))
  return tissues

#############################################################################
def DescribeDf(df):
  buff = io.StringIO()
  #df.info(buf=buff, verbose=bool(verbose>0), null_counts=bool(verbose>0))
  df.info(buf=buff, verbose=True, null_counts=True)
  logging.info(re.sub(re.compile('^', re.M), '\t', buff.getvalue()))

#############################################################################
### Read gene symbols.
#############################################################################
def ReadGenes(ifile):
  logging.info("=== ReadGenes:")
  fin = open(ifile)
  logging.info('Biomart ENSG2NCBI genes datafile: %s'%fin.name)
  genes = pandas.read_csv(fin, sep='\t', usecols=[1,2,3], na_values=[''], dtype={2:str})
  logging.info("Genes dataset nrows: %d ; ncols: %d:"%(genes.shape[0],genes.shape[1]))
  genes.columns = ['ENSG','NCBI','HGNC']
  genes.dropna(inplace=True)
  return genes

#############################################################################
def TAUs(tpms):
  taus = tpms.groupby(['ENSG']).TPM.agg(TAU)
  taus = pandas.DataFrame(taus).rename(columns={'TPM':'TAU'})
  taus = taus.reset_index(drop=False)
  return taus

#############################################################################
def TAUs_SABV(tpms):
  taus_sex = tpms.groupby(['ENSG','SEX']).TPM.agg(TAU)
  taus_sex = pandas.DataFrame(taus_sex).rename(columns={'TPM':'TAU_BYSEX'})
  taus_sex = taus_sex.reset_index(drop=False)
  return taus_sex

#############################################################################
def RankByTissue(tpms):
  """
	Avoid SettingwithCopyWarning with syntax:
	df[(EXPRESSION), 'COLNAME'] = NEWVALS
  """
  tpms['TPM_RANK'] = pandas.Series(dtype=float)
  tpms['TPM_RANK_BYSEX'] = pandas.Series(dtype=float)
  for ensg in tpms.ENSG.unique():
    tpms_this = tpms.loc[(tpms.ENSG==ensg)]
    tpms.loc[(tpms.ENSG==ensg), 'TPM_RANK'] = tpms_this.TPM.rank(ascending=True, pct=True)
    for sex in tpms_this.SEX.unique():
      tpms_this_bysex = tpms.loc[(tpms.SEX==sex) & (tpms.ENSG==ensg)]
      tpms.loc[(tpms.SEX==sex) & (tpms.ENSG==ensg), 'TPM_RANK_BYSEX'] = tpms_this_bysex.TPM.rank(ascending=True, pct=True)
  tpms['TPM_LEVEL'] = tpms.TPM_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')
  tpms['TPM_LEVEL_BYSEX'] = tpms.TPM_RANK_BYSEX.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')
  return tpms

#############################################################################
def SABV_LogFoldChange(tpms):
  lfc = pandas.merge(
	tpms[tpms.SEX=='F'][['ENSG', 'SMTSD','TPM']].rename(columns={'TPM':'TPM_F'}),
	tpms[tpms.SEX=='M'][['ENSG', 'SMTSD','TPM']].rename(columns={'TPM':'TPM_M'}),
	on=['ENSG', 'SMTSD'], how='inner')
  lfc['log2foldchange'] = ((lfc.TPM_F+1) / (lfc.TPM_M+1)).apply(numpy.log2)
  #lfc['log2foldchange_abs'] = lfc.log2foldchange.apply(numpy.abs)
  return lfc

#############################################################################
### Compute tissue specificity index (Yanai et al., 2004).
### $ \tau = \frac{\sum_{i=0}^N (1 - x_i)}{N - 1} $
### * N = number of tissues
### * x = expression profile component normalized by the maximal component value
###
### Validate with example vector from paper.  Should be 0.95.
### logging.info('%.2f'%TAU([0,8,0,0,0,2,0,2,0,0,0,0]))
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
def WilcoxonRankSum(tpms):
  """May be very slow for large datasets."""
  wrs = tpms[['SMTSD', 'ENSG']].copy().drop_duplicates().sort_values(by=['SMTSD', 'ENSG'])
  wrs.reset_index(drop=True, inplace=True)
  wrs['stat'] = pandas.Series(dtype=float)
  wrs['pval'] = pandas.Series(dtype=float)
  for smtsd in wrs.SMTSD.unique():
    tpms_this = tpms[tpms.SMTSD==smtsd]
    for ensg in tpms_this.ENSG.unique():
      vals_f = tpms_this.TPM[(tpms_this.ENSG==ensg) & (tpms_this.SEX=="F")]
      vals_m = tpms_this.TPM[(tpms_this.ENSG==ensg) &(tpms_this.SEX=="M")]
      stat, pval = scipy.stats.ranksums(x=vals_f, y=vals_m)
      #stat, pval = scipy.stats.mannwhitneyu(x=vals_f.rank(), y=vals_m.rank(), use_continuity=True, alternative="two-sided")
      wrs.stat.loc[(wrs.SMTSD==smtsd) & (wrs.ENSG==ensg)] = stat
      wrs.pval.loc[(wrs.SMTSD==smtsd) & (wrs.ENSG==ensg)] = pval
  return wrs

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='GTEx RNAseq Exfiles/SABV preprocessor')
  parser.add_argument("--i", dest="ifile", help="input median TPM, 1-row/gene+tissue+sex (TSV)")
  parser.add_argument("--i_sample", dest="ifile_sample", help="input sample TPM, 1-row/gene+sample (TSV)")
  parser.add_argument("--i_tissue", dest="ifile_tissue", help="input (ordered) tissue file")
  parser.add_argument("--o", dest="ofile", help="output (TSV)")
  parser.add_argument("--decimals", type=int, default=3, help="output decimal places")
  parser.add_argument("-v", "--verbose", action="count")
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  t0 = time.time()

  if args.verbose:
    logging.info('Python: %s\nPandas: %s'%(sys.version, pandas.__version__))

  if args.ifile_tissue:
    tissues = ReadTissues(args.ifile_tissue)
  else:
    tissues = None

  if not args.ifile:
    parser.error('Input TPM file required.')

  tpms = ReadMedianTPMs(args.ifile)

  if args.verbose:
    #DescribeDf(tpms)
    logging.info("SABV TPM median unique counts: genes: %d"%(tpms.ENSG.nunique()))

  logging.info('=== SABV analysis:')
  logging.info('=== Compute TAU (tissue specificity):')
  taus = TAUs(tpms)
  logging.info('=== Compute TAU (tissue specificity)+SABV:')
  taus_sex = TAUs_SABV(tpms)
  tpms = pandas.merge(tpms, taus, on=['ENSG'], how='left')
  tpms = pandas.merge(tpms, taus_sex, on=['ENSG','SEX'], how='left')

  logging.info('=== Compute TPM rank (quantile) among tissues, including by sex:')
  tpms = RankByTissue(tpms)

  logging.info("=== Compute Log fold-change, log of ratio (F/M):")
  ### (One row per gene+tissue, cols for M and F TPM.)
  lfc = SABV_LogFoldChange(tpms)
  tpms = pandas.merge(tpms, lfc, on=['ENSG','SMTSD'], how='left')

  if args.ifile_sample:
    logging.info("=== From SAMPLE-LEVEL TPMs, compute Wilcoxon rank-sum statistic+pval:")
    sampletpms = ReadSampleTPMs(args.ifile_sample)
    wrs = WilcoxonRankSum(sampletpms)
    wrs = wrs.rename(columns={'stat':'WilcoxonRankSum_stat', 'pval':'WilcoxonRankSum_pval'})
    tpms = pandas.merge(tpms, wrs, on=['ENSG','SMTSD'], how='left')
  else:
    logging.info("NOTE: SAMPLE-LEVEL TPMs not provided, so SKIPPING Wilcoxon rank-sum test.")

  if args.ofile:
    logging.info("=== Output SABV file: %s"%args.ofile)
    tpms.round(args.decimals).to_csv(args.ofile, sep='\t', index=False)

  logging.info("Elapsed: %ds"%((time.time()-t0)))
