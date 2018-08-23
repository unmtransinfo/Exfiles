#!/usr/bin/env python3
"""gtex_rnaseq_prep_app.py

GTEx RNAseq Preprocessing, with Sex As Biological Variable (SABV)

Command-line version; see also Jupyter notebook gtex_rnaseq_prep.ipynb

 - Author: Jeremy Yang
 - Based on R code by Oleg Ursu.
 - Required: Python3, Pandas 0.22+

Workflow:
 - READ: GTEx Subjects data, 1-row/subject.
 - READ: GTEx Samples data, 1-row/sample.
 - READ: GTEx RNAseq expression TPM data, 1-row/gene, 1-col/sample.
 - READ: Ensembl gene IDs file ENSG2NCBI, from Ensembl-Biomart (Gene stable ID version, NCBI gene ID, HGNC symbol).
 - REMOVE: samples with Hardy score >2 (prefer healthier).
 - REMOVE: samples with high degree of autolysis (self-digestion).
 - MERGE: Samples and subjects, to 1-row/sample.
 - RESHAPE: RNAseq data from 1-col/sample, to 3 cols: gene, sample, TPM.
 - REMOVE: genes in pseudoautosomal regions (PAR) of chromosome Y.
 - AGGREGATE: samples, computing median TPM by gene+tissue.
 - AGGREGATE: samples, computing median TPM by gene+tissue+sex.

 - COMPUTE: ranks, each gene+tissue, among tissues.
 - COMPUTE: Wilcoxon signed rank test, F vs M, each gene+tissue.
 - COMPUTE: Log fold-change, log of ratio (TPM_F/TPM_M).
 - COMPUTE: TAU, TAU_F, TAU_M, tissue specificity index (Yanai et al., 2004)., by gene and gene+sex.

 - OUTPUT: expression profiles, 1-row/gene+sex.
 - OUTPUT: TAU, TAU_F, TAU_M, 1-row/gene.
 - OUTPUT: SABV metrics, Wilcoxon signed rank and Log fold-change, 1-row/gene+tissue.

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
### Format: one line per tissue name, in preferred order.
#############################################################################
def ReadTissues(ifile, verbose):
  fin = open(ifile)
  print('=== GTEx Tissues datafile: %s'%fin.name, file=sys.stdout)
  tissues = pandas.read_csv(fin, sep=';', index_col=False, header=None, names=['name'])
  tissues = tissues.name.str.strip()
  if verbose: print("n_tissues: %d:"%(tissues.size), file=sys.stdout)
  if verbose: print("DEBUG: tissues:\n%s"%(str(tissues)), file=sys.stderr)
  return tissues

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
### Clean & tidy cols. Remove samples with high degree of autolysis (self-digestion).
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
  i=0
  for name,val in samples.SMTSD.value_counts().sort_index().iteritems():
    i+=1
    print('\t%d. SMTSD "%s": %4d'%(i,name,val), file=sys.stdout)

#############################################################################
### READ GENE TPMs (full or demo subset)
### Top 2 rows, format:
###	#1.2
###	nrow	ncol
### Full file is ~56k rows, 2.6GB uncompressed.  Demo ~1k rows.
### *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
### *   GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz
#############################################################################
def ReadRnaseq(ifile, verbose):
  print("=== ReadRnaseq:", file=sys.stdout)
  fin = open(ifile, "rb")
  print('GTEx RNAseq TPM datafile: %s'%fin.name, file=sys.stdout)
  rnaseq = pandas.read_table(fin, compression='gzip', sep='\t', skiprows=2)
  print("RNAseq dataset nrows: %d ; ncols: %d:"%(rnaseq.shape[0],rnaseq.shape[1]), file=sys.stdout)
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
  print("DEBUG: CleanRnaseq IN: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  print("=== CleanRnaseq:", file=sys.stdout)

  maxtpm_0  = (rnaseq[['ENSG','SMTSD','TPM']].groupby(by=['ENSG','SMTSD'], as_index=True).max()==0).rename(columns={'TPM':'maxtpm_0'})
  #print(maxtpm_0.maxtpm_0.value_counts(), file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, maxtpm_0, left_on=['ENSG', 'SMTSD'], right_index=True)
  rnaseq = rnaseq[~rnaseq['maxtpm_0']]
  rnaseq.drop(columns=['maxtpm_0'], inplace=True)

  sex_count = (rnaseq[['ENSG','SMTSD','SEX']].groupby(by=['ENSG','SMTSD'], as_index=True).nunique()).rename(columns={'SEX':'sex_count'})
  sex_count = sex_count[['sex_count']] #Why needed?
  #print(sex_count.sex_count.value_counts(), file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, sex_count, left_on=['ENSG', 'SMTSD'], right_index=True)

  rnaseq = rnaseq[rnaseq['sex_count'] == 2]
  rnaseq.drop(columns=['sex_count'], inplace=True)
  rnaseq = rnaseq.reset_index(drop=True)

  print("DEBUG: CleanRnaseq OUT: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  return rnaseq

#############################################################################
def SABV_aggregate_median(rnaseq, verbose):
  print("=== SABV_aggregate_median:", file=sys.stdout)
  print("DEBUG: SABV_aggregate_median IN: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)

  ### Compute median TPM by gene+tissue+sex:
  rnaseq = rnaseq[['ENSG', 'SMTSD', 'SEX', 'TPM']].groupby(by=['ENSG','SMTSD','SEX'], as_index=False).median()

  print("DEBUG: SABV_aggregate_median OUT: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  return rnaseq

#############################################################################
def SABV_TAU(rnaseq, verbose):
  print("=== SABV_TAU:", file=sys.stdout)
  print("DEBUG: SABV_TAU IN: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)

  ### Compute TAU by gene:
  rnaseq_nosex_tau = rnaseq.groupby(['ENSG']).TPM.agg(TAU)
  rnaseq_nosex_tau = pandas.DataFrame(rnaseq_nosex_tau).rename(columns={'TPM':'TAU'})
  rnaseq_nosex_tau['SEX'] = 'ALL'


  ### Compute TAU by gene+sex:
  rnaseq_tau_f = rnaseq.loc[rnaseq['SEX']=='female'].groupby(['ENSG']).TPM.agg(TAU)
  rnaseq_tau_f = pandas.DataFrame(rnaseq_tau_f).rename(columns={'TPM':'TAU'})
  rnaseq_tau_f = rnaseq_tau_f.reset_index(drop=False)
  rnaseq_tau_f['SEX'] = 'female'
  print("DEBUG: rnaseq_tau_f.columns = %s"%(str(rnaseq_tau_f.columns)), file=sys.stderr)

  rnaseq_tau_m = rnaseq.loc[rnaseq['SEX']=='male'].groupby(['ENSG']).TPM.agg(TAU)
  rnaseq_tau_m = pandas.DataFrame(rnaseq_tau_m).rename(columns={'TPM':'TAU'})
  rnaseq_tau_m = rnaseq_tau_m.reset_index(drop=False)
  rnaseq_tau_m['SEX'] = 'male'
  print("DEBUG: rnaseq_tau_m.columns = %s"%(str(rnaseq_tau_m.columns)), file=sys.stderr)

  rnaseq_tau = pandas.concat([rnaseq_nosex_tau, rnaseq_tau_f, rnaseq_tau_m])

  print("DEBUG: SABV_TAU OUT: nrows = %d, cols: %s"%(rnaseq_tau.shape[0],str(rnaseq_tau.columns.tolist())), file=sys.stderr)
  return rnaseq_tau


#############################################################################
def SABV_GTRanks(rnaseq, verbose):
  print("=== SABV_GTRanks:", file=sys.stdout)
  print("DEBUG: SABV_GTRanks IN: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  print("=== Assign gene-tissue rank (quantile) among tissues (F):", file=sys.stdout)
  t0 = time.time()
  rnaseq_ranks_f = GTRanks(rnaseq[['ENSG','SMTSD','TPM_F']].copy(), 'TPM_F')
  print("GTRanks (F) elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  rnaseq_ranks_f['LEVEL_F'] = rnaseq_ranks_f.TPM_F_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')
  print("=== Assign gene-tissue rank (quantile) among tissues (M):", file=sys.stdout)
  t0 = time.time()
  rnaseq_ranks_m = GTRanks(rnaseq[['ENSG','SMTSD','TPM_M']].copy(), 'TPM_M')
  print("GTRanks (M) Elapsed: %ds"%(time.time()-t0), file=sys.stderr)
  rnaseq_ranks_m['LEVEL_M'] = rnaseq_ranks_m.TPM_M_RANK.apply(lambda x: 'Not detected' if x==0 else 'Low' if x<.25 else 'Medium' if x<.75 else 'High')
  rnaseq = pandas.merge(rnaseq_ranks_f, rnaseq_ranks_m, on=['ENSG','SMTSD'], how='inner')
  print("DEBUG: SABV_GTRanks OUT: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  return rnaseq

#############################################################################
def WilcoxonSignedRank(rnaseq, rnaseq_ranks, verbose):
  print("=== WilcoxonSignedRank:", file=sys.stdout)
  print("DEBUG: WilcoxonSignedRank IN: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  ### For each gene, compute sex difference via Wilcoxon signed-rank test,
  ### with Wilcox treatment, discarding all zero-differences.

  results = pandas.DataFrame({'ENSG':rnaseq.ENSG.drop_duplicates().sort_values(), 'WilcoxonSignedRank_stat':None, 'WilcoxonSignedRank_pval':None}).reset_index(drop=True)

  for i in range(results.shape[0]):
    tpm_f_this = rnaseq_ranks.TPM_F_RANK[rnaseq_ranks.ENSG==results.ENSG[i]]
    tpm_m_this = rnaseq_ranks.TPM_M_RANK[rnaseq_ranks.ENSG==results.ENSG[i]]

    ### What is best minimum size??
    if tpm_f_this[tpm_f_this>0].size<8 or tpm_m_this[tpm_m_this>0].size<8:
      continue

    try:
      stat, pval = scipy.stats.wilcoxon(x=tpm_f_this, y=tpm_m_this, zero_method='wilcox')
    except Exception, e:
      print("Exception [i=%d; ensg=%s]: %s"%(i+1,results.ENSG[i],str(e)), file=sys.stderr)
      print("DEBUG: tpm_f_this=%s; tpm_m_this=%s"%(str(tpm_f_this),str(tpm_m_this)), file=sys.stderr)
      continue
      
    results.WilcoxonSignedRank_stat.iloc[i] = stat
    results.WilcoxonSignedRank_pval.iloc[i] = pval 

  rnaseq = pandas.merge(rnaseq, results, on=['ENSG'])
  print("DEBUG: WilcoxonSignedRank OUT: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)

  return rnaseq

#############################################################################
### Reshape to one-row-per-gene format.
### From:   ENSG,SMTSD,SEX,TPM,LOG_TPM
### To:	    ENSG,SEX,TPM_1,TPM_2,...TPM_N (N tissues)
### Preserve tissue order.
#############################################################################
def PivotToProfiles(rnaseq, tissues_ordered, verbose):
  print("DEBUG: PivotToProfiles IN: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  tissues = pandas.Series(pandas.unique(rnaseq.SMTSD.sort_values()))
  if type(tissues_ordered)==pandas.core.series.Series:
    if set(tissues) == set(tissues_ordered):
      tissues = tissues_ordered
      print("Note: input tissues (ordered): %s"%(str(set(tissues))), file=sys.stderr)
    else:
      print("Warning: input tissues missing in samples: %s"%(str(set(tissues_ordered) - set(tissues))), file=sys.stderr)
      print("Warning: sample tissues missing in input: %s"%(str(set(tissues) - set(tissues_ordered))), file=sys.stderr)

  # Assure only 1-row per unique (ensg,smtsd) tuple (or pivot will fail).
  #rnaseq = rnaseq.drop_duplicates(subset=['ENSG','SMTSD'], keep='first')

  rnaseq_f = rnaseq[rnaseq.SEX=='female'].drop(columns=['SEX'])
  rnaseq_m = rnaseq[rnaseq.SEX=='male'].drop(columns=['SEX'])

  rnaseq_f = rnaseq_f[['ENSG','SMTSD','TPM']]
  rnaseq_m = rnaseq_m[['ENSG','SMTSD','TPM']]

  exfiles_f = rnaseq_f.pivot(index='ENSG', columns='SMTSD')
  exfiles_f.columns = exfiles_f.columns.get_level_values(1)
  exfiles_f = exfiles_f.reset_index(drop=False)
  exfiles_f['SEX'] = 'female'
  exfiles_m = rnaseq_m.pivot(index='ENSG', columns='SMTSD')
  exfiles_m.columns = exfiles_m.columns.get_level_values(1)
  exfiles_m = exfiles_m.reset_index(drop=False)
  exfiles_m['SEX'] = 'male'
  exfiles = pandas.concat([exfiles_f, exfiles_m])
  cols = ['ENSG','SEX']+tissues.tolist()
  exfiles = exfiles[cols]
  DescribeDf(exfiles,verbose)
  print("DEBUG: PivotToProfiles OUT: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())), file=sys.stderr)
  return exfiles

#############################################################################
def SABV_LogFoldChange(rnaseq, verbose):
  ### Combine rows into one row per gene+tissue, cols TPM_F, TPM_M.
  rnaseq_f = rnaseq.loc[rnaseq['SEX']=='female']
  rnaseq_f = rnaseq_f[['ENSG', 'SMTSD', 'TPM']].rename(columns={'TPM':'TPM_F'})
  rnaseq_m = rnaseq.loc[rnaseq['SEX']=='male']
  rnaseq_m = rnaseq_m[['ENSG','SMTSD','TPM']].rename(columns={'TPM':'TPM_M'})
  rnaseq = pandas.merge(rnaseq_f, rnaseq_m, how='inner', on=['ENSG','SMTSD'])
  rnaseq['Log2FoldChange'] = ((rnaseq.TPM_F+1) / (rnaseq.TPM_M+1)).apply(lambda x: numpy.log2(max(x, 1/x)))
  return rnaseq

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
  parser.add_argument("--i_tissue",dest="ifile_tissue",help="input (ordered) tissue file")
  #parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("--o_sabv",dest="ofile_sabv",help="output SABV (TSV)")
  parser.add_argument("--o_tau",dest="ofile_tau",help="output TAU (TSV)")
  parser.add_argument("--o_tissue",dest="ofile_tissue",help="output tissues (TSV)")
  parser.add_argument("--o_profiles",dest="ofile_profiles",help="output profiles, 1-row/gene (TSV)")
  parser.add_argument("--decimals",type=int,default=3,help="output decimal places")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  PROG=os.path.basename(sys.argv[0])
  t0 = time.time()

  if args.verbose:
    print('Python: %s\nPandas: %s'%(sys.version,pandas.__version__), file=sys.stdout)

  if args.ifile_tissue:
    tissues = ReadTissues(args.ifile_tissue, args.verbose)
  else:
    tissues = None

  if not args.ifile_subject:
    parser.error('Input subject file required.')

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

  if args.ofile_tissue:
    sample_tissues = samples[['SMTS','SMTSD']].reset_index(drop=True)
    sample_tissues = sample_tissues.drop_duplicates().sort_values(['SMTS', 'SMTSD'])
    print("=== Output tissues file: %s"%args.ofile_tissue, file=sys.stdout)
    sample_tissues.round(args.decimals).to_csv(args.ofile_tissue, sep='\t', index=False)

  samples = CleanSamples(samples, args.verbose)

  if not args.ifile_gene:
    parser.error('Input gene file required.')
  genes = ReadGenes(args.ifile_gene, args.verbose)

  if not args.ifile_rnaseq:
    parser.error('Input RNAseq file required.')
  t1 = time.time()
  rnaseq = ReadRnaseq(args.ifile_rnaseq, args.verbose)
  print("ReadRnaseq elapsed: %ds"%(time.time()-t1), file=sys.stdout)

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

  print('=== Non-SABV analysis:', file=sys.stdout)
  print('=== Compute median TPM by gene+tissue:', file=sys.stdout)
  rnaseq_nosex = rnaseq[['ENSG', 'SMTSD', 'TPM']].groupby(by=['ENSG','SMTSD'], as_index=False).median()
  rnaseq_nosex['SEX'] = 'ALL'
  print("RNAseq unique counts: genes: %d ; tissues: %d ; gene-tissue pairs: %d"%(rnaseq_nosex.ENSG.nunique(), rnaseq_nosex.SMTSD.nunique(), rnaseq_nosex.shape[0]), file=sys.stdout)

  ### Needed?
  #print('=== Compute LOG10(TPM+1):', file=sys.stdout)
  #rnaseq_nosex['LOG_TPM'] = rnaseq_nosex.TPM.apply(lambda x: numpy.log10(x+1))
  ### 

  print('=== Compute gene-tissue ranks (GTRanks):', file=sys.stdout)
  rnaseq_nosex = GTRanks(rnaseq_nosex.copy(), 'TPM')
  print("TPM level unique counts: genes: %d"%(rnaseq_nosex.ENSG.nunique()), file=sys.stdout)

  ### Needed?
  #if args.ofile:
  #  print("=== Output file (non-SABV): %s"%args.ofile, file=sys.stdout)
  #  rnaseq_nosex.round(args.decimals).to_csv(args.ofile, sep='\t', index=False)
  ###

  print('=== SABV analysis:', file=sys.stdout)
  print('=== Compute median TPM by gene+tissue+sex:', file=sys.stdout)
  rnaseq = SABV_aggregate_median(rnaseq, args.verbose)

  print("SABV TPM median unique counts: genes: %d"%(rnaseq.ENSG.nunique()), file=sys.stdout)

  print("=== Pivot to one-row-per-gene format: %s"%args.ofile, file=sys.stdout)
  rnaseq_profiles = PivotToProfiles(rnaseq, tissues, args.verbose)
  if args.ofile_profiles:
    print("=== Output profiles file: %s"%args.ofile_profiles, file=sys.stdout)
    rnaseq_profiles.round(args.decimals).to_csv(args.ofile_profiles, sep='\t', index=False)

  print('=== Compute TAU (tissue specificity, Yanai et al., 2004):', file=sys.stdout)
  rnaseq_tau = SABV_TAU(rnaseq, args.verbose)
  if args.ofile_tau:
    print("=== Output TAU file: %s"%args.ofile_tau, file=sys.stdout)
    rnaseq_tau.round(args.decimals).to_csv(args.ofile_tau, sep='\t') #Index is ENSG

  print("=== Compute Log fold-change, log of ratio (F/M):", file=sys.stdout)
  ### (Also combine rows into one row per gene+tissue, cols for M and F TPM.)
  rnaseq = SABV_LogFoldChange(rnaseq, args.verbose)


  ### Ranks needed for WilcoxonSignedRank test.
  rnaseq_ranks = SABV_GTRanks(rnaseq, args.verbose)
  rnaseq = WilcoxonSignedRank(rnaseq, rnaseq_ranks, args.verbose)

  if args.ofile_sabv:
    print("=== Output SABV file: %s"%args.ofile_sabv, file=sys.stdout)
    rnaseq.round(args.decimals).to_csv(args.ofile_sabv, sep='\t', index=False)

  print("%s Elapsed: %ds"%(PROG,(time.time()-t0)), file=sys.stderr)
