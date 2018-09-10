#!/usr/bin/env python3
"""gtex_rnaseq_prep_app.py

GTEx RNAseq Preprocessing, with Sex As Biological Variable (SABV)

Command-line version; see also Jupyter notebook gtex_rnaseq_prep.ipynb

 - Author: Jeremy Yang
 - Based on R code by Oleg Ursu.
 - Required: Python3, Pandas 0.22+

Workflow (prep):
 - READ: GTEx Subjects data, 1-row/subject.
 - READ: GTEx Samples data, 1-row/sample.
 - READ: GTEx RNAseq expression TPM data, 1-row/gene, 1-col/sample.
 - REMOVE: samples with Hardy score >2 (prefer healthier).
 - REMOVE: samples with high degree of autolysis (self-digestion).
 - MERGE: Samples and subjects, to 1-row/sample.
 - RESHAPE: RNAseq data from 1-col/sample, to 3 cols: gene, sample, TPM.
 - REMOVE: genes in pseudoautosomal regions (PAR) of chromosome Y.
 - AGGREGATE: samples, computing median TPM by gene+tissue.
 - AGGREGATE: samples, computing median TPM by gene+tissue+sex.
 - OUTPUT: median TPMs, 1-row/gene+tissue+sex.
 - OUTPUT: expression profiles, 1-row/gene+sex.

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
  subjects = pandas.read_csv(fin, sep='\t')
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
  df.info(buf=buff,verbose=bool(verbose),null_counts=bool(verbose))
  print(re.sub(re.compile('^', re.M), '\t', buff.getvalue()), file=sys.stdout)

#############################################################################
### (GTEx_v7_Annotations_SampleAttributesDS.txt)
#############################################################################
def ReadSamples(ifile, verbose):
  print("=== ReadSamples:", file=sys.stdout)
  fin = open(ifile)
  print('GTEx Samples datafile: %s'%fin.name, file=sys.stdout)
  samples = pandas.read_csv(fin, sep='\t')
  samples = samples[['SAMPID', 'SMATSSCR', 'SMTS', 'SMTSD']]
  print("Samples dataset nrows: %d ; ncols: %d:"%(samples.shape[0],samples.shape[1]), file=sys.stdout)
  ### SUBJID is first two hyphen-delimted fields of SAMPID.
  samples['SUBJID'] = samples.SAMPID.str.extract('^([^-]+-[^-]+)-', expand=True)
  DescribeDf(samples, verbose)
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
    print('\t%d. "%s": %4d'%(i,name,val), file=sys.stdout)

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
  samples = rnaseq.columns[1:]
  print("RNAseq samples count: %d:"%(samples.size), file=sys.stdout)
  print("RNAseq unique samples count: %d:"%(samples.nunique()), file=sys.stdout)
  print("RNAseq genes (ENSG) count: %d:"%(rnaseq.ENSG.size), file=sys.stdout)
  print("RNAseq unique genes (ENSG) count: %d:"%(rnaseq.ENSG.nunique()), file=sys.stdout)
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
### Memory intensive; may need to divide task.
#############################################################################
def CleanRnaseq(rnaseq, verbose):
  print("DEBUG: CleanRnaseq IN: nrows = %d, cols: %s"%(rnaseq.shape[0],str(rnaseq.columns.tolist())), file=sys.stderr)
  print("=== CleanRnaseq:", file=sys.stdout)

  print("Removing data for gene-tissue pairs with all zero expression...", file=sys.stderr)
  maxtpm_0  = (rnaseq[['ENSG','SMTSD','TPM']].groupby(by=['ENSG','SMTSD'], as_index=True).max()==0).rename(columns={'TPM':'maxtpm_0'})
  print("DEBUG: maxtpm_0 value_counts: %s"str(maxtpm_0.maxtpm_0.value_counts()), file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, maxtpm_0, left_on=['ENSG','SMTSD'], right_index=True)
  del(maxtpm_0)
  rnaseq = rnaseq[~rnaseq['maxtpm_0']]
  rnaseq.drop(columns=['maxtpm_0'], inplace=True)

  print("Removing data for gene-tissue pairs not present in both sexes...", file=sys.stderr)
  sex_count = (rnaseq[['ENSG','SMTSD','SEX']].groupby(by=['ENSG','SMTSD'], as_index=True).nunique()).rename(columns={'SEX':'sex_count'})
  sex_count = sex_count[['sex_count']] #Why needed?
  print("DEBUG: sex_count value_counts: %s"str(sex_count.sex_count.value_counts()), file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, sex_count, left_on=['ENSG','SMTSD'], right_index=True)
  del(sex_count)
  rnaseq = rnaseq[rnaseq['sex_count']==2]
  rnaseq.drop(columns=['sex_count'], inplace=True)
  rnaseq = rnaseq.reset_index(drop=True)
  ### This removes most sex-specific tissues, but not breast.
  rnaseq = rnaseq[~rnaseq.SMTSD.str.match('^Breast')]

  rnaseq = rnaseq[['ENSG','SMTSD','SAMPID','SMATSSCR','SEX','AGE','DTHHRDY','TPM']]
  rnaseq = rnaseq.sort_values(by=['ENSG','SMTSD','SAMPID'])
  print("RNAseq unique samples count: %d:"%(rnaseq.SAMPID.nunique()), file=sys.stdout)
  print("RNAseq unique tissues count: %d:"%(rnaseq.SMTSD.nunique()), file=sys.stdout)
  print("RNAseq unique gene count: %d"%(rnaseq.ENSG.nunique()), file=sys.stdout)
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
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='GTEx RNAseq Exfiles/SABV preprocessor')
  parser.add_argument("--i_subject",dest="ifile_subject",help="input subjects file")
  parser.add_argument("--i_sample",dest="ifile_sample",help="input samples file")
  parser.add_argument("--i_rnaseq",dest="ifile_rnaseq",help="input rnaseq file")
  #parser.add_argument("--i_gene",dest="ifile_gene",help="input gene file")
  parser.add_argument("--i_tissue",dest="ifile_tissue",help="input (ordered) tissue file")
  parser.add_argument("--o_median",dest="ofile_median",help="output median TPM, 1-row/gene+tissue+sex (TSV)")
  parser.add_argument("--o_sample",dest="ofile_sample",help="output sample TPM, 1-row/gene+sample (TSV)")
  parser.add_argument("--o_profiles",dest="ofile_profiles",help="output profiles, 1-row/gene+sex (TSV)")
  parser.add_argument("--o_tissue",dest="ofile_tissue",help="output tissues (TSV)")
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
  samples = pandas.merge(samples, subjects, how='inner', on='SUBJID')

  if args.verbose:
    DescribeSamples(samples)

  if args.ofile_tissue:
    sample_tissues = samples[['SMTS','SMTSD']].reset_index(drop=True)
    sample_tissues = sample_tissues.drop_duplicates().sort_values(['SMTS', 'SMTSD'])
    print("=== Output tissues file: %s"%args.ofile_tissue, file=sys.stdout)
    sample_tissues.round(args.decimals).to_csv(args.ofile_tissue, sep='\t', index=False)

  samples = CleanSamples(samples, args.verbose)

#  if not args.ifile_gene:
#    parser.error('Input gene file required.')
#  genes = ReadGenes(args.ifile_gene, args.verbose)

  if not args.ifile_rnaseq:
    parser.error('Input RNAseq file required.')
  t1 = time.time()
  rnaseq = ReadRnaseq(args.ifile_rnaseq, args.verbose)
  print("ReadRnaseq elapsed: %ds"%(time.time()-t1), file=sys.stdout)

  print('=== MELT: One row per ENSG+SAMPID+TPM triplet:', file=sys.stdout)
  ### Easier to handle but ~3x storage.
  rnaseq = rnaseq.melt(id_vars = "ENSG", var_name = "SAMPID", value_name = "TPM")
  DescribeDf(rnaseq,args.verbose)
  print("RNAseq unique gene count (after melt): %d"%(rnaseq.ENSG.nunique()), file=sys.stdout)

  # Note could lose genes via inner join:
#  rnaseq = pandas.merge(genes, rnaseq, on='ENSG', how='inner')
#  print("RNAseq unique gene count (after merge with gene names): %d"%(rnaseq.ENSG.nunique()), file=sys.stdout)

  print('=== Remove genes in pseudoautosomal regions (PAR) of chromosome Y ("ENSGR"):', file=sys.stdout)
  n_ensgr = rnaseq.ENSG.str.startswith('ENSGR').sum()
  print('ENSGR gene TPMs: %d (%.2f%%)'%(n_ensgr,100*n_ensgr/rnaseq.shape[0]), file=sys.stdout)
  rnaseq = rnaseq[~rnaseq.ENSG.str.startswith('ENSGR')]
  print("RNAseq unique gene count (after PAR removal): %d"%(rnaseq.ENSG.nunique()), file=sys.stdout)

  print('=== Merge with samples:', file=sys.stdout)
  rnaseq = pandas.merge(rnaseq, samples, how="inner", on="SAMPID")
  print("RNAseq unique gene count (after merge with samples): %d"%(rnaseq.ENSG.nunique()), file=sys.stdout)

  rnaseq = CleanRnaseq(rnaseq, args.verbose)

  if args.ofile_sample:
    print("=== Output sample TPM file: %s"%args.ofile_sample, file=sys.stdout)
    rnaseq.round(args.decimals).to_csv(args.ofile_sample, sep='\t', index=False)

  print('=== Compute median TPM by gene+tissue+sex:', file=sys.stdout)
  rnaseq = SABV_aggregate_median(rnaseq, args.verbose)

  print("SABV TPM median unique counts: genes: %d"%(rnaseq.ENSG.nunique()), file=sys.stdout)

  if args.ofile_median:
    print("=== Output median (by gene+tissue+sex) TPM file: %s"%args.ofile_median, file=sys.stdout)
    rnaseq.round(args.decimals).to_csv(args.ofile_median, sep='\t', index=False)

  print("=== Pivot to one-row-per-gene format (profiles).", file=sys.stdout)
  rnaseq_profiles = PivotToProfiles(rnaseq, tissues, args.verbose)
  if args.ofile_profiles:
    print("=== Output profiles file: %s"%args.ofile_profiles, file=sys.stdout)
    rnaseq_profiles.round(args.decimals).to_csv(args.ofile_profiles, sep='\t', index=False)

  print("%s Elapsed: %ds"%(PROG,(time.time()-t0)), file=sys.stderr)
