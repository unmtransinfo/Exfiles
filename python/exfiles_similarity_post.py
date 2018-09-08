#!/usr/bin/env python3
"""exfiles_similarity_post.py

Post-process output of exfiles_similarity.py.

Input cols:
	ENSGA,SEXA,ENSGB,SEXB,SpearmanRho,SpearmanP

Output cols:
	ENSGA,ENSGB,Cluster,Rho

where
	Cluster = F (F-F)
	Cluster = M (M-M)
	Cluster = FM (mean)

"""
#############################################################################
import sys,os,io,re,time,argparse
import numpy,scipy,scipy.stats
import pandas



#############################################################################
def ReadCorrfile(ifile, verbose):
  fin = open(ifile)
  print('=== Correlations datafile: %s'%fin.name, file=sys.stdout)
  corrs = pandas.read_csv(fin, sep='\t')
  print("Correlations dataset nrows: %d ; ncols: %d:"%(corrs.shape[0],corrs.shape[1]), file=sys.stdout)
  print("Correlations cols: %s:"%(str(corrs.columns.tolist())), file=sys.stdout)
  return corrs

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
def GroupCorrs(corrs, genes, verbose):
  corrs = corrs[['ENSGA', 'SEXA', 'ENSGB', 'SEXB', 'SpearmanRho']]
  corrs = pandas.merge(corrs,
	genes.rename(columns={'ENSG':'ENSGA','HGNC':'Ga'}),
	on=['ENSGA'], how='inner')
  corrs = pandas.merge(corrs,
	genes.rename(columns={'ENSG':'ENSGB','HGNC':'Gb'}),
	on=['ENSGB'], how='inner')
  corrs = corrs[['Ga', 'Gb', 'SEXA', 'SEXB', 'SpearmanRho']]

  corrs = corrs.drop_duplicates()

  corrs_grouped_f = corrs[
	(corrs.SEXA=='female') & (corrs.SEXB=='female')
	].drop(columns=['SEXA','SEXB'])
  corrs_grouped_f['Cluster'] = 'F'
  corrs_grouped_m = corrs[
	(corrs.SEXA=='male') & (corrs.SEXB=='male')
	].drop(columns=['SEXA','SEXB'])
  corrs_grouped_m['Cluster'] = 'M'

  corrs_grouped = pandas.concat([corrs_grouped_f,corrs_grouped_m])

  ## Now handle MF; aggregate on Ga, Gb in alpha order.
  corrs_grouped_mf = corrs[
	((corrs.SEXA=='female') & (corrs.SEXB=='male'))
	| ((corrs.SEXA=='male') & (corrs.SEXB=='female'))
	].drop(columns=['SEXA','SEXB'])
  corrs_grouped_mf['Gmin'] = corrs_grouped_mf[['Ga','Gb']].min(axis=1)
  corrs_grouped_mf['Gmax'] = corrs_grouped_mf[['Ga','Gb']].max(axis=1)
  corrs_grouped_mf['Ga'] = corrs_grouped_mf['Gmin']
  corrs_grouped_mf['Gb'] = corrs_grouped_mf['Gmax']
  corrs_grouped_mf = corrs_grouped_mf.drop(columns=['Gmin','Gmax'])
  corrs_grouped_mf = corrs_grouped_mf.groupby(by=['Ga','Gb'], as_index=False).mean()
  corrs_grouped_mf['Cluster'] = 'MF'
  corrs_grouped = pandas.concat([corrs_grouped,corrs_grouped_mf])
  corrs_grouped = corrs_grouped[['Ga','Gb','Cluster','SpearmanRho']]
  return corrs_grouped


#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Exfiles similarity')
  parser.add_argument("--i",dest="ifile",help="input gene-gene correlations (TSV)")
  parser.add_argument("--i_gene",dest="ifile_gene",help="input gene IDs (TSV)")
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("--decimals",type=int,default=3,help="output decimal places")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  PROG=os.path.basename(sys.argv[0])
  t0 = time.time()

  if args.verbose:
    print('Python: %s; pandas: %s; numpy: %s; scipy: %s'%(sys.version.split()[0],pandas.__version__, numpy.__version__, scipy.__version__), file=sys.stdout)

  if not args.ifile:
    parser.error('Input file required.')

  if not args.ifile_gene:
    parser.error('Input gene file required.')
  genes = ReadGenes(args.ifile_gene, args.verbose)


  corrs = ReadCorrfile(args.ifile, args.verbose)

  corrs_grouped = GroupCorrs(corrs, genes, args.verbose)

  if args.ofile:
    print("=== Output file: %s"%args.ofile, file=sys.stdout)
    corrs_grouped.round(args.decimals).to_csv(args.ofile, sep='\t', index=False)

  print("%s Elapsed: %ds"%(PROG,(time.time()-t0)), file=sys.stderr)
