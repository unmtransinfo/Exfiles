#!/usr/bin/env python3
"""exfiles_similarity_post.py

Post-process output of (1) exfiles_similarity.py, and (2) exfiles_similarity_wcorr.R

Input cols:
	ENSGA,SEXA,ENSGB,SEXB,Ruzicka
	ENSGA,SEXA,ENSGB,SEXB,wRho

Output cols:
	ENSGA,ENSGB,Cluster,wRho,Ruzicka

where
	Cluster = F (F-F)
	Cluster = M (M-M)
	Cluster = FM (mean)

"""
#############################################################################
import sys,os,io,re,time,argparse
import pandas,numpy,scipy,scipy.stats

#############################################################################
def ReadCorrfile(ifile, verbose):
  fin = open(ifile)
  print('=== Correlation datafile: %s'%fin.name, file=sys.stdout)
  cors = pandas.read_csv(fin, sep='\t', na_values=['','NA','NaN'])
  cors.dropna(how='any', inplace=True)
  print("Correlation dataset nrows: %d ; ncols: %d:"%(cors.shape[0],cors.shape[1]), file=sys.stdout)
  print("Correlation cols: %s:"%(str(cors.columns.tolist())), file=sys.stdout)
  return cors

#############################################################################
def ReadSimfile(ifile, verbose):
  fin = open(ifile)
  print('=== Similarity datafile: %s'%fin.name, file=sys.stdout)
  sims = pandas.read_csv(fin, sep='\t', na_values=['','NA','NaN'])
  sims.dropna(how='any', inplace=True)
  print("Similarity dataset nrows: %d ; ncols: %d:"%(sims.shape[0],sims.shape[1]), file=sys.stdout)
  print("Similarity cols: %s:"%(str(sims.columns.tolist())), file=sys.stdout)
  return sims

#############################################################################
### Read gene IDs, etc.: ENSG,NCBI,HGNCID,symbol,name
#############################################################################
def ReadGenes(ifile, verbose):
  LOG("=== ReadGenes:")
  fin = open(ifile)
  LOG('GTEx/Ensembl/HGNC genes datafile: %s'%fin.name)
  genes = pandas.read_csv(fin, sep='\t', na_values=[''], dtype={2:str})
  LOG("Genes dataset nrows: %d ; ncols: %d:"%(genes.shape[0],genes.shape[1]))
  return genes

#############################################################################
def GroupComparisons(cors, sims, genes, verbose):
  cors = cors[['ENSGA','SEXA','ENSGB','SEXB','wRho']]
  cors = pandas.merge(cors,genes.rename(columns={'ENSG':'ENSGA','symbol':'Ga'}), on=['ENSGA'], how='inner')
  cors = pandas.merge(cors, genes.rename(columns={'ENSG':'ENSGB','symbol':'Gb'}), on=['ENSGB'], how='inner')
  cors = cors[['Ga','Gb','SEXA','SEXB','wRho']]
  cors = cors.drop_duplicates()
  #LOG("DEBUG: cors nrows: %d ; ncols: %d:"%(cors.shape[0],cors.shape[1]))
  cors_ff = cors[(cors.SEXA=='F')&(cors.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  cors_ff['Cluster']='F'
  cors_mm = cors[(cors.SEXA=='M')&(cors.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  cors_mm['Cluster']='M'
  ## MF: aggregate on Ga, Gb in alpha order.
  cors_mf = cors[((cors.SEXA=='F')&(cors.SEXB=='M'))|((cors.SEXA=='M')&(cors.SEXB=='F'))].drop(columns=['SEXA','SEXB'])
  cors_mf['Gmin'] = cors_mf[['Ga','Gb']].min(axis=1)
  cors_mf['Gmax'] = cors_mf[['Ga','Gb']].max(axis=1)
  cors_mf['Ga'] = cors_mf['Gmin']
  cors_mf['Gb'] = cors_mf['Gmax']
  cors_mf = cors_mf.drop(columns=['Gmin','Gmax'])
  cors_mf = cors_mf.groupby(by=['Ga','Gb'], as_index=False).mean()
  cors_mf['Cluster'] = 'FM'
  cors_grouped = pandas.concat([cors_ff,cors_mm,cors_mf])
  cors_grouped = cors_grouped[['Ga','Gb','Cluster','wRho']]
  #LOG("DEBUG: cors_grouped nrows: %d ; ncols: %d:"%(cors_grouped.shape[0],cors_grouped.shape[1]))
  #
  sims = sims[['ENSGA','SEXA','ENSGB','SEXB','Ruzicka']]
  sims = pandas.merge(sims,genes.rename(columns={'ENSG':'ENSGA','symbol':'Ga'}), on=['ENSGA'], how='inner')
  sims = pandas.merge(sims, genes.rename(columns={'ENSG':'ENSGB','symbol':'Gb'}), on=['ENSGB'], how='inner')
  sims = sims[['Ga','Gb','SEXA','SEXB','Ruzicka']]
  sims = sims.drop_duplicates()
  #LOG("DEBUG: sims nrows: %d ; ncols: %d:"%(sims.shape[0],sims.shape[1]))
  sims_ff = sims[(sims.SEXA=='F')&(sims.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  sims_ff['Cluster']='F'
  sims_mm = sims[(sims.SEXA=='M')&(sims.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  sims_mm['Cluster']='M'
  ## MF: aggregate on Ga, Gb in alpha order.
  sims_mf = sims[((sims.SEXA=='F')&(sims.SEXB=='M'))|((sims.SEXA=='M')&(sims.SEXB=='F'))].drop(columns=['SEXA','SEXB'])
  sims_mf['Gmin'] = sims_mf[['Ga','Gb']].min(axis=1)
  sims_mf['Gmax'] = sims_mf[['Ga','Gb']].max(axis=1)
  sims_mf['Ga'] = sims_mf['Gmin']
  sims_mf['Gb'] = sims_mf['Gmax']
  sims_mf = sims_mf.drop(columns=['Gmin','Gmax'])
  sims_mf = sims_mf.groupby(by=['Ga','Gb'], as_index=False).mean()
  sims_mf['Cluster'] = 'FM'
  sims_grouped = pandas.concat([sims_ff,sims_mm,sims_mf])
  sims_grouped = sims_grouped[['Ga','Gb','Cluster','Ruzicka']]
  LOG("DEBUG: sims_grouped nrows: %d ; ncols: %d:"%(sims_grouped.shape[0],sims_grouped.shape[1]))
  #
  cmps = pandas.merge(cors_grouped, sims_grouped, on=['Ga','Gb','Cluster'])
  cmps = cmps[['Ga','Gb','Cluster','wRho','Ruzicka']]
  #
  return cmps

#############################################################################
def LOG(msg, file=sys.stdout):
  print(msg, file=file, flush=True)

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Exfiles similarity')
  parser.add_argument("--i_cor",dest="ifile_cor",help="input gene-gene correlation (TSV)")
  parser.add_argument("--i_sim",dest="ifile_sim",help="input gene-gene similarity (TSV)")
  parser.add_argument("--i_gene",dest="ifile_gene",help="input gene IDs (TSV)")
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("--decimals",type=int,default=3,help="output decimal places")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  PROG=os.path.basename(sys.argv[0])
  t0 = time.time()

  if args.verbose:
    print('Python: %s; pandas: %s; numpy: %s; scipy: %s'%(sys.version.split()[0],pandas.__version__, numpy.__version__, scipy.__version__), file=sys.stdout)

  if not args.ifile_gene: parser.error('Input gene file required.')
  if not args.ifile_cor: parser.error('Input cor file required.')
  if not args.ifile_sim: parser.error('Input sim file required.')

  genes = ReadGenes(args.ifile_gene, args.verbose)
  cors = ReadCorrfile(args.ifile_cor, args.verbose)
  sims = ReadSimfile(args.ifile_sim, args.verbose)

  cmps = GroupComparisons(cors, sims, genes, args.verbose)

  if args.ofile:
    print("=== Output file: %s"%args.ofile, file=sys.stdout)
    LOG("Output nrows: %d ; ncols: %d:"%(cmps.shape[0],cmps.shape[1]))
    cmps.round(args.decimals).to_csv(args.ofile, sep='\t', index=False)

  print("%s Elapsed: %ds"%(PROG,(time.time()-t0)), file=sys.stderr)
