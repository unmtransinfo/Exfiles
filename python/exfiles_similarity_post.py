#!/usr/bin/env python3
"""exfiles_similarity_post.py

Post-process output from comparison calculations.

Input cols:
	ENSGA,SEXA,ENSGB,SEXB,Ruzicka
	ENSGA,SEXA,ENSGB,SEXB,wRho

Output cols:
	ENSGA,ENSGB,Group,wRho,Ruzicka

where
	Group = F|M|C

F = F vs F
M = M vs M
C = Combined-sexes, profiles = (F+M)/2

(No longer filtering in this code. MIN_KEEP flawed idea since it helps 
with Search mode but loses many comparisons for Compare mode. Cutoffs
now in similarity and correlation calculation code.)

"""
#############################################################################
import sys,os,io,re,time,argparse,logging
import pandas,numpy,scipy,scipy.stats

#############################################################################
def ReadCorrfile(ifile):
  logging.info('=== Correlation datafile: %s'%ifile)
  cors = pandas.read_csv(ifile, sep='\t', na_values=['','NA','NaN'], compression='infer')
  cors.dropna(how='any', inplace=True)
  logging.info("Correlation dataset nrows: %d ; ncols: %d:"%(cors.shape[0],cors.shape[1]))
  logging.info("Correlation cols: %s:"%(str(cors.columns.tolist())))
  return cors

#############################################################################
def ReadSimfile(ifile):
  logging.info('=== Similarity datafile: %s'%ifile)
  sims = pandas.read_csv(ifile, sep='\t', na_values=['','NA','NaN'], compression='infer')
  sims.dropna(how='any', inplace=True)
  logging.info("Similarity dataset nrows: %d ; ncols: %d:"%(sims.shape[0],sims.shape[1]))
  logging.info("Similarity cols: %s:"%(str(sims.columns.tolist())))
  return sims

#############################################################################
def GroupComparisons(cors, sims):
  cors = cors[['ENSGA','SEXA','ENSGB','SEXB','wRho']]
  logging.debug("cors nrows: %d ; ncols: %d:"%(cors.shape[0],cors.shape[1]))
  cors_f = cors[(cors.SEXA=='F')&(cors.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  cors_f['Group']='F'
  cors_m = cors[(cors.SEXA=='M')&(cors.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  cors_m['Group']='M'
  cors_n = cors[(cors.SEXA=='C')&(cors.SEXB=='C')].drop(columns=['SEXA','SEXB'])
  cors_n['Group'] = 'C'
  #
  cors_grouped = pandas.concat([cors_f,cors_m,cors_n])
  cors_grouped = cors_grouped[['ENSGA','ENSGB','Group','wRho']]
  logging.debug("cors_grouped nrows: %d ; ncols: %d:"%(cors_grouped.shape[0],cors_grouped.shape[1]))
  logging.debug("cors_grouped.Group.value_counts():\n%s"%(str(cors_grouped.Group.value_counts())))
  ###
  sims = sims[['ENSGA','SEXA','ENSGB','SEXB','Ruzicka']]
  logging.debug("sims nrows: %d ; ncols: %d:"%(sims.shape[0],sims.shape[1]))
  sims_f = sims[(sims.SEXA=='F')&(sims.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  sims_f['Group']='F'
  sims_m = sims[(sims.SEXA=='M')&(sims.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  sims_m['Group']='M'
  sims_n = sims[(sims.SEXA=='C')&(sims.SEXB=='C')].drop(columns=['SEXA','SEXB'])
  sims_n['Group'] = 'C'
  #
  sims_grouped = pandas.concat([sims_f,sims_m,sims_n])
  sims_grouped = sims_grouped[['ENSGA','ENSGB','Group','Ruzicka']]
  logging.debug("sims_grouped nrows: %d ; ncols: %d:"%(sims_grouped.shape[0],sims_grouped.shape[1]))
  logging.debug("sims_grouped.Group.value_counts():\n%s"%(str(sims_grouped.Group.value_counts())))
  #
  cmps = pandas.merge(cors_grouped, sims_grouped, on=['ENSGA','ENSGB','Group'])
  cmps = cmps[['ENSGA','ENSGB','Group','wRho','Ruzicka']]
  logging.debug("cmps.Group.value_counts():\n%s"%(str(cmps.Group.value_counts())))
  #
  return cmps

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Exfiles similarity')
  parser.add_argument("--i_cor",dest="ifile_cor",help="input gene-gene correlation (TSV)")
  parser.add_argument("--i_sim",dest="ifile_sim",help="input gene-gene similarity (TSV)")
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("--min_sim",type=float,default=.7,help="min similarity")
  parser.add_argument("--min_cor",type=float,default=.7,help="min correlation")
  parser.add_argument("--max_anticor",type=float,default=-.7,help="max anti-correlation")
  parser.add_argument("--decimals",type=int,default=3,help="output decimal places")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  t0 = time.time()

  if args.verbose:
    logging.info('Python: %s; pandas: %s; numpy: %s; scipy: %s'%(sys.version.split()[0],pandas.__version__, numpy.__version__, scipy.__version__))

  if not args.ifile_cor: parser.error('Input cor file required.')
  if not args.ifile_sim: parser.error('Input sim file required.')

  cors = ReadCorrfile(args.ifile_cor)
  sims = ReadSimfile(args.ifile_sim)
  cmps = GroupComparisons(cors, sims)

  if args.ofile:
    logging.info("=== Output file: %s"%args.ofile)
    cmps.round(args.decimals).to_csv(args.ofile, sep='\t', index=False)

  logging.info("Elapsed: %ds"%((time.time()-t0)))
