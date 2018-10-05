#!/usr/bin/env python3
"""exfiles_similarity_post.py

Post-process output from comparison calculations.

Input cols:
	ENSGA,SEXA,ENSGB,SEXB,Ruzicka
	ENSGA,SEXA,ENSGB,SEXB,wRho

Output cols:
	ENSGA,ENSGB,Group,wRho,Ruzicka

where
	Group = F|M|N

F = F vs F
M = M vs M
N = Non-sexed, profiles = (F+M)/2

(No longer filtering in this code. MIN_KEEP flawed idea since it helps 
with Search mode but loses many comparisons for Compare mode. Cutoffs
now in similarity and correlation calculation code.)

"""
#############################################################################
import sys,os,io,re,time,argparse
import pandas,numpy,scipy,scipy.stats

#############################################################################
def ReadCorrfile(ifile, verbose):
  print('=== Correlation datafile: %s'%ifile, file=sys.stdout)
  cors = pandas.read_csv(ifile, sep='\t', na_values=['','NA','NaN'], compression='infer')
  cors.dropna(how='any', inplace=True)
  print("Correlation dataset nrows: %d ; ncols: %d:"%(cors.shape[0],cors.shape[1]), file=sys.stdout)
  print("Correlation cols: %s:"%(str(cors.columns.tolist())), file=sys.stdout)
  return cors

#############################################################################
def ReadSimfile(ifile, verbose):
  print('=== Similarity datafile: %s'%ifile, file=sys.stdout)
  sims = pandas.read_csv(ifile, sep='\t', na_values=['','NA','NaN'], compression='infer')
  sims.dropna(how='any', inplace=True)
  print("Similarity dataset nrows: %d ; ncols: %d:"%(sims.shape[0],sims.shape[1]), file=sys.stdout)
  print("Similarity cols: %s:"%(str(sims.columns.tolist())), file=sys.stdout)
  return sims

#############################################################################
def GroupComparisons(cors, sims, verbose):
  cors = cors[['ENSGA','SEXA','ENSGB','SEXB','wRho']]
  #LOG("DEBUG: cors nrows: %d ; ncols: %d:"%(cors.shape[0],cors.shape[1]))
  cors_f = cors[(cors.SEXA=='F')&(cors.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  cors_f['Group']='F'
  cors_m = cors[(cors.SEXA=='M')&(cors.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  cors_m['Group']='M'
  cors_n = cors[(cors.SEXA=='N')&(cors.SEXB=='N')].drop(columns=['SEXA','SEXB'])
  cors_n['Group'] = 'N'
  #
  cors_grouped = pandas.concat([cors_f,cors_m,cors_n])
  cors_grouped = cors_grouped[['ENSGA','ENSGB','Group','wRho']]
  LOG("DEBUG: cors_grouped nrows: %d ; ncols: %d:"%(cors_grouped.shape[0],cors_grouped.shape[1]))
  LOG("DEBUG: cors_grouped.Group.value_counts():\n%s"%(str(cors_grouped.Group.value_counts())))
  ###
  sims = sims[['ENSGA','SEXA','ENSGB','SEXB','Ruzicka']]
  #LOG("DEBUG: sims nrows: %d ; ncols: %d:"%(sims.shape[0],sims.shape[1]))
  sims_f = sims[(sims.SEXA=='F')&(sims.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  sims_f['Group']='F'
  sims_m = sims[(sims.SEXA=='M')&(sims.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  sims_m['Group']='M'
  sims_n = sims[(sims.SEXA=='N')&(sims.SEXB=='N')].drop(columns=['SEXA','SEXB'])
  sims_n['Group'] = 'N'
  #
  sims_grouped = pandas.concat([sims_f,sims_m,sims_n])
  sims_grouped = sims_grouped[['ENSGA','ENSGB','Group','Ruzicka']]
  LOG("DEBUG: sims_grouped nrows: %d ; ncols: %d:"%(sims_grouped.shape[0],sims_grouped.shape[1]))
  LOG("DEBUG: sims_grouped.Group.value_counts():\n%s"%(str(sims_grouped.Group.value_counts())))
  #
  cmps = pandas.merge(cors_grouped, sims_grouped, on=['ENSGA','ENSGB','Group'])
  cmps = cmps[['ENSGA','ENSGB','Group','wRho','Ruzicka']]
  LOG("DEBUG: cmps.Group.value_counts():\n%s"%(str(cmps.Group.value_counts())))
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
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("--min_sim",type=float,default=.7,help="min similarity")
  parser.add_argument("--min_cor",type=float,default=.7,help="min correlation")
  parser.add_argument("--max_anticor",type=float,default=-.7,help="max anti-correlation")
  parser.add_argument("--decimals",type=int,default=3,help="output decimal places")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  PROG=os.path.basename(sys.argv[0])
  t0 = time.time()

  if args.verbose:
    print('Python: %s; pandas: %s; numpy: %s; scipy: %s'%(sys.version.split()[0],pandas.__version__, numpy.__version__, scipy.__version__), file=sys.stdout)

  if not args.ifile_cor: parser.error('Input cor file required.')
  if not args.ifile_sim: parser.error('Input sim file required.')

  cors = ReadCorrfile(args.ifile_cor, args.verbose)
  sims = ReadSimfile(args.ifile_sim, args.verbose)

  cmps = GroupComparisons(cors, sims, args.verbose)

  if args.ofile:
    LOG("=== Output file: %s"%args.ofile)
    cmps.round(args.decimals).to_csv(args.ofile, sep='\t', index=False)

  print("%s Elapsed: %ds"%(PROG,(time.time()-t0)), file=sys.stderr)
