#!/usr/bin/env python3
"""exfiles_similarity_post.py

Post-process output from exfiles_similarity*.py|R.

Input cols:
	ENSGA,SEXA,ENSGB,SEXB,Ruzicka
	ENSGA,SEXA,ENSGB,SEXB,wRho

Output cols:
	ENSGA,ENSGB,Cluster,wRho,Ruzicka

where
	Cluster = F (F-F)
	Cluster = M (M-M)
	Cluster = FM (mean)

and FILTER:

	For all genes keep top 10 by combo score (sim*cor).
	Among rest, delete sim data below min cutoff.
	Among rest, delete cor data: (< min cutoff) & (> max anti-cor cutoff)

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

  cors_ff = cors[(cors.SEXA=='F')&(cors.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  cors_ff['Cluster']='F'
  cors_mm = cors[(cors.SEXA=='M')&(cors.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  cors_mm['Cluster']='M'

  ## MF: aggregate on ENSGA, ENSGB in alpha order.
  cors_mf = cors[((cors.SEXA=='F')&(cors.SEXB=='M'))|((cors.SEXA=='M')&(cors.SEXB=='F'))].drop(columns=['SEXA','SEXB'])
  cors_mf['ENSGmin'] = cors_mf[['ENSGA','ENSGB']].min(axis=1)
  cors_mf['ENSGmax'] = cors_mf[['ENSGA','ENSGB']].max(axis=1)
  cors_mf['ENSGA'] = cors_mf['ENSGmin']
  cors_mf['ENSGB'] = cors_mf['ENSGmax']
  cors_mf = cors_mf.drop(columns=['ENSGmin','ENSGmax'])
  cors_mf = cors_mf.groupby(by=['ENSGA','ENSGB'], as_index=False).mean()
  cors_mf['Cluster'] = 'FM'
  cors_grouped = pandas.concat([cors_ff,cors_mm,cors_mf])
  cors_grouped = cors_grouped[['ENSGA','ENSGB','Cluster','wRho']]
  LOG("DEBUG: cors_grouped nrows: %d ; ncols: %d:"%(cors_grouped.shape[0],cors_grouped.shape[1]))
  LOG("DEBUG: cors_grouped.Cluster.value_counts():\n%s"%(str(cors_grouped.Cluster.value_counts())))

  ###
  sims = sims[['ENSGA','SEXA','ENSGB','SEXB','Ruzicka']]
  #LOG("DEBUG: sims nrows: %d ; ncols: %d:"%(sims.shape[0],sims.shape[1]))
  sims_ff = sims[(sims.SEXA=='F')&(sims.SEXB=='F')].drop(columns=['SEXA','SEXB'])
  sims_ff['Cluster']='F'
  sims_mm = sims[(sims.SEXA=='M')&(sims.SEXB=='M')].drop(columns=['SEXA','SEXB'])
  sims_mm['Cluster']='M'
  ## MF: aggregate on ENSGA, ENSGB in alpha order.
  sims_mf = sims[((sims.SEXA=='F')&(sims.SEXB=='M'))|((sims.SEXA=='M')&(sims.SEXB=='F'))].drop(columns=['SEXA','SEXB'])
  sims_mf['ENSGmin'] = sims_mf[['ENSGA','ENSGB']].min(axis=1)
  sims_mf['ENSGmax'] = sims_mf[['ENSGA','ENSGB']].max(axis=1)
  sims_mf['ENSGA'] = sims_mf['ENSGmin']
  sims_mf['ENSGB'] = sims_mf['ENSGmax']
  sims_mf = sims_mf.drop(columns=['ENSGmin','ENSGmax'])
  sims_mf = sims_mf.groupby(by=['ENSGA','ENSGB'], as_index=False).mean()
  sims_mf['Cluster'] = 'FM'
  sims_grouped = pandas.concat([sims_ff,sims_mm,sims_mf])
  sims_grouped = sims_grouped[['ENSGA','ENSGB','Cluster','Ruzicka']]
  LOG("DEBUG: sims_grouped nrows: %d ; ncols: %d:"%(sims_grouped.shape[0],sims_grouped.shape[1]))
  LOG("DEBUG: sims_grouped.Cluster.value_counts():\n%s"%(str(sims_grouped.Cluster.value_counts())))
  #
  cmps = pandas.merge(cors_grouped, sims_grouped, on=['ENSGA','ENSGB','Cluster'])
  cmps = cmps[['ENSGA','ENSGB','Cluster','wRho','Ruzicka']]
  LOG("DEBUG: cmps.Cluster.value_counts():\n%s"%(str(cmps.Cluster.value_counts())))
  #
  return cmps

#############################################################################
def FilterComparisons2File(cmps, min_keep, min_sim, min_cor, max_anticor, decimals, ofile, verbose):
  fout = open(ofile, "w")
  tags = cmps.columns.tolist()
  fout.write('\t'.join(tags)+'\n')
  nrow_i=cmps.shape[0]
  LOG("FilterComparisons input: nrows: %d:"%(cmps.shape[0]))
  LOG('Min keep: %d'%min_keep)
  LOG('Min similarity: %f'%min_sim)
  LOG('Min correlation: %f'%min_cor)
  LOG('Max anticorrelation: %f'%max_anticor)

  cmps.dropna(how='any', inplace=True)
  cmps = cmps.round(decimals)
  cmps = cmps.set_index(['ENSGA'])

  ### Compute combo score and mark for keeping top-10 for each gene:
# Maybe first pivot into matrix?
  ### 
  cmps['Combo'] = cmps.wRho.abs() * cmps.Ruzicka
  cmps['out'] = False #flag: written
  n_out=0; n_genes_sub=0;
  for ensgA in cmps.index.get_level_values(0):
    cmps_this = cmps.loc[ensgA]
    if cmps_this.shape[0]>min_keep:
      combo_min = cmps_this.Combo.sort_values(ascending=False).tolist()[min_keep-1]
    else:
      n_genes_sub+=1
      combo_min = 0
    cmps.reset_index(drop=False)
    for row in cmps_this[cmps.Combo>=combo_min].itertuples():
      if not row.out:
        fout.write('\t'.join([str(getattr(row,tag)) for tag in tags])+'\n')
        n_out+=1
        if n_out%1000==0:
          LOG("Progress (top hits step): n_out: %d ; genes %d / %d (%.1f%%)"%(n_out,i_gene,ensgs.size,100*i_gene/ensgs.size), file=sys.stderr)
    cmps.out[((cmps.ENSGA==ensg)|(cmps.ENSGB==ensg))&(cmps.combo>=combo_min)] = True

  LOG("FilterComparisons genes with <%d similars: %d / %d: (%.1f%%)"%(min_keep,n_genes_sub,ensgs.size,100*n_genes_sub/ensgs.size))

  i_row=0
  for row in cmps.itertuples():
    i_row+=1
    if row.out:
      continue
    elif row.Ruzicka>min_sim or row.wRho>=min_cor or row.wRho<=max_anticor:
      fout.write('\t'.join([str(getattr(row,tag)) for tag in tags])+'\n')
      n_out+=1
      if n_out%1000==0:
        LOG("Progress (sim/cor filter step): n_out: %d ; rows %d / %d (%.1f%%)"%(n_out,i_row,cmps.shape[0],100*i_row/cmps.shape[0]), file=sys.stderr)

  LOG("FilterComparisons output: nrows: %d: (%.1f%%)"%(n_out,100*n_out/nrow_i))

  fout.close()

#############################################################################
def LOG(msg, file=sys.stdout):
  print(msg, file=file, flush=True)

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Exfiles similarity')
  parser.add_argument("--i_cor",dest="ifile_cor",help="input gene-gene correlation (TSV)")
  parser.add_argument("--i_sim",dest="ifile_sim",help="input gene-gene similarity (TSV)")
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("--o_unfiltered",dest="ofile_unfiltered",help="output (TSV)")
  parser.add_argument("--min_keep",type=int,default=10,help="min count, sim-genes to keep per gene")
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

  if args.ofile_unfiltered:
    LOG("=== Output unfiltered file: %s"%args.ofile_unfiltered)
    cmps.round(args.decimals).to_csv(args.ofile_unfiltered, sep='\t', index=False)

  if args.ofile:
    ### Directly to file saves memory.
    FilterComparisons2File(cmps, args.min_keep, args.min_sim, args.min_cor, args.max_anticor, args.decimals, args.ofile, args.verbose)

  print("%s Elapsed: %ds"%(PROG,(time.time()-t0)), file=sys.stderr)
