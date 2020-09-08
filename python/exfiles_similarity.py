#!/usr/bin/env python3
"""exfiles_similarity.py

Expression profiles similarity computation.

 - Author: Jeremy Yang
 - Required: Python3, Pandas 0.22+

 - Input expression profiles format expected: TSV, 2 columns of identifiers (ENSG, SEX) followed by multiple columns of expression values.
 - File may be produced by gtex_rnaseq_prep_app.py

 - Parallelize? import multiprocessing?

"""
#############################################################################
import sys,os,io,re,time,argparse,logging,tqdm
import pandas,numpy,scipy,scipy.stats

#############################################################################
### Pearson correlation coefficient.
#############################################################################
def Pearson(exfiles, idcols, datacols, minval, ofile):
  idcoltags = exfiles.columns[idcols]
  logging.debug("Pearson IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())))
  logging.debug("idcols = %s ; datacols = %s"%(str(idcols),str(datacols)))
  logging.debug("idcoltags = %s"%(str(idcoltags)))
  fout = open(ofile, 'w')
  fout.write('%s\tPearson\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  n_out=0; n_nan=0; n_submin=0;
  for iA in range(exfiles.shape[0]):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      c = numpy.corrcoef(numpy.array([A.tolist(),B.tolist()]))[0][1]
      if numpy.isnan(c):
        n_nan+=1
        continue
      elif c<minval:
        n_submin+=1
        continue
      n_out+=1
      fout.write('%s\t%f\n'%('\t'.join(exfiles.iloc[iA,idcols].tolist()+exfiles.iloc[iB,idcols].tolist()),c))
  logging.info("n_out: %d; n_nan: %d; n_submin: %d"%(n_out, n_nan, n_submin))

#############################################################################
def Spearman(exfiles, idcols, datacols, minval, ofile):
  idcoltags = exfiles.columns[idcols]
  logging.debug("Spearman IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())))
  logging.debug("idcols = %s ; datacols = %s"%(str(idcols),str(datacols)))
  logging.debug("idcoltags = %s"%(str(idcoltags)))
  fout = open(ofile, 'w')
  fout.write('%s\tSpearmanRho\tSpearmanP\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  n_out=0; n_nan=0; n_submin=0; n_err=0;
  for iA in range(exfiles.shape[0]):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      try:
        rho,pval = scipy.stats.spearmanr(A,B)
      except Exception as e:
        n_err+=1
        continue
      if numpy.isnan(rho):
        n_nan+=1
        continue
      elif rho<minval:
        n_submin+=1
        continue
      n_out+=1
      fout.write('%s\t%f\t%f\n'%('\t'.join(exfiles.iloc[iA,idcols].tolist()+exfiles.iloc[iB,idcols].tolist()),rho,pval))
  logging.info("n_out: %d; n_nan: %d; n_submin: %d"%(n_out, n_nan, n_submin))

#############################################################################
def ABC(A,B):
  abc = 0.0
  for i in range(len(A)-1):
    Amid = numpy.mean(A[i:i+1])
    Bmid = numpy.mean(B[i:i+1])
    if (A[i]>=B[i]):
      if (A[i+1]>=B[i+1]):
        abc = abc + (AULS(A[i], A[i+1], 1) - AULS(B[i], B[i+1], 1))
      else:
        abc = abc + (AULS(A[i], Amid, .5) - AULS(B[i], Bmid , .5))
        abc = abc + (AULS(Bmid, B[i+1], .5) - AULS(Amid, A[i+1], .5))
    else:
      if (A[i+1]<B[i+1]):
        abc = abc + (AULS(B[i], B[i+1], 1) - AULS(A[i], A[i+1], 1))
      else:
        abc = abc + (AULS(B[i], Bmid, .5) - AULS(A[i], Amid , .5))
        abc = abc + (AULS(Amid, A[i+1], .5) - AULS(Bmid, B[i+1], .5))
  return(abc)
#
#############################################################################
### Area Under Line Segment
def AULS(y1, y2, w):
  a = min(y1,y2) * w
  a = a + 0.5 * w * abs(y1-y2)
  return(a)

#############################################################################
def ABC(exfiles, idcols, datacols, minval, ofile):
  idcoltags = exfiles.columns[idcols]
  logging.debug("ABC IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())))
  logging.debug("idcols = %s ; datacols = %s"%(str(idcols),str(datacols)))
  logging.debug("idcoltags = %s"%(str(idcoltags)))
  fout = open(ofile, 'w')
  fout.write('%s\tABC\tABC_sim\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  n_out=0; n_nan=0; n_submin=0;
  for iA in range(exfiles.shape[0]):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      abc = ABC(A,B)
      abc_sim = (1 / (1 + abc/len(A)))
      if numpy.isnan(abc_sim):
        n_nan+=1
        continue
      elif abc_sim<minval:
        n_submin+=1
        continue
      n_out+=1
      fout.write('%s\t%f\t%f\n'%('\t'.join(exfiles.iloc[iA,idcols].tolist()+exfiles.iloc[iB,idcols].tolist()),abc,abc_sim))
  logging.info("n_out: %d; n_nan: %d; n_submin: %d"%(n_out, n_nan, n_submin))

#############################################################################
def Cosine(exfiles, idcols, datacols, minval, ofile):
  idcoltags = exfiles.columns[idcols]
  logging.debug("Cosine IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())))
  logging.debug("idcols = %s ; datacols = %s"%(str(idcols),str(datacols)))
  logging.debug("idcoltags = %s"%(str(idcoltags)))
  fout = open(ofile, 'w')
  fout.write('%s\tCosine\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  #First compute |V| for each vector.
  Vlen = numpy.ndarray(shape=(exfiles.shape[0],1), dtype=float)
  for i in range(exfiles.shape[0]):
    V = exfiles.iloc[i,datacols]
    Vlen[i] = numpy.sqrt(numpy.dot(V,V))
  n_out=0; n_nan=0; n_submin=0;
  for iA in range(exfiles.shape[0]):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      AB = numpy.dot(A,B)
      den = Vlen[iA]*Vlen[iB]
      c = (AB / den) if den>0 else numpy.nan
      if numpy.isnan(c):
        n_nan+=1
        continue
      elif c<minval:
        n_submin+=1
        continue
      idvalsA = exfiles.iloc[iA,idcols].tolist()
      idvalsB = exfiles.iloc[iB,idcols].tolist()
      fout.write('%s\t%f\n'%('\t'.join(exfiles.iloc[iA,idcols].tolist()+exfiles.iloc[iB,idcols].tolist()),c))
      n_out+=1
  logging.info("n_out: %d; n_nan: %d; n_submin: %d"%(n_out, n_nan, n_submin))

#############################################################################
###  Tanimoto(A,B) = (A %*% B) / (A %*% A + B %*% B - A %*% B)
###  where "%*%" denotes vector dot product
### Vectorize: use matrix form of numpy.dot(), or maybe numpy.matmul()
### ~14min for ~1000 input profiles, ~500k calculations.
#############################################################################
def Tanimoto(exfiles, idcols, datacols, minval, ofile):
  idcoltags = exfiles.columns[idcols]
  logging.debug("Tanimoto IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())))
  logging.debug("idcols = %s ; datacols = %s"%(str(idcols),str(datacols)))
  logging.debug("idcoltags = %s"%(str(idcoltags)))
  fout = open(ofile, 'w')
  fout.write('%s\tTanimoto\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  #First compute |V|^2 for each vector.
  VV = numpy.ndarray(shape=(exfiles.shape[0],1), dtype=float)
  for i in range(exfiles.shape[0]):
    V = exfiles.iloc[i,datacols]
    VV[i] = numpy.dot(V,V)
  n_out=0; n_nan=0; n_submin=0;
  for iA in range(exfiles.shape[0]):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      AA = VV[iA]
      BB = VV[iB]
      AB = numpy.dot(A,B)
      den = AA + BB - AB
      t = (AB / den) if den>0 else numpy.nan
      if numpy.isnan(t):
        n_nan+=1
        continue
      elif t<minval:
        n_submin+=1
        continue
      n_out+=1
      fout.write('%s\t%f\n'%('\t'.join(exfiles.iloc[iA,idcols].tolist()+exfiles.iloc[iB,idcols].tolist()),t))
  logging.info("n_out: %d; n_nan: %d; n_submin: %d"%(n_out, n_nan, n_submin))

#############################################################################
def Ruzicka(exfiles, idcols, datacols, minval, ofile):
  tq = None;
  t0 = time.time()
  idcoltags = exfiles.columns[idcols]
  logging.debug("Ruzicka IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())))
  logging.debug("idcols = %s ; datacols = %s"%(str(idcols),str(datacols)))
  logging.debug("idcoltags = %s"%(str(idcoltags)))
  fout = open(ofile, 'w')
  fout.write('%s\tRuzicka\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  n_out=0; n_nan=0; n_submin=0; n_calc=0;
  n_calc_total = exfiles.shape[0]*(exfiles.shape[0]-1)/2
  for iA in range(exfiles.shape[0]):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      s = sum(numpy.fmin(A,B))/sum(numpy.fmax(A,B))
      n_calc+=1
      if tq is not None: tq.update()
      if numpy.isnan(s):
        n_nan+=1
        continue
      elif s<minval:
        n_submin+=1
        continue
      idvalsA = exfiles.iloc[iA,idcols].tolist()
      idvalsB = exfiles.iloc[iB,idcols].tolist()
      n_out+=1
      fout.write('%s\t%f\n'%('\t'.join(exfiles.iloc[iA,idcols].tolist()+exfiles.iloc[iB,idcols].tolist()),s))
    if not tq: tq = tqdm.tqdm(total=n_calc_total, unit="calcs")
    #logging.info('Progress: %d / %d (%.1f%%) ; elapsed: %s'%(n_calc, n_calc_total, 100*n_calc/n_calc_total, time.strftime("%H:%M:%S",time.gmtime(time.time()-t0))))
  logging.info("n_out: %d; n_nan: %d; n_submin: %d"%(n_out, n_nan, n_submin))

#############################################################################
def ReadExfiles(ifile):
  fin = open(ifile)
  logging.info('=== Expression profiles datafile: %s'%fin.name)
  exfiles = pandas.read_csv(fin, sep='\t')
  logging.info("Exfiles dataset nrows: %d ; ncols: %d:"%(exfiles.shape[0],exfiles.shape[1]))
  for name,val in exfiles.SEX.value_counts().sort_index().iteritems():
    logging.info('\tExfiles (SEX=%s): %5d'%(name,val))
  return exfiles

#############################################################################
### Replace missing (nan) values with 0.
#############################################################################
def CleanExfiles(exfiles):
  exfiles = exfiles.fillna(0)
  return exfiles

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Exfiles similarity')
  parser.add_argument("--i",dest="ifile",help="input profiles, 1-row/gene (TSV)")
  parser.add_argument("--o_ruzicka",dest="ofile_ruzicka",help="output (TSV)")
  parser.add_argument("--ruzicka_min",type=float,default=0,help="minimum values output")
  parser.add_argument("--o_cosine",dest="ofile_cosine",help="output (TSV)")
  parser.add_argument("--cosine_min",type=float,default=0,help="minimum values output")
  parser.add_argument("--o_tanimoto",dest="ofile_tanimoto",help="output (TSV)")
  parser.add_argument("--tanimoto_min",type=float,default=0,help="minimum values output")
  parser.add_argument("--o_abc",dest="ofile_abc",help="output (TSV)")
  parser.add_argument("--abc_min",type=float,default=0,help="maximum values output")
  parser.add_argument("--o_pearson",dest="ofile_pearson",help="output (TSV)")
  parser.add_argument("--pearson_min",type=float,default=-1,help="minimum values output")
  parser.add_argument("--o_spearman",dest="ofile_spearman",help="output (TSV)")
  parser.add_argument("--spearman_min",type=float,default=-1,help="minimum values output")
  parser.add_argument("--n_idcols",type=int,default=2,help="number of ID cols preceding TPMs")
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  t0 = time.time()

  if args.verbose:
    logging.info('Python: %s; pandas: %s; numpy: %s; scipy: %s'%(sys.version.split()[0],pandas.__version__, numpy.__version__, scipy.__version__))

  if not args.ifile:
    parser.error('Input file required.')

  exfiles = ReadExfiles(args.ifile)

  exfiles = CleanExfiles(exfiles)

  if args.ofile_pearson:
    Pearson(exfiles, list(range(args.n_idcols)), list(range(args.n_idcols,exfiles.shape[1])), args.pearson_min, args.ofile_pearson)

  if args.ofile_spearman:
    Spearman(exfiles, list(range(args.n_idcols)), list(range(args.n_idcols,exfiles.shape[1])), args.spearman_min, args.ofile_spearman)

  if args.ofile_cosine:
    Cosine(exfiles, list(range(args.n_idcols)), list(range(args.n_idcols,exfiles.shape[1])), args.cosine_min, args.ofile_cosine)

  if args.ofile_ruzicka:
    Ruzicka(exfiles, list(range(args.n_idcols)), list(range(args.n_idcols,exfiles.shape[1])), args.ruzicka_min, args.ofile_ruzicka)

  if args.ofile_tanimoto:
    Tanimoto(exfiles, list(range(args.n_idcols)), list(range(args.n_idcols,exfiles.shape[1])), args.tanimoto_min, args.ofile_tanimoto)

  if args.ofile_abc:
    ABC(exfiles, list(range(args.n_idcols)), list(range(args.n_idcols,exfiles.shape[1])), args.abc_min, args.ofile_abc)

  logging.info("Elapsed: %ds"%((time.time()-t0)))
