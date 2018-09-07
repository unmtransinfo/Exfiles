#!/usr/bin/env python3
"""exfiles_similarity.py

Expression profiles similarity computation.

 - Author: Jeremy Yang
 - Required: Python3, Pandas 0.22+

 - Input expression profiles format expected: TSV, 2 columns of identifiers (ENSG, SEX) followed by multiple columns of expression values.
 - File may be produced by gtex_rnaseq_prep_app.py

"""
#############################################################################
import sys,os,io,re,time,argparse
import numpy,scipy,scipy.stats
import pandas


#############################################################################
### Pearson correlation coefficient.
#############################################################################
def Pearson_NxN(exfiles, idcols, datacols, minval, ofile, verbose):
  idcoltags = exfiles.columns[idcols]
  print("DEBUG: Pearson_NxN IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())), file=sys.stderr)
  print("DEBUG: idcols = %s ; datacols = %s"%(str(idcols),str(datacols)), file=sys.stderr)
  print("DEBUG: idcoltags = %s"%(str(idcoltags)), file=sys.stderr)
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
  print("n_out: %d"%(n_out), file=sys.stdout)
  print("n_nan: %d"%(n_nan), file=sys.stdout)
  print("n_submin: %d"%(n_submin), file=sys.stdout)

#############################################################################
def Spearman_NxN(exfiles, idcols, datacols, minval, ofile, verbose):
  idcoltags = exfiles.columns[idcols]
  print("DEBUG: Spearman_NxN IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())), file=sys.stderr)
  print("DEBUG: idcols = %s ; datacols = %s"%(str(idcols),str(datacols)), file=sys.stderr)
  print("DEBUG: idcoltags = %s"%(str(idcoltags)), file=sys.stderr)
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
  print("n_out: %d"%(n_out), file=sys.stdout)
  print("n_nan: %d"%(n_nan), file=sys.stdout)
  print("n_submin: %d"%(n_submin), file=sys.stdout)

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
def ABC_NxN(exfiles, idcols, datacols, minval, ofile, verbose):
  idcoltags = exfiles.columns[idcols]
  print("DEBUG: ABC_NxN IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())), file=sys.stderr)
  print("DEBUG: idcols = %s ; datacols = %s"%(str(idcols),str(datacols)), file=sys.stderr)
  print("DEBUG: idcoltags = %s"%(str(idcoltags)), file=sys.stderr)
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
  print("n_out: %d"%(n_out), file=sys.stdout)
  print("n_nan: %d"%(n_nan), file=sys.stdout)
  print("n_submin: %d"%(n_submin), file=sys.stdout)

#############################################################################
def Cosine_NxN(exfiles, idcols, datacols, minval, ofile, verbose):
  idcoltags = exfiles.columns[idcols]
  print("DEBUG: Cosine_NxN IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())), file=sys.stderr)
  print("DEBUG: idcols = %s ; datacols = %s"%(str(idcols),str(datacols)), file=sys.stderr)
  print("DEBUG: idcoltags = %s"%(str(idcoltags)), file=sys.stderr)
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
  print("n_out: %d"%(n_out), file=sys.stdout)
  print("n_nan: %d"%(n_nan), file=sys.stdout)
  print("n_submin: %d"%(n_submin), file=sys.stdout)

#############################################################################
###  Tanimoto(A,B) = (A %*% B) / (A %*% A + B %*% B - A %*% B)
###  where "%*%" denotes vector dot product
### Vectorize: use matrix form of numpy.dot(), or maybe numpy.matmul()
### ~14min for ~1000 input profiles, ~500k calculations.
#############################################################################
def Tanimoto_NxN(exfiles, idcols, datacols, minval, ofile, verbose):
  idcoltags = exfiles.columns[idcols]
  print("DEBUG: Tanimoto_NxN IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())), file=sys.stderr)
  print("DEBUG: idcols = %s ; datacols = %s"%(str(idcols),str(datacols)), file=sys.stderr)
  print("DEBUG: idcoltags = %s"%(str(idcoltags)), file=sys.stderr)
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
  print("n_out: %d"%(n_out), file=sys.stdout)
  print("n_nan: %d"%(n_nan), file=sys.stdout)
  print("n_submin: %d"%(n_submin), file=sys.stdout)

#############################################################################
def Ruzicka_NxN(exfiles, idcols, datacols, minval, ofile, verbose):
  idcoltags = exfiles.columns[idcols]
  print("DEBUG: Ruzicka_NxN IN: nrows = %d, cols: %s"%(exfiles.shape[0],str(exfiles.columns.tolist())), file=sys.stderr)
  print("DEBUG: idcols = %s ; datacols = %s"%(str(idcols),str(datacols)), file=sys.stderr)
  print("DEBUG: idcoltags = %s"%(str(idcoltags)), file=sys.stderr)
  fout = open(ofile, 'w')
  fout.write('%s\tRuzicka\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  n_out=0; n_nan=0; n_submin=0;
  for iA in range(exfiles.shape[0]):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      s = sum(numpy.fmin(A,B))/sum(numpy.fmax(A,B))
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
  print("n_out: %d"%(n_out), file=sys.stdout)
  print("n_nan: %d"%(n_nan), file=sys.stdout)
  print("n_submin: %d"%(n_submin), file=sys.stdout)

#############################################################################
def ReadExfiles(ifile, verbose):
  fin = open(ifile)
  print('=== Expression profiles datafile: %s'%fin.name, file=sys.stdout)
  exfiles = pandas.read_csv(fin, sep='\t')
  print("Exfiles dataset nrows: %d ; ncols: %d:"%(exfiles.shape[0],exfiles.shape[1]), file=sys.stdout)
  for name,val in exfiles.SEX.value_counts().sort_index().iteritems():
    print('\tExfiles (SEX=%s): %5d'%(name,val), file=sys.stdout)
  return exfiles

#############################################################################
### Replace missing (nan) values with 0.
#############################################################################
def CleanExfiles(exfiles, verbose):
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
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("-v","--verbose",action="count")
  args = parser.parse_args()

  PROG=os.path.basename(sys.argv[0])
  t0 = time.time()

  if args.verbose:
    print('Python: %s; pandas: %s; numpy: %s; scipy: %s'%(sys.version.split()[0],pandas.__version__, numpy.__version__, scipy.__version__), file=sys.stdout)

  if not args.ifile:
    parser.error('Input file required.')

  exfiles = ReadExfiles(args.ifile, args.verbose)

  exfiles = CleanExfiles(exfiles, args.verbose)

  if args.ofile_pearson:
    Pearson_NxN(exfiles, [0,1], list(range(2,exfiles.shape[1])), args.pearson_min, args.ofile_pearson, args.verbose)

  if args.ofile_spearman:
    Spearman_NxN(exfiles, [0,1], list(range(2,exfiles.shape[1])), args.spearman_min, args.ofile_spearman, args.verbose)

  if args.ofile_cosine:
    Cosine_NxN(exfiles, [0,1], list(range(2,exfiles.shape[1])), args.cosine_min, args.ofile_cosine, args.verbose)

  if args.ofile_ruzicka:
    Ruzicka_NxN(exfiles, [0,1], list(range(2,exfiles.shape[1])), args.ruzicka_min, args.ofile_ruzicka, args.verbose)

  if args.ofile_tanimoto:
    Tanimoto_NxN(exfiles, [0,1], list(range(2,exfiles.shape[1])), args.tanimoto_min, args.ofile_tanimoto, args.verbose)

  if args.ofile_abc:
    ABC_NxN(exfiles, [0,1], list(range(2,exfiles.shape[1])), args.abc_min, args.ofile_abc, args.verbose)

  print("%s Elapsed: %ds"%(PROG,(time.time()-t0)), file=sys.stderr)
