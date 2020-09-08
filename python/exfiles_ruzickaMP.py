#!/usr/bin/env python3
"""exfiles_similarity.py

Expression profiles similarity computation.

 - Author: Jeremy Yang
 - Required: Python3, Pandas 0.22+

 - Input expression profiles format expected: TSV, 2 columns of identifiers (ENSG, SEX) followed by multiple columns of expression values.
 - File may be produced by gtex_rnaseq_prep_app.py

"""
#############################################################################
import sys,os,io,re,time,argparse,logging
import pandas,numpy,scipy,scipy.stats
import multiprocessing as mp

#############################################################################
def Ruzicka_MPWorker(iotuple):
  exfiles, i_start, i_end, n_idcols, minval, ofile, i_proc = iotuple
  idcols = list(range(n_idcols))
  idcoltags = exfiles.columns[idcols]
  datacols = list(range(n_idcols,exfiles.shape[1]))
  fout = open('%s_%02d'%(ofile,i_proc), 'w')
  fout.write('%s\tRuzicka\n'%('\t'.join([tag+'A' for tag in idcoltags]+[tag+'B' for tag in idcoltags])))
  for iA in range(i_start,i_end):
    A = exfiles.iloc[iA,datacols]
    for iB in range(iA+1, exfiles.shape[0]):
      B = exfiles.iloc[iB,datacols]
      s = sum(numpy.fmin(A,B))/sum(numpy.fmax(A,B))
      if numpy.isnan(s):
        continue
      elif s<minval:
        continue
      idvalsA = exfiles.iloc[iA,idcols].tolist()
      idvalsB = exfiles.iloc[iB,idcols].tolist()
      fout.write('%s\t%f\n'%('\t'.join(exfiles.iloc[iA,idcols].tolist()+exfiles.iloc[iB,idcols].tolist()),s))
      fout.flush()
  fout.close()

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
def CleanExfiles(exfiles):
  #exfiles = exfiles.fillna(0)
  exfiles = exfiles.dropna(how="any")
  return exfiles

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Exfiles similarity')
  parser.add_argument("--i", dest="ifile", help="input profiles, 1-row/gene (TSV)")
  parser.add_argument("--ruzicka_min", type=float, default=0, help="minimum values output")
  parser.add_argument("--n_idcols", type=int, default=2, help="number of ID cols preceding TPMs")
  parser.add_argument("--nproc", type=int, default=2, help="number of parallel processes")
  parser.add_argument("--o", dest="ofile", help="output (TSV)")
  parser.add_argument("-v", "--verbose", action="count")
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  t0 = time.time()

  if args.verbose:
    logging.info('Python: %s; pandas: %s; numpy: %s; scipy: %s'%(sys.version.split()[0],pandas.__version__, numpy.__version__, scipy.__version__))

  if not args.ifile:
    parser.error('Input file required.')

  exfiles = ReadExfiles(args.ifile)

  exfiles = CleanExfiles(exfiles)

  if not args.ofile:
    parser.error('Output file required.')

  t0 = time.time()
  idcols = list(range(args.n_idcols))
  idcoltags = exfiles.columns[idcols]

  iotuples = []
  i_proc=0;
  for i in range(args.nproc):
    i_proc+=1
    i_start = int((exfiles.shape[0])*(i/args.nproc)) 
    i_end = int((exfiles.shape[0])*((i+1)/args.nproc)) 
    iotuples.append((exfiles,i_start,i_end,args.n_idcols,args.ruzicka_min,args.ofile,i_proc))

  pool = mp.Pool(processes=args.nproc)
  i=0;
  for proc in pool.imap_unordered(Ruzicka_MPWorker, iotuples, chunksize=1):
    exfiles, i_start, i_end, n_idcols, ruzicka_min, ofile, i_proc  = iotuples[i]
    i+=1

  logging.info("Elapsed: %ds"%((time.time()-t0)))
