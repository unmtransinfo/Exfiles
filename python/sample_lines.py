#!/usr/bin/env python
#############################################################################
### sample_lines.py
### Jeremy Yang
### 23 Oct 2012
#############################################################################
import sys,os,getopt
import random

PROG=os.path.basename(sys.argv[0])


#############################################################################
def SampleKofN(ifile, k, n, verbose):
  fin=open(ifile)

  if not n:
    n=0
    while True:
      line=fin.readline()
      if not line: break
      n+=1
    fin.close()
    fin=open(ifile)
  print >>sys.stderr, "n = %d"%(n)

  i_sams=random.sample(range(n),k)
  i_sams.sort()
  i_sams_hash={}
  for i in i_sams:
    i_sams_hash[i]=True

  n_in=0; n_out=0;
  while True:
    line=fin.readline()
    if not line: break
    if i_sams_hash.has_key(n_in):
      sys.stdout.write(line)
      n_out+=1
      if n_out>=k: break
    n_in+=1

  print >>sys.stderr, "n_in: %d"%(n_in)
  print >>sys.stderr, "n_out: %d"%(n_out)

#############################################################################
def SampleP(ifile, p, verbose):
  if ifile=='-':
    fin = sys.stdin
  else:
    fin=open(ifile)

  n_in=0; n_out=0;
  while True:
    line=fin.readline()
    if not line: break
    if random.random()<p:
      sys.stdout.write(line)
      n_out+=1
    n_in+=1

  print >>sys.stderr, "n_in: %d"%(n_in)
  print >>sys.stderr, "n_out: %d"%(n_out)
  print >>sys.stderr, "p: %.2f ; p_sample: %.2f"%(p,float(n_out)/n_in)

#############################################################################
if __name__=='__main__':

  usage='''
  %(PROG)s - randomly sample lines from file

  required:
  --i INFILE ............... input file (or "-" for stdin)

  either:
  --k K .................... choose K lines
  or
  --p P .................... probability each line (0-1)

  options:
  --n N .................... population size [default=all-lines] (disallowed for stdin)
  --v ...................... verbose
  --h ...................... this help
'''%{'PROG':PROG}

  def ErrorExit(msg):
    print >>sys.stderr,msg
    sys.exit(1)

  verbose=0;
  k=None; n=None; ifile=None;
  opts,pargs = getopt.getopt(sys.argv[1:],'',['h','v','vv','i=','n=','k=', 'p='])
  if not opts: ErrorExit(usage)
  for (opt,val) in opts:
    if opt=='--h': ErrorExit(usage)
    elif opt=='--i': ifile=val
    elif opt=='--k': k=int(val)
    elif opt=='--p': p=float(val)
    elif opt=='--n': n=int(val)
    elif opt=='--v': verbose=1
    elif opt=='--vv': verbose=2
    else: ErrorExit('Illegal option: %s'%val)

  if ifile=='-' and k:
    ErrorExit('ERROR: k parameter disallowed with stdin.')

  if k:
    SampleKofN(ifile, k, n, verbose)

  elif p:
    SampleP(ifile, p, verbose)

  else:
    ErrorExit('ERROR: either --k or --p required.')

