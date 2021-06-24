#!/usr/bin/env python3
###

import sys,os,argparse,logging
import numpy as np
import h5py

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='H5 file operations', epilog="")
  OPS = ['summary']
  parser.add_argument("op", choices=OPS, help='OPERATION')
  parser.add_argument("--i", dest="ifile", required=True, help="input file")
  parser.add_argument("--o", dest="ofile", help="output (TSV)")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  f = h5py.File(args.ifile, 'r')

  logging.debug(list(f.keys()))

  for key in list(f.keys()):
    dset = f[key]
    logging.debug(f"type(f[{key}]): {type(f[key])}")
    #logging.debug(f"Dataset {key} type: {dset.dtype}")
    #logging.debug(f"Dataset {key} shape: {dset.shape}")


