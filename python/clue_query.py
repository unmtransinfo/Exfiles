#!/usr/bin/env python3
#############################################################################
### Clue CMap L1000 API client.
### https://clue.io/api
### https://clue.io/connectopedia/query_api_tutorial
### l1000_type: "landmark"|"inferred"|"best inferred"|"not inferred"
#############################################################################
import sys,os,argparse,json
import urllib.parse

import rest_utils_py3 as rest_utils

API_HOST="api.clue.io"
API_BASE_URL="/api"
API_KEY="e404c56be1d299b6a0b4253808e88fdf" # user_key (JJYang)

N_CHUNK=100

#############################################################################
def ListDatatypes(verbose):
  headers={"Accept": "application/json", "user_key": API_KEY}
  url=('https://'+API_HOST+API_BASE_URL+'/dataTypes')
  response = rest_utils.PostURL(url, headers=headers, parse_json=True, verbose=verbose)
  print(json.dumps(response, indent=2), file=sys.stderr)

#############################################################################
def ListDatasets(verbose):
  headers={"Accept": "application/json", "user_key": API_KEY}
  url=('https://'+API_HOST+API_BASE_URL+'/datasets')
  response = rest_utils.PostURL(url, headers=headers, parse_json=True, verbose=verbose)
  print(json.dumps(response, indent=2), file=sys.stderr)

#############################################################################
def GetGenes(ids, id_type, fout, verbose):
  url_base=('https://'+API_HOST+API_BASE_URL+'/genes?user_key='+API_KEY)
  tags = None
  n_gene=0; i_chunk=0;
  for id_this in ids:
    while True:
      qry=('{"where":{"%s":"%s"},"skip":%d,"limit":%d}'%(
	id_type, urllib.parse.quote(id_this), i_chunk*N_CHUNK, N_CHUNK))
      url=url_base+('&filter=%s'%(qry))
      genes = rest_utils.GetURL(url, parse_json=True, verbose=verbose)
      if not genes:
        break
      for gene in genes:
        n_gene+=1
        if not tags:
          tags = gene.keys()
          fout.write('\t'.join(tags)+'\n')
        vals = []
        for tag in tags:
          if tag in gene:
            if type(gene[tag]) in (list, tuple):
              vals.append(';'.join([str(x) for x in gene[tag]]))
            else:
              vals.append(str(gene[tag]))
        fout.write('\t'.join(vals)+'\n')
      i_chunk+=1
  print('genes: %d'%n_gene, file=sys.stderr)

#############################################################################
def GetGenes_Landmark(fout, verbose):
  GetGenes(['landmark'], 'l1000_type', fout, verbose)

#############################################################################
def GetGenes_All(fout, verbose):
  GetGenes(['landmark', 'inferred', 'best inferred', 'not inferred'], 'l1000_type', fout, verbose)

#############################################################################
if __name__=="__main__":
  PROG=os.path.basename(sys.argv[0])
  parser = argparse.ArgumentParser(description='CLUE.IO REST API client utility')
  ops = ['getGenes_landmark', 'getGenes', 'getGenes_all', 'listDatasets', 'listDatatypes']
  id_types=['gene_symbol', 'entrez_id', 'l1000_type']
  parser.add_argument("op", choices=ops, help='operation')
  parser.add_argument("--id_type", choices=id_types, default='gene_symbol', help='gene ID type')
  parser.add_argument("--id", help="ID")
  parser.add_argument("--ifile", help="input file, IDs")
  parser.add_argument("--nmax", type=int, help="max results")
  parser.add_argument("--skip", type=int, help="skip results", default=0)
  parser.add_argument("--o", dest="ofile", help="output (CSV)")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0)

  args = parser.parse_args()

  if args.ofile:
    fout=open(args.ofile,"w+")
    if not fout: parser.error('ERROR: cannot open outfile: %s'%args.ofile)
  else:
    fout=sys.stdout

  ids=[];
  if args.ifile:
    fin=open(args.ifile)
    if not fin: parser.error('ERROR: cannot open: %s'%args.ifile)
    while True:
      line=fin.readline()
      if not line: break
      ids.append(line.rstrip())
    if args.verbose:
      print('%s: input queries: %d'%(PROG,len(ids)),file=sys.stderr)
    fin.close()
  elif args.id:
    ids.append(args.id)

  if args.op=='getGenes':
    if not ids: parser.error('--id or --ifile required.')
    GetGenes(ids, args.id_type, fout, args.verbose)

  elif args.op=='getGenes_landmark':
    GetGenes_Landmark(fout, args.verbose)

  elif args.op=='getGenes_all':
    GetGenes_All(fout, args.verbose)

  elif args.op=='listDatasets':
    ListDatasets(args.verbose)

  elif args.op=='listDatatypes':
    ListDatatypes(args.verbose)

  else:
    parser.print_help()
