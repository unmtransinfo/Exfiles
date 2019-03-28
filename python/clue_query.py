#!/usr/bin/env python3
#############################################################################
### Clue CMap API client.
### https://clue.io/api
### https://clue.io/connectopedia/query_api_tutorial
### https://clue.io/connectopedia/perturbagen_types_and_controls
### l1000_type: "landmark"|"inferred"|"best inferred"|"not inferred"
### API services: cells, genes, perts, pcls, plates, profiles, sigs,
### probeset_to_entrez_id, rep_fda_exclusivity
#############################################################################
### Cell app: https://clue.io/cell-app
#############################################################################
### curl -X GET --header "user_key: e404c56be1d299b6a0b4253808e88fdf" 'https://api.clue.io/api/cells?filter=\{"where":\{"provider_name":"ATCC"\}\}' |python3 -m json.tool
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
  response = rest_utils.GetURL(url, headers=headers, parse_json=True, verbose=verbose)
  print(json.dumps(response, indent=2), file=sys.stderr)

#############################################################################
def ListDatasets(verbose):
  headers={"Accept": "application/json", "user_key": API_KEY}
  url=('https://'+API_HOST+API_BASE_URL+'/datasets')
  response = rest_utils.GetURL(url, headers=headers, parse_json=True, verbose=verbose)
  print(json.dumps(response, indent=2), file=sys.stderr)

#############################################################################
def ListPerturbagenClasses(fout, verbose):
  headers={"user_key": API_KEY}
  url=('https://'+API_HOST+API_BASE_URL+'/pcls')
  pcls = rest_utils.GetURL(url, headers=headers, parse_json=True, verbose=verbose)
  tags=None; n_pcl=0;
  for pcl in pcls:
    n_pcl+=1
    if not tags:
      tags=pcl.keys()
      fout.write('\t'.join(tags)+'\n')
    vals=[];
    for tag in tags:
      if tag not in pcl:
        vals.append('')
      else:
        vals.append(str(pcl[tag]))
    fout.write('\t'.join(vals)+'\n')
  print('pcls: %d'%n_pcl, file=sys.stderr)

#############################################################################
def GetGenes(ids, id_type, fout, verbose):
  url_base=('https://'+API_HOST+API_BASE_URL+'/genes?user_key='+API_KEY)
  tags = None
  n_gene=0;
  for id_this in ids:
    n_gene_this=0;
    if verbose:
      print('id: %s'%id_this, file=sys.stderr)
    i_chunk=0;
    while True:
      qry=('{"where":{"%s":"%s"},"skip":%d,"limit":%d}'%(
	id_type, urllib.parse.quote(id_this), i_chunk*N_CHUNK, N_CHUNK))
      url=url_base+('&filter=%s'%(qry))
      try:
        genes = rest_utils.GetURL(url, parse_json=True, verbose=verbose)
      except:
        break
      if not genes:
        break
      for gene in genes:
        n_gene_this+=1
        if not tags:
          tags = gene.keys()
          fout.write('\t'.join(tags)+'\n')
        vals = []
        for tag in tags:
          if tag not in gene:
            vals.append('')
          elif type(gene[tag]) in (list, tuple):
            vals.append(';'.join([str(x) for x in gene[tag]]))
          else:
            vals.append(str(gene[tag]))
        fout.write('\t'.join(vals)+'\n')
      i_chunk+=1
    if verbose:
      print('\tgenes: %d'%n_gene_this, file=sys.stderr)
    n_gene+=n_gene_this
  print('genes: %d'%n_gene, file=sys.stderr)

#############################################################################
def GetGenes_Landmark(fout, verbose):
  GetGenes(['landmark'], 'l1000_type', fout, verbose)

#############################################################################
def GetGenes_All(fout, verbose):
  GetGenes(['landmark', 'inferred', 'best inferred', 'not inferred'], 'l1000_type', fout, verbose)

#############################################################################
### pert_type:
###	trt_cp - Compound
###	trt_sh - shRNA for loss of function (LoF) of gene
###	trt_lig - Peptides and other biological agents (e.g. cytokine)
###	trt_sh.cgs - Consensus signature from shRNAs targeting the same gene
#############################################################################
def GetPerturbagens(ids, id_type, fout, verbose):
  headers={"Accept": "application/json", "user_key": API_KEY}
  url_base=('https://'+API_HOST+API_BASE_URL+'/perts')
  tags = None
  fields = ['pert_id', 'pert_iname', 'pert_type', 'pert_vendor', 'pert_url', 'id', 'pubchem_cid', 'entrez_geneId', 'vector_id', 'clone_name', 'oligo_seq', 'description', 'target', 'structure_url', 'moa', 'pcl_membership', 'tas', 'num_sig', 'status']
  n_pert=0;
  for id_this in ids:
    n_pert_this=0;
    if verbose:
      print('id: %s'%id_this, file=sys.stderr)
    i_chunk=0;
    while True:
      qry=('{"where":{"%s":"%s"},"fields":[%s],"skip":%d,"limit":%d}'%(
	id_type, urllib.parse.quote(id_this),
	(','.join(['"%s"'%f for f in fields])),
	i_chunk*N_CHUNK, N_CHUNK))
      url=url_base+('?filter=%s'%(qry))
      try:
        perts = rest_utils.GetURL(url, headers=headers, parse_json=True, verbose=verbose)
      except:
        continue
      if not perts:
        break
      #print(json.dumps(perts, indent=2), file=sys.stderr) #DEBUG
      for pert in perts:
        n_pert_this+=1
        if not tags:
          tags = pert.keys()
          fout.write('\t'.join(tags)+'\n')
        vals = []
        for tag in tags:
          if tag not in pert:
            vals.append('')
          elif type(pert[tag]) in (list, tuple):
            vals.append(';'.join([str(x) for x in pert[tag]]))
          else:
            vals.append(str(pert[tag]))
        fout.write('\t'.join(vals)+'\n')
      i_chunk+=1
    if verbose:
      print('\tperturbagens: %d'%n_pert_this, file=sys.stderr)
    n_pert+=n_pert_this
  print('perturbagens: %d'%n_pert, file=sys.stderr)

#############################################################################
def GetPerturbagens_All(fout, verbose):
  pert_types=['trt_cp', 'trt_lig', 'trt_sh', 'trt_sh.cgs', 'trt_oe', 'trt_oe.mut', 'trt_xpr', 'trt_sh.css', 'ctl_vehicle.cns', 'ctl_vehicle', 'ctl_vector', 'ctl_vector.cns', 'ctl_untrt.cns', 'ctl_untrt']
  GetPerturbagens(pert_types, 'pert_type', fout, verbose)

#############################################################################
def GetCells(ids, id_type, fout, verbose):
  headers={"Accept": "application/json", "user_key": API_KEY}
  url_base=('https://'+API_HOST+API_BASE_URL+'/cells')
  tags = None
  n_cell=0;
  for id_this in ids:
    n_cell_this=0;
    if verbose:
      print('id: %s'%id_this, file=sys.stderr)
    i_chunk=0;
    while True:
      qry=('{"where":{"%s":"%s"},"skip":%d,"limit":%d}'%(
	id_type, urllib.parse.quote(id_this),
	i_chunk*N_CHUNK, N_CHUNK))
      url=url_base+('?filter=%s'%(qry))
      try:
        cells = rest_utils.GetURL(url, headers=headers, parse_json=True, verbose=verbose)
      except:
        break
      if not cells:
        break
      #print(json.dumps(cells, indent=2), file=sys.stderr) #DEBUG
      for cell in cells:
        n_cell_this+=1
        if not tags:
          tags = cell.keys()
          fout.write('\t'.join(tags)+'\n')
        vals = []
        for tag in tags:
          if tag not in cell:
            vals.append('')
          elif type(cell[tag]) in (list, tuple):
            vals.append(';'.join([str(x) for x in cell[tag]]))
          else:
            vals.append(str(cell[tag]))
        fout.write('\t'.join(vals)+'\n')
      i_chunk+=1
    if verbose:
      print('\tcells: %d'%n_cell_this, file=sys.stderr)
    n_cell+=n_cell_this
  print('cells: %d'%n_cell, file=sys.stderr)

#############################################################################
# 2570 (3/28/2019)
def GetCells_All(fout, verbose):
  headers={"Accept": "application/json", "user_key": API_KEY}
  url_base=('https://'+API_HOST+API_BASE_URL+'/cells')
  tags = None
  n_cell=0; i_chunk=0;
  while True:
    qry=('{"skip":%d,"limit":%d}'%(i_chunk*N_CHUNK, N_CHUNK))
    url=url_base+('?filter=%s'%(qry))
    try:
      cells = rest_utils.GetURL(url, headers=headers, parse_json=True, verbose=verbose)
    except Exception as e:
      print('Exception: %s'%e, file=sys.stderr)
      continue
    if not cells:
      break
    for cell in cells:
      n_cell+=1
      if not tags:
        tags = cell.keys()
        fout.write('\t'.join(tags)+'\n')
      vals = []
      for tag in tags:
        if tag not in cell:
          vals.append('')
        elif type(cell[tag]) in (list, tuple):
          vals.append(';'.join([str(x) for x in cell[tag]]))
        else:
          vals.append(str(cell[tag]))
      fout.write('\t'.join(vals)+'\n')
    i_chunk+=1
  print('cells: %d'%n_cell, file=sys.stderr)

#############################################################################
def GetSignatures(ids, id_type, fout, verbose):
  headers={"Accept": "application/json", "user_key": API_KEY}
  url_base=('https://'+API_HOST+API_BASE_URL+'/sigs')
  tags = None
  n_sig=0;
  for id_this in ids:
    n_sig_this=0;
    if verbose:
      print('id: %s'%id_this, file=sys.stderr)
    i_chunk=0;
    while True:
      qry=('{"where":{"%s":"%s"},"skip":%d,"limit":%d}'%(
	id_type, urllib.parse.quote(id_this),
	i_chunk*N_CHUNK, N_CHUNK))
      url=url_base+('?filter=%s'%(qry))
      try:
        sigs = rest_utils.GetURL(url, headers=headers, parse_json=True, verbose=verbose)
      except:
        continue
      if not sigs:
        break
      #print(json.dumps(sigs, indent=2), file=sys.stderr) #DEBUG
      for sig in sigs:
        n_sig_this+=1
        if not tags:
          tags = sig.keys()
          fout.write('\t'.join(tags)+'\n')
        vals = []
        for tag in tags:
          if tag not in sig:
            vals.append('')
          elif type(sig[tag]) in (list, tuple):
            vals.append(';'.join([str(x) for x in sig[tag]]))
          else:
            vals.append(str(sig[tag]))
        fout.write('\t'.join(vals)+'\n')
      i_chunk+=1
    if verbose:
      print('\tsigs: %d'%n_sig_this, file=sys.stderr)
    n_sig+=n_sig_this
  print('sigs: %d'%n_sig, file=sys.stderr)

#############################################################################
def GetThings(service, ids, id_type, fout, verbose):
  url_base=('https://'+API_HOST+API_BASE_URL+'/'+service+'?user_key='+API_KEY)
  tags = None
  n_thing=0;
  for id_this in ids:
    n_thing_this=0;
    if verbose:
      print('id: %s'%id_this, file=sys.stderr)
    i_chunk=0;
    while True:
      qry=('{"where":{"%s":"%s"},"skip":%d,"limit":%d}'%(
	id_type, urllib.parse.quote(id_this), i_chunk*N_CHUNK, N_CHUNK))
      url=url_base+('&filter=%s'%(qry))
      try:
        things = rest_utils.GetURL(url, parse_json=True, verbose=verbose)
      except:
        break
      if not things:
        break
      for thing in things:
        n_thing_this+=1
        if not tags:
          tags = thing.keys()
          fout.write('\t'.join(tags)+'\n')
        vals = []
        for tag in tags:
          if tag not in thing:
            vals.append('')
          elif type(thing[tag]) in (list, tuple):
            vals.append(';'.join([str(x) for x in thing[tag]]))
          else:
            vals.append(str(thing[tag]))
        fout.write('\t'.join(vals)+'\n')
      i_chunk+=1
    if verbose:
      print('\t%ss: %d'%(service, n_thing_this), file=sys.stderr)
    n_thing+=n_thing_this
  print('%ss: %d'%(service, n_thing), file=sys.stderr)

#############################################################################
if __name__=="__main__":
  PROG=os.path.basename(sys.argv[0])
  parser = argparse.ArgumentParser(description='CLUE.IO REST API client utility')
  ops = ['getGenes', 'getGenes_all', 'getGenes_landmark', 
	'getPerturbagens', 'getPerturbagens_all',
	'getSignatures',
	'getProfiles',
	'getCells', 'getCells_all',
	'listPerturbagenClasses',
	'listDatasets', 'listDatatypes']
  #id_types=['gene_symbol', 'entrez_id', 'l1000_type', 'pert_type', 'pert_id', 'pert_iname', 'pert_desc', 'cell_type', 'cell_iname']
  parser.add_argument("op", choices=ops, help='operation')
  parser.add_argument("--id_type", help='query ID or field type, e.g. gene_symbol')
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
  elif args.op=='getGenes_all':
    GetGenes_All(fout, args.verbose)
  elif args.op=='getGenes_landmark':
    GetGenes_Landmark(fout, args.verbose)

  elif args.op=='getPerturbagens':
    if not ids: parser.error('--id or --ifile required.')
    GetPerturbagens(ids, args.id_type, fout, args.verbose)
  elif args.op=='getPerturbagens_all':
    GetPerturbagens_All(fout, args.verbose)

  elif args.op=='getCells':
    if not ids: parser.error('--id or --ifile required.')
    GetCells(ids, args.id_type, fout, args.verbose)
  elif args.op=='getCells_all':
    GetCells_All(fout, args.verbose)

  elif args.op=='getSignatures':
    if not ids: parser.error('--id or --ifile required.')
    GetSignatures(ids, args.id_type, fout, args.verbose)

  elif args.op=='getProfiles':
    if not ids: parser.error('--id or --ifile required.')
    GetThings('profiles', ids, args.id_type, fout, args.verbose)

  elif args.op=='listDatasets':
    ListDatasets(args.verbose)

  elif args.op=='listDatatypes':
    ListDatatypes(args.verbose)

  elif args.op=='listPerturbagenClasses':
    ListPerturbagenClasses(fout, args.verbose)

  else:
    parser.print_help()
