#!/usr/bin/env python
"""
"""
import sys
import os.path
import re
import optparse
import urllib
import urllib2
from optparse import OptionParser

mhci_methods = ['recommended','consensus','netmhcpan','ann','smmpmbec','smm','comblib_sidney2008','netmhccons','pickpocket']
mhcii_methods = ['recommended','consensus3','NetMHCIIpan','nn_align','smm_align','comblib','tepitope']
processing_methods = ['recommended','consensus','netmhcpan','ann','smmpmbec','smm','comblib_sidney2008']
mhcnp_methods = ['mhcnp']
bcell_methods = ['Bepipred','Chou-FasmanEmini','Karplus-Schulz','Kolaskar-Tongaonkar','Parker']
prediction_methods = {'mhci':mhci_methods,'mhcii':mhcii_methods,'processing':processing_methods,'mhcnp':mhcnp_methods,'bcell':bcell_methods}

def warn_err(msg,exit_code=1):
  sys.stderr.write(msg)
  if exit_code:
    sys.exit(exit_code)


def __main__():
  #Parse Command Line
  parser = optparse.OptionParser()
  parser.add_option( '-p', '--prediction', dest='prediction', default='mhci', choices=['mhci','mhcii','processing','mhcnp','bcell'], help='IEDB API prediction service' )
  parser.add_option( '-s', '--sequence', dest='sequence', action="append", default=None, help='Peptide Sequence' )
  parser.add_option( '-m', '--method', dest='method', default='recommended', choices=['recommended','consensus','netmhcpan','ann','smmpmbec','smm','comblib_sidney2008','netmhccons','pickpocket' ], help='prediction method' )
  parser.add_option( '-a', '--allele', dest='allele', action="append", default=[], help='Alleles for which to make predictions' )
  parser.add_option( '-l', '--length', dest='length', action="append", default=[], choices=['8', '9', '10', '11', '12', '13', '14', '15'], help='lengths for which to make predictions, 1 per allele' )
  parser.add_option( '-i', '--input', dest='input', default=None, help='Input file for peptide sequences (fasta or tabular)' )
  parser.add_option( '-c', '--column', dest='column', default=None, help='Peptide Column in a tabular input file' )
  parser.add_option( '-C', '--id_column', dest='id_column', default=None, help='ID Column in a tabular input file' )
  parser.add_option( '-o', '--output', dest='output', default=None, help='Output file for query results' )
  parser.add_option( '-d', '--debug', dest='debug', action='store_true', default=False, help='Turn on wrapper debugging to stderr'  )
  (options, args) = parser.parse_args()

  aapat = '^[ABCDEFGHIKLMNPQRSTVWY]+$'          

  if not options.allele and options.prediction != 'bcell':
    warn_err('-a allele required\n', exit_code=1)

  if not (options.sequence or options.input): 
    warn_err('NO Sequences given: either -s sequence or -i input_file is required\n', exit_code=1)

  if options.output != None:
    try:
      outputPath = os.path.abspath(options.output)
      outputFile = open(outputPath, 'w')
    except Exception, e:
      warn_err("Unable to open output file: %s\n" % e, exit_code=1)
  else:
    outputFile = sys.stdout

  url = 'http://tools-api.iedb.org/tools_api/%s/' % options.prediction

  #TODO parse alleles from the options.alleles file
  alleles = ','.join(options.allele)
  lengths = ','.join(options.length)
  method = options.method

  results = []
  global header
  header = None

  sequence_text = []
  def add_seq(seqid,seq):
    sequence_text.append(">%s\n%s" % (seqid if seqid else "peptide%d" % len(sequence_text),seq))

  def query(url,seq,allele,length,seqid=None,method='recommended'):
    global header
    params = dict()
    if method:
      params['method'] = method
    params['sequence_text'] = seq
    params['allele'] = allele
    params['length'] = length
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    if options.debug:
      print >> sys.stderr, "url %s %s %s" % (request.get_full_url(), seqid if seqid else "None", seq)
    response = None
    response = urllib2.urlopen(request)
    if response and response.getcode() == 200:
      resp_data = response.readlines()
      for line in resp_data:
        if line.find('eptide') > 0:
          header = "#%s%s" % ("ID\t" if seqid else "", line)
          continue
        if seqid:
          results.append("%s\t%s" % (seqid,line))
        else:
          results.append(line)
    elif not response:
      warn_err("NO response from IEDB server\n",  exit_code=3)
    else:
      warn_err("Error connecting to IEDB server\n",  exit_code=response.getcode())

  if options.sequence:
    for i,seq in enumerate(options.sequence):
      query(url,seq,alleles,lengths,seqid=None,method=method)
  if options.input:
    try:
      fh = open(options.input,'r')
      if options.column: ## tabular
        col = int(options.column)
        idcol = int(options.id_column) if options.id_column else None
        for i,line in enumerate(fh):
          fields = line.split('\t')
          if len(fields) > col:
            seq = re.sub('[_*]','',fields[col])
            if re.match(aapat,seq):
              seqid = fields[idcol] if idcol != None and idcol < len(fields) else None
              query(url,seq,alleles,lengths,seqid=seqid,method=method)
            else:
              warn_err('Line %d, Not a peptide: %s\n' % (i,seq),exit_code=None)
      else:  ## fasta
        seqid = None
        seq = ''
        for i,line in enumerate(fh):
          if line.startswith('>'):
            if seqid and len(seq) > 0:
              query(url,seq,alleles,lengths,seqid=seqid,method=method)
            seqid = line[1:].strip()
            seq = ''
          else:
            seq += line.strip()
        if seqid and len(seq) > 0:
          query(url,seq,alleles,lengths,seqid=seqid,method=method)
      fh.close()
    except Exception, e:
      warn_err("Unable to open input file: %s\n" % e, exit_code=1)

  if header:
    outputFile.write(header)  
  for line in results:
    outputFile.write(line)  

if __name__ == "__main__": __main__()

