#!/usr/bin/env python
"""
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#         Copyright 2016, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  James E Johnson
#
#------------------------------------------------------------------------------
"""

"""
For each line in a file:
  match regex or col
  get value of match
  lookup repalcement for match from dict or SQLite
  replace and write
  
"""

import sys,re,os.path
import sqlite3 as sqlite
import tempfile
import optparse
from optparse import OptionParser

def __main__():
  #Parse Command Line
  parser = optparse.OptionParser()
  parser.add_option( '-i', '--input', dest='input', default=None, help='Tabular input file' )
  parser.add_option( '-o', '--output', dest='output', default=None, help='Tabular output file' )
  parser.add_option( '-c', '--column', type='int', dest='column', default=1, help='column ordinal in input to replace' )
  parser.add_option( '-p', '--pattern', dest='pattern', default=None, help='Regex search pattern, there must be one group for lookup key' )
  parser.add_option( '-m', '--match_group', type='int', dest='match_group', default=0, help='Regex match group for lookup key' )
  parser.add_option( '-r', '--replace', dest='replace', default=None, help='replace pattern, there must be one %s' )
  parser.add_option( '-s', '--sqlitedb', dest='sqlitedb', default=None, help='SQLite DB for lookup' )
  parser.add_option( '-q', '--query', dest='query', default=None, help='SQLite DB query for lookup' )
  parser.add_option( '-t', '--table', dest='table', default=None, help='Tabular lookup table' )
  parser.add_option( '-k', '--key_column', type='int', dest='key_column', default=1, help='column ordinal  for lookup key' )
  parser.add_option( '-v', '--value_column', type='int', dest='value_column', default=2, help='column ordinal  for lookup value' )
  parser.add_option( '-d', '--debug', dest='debug', action='store_true', default=False, help='Turn on wrapper debugging to stderr'  )
  (options, args) = parser.parse_args()
  # Input file
  if options.input != None:
    try:
      inputPath = os.path.abspath(options.input)
      inputFile = open(inputPath, 'r')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(2)
  else:
    inputFile = sys.stdin
  # Output file
  if options.output != None:
    try:
      outputPath = os.path.abspath(options.output)
      outputFile = open(outputPath, 'w')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(3)
  else:
    outputFile = sys.stdout
  lookup_dict = None
  if options.table:
    lookup_dict = dict()
    kcol = options.key_column - 1
    vcol = options.value_column - 1
    mcol = max(kcol,vcol)
    with open(options.table,'r') as fh:
      for i,line in enumerate(fh):
        fields = line.rstrip().split('\t')
        if len(fields) > mcol:
          lookup_dict[fields[kcol]] = fields[vcol]
  search_pat = None
  if options.pattern:
    search_pat = re.compile(options.pattern)
  replace = '%s'
  if options.replace:
    replace = options.replace
  conn = sqlite.connect(options.sqlitedb) if options.sqlitedb and options.query else None
  def lookup(k):
    if conn:
      query = options.query % k
      cur = conn.cursor()
      cur.execute(query)
      v = cur.fetchone()
      return v
    if lookup_dict:
      return lookup_dict[k] if k and k in lookup_dict else None
  try:
    for i, line in enumerate( inputFile ):
      if search_pat:
        m = search_pat.search(line)
        if m:
          if len(m.groups()) > 0:
            k = m.groups()[options.match_group]
            v = lookup(k)
            if v:
              repval = replace % v if replace and replace.find('%s') >= 0 else v
              line = search_pat.sub(repval,line)
      outputFile.write(line)
  except Exception, e:
    print >> sys.stderr, "failed: Error reading %s - %s" % (options.input if options.input else 'stdin',e)

if __name__ == "__main__" : __main__()

