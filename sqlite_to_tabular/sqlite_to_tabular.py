#!/usr/bin/env python

import sys
import os.path
import optparse
import sqlite3 as sqlite

def __main__():
  #Parse Command Line
  parser = optparse.OptionParser()
  parser.add_option( '-s', '--sqlitedb', dest='sqlitedb', default=None, help='The SQLite Database' )
  parser.add_option( '-q', '--query', dest='query', default=None, help='SQL query' )
  parser.add_option( '-Q', '--query_file', dest='query_file', default=None, help='SQL query file' )
  parser.add_option( '-n', '--no_header', dest='no_header', action='store_true', default=False, help='Include a column headers line' )
  parser.add_option( '-o', '--output', dest='output', default=None, help='Output file for query results' )
  (options, args) = parser.parse_args()

  # determine output destination
  if options.output != None:
    try:
      outputPath = os.path.abspath(options.output)
      outputFile = open(outputPath, 'w')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(3)
  else:
    outputFile = sys.stdout

  query = None
  if (options.query_file != None):
    with open(options.query_file,'r') as fh:
      query = ''
      for line in fh:
        query += line
  elif (options.query != None):
    query = options.query

  if (query is None):
    try:
      conn = sqlite.connect(options.sqlitedb)
      c = conn.cursor()
      tables_query = "SELECT name,sql FROM sqlite_master WHERE type='table' ORDER BY name"
      rslt = c.execute(tables_query).fetchall()
      for table,sql in rslt:
        print >> sys.stderr, "Table %s:" % table
        try:
          col_query = 'SELECT * FROM %s LIMIT 0' % table
          cur = conn.cursor().execute(col_query)
          cols = [col[0] for col in cur.description]
          print >> sys.stderr, " Columns: %s" % cols
        except Exception, exc:
          print >> sys.stderr, "Error: %s" % exc
    except Exception, exc:
      print >> sys.stderr, "Error: %s" % exc
    exit(0)
  #if not sqlite.is_read_only_query(query):
  #  print >> sys.stderr, "Error: Must be a read only query"
  #  exit(2)
  try:
    conn = sqlite.connect(options.sqlitedb)
    cur = conn.cursor()
    results = cur.execute(query)
    if not options.no_header:
      outputFile.write("#%s\n" % '\t'.join([str(col[0]) for col in cur.description]))
        # yield [col[0] for col in cur.description]
    for i,row in enumerate(results):
        # yield [val for val in row]
      outputFile.write("%s\n" % '\t'.join([str(val) for val in row]))
  except Exception, exc:
    print >> sys.stderr, "Error: %s" % exc
    exit(1)

if __name__ == "__main__": __main__()

