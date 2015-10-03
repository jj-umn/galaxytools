#!/usr/bin/env python
"""
"""
import sys
import os.path
import sqlite3 as sqlite
import optparse
from optparse import OptionParser

"""
TODO: 
- could add some transformations on tabular columns, e.g. a regex to format date/time strings
- allow multiple queries and outputs
- add a --json input for table definitions (or yaml)
JSON config:
{ tables : [
    { file_path : '/home/galaxy/dataset_101.dat', 
      table_name : 't1',
      column_names : ['c1','c2','c3']
    },
    { file_path : '/home/galaxy/dataset_102.dat', 
      table_name : 't2',
      column_names : ['c1','c2','c3']
    }
  ]
}
"""

def getValueType(val):
  if val or 0. == val:
    try:
      int(val)
      return 'INTEGER'
    except:
      try:
        float(val)
        return 'REAL'
      except:
        return 'TEXT'
  return None
    

def get_column_def(file_path,table_name,skip=0,comment_char='#',column_names=None,max_lines=100):
  col_pref = ['TEXT','REAL','INTEGER',None]
  col_types = []
  data_lines = 0
  try:
    with open(file_path,"r") as fh:
      for linenum,line in enumerate(fh):
        if linenum < skip:
          continue
        if line.startswith(comment_char):
          continue
        data_lines += 1
        try:
          fields = line.split('\t')
          while len(col_types) < len(fields):
            col_types.append(None)
          for i,val in enumerate(fields):
            colType = getValueType(val)
            if col_pref.index(colType) < col_pref.index(col_types[i]):
              col_types[i] = colType
        except Exception, e:
          print >> sys.stderr, 'Failed at line: %d err: %s' % (linenum,e)
  except Exception, e:
    print >> sys.stderr, 'Failed: %s' % (e)
  for i,col_type in enumerate(col_types):
    if not col_type:
      col_types[i] = 'TEXT'
  col_names = ['c%d' % i for i in range(1,len(col_types) + 1)]
  if column_names:
    for i,cname in enumerate([cn.strip() for cn in column_names.split(',')]):
      if cname and i < len(col_names):
        col_names[i] = cname
  col_def = []
  for i,col_name in enumerate(col_names):
    col_def.append('%s %s' % (col_names[i],col_types[i]))
  return col_names,col_types,col_def
  
def create_table(conn,file_path,table_name,skip=0,comment_char='#',column_names=None):
  col_names,col_types,col_def = get_column_def(file_path,table_name,skip=skip,comment_char=comment_char,column_names=column_names)
  col_func = [float if t == 'REAL' else int if t == 'INTEGER' else str  for t in col_types]
  table_def = 'CREATE TABLE %s (\n  %s\n);' % (table_name,',\n  '.join(col_def))
  # print >> sys.stdout, table_def
  insert_stmt = 'INSERT INTO %s(%s) VALUES(%s)' % (table_name,','.join(col_names),','.join([ "?" for x in col_names]))
  # print >> sys.stdout, insert_stmt
  data_lines = 0
  try:
    c = conn.cursor()
    c.execute(table_def)
    with open(file_path,"r") as fh:
      for linenum,line in enumerate(fh):
        if linenum < skip or line.startswith(comment_char):
          continue
        data_lines += 1
        try:
          fields = line.split('\t')
          vals = [col_func[i](x) if x else None for i,x in enumerate(fields)]
          c.execute(insert_stmt,vals)
        except Exception, e:
          print >> sys.stderr, 'Failed at line: %d err: %s' % (linenum,e)
    conn.commit()
    c.close()
  except Exception, e:
    print >> sys.stderr, 'Failed: %s' % (e)
    exit(1)

def __main__():
  #Parse Command Line
  parser = optparse.OptionParser()
  parser.add_option( '-s', '--sqlitedb', dest='sqlitedb', default=None, help='The SQLite Database' )
  parser.add_option( '-t', '--table', dest='tables', action="append", default=[], help='Tabular file: file_path[=table_name[:column_name,...]' )
  parser.add_option( '-q', '--query', dest='query', default=None, help='SQL query' )
  parser.add_option( '-n', '--no_header', dest='no_header', action='store_true', default=False, help='Include a column headers line' )
  parser.add_option( '-o', '--output', dest='output', default=None, help='Output file for query results' )
  (options, args) = parser.parse_args()

  # oprn sqlite connection
  conn = sqlite.connect(options.sqlitedb)
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

  # determine output destination
  for ti,table in enumerate(options.tables):
    table_name = 't%d' % (ti + 1)
    column_names = None
    fields = table.split('=') 
    path = fields[0]
    if len(fields) > 1:
      names = fields[1].split(':')
      table_name = names[0] if names[0] else table_name
      if len(names) > 1:
        column_names = names[1]
    # print >> sys.stdout, '%s %s' % (table_name, path)
    create_table(conn,path,table_name,column_names=column_names)
  conn.close()

  if (options.query is None):
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
  #if not sqlite.is_read_only_query(options.query):
  #  print >> sys.stderr, "Error: Must be a read only query"
  #  exit(2)
  try:
    conn = sqlite.connect(options.sqlitedb)
    cur = conn.cursor()
    results = cur.execute(options.query)
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


