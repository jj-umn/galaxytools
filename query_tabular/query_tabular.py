#!/usr/bin/env python
"""
"""
import sys
import re
import os.path
import json
import sqlite3 as sqlite
import optparse
from optparse import OptionParser

"""
TODO:
- could read column names from comment lines, but issues with legal names
- could add some transformations on tabular columns,
  e.g. a regex to format date/time strings
    index: ['c2','c4,c5']
    unique: ['c1']
    format: {
      c2 : re.sub('pat', 'sub', c2)
      c3 : len(c3)
   }
   def format(colname,val, expr):
     
- allow optional autoincrement id column - user supplied name?
    autoincrement : 'id'
- column_defs dict of columns to create from tabular input
    column_defs : { 'name1' : 'expr', 'name2' : 'expr'}
- allow multiple queries and outputs
- add a --json input for table definitions (or yaml)
JSON config:
{ tables : [
    { file_path : '/home/galaxy/dataset_101.dat',
            table_name : 't1',
            column_names : ['c1', 'c2', 'c3'],
            comment_lines : 1
    },
    { file_path : '/home/galaxy/dataset_102.dat',
            table_name : 'gff',
            column_names : ['seqname',,,'start','end']
            comment_lines : 1
            load_named_columns : True
    },
    { file_path : '/home/galaxy/dataset_103.dat',
            table_name : 'test',
            column_names : ['c1', 'c2', 'c3']
    }
    ]
}
"""

tables_query = \
    "SELECT name, sql FROM sqlite_master WHERE type='table' ORDER BY name"


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


def get_column_def(file_path, table_name, skip=0, comment_char='#',
                   column_names=None, max_lines=100,load_named_columns=False):
    col_pref = ['TEXT', 'REAL', 'INTEGER', None]
    col_types = []
    col_idx = None
    data_lines = 0

    try:
        with open(file_path, "r") as fh:
            for linenum, line in enumerate(fh):
                if linenum < skip:
                    continue
                if line.startswith(comment_char):
                    continue
                data_lines += 1
                try:
                    fields = line.split('\t')
                    while len(col_types) < len(fields):
                        col_types.append(None)
                    for i, val in enumerate(fields):
                        colType = getValueType(val)
                        if col_pref.index(colType) < col_pref.index(col_types[i]):
                            col_types[i] = colType
                except Exception, e:
                    print >> sys.stderr, 'Failed at line: %d err: %s' % (linenum, e)
    except Exception, e:
        print >> sys.stderr, 'Failed: %s' % (e)
    for i,col_type in enumerate(col_types):
        if not col_type:
            col_types[i] = 'TEXT'
    if column_names: 
        col_names = []
        if load_named_columns:
            col_idx = []
            for i, cname in enumerate([cn.strip() for cn in column_names.split(',')]):
                if cname != '':
                    col_idx.append(i)
                    col_names.append(cname)                
            col_types = [col_types[i] for i in col_idx]
        else:
            col_names = ['c%d' % i for i in range(1, len(col_types) + 1)]
            for i, cname in enumerate([cn.strip() for cn in column_names.split(',')]):
                if cname and i < len(col_names):
                    col_names[i] = cname
    else:
        col_names = ['c%d' % i for i in range(1, len(col_types) + 1)]
    col_def = []
    for i, col_name in enumerate(col_names):
        col_def.append('%s %s' % (col_names[i], col_types[i]))
    return col_names, col_types, col_def, col_idx


def create_table(conn, file_path, table_name, skip=0, comment_char='#', column_names=None,load_named_columns=False,unique_indexes=[],indexes=[]):
    col_names, col_types, col_def, col_idx = get_column_def(file_path, table_name, skip=skip, comment_char=comment_char, column_names=column_names,load_named_columns=load_named_columns)
    col_func = [float if t == 'REAL' else int if t == 'INTEGER' else str for t in col_types]
    table_def = 'CREATE TABLE %s (\n    %s\n);' % (table_name, ', \n    '.join(col_def))
    # print >> sys.stdout, table_def
    insert_stmt = 'INSERT INTO %s(%s) VALUES(%s)' % (table_name, ','.join(col_names), ','.join(["?" for x in col_names]))
    # print >> sys.stdout, insert_stmt
    data_lines = 0
    try:
        c = conn.cursor()
        c.execute(table_def)
        conn.commit()
        c.close()
        for i,index in enumerate(unique_indexes):
            index_name='idx_uniq_%s_%d' % (table_name,i)
            index_columns = index.split(',')
            create_index(conn, table_name, index_name, index_columns,unique=True)
        for i,index in enumerate(indexes):
            index_name='idx_%s_%d' % (table_name,i)
            index_columns = index.split(',')
            create_index(conn, table_name, index_name, index_columns)
        c = conn.cursor()
        with open(file_path, "r") as fh:
            for linenum, line in enumerate(fh):
                if linenum < skip or line.startswith(comment_char):
                    continue
                data_lines += 1
                try:
                    fields = line.rstrip('\r\n').split('\t')
                    if col_idx:
                        fields = [fields[i] for i in col_idx]
                    vals = [col_func[i](x) if x else None for i, x in enumerate(fields)]
                    c.execute(insert_stmt, vals)
                except Exception, e:
                    print >> sys.stderr, 'Failed at line: %d err: %s' % (linenum, e)
        conn.commit()
        c.close()
    except Exception, e:
        print >> sys.stderr, 'Failed: %s' % (e)
        exit(1)

def create_index(conn, table_name, index_name, index_columns, unique=False):
    index_def = "CREATE %s INDEX %s on %s(%s)" % ('UNIQUE' if unique else '', index_name, table_name, ','.join(index_columns))
    c = conn.cursor()
    c.execute(index_def)
    conn.commit()
    c.close()

def regex_match(expr, item):
    return re.match(expr, item) is not None


def regex_search(expr, item):
    return re.search(expr, item) is not None


def regex_sub(expr, replace, item):
    return re.sub(expr, replace, item)


def get_connection(sqlitedb_path, addfunctions=False):
    conn = sqlite.connect(sqlitedb_path)
    if addfunctions:
        conn.create_function("re_match", 2, regex_match)
        conn.create_function("re_search", 2, regex_search)
        conn.create_function("re_sub", 3, regex_sub)
    return conn


def __main__():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('-s', '--sqlitedb', dest='sqlitedb', default=None, help='The SQLite Database')
    parser.add_option('-t', '--table', dest='tables', action="append", default=[], help='Tabular file: file_path[=table_name[:column_name, ...]')
    parser.add_option('-j', '--jsonfile', dest='jsonfile', default=None, help='Tabular file: file_path[=table_name[:column_name, ...]')
    parser.add_option('-q', '--query', dest='query', default=None, help='SQL query')
    parser.add_option('-Q', '--query_file', dest='query_file', default=None, help='SQL query file')
    parser.add_option('-n', '--no_header', dest='no_header', action='store_true', default=False, help='Include a column headers line')
    parser.add_option('-o', '--output', dest='output', default=None, help='Output file for query results')
    (options, args) = parser.parse_args()

    # open sqlite connection
    conn = get_connection(options.sqlitedb)
    # determine output destination
    if options.output is not None:
        try:
            outputPath = os.path.abspath(options.output)
            outputFile = open(outputPath, 'w')
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(3)
    else:
        outputFile = sys.stdout

    # get table defs
    if options.tables:
        for ti, table in enumerate(options.tables):
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
            create_table(conn, path, table_name, column_names=column_names)
    if options.jsonfile:
        try:
            fh = open(options.jsonfile)
            tdef = json.load(fh)
            if 'tables' in tdef:
                for ti, table in enumerate(tdef['tables']):
                    path = table['file_path']
                    table_name = table['table_name'] if 'table_name' in table else 't%d' % (ti + 1)
                    comment_lines = table['comment_lines'] if 'comment_lines' in table else 0
                    column_names = table['column_names'] if 'column_names' in table else None
                    if column_names:
                        load_named_columns = table['load_named_columns'] if 'load_named_columns' in table else False
                    else:   
                        load_named_columns = False
                    unique_indexes = table['unique'] if 'unique' in table else []
                    indexes = table['index'] if 'index' in table else []
                    create_table(conn, path, table_name, column_names=column_names, 
                                 skip=comment_lines, load_named_columns=load_named_columns, 
                                 unique_indexes=unique_indexes, indexes=indexes)
        except Exception, exc:
            print >> sys.stderr, "Error: %s" % exc
    conn.close()

    query = None
    if (options.query_file is not None):
        with open(options.query_file, 'r') as fh:
            query = ''
            for line in fh:
                query += line
    elif (options.query is not None):
        query = options.query

    if (query is None):
        try:
            conn = get_connection(options.sqlitedb)
            c = conn.cursor()
            rslt = c.execute(tables_query).fetchall()
            for table, sql in rslt:
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
    # if not sqlite.is_read_only_query(query):
    #    print >> sys.stderr, "Error: Must be a read only query"
    #    exit(2)
    try:
        conn = get_connection(options.sqlitedb, addfunctions=True)
        cur = conn.cursor()
        results = cur.execute(query)
        if not options.no_header:
            outputFile.write("#%s\n" % '\t'.join([str(col[0]) for col in cur.description]))
            # yield [col[0] for col in cur.description]
        for i, row in enumerate(results):
            # yield [val for val in row]
            outputFile.write("%s\n" % '\t'.join([str(val) for val in row]))
    except Exception, exc:
        print >> sys.stderr, "Error: %s" % exc
        exit(1)

if __name__ == "__main__":
    __main__()
