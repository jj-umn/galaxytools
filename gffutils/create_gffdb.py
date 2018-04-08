#!/usr/bin/env python

from __future__ import print_function

import argparse
import json
import sys

import gffutils


def __main__():
    parser = argparse.ArgumentParser(
        description='Convert mga output to bed and tsv')
    parser.add_argument(
        'input_gff',
        help="gff3 or gtf file to load")
    parser.add_argument(
        'gff_sqlite',
        help="sqlite file")
    args = parser.parse_args()

    db = gffutils.create_db(args.input_gff, args.gff_sqlite)
    version = db.execute('SELECT version FROM meta').fetchone()[0]
    print('version:%s' % version,file=sys.stdout)
    dialect = db.execute('SELECT dialect FROM meta').fetchone()[0]
    info = json.loads(dialect)
    if 'fmt' in info:
        print('format:%s' % info['fmt'],file=sys.stdout)
    seqids = [r[0] for r in db.execute('SELECT distinct seqid  FROM features').fetchall()]
    print('seqids:%s' % ','.join(seqids),file=sys.stdout)


if __name__ == "__main__":
    __main__()
