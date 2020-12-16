#!/usr/bin/env python

import argparse
import json
import pandas as pd
import sys
from json.decoder import JSONDecodeError


def __main__():
    p = argparse.ArgumentParser()
    p.add_argument(
        '-i', '--input',
        type=argparse.FileType('r'),
        required=True,
        help='Tabular input file to pivot'
    )
    p.add_argument(
        '-o', '--output',
        type=argparse.FileType('w'),
        required=True,
        help='Output file'
    )
    p.add_argument(
        '-S', '--skiprows',
        type=int,
        default=0,
        help='Input column names'
    )
    p.add_argument(
        '-H', '--header',
        default=None,
        help='Input column names'
    )
    p.add_argument(
        '-P', '--prefix',
        default=None,
        help='Prefix for input column names'
    )
    p.add_argument(
        '-I', '--index',
        help='index columns'
    )
    p.add_argument(
        '-C', '--columns',
        help='columns values which are returned as columns'
    )
    p.add_argument(
        '-V', '--values',
        help='values'
    )
    p.add_argument(
        '-F', '--aggfunc',
        help='aggregate functions on the values'
    )
    p.add_argument(
        '-N', '--fill_value',
        default=None,
        help='fill value for missing values'
    )
    args = p.parse_args()

    def getValueType(val):
        if val or 0. == val:
            try:
                return int(val)
            except ValueError:
                try:
                    return float(val)
                except ValueError:
                    return val
        return None

    def getColumn(name, dfcols):
        if name in dfcols:
            return name
        else:
            try:
                i = int(name)
                return dfcols[i]
            except Exception:
                print('%s not a column in %s' % (name, dfcols),
                      file=sys.stderr)
                exit(1)

    def getColumns(val, dfcols):
        fields = [v.strip() for v in val.split(',')]
        cols = []
        for name in fields:
            cols.append(getColumn(name, dfcols))
        return cols

    def getAggFunc(funcStr, dfcols):
        af = funcStr
        try:
            af = json.loads(funcStr)
        except JSONDecodeError as de:
            print('"%s" is not a json string: ' % funcStr, de.msg,
                  file=sys.stderr)
            exit(1)
        if isinstance(af, dict):
            aggfunc = {getColumn(k, dfcols): v for k, v in af.items()}
        elif isinstance(af, list):
            aggfunc = af
        else:
            aggfunc = af
        return aggfunc

    if args.prefix:
        df = pd.read_table(args.input,
                           skiprows=args.skiprows,
                           header=None,
                           prefix=args.prefix)
    elif args.header:
        df = pd.read_table(args.input,
                           skiprows=args.skiprows,
                           header=args.header)
    else:
        df = pd.read_table(args.input, skiprows=args.skiprows)
    df_columns = df.columns.tolist()
    index = getColumns(args.index, df_columns)
    columns = getColumns(args.columns, df_columns)
    values = getColumns(args.values, df_columns)
    fill_value = getValueType(args.fill_value)
    aggfunc = getAggFunc(args.aggfunc, values)
    pdf = df.pivot_table(index=index, columns=columns,
                         values=values, aggfunc=aggfunc,
                         fill_value=fill_value)
    pdf_cols = ['_'.join(reversed(p)) if isinstance(p, tuple) else p
                for p in pdf.columns.tolist()]
    pdf.to_csv(args.output, sep='\t', float_format='%0.6f', header=pdf_cols)


if __name__ == "__main__":
    __main__()
