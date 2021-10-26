#!/usr/bin/env python
from __future__ import print_function

import argparse
import multiprocessing as mp
import subprocess
from operator import itemgetter


gresults = []


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd',
                                       'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]


def get_count(pat, fastq):
    cmd = 'HITS=`grep -E \'' + pat + '\' ' + fastq + ' | wc -l` && echo $HITS'
    process = subprocess.Popen(['bash', '-c', cmd],
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    output = process.stdout.read().strip()
    return int(output) if output else 0


def get_seq_count(path, name, ln, seq, label, strand):
    # print('%s\t%d\t%s\t%s\t' % (path, ln, label, strand), file=sys.stderr)
    cmd = 'HITS=`grep -E \'' + seq + '\' ' + path + ' | wc -l` && echo $HITS'
    process = subprocess.Popen(['bash', '-c', cmd],
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    output = process.stdout.read().strip()
    n = int(output) if output else 0
    return (path, name, ln, seq, label, strand, n)


def process_seq_count(result):
    global results
    global csh
    gresults.append(result)
    """
    if csh:
        (path, name, ln, seq, label, strand, n) = result
        print('%s\t%d\t%s\t%s\t%d' % (path, ln, label, strand, n), file=cfh)
    """


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path_file', required=True,
                        help='file files paths and names')
    parser.add_argument('-i', '--query_file', required=True,
                        help='queries file')
    parser.add_argument('-s', '--summary', required=True,
                        help='summary report file')
    parser.add_argument('-n', '--counts', required=False,
                        help='counts')
    parser.add_argument('-I', '--id_col', required=False, default=None,
                        help='identifier column indices')
    parser.add_argument('-q', '--q_col', type=int, required=True,
                        help='column in queries for search sequence')
    parser.add_argument('-Q', '--q_label', required=False, default='mutant',
                        help='label for search sequence column')
    parser.add_argument('-c', '--c_col', type=int,
                        default=None, required=False,
                        help='column in queries for comparison sequence')
    parser.add_argument('-C', '--c_label', required=False, default='normal',
                        help='label for comparison sequence column')
    parser.add_argument('-r', '--reverse_complement',
                        action='store_true', default=False,
                        help='Also search for reverse complements')
    parser.add_argument('-P', '--per_file', action='store_true', default=False,
                        help='report per file')
    parser.add_argument('-T', '--threads', type=int, default=4, required=False,
                        help='threads')
    args = parser.parse_args()

    id_col_list = []
    if args.id_col:
        id_col_list = [int(x) for x in str(args.id_col).split(',')]
    strands = ['+', '-'] if args.reverse_complement else ['+']
    labels = [args.q_label]
    if args.c_col is not None:
        labels.append(args.c_label)
    ids = dict()
    qseqs = dict()
    cseqs = dict()
    files = dict()
    with open(args.query_file, 'r') as fh:
        for ln, line in enumerate(fh):
            qnum = ln + 1
            fields = str(line).rstrip().split('\t')
            id = '\t'.join([str(qnum)] +
                           [fields[i] for i in id_col_list])
            ids[ln] = id
            qseqs[ln] = fields[args.q_col]
            if args.c_col is not None:
                cseqs[ln] = fields[args.c_col]

    with open(str(args.path_file), 'r') as fh:
        for ln, line in enumerate(fh):
            path, name = str(line).rstrip().split('\t')
            files[name] = path

    queries = []
    for name, path in files.items():
        for i in range(len(qseqs)):
            ln = i + 1
            seq = qseqs[i]
            label = args.q_label
            strand = '+'
            queries.append((path, name, ln, seq, label, strand))
            if args.reverse_complement:
                strand = '-'
                queries.append((path, name, ln, revcomp(seq), label, strand))
            if i in cseqs:
                seq = cseqs[i]
                label = args.c_label
                strand = '+'
                queries.append((path, name, ln, seq, label, strand))
                if args.reverse_complement:
                    strand = '-'
                    queries.append((path, name, ln, revcomp(seq),
                                    label, strand))

    pool = mp.Pool(args.threads)
    for i, query in enumerate(queries):
        pool.apply_async(get_seq_count,
                         args=(query),
                         callback=process_seq_count)
    pool.close()
    pool.join()

    if args.counts:
        with open(args.counts, 'w') as cfh:
            for result in gresults:
                print('\t'.join([str(x) for x in result[1:]]), file=cfh)

    # count ln name label
    counts = dict()
    for i in range(len(ids)):
        counts[i] = dict()
        for name, path in files.items():
            counts[i][name] = dict()
            for j, label in enumerate(labels):
                counts[i][name][label] = dict()
                for strand in strands:
                    counts[i][name][label][strand] = 0

    results = sorted(gresults, key=itemgetter(2, 1, 4, 5))
    for i, result in enumerate(results):
        (path, name, ln, seq, label, strand, n) = result
        counts[ln-1][name][label][strand] = n

    with open(args.summary, 'w') as ofh:
        labels = [args.q_label]
        if args.c_col is not None:
            labels.append(args.c_label)
        for i in range(len(ids)):
            tcnts = [0, 0]
            print(ids[i], end='\t', file=ofh)
            for name, path in files.items():
                frac = 1.
                cnts = [0, 0]
                for j, label in enumerate(labels):
                    for strand in strands:
                        cnts[j] += counts[i][name][label][strand]
                    if args.per_file:
                        print(cnts[j], end='\t', file=ofh)
                if args.per_file:
                    cnt = sum(cnts)
                    frac = float(cnts[0])/cnt if cnt else 1.0
                    print(frac, end='\t', file=ofh)
                tcnts[0] += cnts[0]
                tcnts[1] += cnts[1]
            for k, label in enumerate(labels):
                print(tcnts[k], end='\t', file=ofh)
            tcnt = sum(tcnts)
            frac = float(tcnts[0])/tcnt if tcnt else 1.0
            print(frac, end='\n', file=ofh)


if __name__ == "__main__":
    main()
