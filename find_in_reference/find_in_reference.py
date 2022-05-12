#!/usr/bin/env python3


import os.path
import sys
import optparse


"""
#
#------------------------------------------------------------------------------
#                                                 University of Minnesota
#                 Copyright 2013, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#    James E Johnson
#
#------------------------------------------------------------------------------
"""

"""
Takes 2 tabular files as input:
    1. The file to be filtered
    2. The reference file

The string value of selected column of the input file is searched for
in the string values of the selected column of the reference file.

The intended purpose is to filter a peptide fasta file in tabular format
by whether those peptide sequences are found in a reference fasta file.

"""


def __main__():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest='input', help='The input file to filter. (Otherwise read from stdin)')
    parser.add_option('-r', '--reference', dest='reference', help='The reference file to filter against')
    parser.add_option('-o', '--output', dest='output', help='The output file for input lines filtered by reference')
    parser.add_option('-f', '--filtered', dest='filtered', help='The output file for input lines not in the output')
    parser.add_option('-c', '--input_column', dest='input_column', type="int", default=None, help='The column for the value in the input file. (first column = 1, default to last column)')
    parser.add_option('-C', '--reference_column', dest='reference_column', type="int", default=None, help='The column for the value in the reference file. (first column = 1, default to last column)')
    parser.add_option('-I', '--case_insensitive', dest='ignore_case', action="store_true", default=False, help='case insensitive')
    parser.add_option('-R', '--reverse_find', dest='reverse_find', action="store_true", default=False, help='find the reference string in the input string')
    parser.add_option('-B', '--test_reverse', dest='test_reverse', action="store_true", default=False, help='Also search for reversed input string in reference')
    parser.add_option('-D', '--test_dna_reverse_complement', dest='test_reverse_comp', action="store_true", default=False, help='Also search for the DNA reverse complement of input string')
    parser.add_option('-k', '--keep', dest='keep', action="store_true", default=False, help='')
    parser.add_option('-a', '--annotation_columns', dest='annotation_columns', default=None, help='If string is found, add these columns from reference')
    parser.add_option('-s', '--annotation_separator', dest='annotation_separator', default=';', help='separator character between annotations from different lines')
    parser.add_option('-S', '--annotation_col_sep', dest='annotation_col_sep', default=', ', help='separator character between annotation column from the same line')
    parser.add_option('-d', '--debug', dest='debug', action='store_true', default=False, help='Turn on wrapper debugging to stdout')
    (options, args) = parser.parse_args()

    # revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'N': 'N', 'n': 'n'}[B] for B in x][: : -1])

    COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'N': 'N', 'n': 'n'}

    def revcompl(seq):
        return ''.join([COMP[B] for B in seq][::-1])

    def test_rcomplement(seq, target):
        if options.test_reverse_comp:
            try:
                comp = revcompl(seq)
                return comp in target
            except Exception:
                pass
        return False

    def test_reverse(seq, target):
        return options.test_reverse and seq and seq[::-1] in target

    # Input files
    if options.input is not None:
        try:
            inputPath = os.path.abspath(options.input)
            inputFile = open(inputPath, 'r')
        except Exception as e:
            print("failed: %s" % e, file=sys.stderr)
            exit(2)
    else:
        inputFile = sys.stdin
    # Reference
    if options.reference is None:
        print("failed: reference file is required", file=sys.stderr)
        exit(2)
    # Output files
    outFile = None
    filteredFile = None
    if options.filtered is None and options.output is None:
        # write to stdout
        outFile = sys.stdout
    else:
        if options.output is not None:
            try:
                outPath = os.path.abspath(options.output)
                outFile = open(outPath, 'w')
            except Exception as e:
                print("failed: %s" % e, file=sys.stderr)
                exit(3)
        if options.filtered is not None:
            try:
                filteredPath = os.path.abspath(options.filtered)
                filteredFile = open(filteredPath, 'w')
            except Exception as e:
                print("failed: %s" % e, file=sys.stderr)
                exit(3)
    incol = -1
    if options.input_column and options.input_column > 0:
        incol = int(options.input_column)-1
    refcol = -1
    if options.reference_column and options.reference_column > 0:
        refcol = int(options.reference_column)-1
    if options.annotation_columns:
        annotate = True
        annotation_columns = [int(x) - 1 for x in options.annotation_columns.split(', ')]
    else:
        annotate = False
    refFile = None
    num_found = 0
    num_novel = 0
    for ln, line in enumerate(inputFile):
        annotations = []
        try:
            found = False
            search_string = line.split('\t')[incol].rstrip('\r\n')
            if options.ignore_case:
                search_string = search_string.upper()
            if options.debug:
                print("search: %s" % (search_string), file=sys.stderr)
            refFile = open(options.reference, 'r')
            for tn, fline in enumerate(refFile):
                fields = fline.split('\t')
                target_string = fields[refcol].rstrip('\r\n')
                if options.ignore_case:
                    target_string = target_string.upper()
                search = search_string if not options.reverse_find else target_string
                target = target_string if not options.reverse_find else search_string
                if options.debug:
                    print("in: %s %s %s" % (search, search in target, target), file=sys.stderr)
                if search in target or test_reverse(search, target) or test_rcomplement(search, target):
                    found = True
                    if annotate:
                        annotation = options.annotation_col_sep.join([fields[i] for i in annotation_columns])
                        annotations.append(annotation)
                    else:
                        break
            if found:
                num_found += 1
                if annotate:
                    line = '%s\t%s\n' % (line.rstrip('\r\n'), options.annotation_separator.join(annotations))
                if options.keep is True:
                    if outFile:
                        outFile.write(line)
                else:
                    if filteredFile:
                        filteredFile.write(line)
            else:
                num_novel += 1
                if options.keep is True:
                    if filteredFile:
                        filteredFile.write(line)
                else:
                    if outFile:
                        outFile.write(line)
        except Exception as e:
            print("failed: Error reading %s - %s" % (options.reference, e), file=sys.stderr)
        finally:
            if refFile:
                refFile.close()
    print("found: %d novel: %d" % (num_found, num_novel), file=sys.stdout)


if __name__ == "__main__":
    __main__()
