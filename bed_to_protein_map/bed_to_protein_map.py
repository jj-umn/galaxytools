#!/usr/bin/env python
"""
"""
import sys
import os.path
import optparse

"""
X	276352	291629	ENST00000430923	20	+	284187	291629	80,80,80	5	42,148,137,129,131	0,7814,12380,14295,15146
X	304749	318819	ENST00000326153	20	-	305073	318787	80,80,80	10	448,153,149,209,159,68,131,71,138,381	0,2610,2982,6669,8016,9400,10140,10479,12164,13689
grep ENST bed_to_protein_map.py | grep -v grep | ./bed_to_protein_map.py
"""

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    #I/O
    parser.add_option( '-i', '--input', dest='input', default=None, help='Tabular file with peptide_sequence column' )
    parser.add_option( '-c', '--column', type='int', dest='column', default=1, help='column ordinal with Ensembl transcript ID' )
    parser.add_option( '-g', '--gtf', dest='gtf', default=None, help='GTF gene model file. Used to annotate NSJ peptide entries.')
    parser.add_option( '-2', '--twobit', dest='twobit', default=None, help='Reference genome in UCSC twobit format')
    parser.add_option( '-C', '--coding', dest='coding', action='store_true', default=False, help='Output coding BED line')
    parser.add_option( '-o', '--output', dest='output', default=None, help='The output bed (else write to stdout)' )
    parser.add_option( '-s', '--sqlitedb', dest='sqlitedb', default=None, help='The sqlitedb' )
    parser.add_option( '--debug', dest='debug', action='store_true', default=False, help='Print debugging messages')
    (options, args) = parser.parse_args()
    ##INPUTS##
    if options.input != None:
        try:
            inputPath = os.path.abspath(options.input)
            inputFile = open(inputPath, 'r')
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(2)
    else:
        inputFile = sys.stdin

    if options.output != None:
        try:
            outputPath = os.path.abspath(options.output)
            outputFile = open(outputPath, 'w')
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(3)
    else:
        outputFile = sys.stdout

    dbFile = None
    if options.sqlitedb != None:
        try:
            dbPath = os.path.abspath(options.sqlitedb)
            dbFile = open(dbPath, 'w')
        except Exception, e:
            print >> sys.stderr, "failed: %s" % e
            exit(3)


    try:
        for linenum,line in enumerate(inputFile):
            if options.debug:
                print >> sys.stderr, "%d: %s\n" % (linenum,line)
            if line.startswith('#'):
                continue
            if line.strip() == '':
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 12:
                print >> sys.stderr, "%d: %s\n" % (linenum,line)
                continue
            (chrom,_chromStart,_chromEnd,name,score,strand,_thickStart,_thickEnd,itemRgb,_blockCount,blockSizes,blockStarts) = fields[0:12]
            chromStart = int(_chromStart)
            chromEnd = int(_chromEnd)
            thickStart = int(_thickStart)
            thickEnd = int(_thickEnd)
            blockCount = int(_blockCount)
            blockSizes = [int(x) for x in blockSizes.split(',')]
            blockStarts = [int(x) for x in blockStarts.split(',')]
            if strand == '+':
                cds_start = 0
                cds_end = 0
                for i in range(blockCount):
                    start = chromStart + blockStarts[i]
                    end = start + blockSizes[i]
                    if end < thickStart:
                        continue
                    if start > thickEnd:
                        break
                    if start < thickStart:  
                        start = thickStart
                    if end > thickEnd:
                        end = thickEnd
                    cds_end = cds_start + (end - start)
                    ##dbFile.write('%s\t%s\t%d\t%d\t%s\t%d\t%d\n' % (name,chrom,start,end,strand,cds_start,cds_end))
                    outputFile.write('%s\t%s\t%d\t%d\t%s\t%d\t%d\n' % (name,chrom,start,end,strand,cds_start,cds_end))
                    cds_start = cds_end
            elif strand == '-':
                cds_start = 0
                cds_end = 0
                for i in reversed(range(blockCount)):
                    start = chromStart + blockStarts[i]
                    end = start + blockSizes[i]
                    if end < thickStart:
                        break
                    if start > thickEnd:
                        continue
                    if start < thickStart:  
                        start = thickStart
                    if end > thickEnd:
                        end = thickEnd
                    cds_end = cds_start + (end - start)
                    outputFile.write('%s\t%s\t%d\t%d\t%s\t%d\t%d\n' % (name,chrom,start,end,strand,cds_start,cds_end))
                    cds_start = cds_end
                pass
    except Exception, e:
        print >> sys.stderr, "failed: %s" % e
        exit(1)


if __name__ == "__main__" : __main__()

