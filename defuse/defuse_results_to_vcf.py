#!/usr/bin/env python
"""
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#         Copyright 2012, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  James E Johnson
#  Jesse Erdmann
#
#------------------------------------------------------------------------------
"""


"""
This tool takes the defuse results.tsv  tab-delimited file as input and creates a Variant Call Format file as output.
"""

import sys,re,os.path
import optparse
from optparse import OptionParser

"""
http://www.1000genomes.org/wiki/analysis/variant-call-format/vcf-variant-call-format-version-42

5. INFO keys used for structural variants
When the INFO keys reserved for encoding structural variants are used for imprecise variants, the values should be best estimates. When a key reflects a property of a single alt allele (e.g. SVLEN), then when there are multiple alt alleles there will be multiple values for the key corresponding to each alelle (e.g. SVLEN=-100,-110 for a deletion with two distinct alt alleles).
The following INFO keys are reserved for encoding structural variants. In general, when these keys are used by imprecise variants, the values should be best estimates. When a key reflects a property of a single alt allele (e.g. SVLEN), then when there are multiple alt alleles there will be multiple values for the key corresponding to each alelle (e.g. SVLEN=-100,-110 for a deletion with two distinct alt alleles).
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
For precise variants, END is POS + length of REF allele - 1, and the for imprecise variants the corresponding best estimate.
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
Value should be one of DEL, INS, DUP, INV, CNV, BND. This key can be derived from the REF/ALT fields but is useful for filtering.
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
One value for each ALT allele. Longer ALT alleles (e.g. insertions) have positive values, shorter ALT alleles (e.g. deletions) have negative values.
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">
For precise variants, the consensus sequence the alternate allele assembly is derivable from the REF and ALT fields. However, the alternate allele assembly file may contain additional information about the characteristics of the alt allele contigs.
##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
##INFO=<ID=METRANS,Number=4,Type=String,Description="Mobile element transduction info of the form CHR,START,END,POLARITY">
##INFO=<ID=DGVID,Number=1,Type=String,Description="ID of this element in Database of Genomic Variation">
##INFO=<ID=DBVARID,Number=1,Type=String,Description="ID of this element in DBVAR">
##INFO=<ID=DBRIPID,Number=1,Type=String,Description="ID of this element in DBRIP">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around the length of the inserted material between breakends">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth of segment containing breakend">
##INFO=<ID=DPADJ,Number=.,Type=Integer,Description="Read Depth of adjacency">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number of segment containing breakend">
##INFO=<ID=CNADJ,Number=.,Type=Integer,Description="Copy number of adjacency">
##INFO=<ID=CICN,Number=2,Type=Integer,Description="Confidence interval around copy number for the segment">
##INFO=<ID=CICNADJ,Number=.,Type=Integer,Description="Confidence interval around copy number for the adjacency">
6. FORMAT keys used for structural variants
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description="Copy number genotype likelihood for imprecise events">
##FORMAT=<ID=NQ,Number=1,Type=Integer,Description="Phred style probability score that the variant is novel with respect to the genome's ancestor">
##FORMAT=<ID=HAP,Number=1,Type=Integer,Description="Unique haplotype identifier">
##FORMAT=<ID=AHAP,Number=1,Type=Integer,Description="Unique identifier of ancestral haplotype">
These keys are analogous to GT/GQ/GL and are provided for genotyping imprecise events by copy number (either because there is an unknown number of alternate alleles or because the haplotypes cannot be determined). CN specifies the integer copy number of the variant in this sample. CNQ is encoded as a phred quality -10log_10p(copy number genotype call is wrong). CNL specifies a list of log10 likelihoods for each potential copy number, starting from zero. When possible, GT/GQ/GL should be used instead of (or in addition to) these keys.

Specifying Complex Rearrangements with Breakends
An arbitrary rearrangement event can be summarized as a set of novel adjacencies.
Each adjacency ties together 2 breakends. The two breakends at either end of a novel adjacency are called mates.
There is one line of VCF (i.e. one record) for each of the two breakends in a novel adjacency. A breakend record is identified with the tag SYTYPE=BND" in the INFO field. The REF field of a breakend record indicates a base or sequence s of bases beginning at position POS, as in all VCF records. The ALT field of a breakend record indicates a replacement for s. This "breakend replacement" has three parts:
the string t that replaces places s. The string t may be an extended version of s if some novel bases are inserted during the formation of the novel adjacency.
The position p of the mate breakend, indicated by a string of the form "chr:pos". This is the location of the first mapped base in the piece being joined at this novel adjacency.
The direction that the joined sequence continues in, starting from p. This is indicated by the orientation of square brackets surrounding p.
These 3 elements are combined in 4 possible ways to create the ALT. In each of the 4 cases, the assertion is that s is replaced with t, and then some piece starting at position p is joined to t. The cases are:
REF   ALT    Meaning
s     t[p[   piece extending to the right of p is joined after t
s     t]p]   reverse comp piece extending left of p is joined after t
s     ]p]t   piece extending to the left of p is joined before t
s     [p[t   reverse comp piece extending right of p is joined before t

Examples:
#CHROM POS    ID     REF ALT           QUAL FILT INFO
2      321681 bnd_W  G   G]17:198982]  6    PASS SVTYPE=BND;MATEID=bnd_Y
2      321682 bnd_V  T   ]13:123456]T  6    PASS SVTYPE=BND;MATEID=bnd_U
13     123456 bnd_U  C   C[2:321682[   6    PASS SVTYPE=BND;MATEID=bnd_V
13     123457 bnd_X  A   [17:198983[A  6    PASS SVTYPE=BND;MATEID=bnd_Z
17     198982 bnd_Y  A   A]2:321681]   6    PASS SVTYPE=BND;MATEID=bnd_W
17     198983 bnd_Z  C   [13:123457[C  6    PASS SVTYPE=BND;MATEID=bnd_X
"""

vcf_header =  """\
##fileformat=VCFv4.1
##source=defuse
##reference=%s
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of the BND mate">
##INFO=<ID=MATELOC,Number=1,Type=String,Description="The chrom:position of the BND mate">
##INFO=<ID=GENESTRAND,Number=2,Type=String,Description="Strands">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth of segment containing breakend">
##INFO=<ID=SPLITCNT,Number=1,Type=Integer,Description="number of split reads supporting the prediction">
##INFO=<ID=SPANCNT,Number=1,Type=Integer,Description="number of spanning reads supporting the fusion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=SPLICESCORE,Number=1,Type=Integer,Description="number of nucleotides similar to GTAG at fusion splice">
##INFO=<ID=GENE,Number=2,Type=String,Description="Gene Names at each breakend">
##INFO=<ID=GENEID,Number=2,Type=String,Description="Gene IDs at each breakend">
##INFO=<ID=GENELOC,Number=2,Type=String,Description="location of breakpoint releative to genes">
##INFO=<ID=EXPR,Number=2,Type=Integer,Description="expression of genes as number of concordant pairs aligned to exons">
##INFO=<ID=ORF,Number=0,Type=Flag,Description="fusion combines genes in a way that preserves a reading frame">
##INFO=<ID=EXONBND,Number=0,Type=Flag,Description="fusion splice at exon boundaries">
##INFO=<ID=INTERCHROM,Number=0,Type=Flag,Description="fusion produced by an interchromosomal translocation">
##INFO=<ID=READTHROUGH,Number=0,Type=Flag,Description="fusion involving adjacent potentially resulting from co-transcription rather than genome rearrangement">
##INFO=<ID=ADJACENT,Number=0,Type=Flag,Description="fusion between adjacent genes">
##INFO=<ID=ALTSPLICE,Number=0,Type=Flag,Description="fusion likely the product of alternative splicing between adjacent genes">
##INFO=<ID=DELETION,Number=0,Type=Flag,Description="fusion produced by a genomic deletion">
##INFO=<ID=EVERSION,Number=0,Type=Flag,Description="fusion produced by a genomic eversion">
##INFO=<ID=INVERSION,Number=0,Type=Flag,Description="fusion produced by a genomic inversion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\
"""

def cmp_alphanumeric(s1,s2):
  if s1 == s2:
    return 0
  a1 = re.findall("\d+|[a-zA-Z]+",s1)
  a2 = re.findall("\d+|[a-zA-Z]+",s2)
  for i in range(min(len(a1),len(a2))):
    if a1[i] == a2[i]:
      continue
    if a1[i].isdigit() and a2[i].isdigit():
      return int(a1[i]) - int(a2[i])
    return 1 if a1[i] >  a2[i] else -1
  return len(a1) - len(a2)

def __main__():
  # VCF functions
  chr_dict = dict()
  def add_vcf_line(chr,pos,id,line):
    if chr not in chr_dict:
      pos_dict = dict()
      chr_dict[chr] = pos_dict
    if pos not in chr_dict[chr]:
      id_dict = dict()
      chr_dict[chr][pos] = id_dict
    chr_dict[chr][pos][id] = line

  def write_vcf():
    print >> outputFile, vcf_header % (refname)
    for chr in sorted(chr_dict.keys(),cmp=cmp_alphanumeric):
      for pos in sorted(chr_dict[chr].keys()):
        for id in chr_dict[chr][pos]:
          print >> outputFile, chr_dict[chr][pos][id]
  #Parse Command Line
  parser = optparse.OptionParser()
  # files
  parser.add_option( '-i', '--input', dest='input', help='The input defuse results.tsv file (else read from stdin)' )
  parser.add_option( '-o', '--output', dest='output', help='The output vcf file (else write to stdout)' )
  parser.add_option( '-r', '--reference', dest='reference', default=None, help='The genomic reference id' )
  (options, args) = parser.parse_args()

  # results.tsv input 
  if options.input != None:
    try:
      inputPath = os.path.abspath(options.input)
      inputFile = open(inputPath, 'r')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(2)
  else:
    inputFile = sys.stdin
  # vcf output 
  if options.output != None:
    try:
      outputPath = os.path.abspath(options.output)
      outputFile = open(outputPath, 'w')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(3)
  else:
    outputFile = sys.stdout

  refname = options.reference if options.reference else 'unknown'

  svtype = 'SVTYPE=BND'
  filt = 'PASS'
  columns = []
  try:
    for linenum,line in enumerate(inputFile):
      ## print >> sys.stderr, "%d: %s\n" % (linenum,line)
      fields = line.strip().split('\t')
      if line.startswith('cluster_id'):
        columns = fields
        ## print >> sys.stderr, "columns: %s\n" % columns
        continue
      cluster_id = fields[columns.index('cluster_id')]
      gene_chromosome1 = fields[columns.index('gene_chromosome1')]
      gene_chromosome2 = fields[columns.index('gene_chromosome2')]
      genomic_strand1 = fields[columns.index('genomic_strand1')]
      genomic_strand2 = fields[columns.index('genomic_strand2')]
      gene1 = fields[columns.index('gene1')]
      gene2 = fields[columns.index('gene2')]
      gene_info = 'GENEID=%s,%s' % (gene1,gene2)
      gene_name1 = fields[columns.index('gene_name1')]
      gene_name2 = fields[columns.index('gene_name2')]
      gene_name_info = 'GENE=%s,%s' % (gene_name1,gene_name2)
      gene_location1 = fields[columns.index('gene_location1')]
      gene_location2 = fields[columns.index('gene_location2')]
      gene_loc = 'GENELOC=%s,%s' % (gene_location1,gene_location2)
      expression1 = int(fields[columns.index('expression1')])
      expression2 = int(fields[columns.index('expression2')])
      expr = 'EXPR=%d,%d' % (expression1,expression2)
      genomic_break_pos1 = int(fields[columns.index('genomic_break_pos1')])
      genomic_break_pos2 = int(fields[columns.index('genomic_break_pos2')])
      breakpoint_homology = int(fields[columns.index('breakpoint_homology')])
      homlen = 'HOMLEN=%s' % breakpoint_homology
      orf = fields[columns.index('orf')] == 'Y'
      exonboundaries = fields[columns.index('exonboundaries')] == 'Y'
      read_through = fields[columns.index('read_through')] == 'Y'
      interchromosomal = fields[columns.index('interchromosomal')] == 'Y'
      adjacent = fields[columns.index('adjacent')] == 'Y'
      altsplice = fields[columns.index('altsplice')] == 'Y'
      deletion = fields[columns.index('deletion')] == 'Y'
      eversion = fields[columns.index('eversion')] == 'Y'
      inversion = fields[columns.index('inversion')] == 'Y'
      span_count = int(fields[columns.index('span_count')])
      splitr_count = int(fields[columns.index('splitr_count')])
      splice_score = int(fields[columns.index('splice_score')])
      probability = fields[columns.index('probability')] if columns.index('probability') else '.'
      splitr_sequence = fields[columns.index('splitr_sequence')]
      split_seqs = splitr_sequence.split('|')
      mate_id1 = "bnd_%s_1" % cluster_id
      mate_id2 = "bnd_%s_2" % cluster_id
      ref1 = split_seqs[0][-1]
      ref2 = split_seqs[1][0]
      b1 = '[' if genomic_strand1 == '+' else ']'
      b2 = '[' if genomic_strand2 == '+' else ']'
      alt1 = "%s%s%s:%d%s" %  (ref1,b2,gene_chromosome2,genomic_break_pos2,b2) 
      alt2 = "%s%s:%d%s%s" %  (b1,gene_chromosome1,genomic_break_pos1,b1,ref2) 
      #TODO evaluate what should be included in the INFO field
      info = ['DP=%d' % (span_count + splitr_count),'SPLITCNT=%d' % splitr_count,'SPANCNT=%d' % span_count,gene_name_info,gene_info,gene_loc,expr,homlen,'SPLICESCORE=%d' % splice_score]
      if orf:
        info.append('ORF')
      if exonboundaries:
        info.append('EXONBND')
      if interchromosomal:
        info.append('INTERCHROM')
      if read_through:
        info.append('READTHROUGH')
      if adjacent:
        info.append('ADJACENT')
      if altsplice:
        info.append('ALTSPLICE')
      if deletion:
        info.append('DELETION')
      if eversion:
        info.append('EVERSION')
      if inversion:
        info.append('INVERSION')
      info1 = [svtype,'MATEID=%s;MATELOC=%s:%d' % (mate_id2,gene_chromosome2,genomic_break_pos2)] + info
      info2 = [svtype,'MATEID=%s;MATELOC=%s:%d' % (mate_id1,gene_chromosome1,genomic_break_pos1)] + info
      qual = int(float(fields[columns.index('probability')]) * 255) if columns.index('probability') else '.'
      vcf1 = '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s'% (gene_chromosome1,genomic_break_pos1, mate_id1, ref1, alt1, qual, filt, ';'.join(info1) )
      vcf2 = '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s'% (gene_chromosome2,genomic_break_pos2, mate_id2, ref2, alt2, qual, filt, ';'.join(info2) )
      add_vcf_line(gene_chromosome1,genomic_break_pos1,mate_id1,vcf1)
      add_vcf_line(gene_chromosome2,genomic_break_pos2,mate_id2,vcf2)
    write_vcf()
  except Exception, e:
    print >> sys.stderr, "failed: %s" % e
    sys.exit(1)

if __name__ == "__main__" : __main__()

