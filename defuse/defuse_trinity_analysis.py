#!/usr/bin/env python
"""
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#         Copyright 2014, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  James E Johnson
#
#------------------------------------------------------------------------------
"""


"""
This tool takes the defuse results.tsv  tab-delimited file, trinity 
and creates a tabular report

Would it be possible to create 2 additional files from the deFuse-Trinity comparison program.  
One containing all the Trinity records matched to deFuse records (with the deFuse ID number), 
and the other with the ORFs records matching back to the Trinity records in the first files?

M045_Report.csv
"","deFuse_subset.count","deFuse.gene_name1","deFuse.gene_name2","deFuse.span_count","deFuse.probability","deFuse.gene_chromosome1","deFuse.gene_location1","deFuse.gene_chromosome2","deFuse.gene_location2","deFuse_subset.type"
"1",1,"Rps6","Dennd4c",7,0.814853504,"4","coding","4","coding","TIC  "



OS03_Matched_Rev.csv
"count","gene1","gene2","breakpoint","fusion","Trinity_transcript_ID","Trinity_transcript","ID1","protein"

"","deFuse.splitr_sequence","deFuse.gene_chromosome1","deFuse.gene_chromosome2","deFuse.gene_location1","deFuse.gene_location2","deFuse.gene_name1","deFuse.gene_name2","deFuse.span_count","deFuse.probability","word1","word2","fusion_part_1","fusion_part_2","fusion_point","fusion_point_rc","count","transcript"

"""

import sys,re,os.path,math
import textwrap
import optparse
from optparse import OptionParser

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a','N':'N','n':'n'}[B] for B in x][::-1])

codon_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

def translate(seq) :
  rna = seq.upper().replace('T','U')
  aa = []
  for i in range(0,len(rna) - 2, 3):
    codon = rna[i:i+3]
    aa.append(codon_map[codon] if codon in codon_map else 'X')
  return ''.join(aa)

def get_stop_codons(seq) :
  rna = seq.upper().replace('T','U')
  stop_codons = []
  for i in range(0,len(rna) - 2, 3):
    codon = rna[i:i+3]
    aa = codon_map[codon] if codon in codon_map else 'X'
    if aa == '*':
      stop_codons.append(codon)
  return stop_codons

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def test_rcomplement(seq, target):
  try:
    comp = revcompl(seq)
    return comp in target
  except:
    pass
  return False

def test_reverse(seq,target):
  return options.test_reverse and seq and seq[::-1] in target

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

def parse_defuse_results(inputFile): 
  defuse_results = []
  columns = []
  coltype_int = ['expression1', 'expression2', 'gene_start1', 'gene_start2', 'gene_end1', 'gene_end2', 'genomic_break_pos1', 'genomic_break_pos2', 'breakpoint_homology', 'span_count', 'splitr_count', 'splice_score']
  coltype_float = ['probability']
  coltype_yn = [ 'orf', 'exonboundaries', 'read_through', 'interchromosomal', 'adjacent', 'altsplice', 'deletion', 'eversion', 'inversion']
  try:
    for linenum,line in enumerate(inputFile):
      ## print >> sys.stderr, "%d: %s\n" % (linenum,line)
      fields = line.strip().split('\t')
      if line.startswith('cluster_id'):
        columns = fields
        ## print >> sys.stderr, "columns: %s\n" % columns
        continue
      elif fields and len(fields) == len(columns):
        cluster_id = fields[columns.index('cluster_id')]
        cluster = dict()
        flags = []
        defuse_results.append(cluster)
        for i,v in enumerate(columns):
          if v in coltype_int:
            cluster[v] = int(fields[i])
          elif v in coltype_float:
            cluster[v] = float(fields[i])
          elif v in coltype_yn:
            cluster[v] = fields[i] == 'Y'
            if cluster[v]:
              flags.append(columns[i]) 
          else:
            cluster[v] = fields[i]
        cluster['flags'] = ','.join(flags)
  except Exception, e:
    print >> sys.stderr, "failed to read cluster_dict: %s" % e
    exit(1)
  return defuse_results

## deFuse params to the mapping application?

def __main__():
  #Parse Command Line
  parser = optparse.OptionParser()
  # files
  parser.add_option( '-i', '--input', dest='input', default=None, help='The input defuse results.tsv file (else read from stdin)' )
  parser.add_option( '-t', '--transcripts', dest='transcripts', default=None, help='Trinity transcripts' )
  parser.add_option( '-p', '--peptides', dest='peptides', default=None, help='Trinity ORFs' )
  parser.add_option( '-o', '--output', dest='output', default=None, help='The output report (else write to stdout)' )
  parser.add_option( '-m', '--matched', dest='matched', default=None, help='The output matched report' )
  parser.add_option( '-a', '--transcript_alignment', dest='transcript_alignment', default=None, help='The output alignment file' )
  parser.add_option( '-A', '--orf_alignment', dest='orf_alignment', default=None, help='The output ORF alignment file' )
  parser.add_option( '-N', '--nbases', dest='nbases', type='int', default=12, help='Number of bases on either side of the fusion to compare' )
  parser.add_option( '-L', '--min_pep_len', dest='min_pep_len', type='int', default=100, help='Minimum length of peptide to report' )
  parser.add_option( '-T', '--ticdist', dest='ticdist', type='int', default=1000000, help='Maximum intrachromosomal distance to be classified a Transcription-induced chimera (TIC)' )
  parser.add_option( '-P', '--prior_aa', dest='prior_aa', type='int', default=11, help='Number of protein AAs to show preceeding fusion point' )
  parser.add_option( '-I', '--incomplete_orfs', dest='incomplete_orfs', action='store_true', default=False, help='Count incomplete ORFs'  )
  parser.add_option( '-O', '--orf_type', dest='orf_type', action='append', default=['complete','5prime_partial'], choices=['complete','5prime_partial','3prime_partial','internal'], help='ORF types to report'  )
  parser.add_option( '-r', '--readthrough', dest='readthrough', type='int', default=3, help='Number of stop_codons to read through' )
  # min_orf_len
  # split_na_len
  # tic_len = 1000000
  # prior
  # deFuse direction reversed 
  # in frame ?
  # contain known protein elements
  # what protein change
  # trinity provides full transctipt, defuse doesn't show full
  #parser.add_option( '-r', '--reference', dest='reference', default=None, help='The genomic reference fasta' )
  #parser.add_option( '-g', '--gtf', dest='gtf', default=None, help='The genomic reference gtf feature file')
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
  outputTxFile = None
  outputOrfFile = None
  if options.transcript_alignment: 
    try:
      outputTxFile = open(options.transcript_alignment,'w')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(3)
  if options.orf_alignment: 
    try:
      outputOrfFile = open(options.orf_alignment,'w')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(3)
  # Add percent match after transcript
  report_fields = ['gene_name1','gene_name2','span_count','probability','gene_chromosome1','gene_location1','gene_chromosome2','gene_location2','fusion_type','Transcript','coverage','Protein','flags','alignments1','alignments2']
  report_fields = ['cluster_id','gene_name1','gene_name2','span_count','probability','genomic_bkpt1','gene_location1','genomic_bkpt2','gene_location2','fusion_type','Transcript','coverage','Protein','flags','alignments1','alignments2']
  report_colnames = {'gene_name1':'Gene 1','gene_name2':'Gene 2','span_count':'Span cnt','probability':'Probability','gene_chromosome1':'From Chr','gene_location1':'Fusion point','gene_chromosome2':'To Chr','gene_location2':'Fusion point', 'cluster_id':'cluster_id', 'splitr_sequence':'splitr_sequence', 'splitr_count':'splitr_count', 'splitr_span_pvalue':'splitr_span_pvalue', 'splitr_pos_pvalue':'splitr_pos_pvalue', 'splitr_min_pvalue':'splitr_min_pvalue', 'adjacent':'adjacent', 'altsplice':'altsplice', 'break_adj_entropy1':'break_adj_entropy1', 'break_adj_entropy2':'break_adj_entropy2', 'break_adj_entropy_min':'break_adj_entropy_min', 'breakpoint_homology':'breakpoint_homology', 'breakseqs_estislands_percident':'breakseqs_estislands_percident', 'cdna_breakseqs_percident':'cdna_breakseqs_percident', 'deletion':'deletion', 'est_breakseqs_percident':'est_breakseqs_percident', 'eversion':'eversion', 'exonboundaries':'exonboundaries', 'expression1':'expression1', 'expression2':'expression2', 'gene1':'gene1', 'gene2':'gene2', 'gene_align_strand1':'gene_align_strand1', 'gene_align_strand2':'gene_align_strand2', 'gene_end1':'gene_end1', 'gene_end2':'gene_end2', 'gene_start1':'gene_start1', 'gene_start2':'gene_start2', 'gene_strand1':'gene_strand1', 'gene_strand2':'gene_strand2', 'genome_breakseqs_percident':'genome_breakseqs_percident', 'genomic_break_pos1':'genomic_break_pos1', 'genomic_break_pos2':'genomic_break_pos2', 'genomic_strand1':'genomic_strand1', 'genomic_strand2':'genomic_strand2', 'interchromosomal':'interchromosomal', 'interrupted_index1':'interrupted_index1', 'interrupted_index2':'interrupted_index2', 'inversion':'inversion', 'library_name':'library_name', 'max_map_count':'max_map_count', 'max_repeat_proportion':'max_repeat_proportion', 'mean_map_count':'mean_map_count', 'min_map_count':'min_map_count', 'num_multi_map':'num_multi_map', 'num_splice_variants':'num_splice_variants', 'orf':'orf', 'read_through':'read_through', 'repeat_proportion1':'repeat_proportion1', 'repeat_proportion2':'repeat_proportion2', 'span_coverage1':'span_coverage1', 'span_coverage2':'span_coverage2', 'span_coverage_max':'span_coverage_max', 'span_coverage_min':'span_coverage_min', 'splice_score':'splice_score', 'splicing_index1':'splicing_index1', 'splicing_index2':'splicing_index2', 'fusion_type':'Type', 'coverage':'fusion%','Transcript':'Transcript?','Protein':'Protein?','flags':'descriptions','fwd_seq':'fusion','alignments1':'alignments1','alignments2':'alignments2','genomic_bkpt1':'From Chr', 'genomic_bkpt2':'To Chr'}

  ## Read defuse results
  fusions = parse_defuse_results(inputFile)
  ## Create a field with the 12 nt before and after the fusion point. 
  ## Create a field with the reverse complement of the 24 nt fusion point field.
  ## Add fusion type filed (INTER, INTRA, TIC)
  for i,fusion in enumerate(fusions):
      fusion['ordinal'] = i + 1
      fusion['genomic_bkpt1'] = "%s:%d" % (fusion['gene_chromosome1'], fusion['genomic_break_pos1'])
      fusion['genomic_bkpt2'] = "%s:%d" % (fusion['gene_chromosome2'], fusion['genomic_break_pos2'])
      fusion['alignments1'] = "%s%s%s" % (fusion['genomic_strand1'], fusion['gene_strand1'], fusion['gene_align_strand1'])
      fusion['alignments2'] = "%s%s%s" % (fusion['genomic_strand2'], fusion['gene_strand2'], fusion['gene_align_strand2'])
      split_seqs = fusion['splitr_sequence'].split('|')
      fusion['split_seqs'] = split_seqs
      fusion['split_seqs'] = split_seqs
      fusion['split_seq_lens'] = [len(split_seqs[0]),len(split_seqs[1])]
      fusion['split_max_lens'] = [len(split_seqs[0]),len(split_seqs[1])]
      fwd_off = min(abs(options.nbases),len(split_seqs[0]))
      rev_off = min(abs(options.nbases),len(split_seqs[1]))
      fusion['fwd_off'] = fwd_off
      fusion['rev_off'] = rev_off
      fwd_seq = split_seqs[0][-fwd_off:] + split_seqs[1][:rev_off]
      rev_seq =  revcompl(fwd_seq)
      fusion['fwd_seq'] = fwd_seq
      fusion['rev_seq'] = rev_seq
      fusion_type = 'inter' if fusion['gene_chromosome1'] != fusion['gene_chromosome2'] else 'intra' if abs(fusion['genomic_break_pos1'] - fusion['genomic_break_pos2']) > options.ticdist else 'TIC'
      fusion['fusion_type'] = fusion_type
      fusion['transcripts'] = dict()
      fusion['Transcript'] = 'No'
      fusion['coverage'] = 0
      fusion['Protein'] = 'No'
      # print >> sys.stdout, "%4d\t%6s\t%s\t%s\t%s\t%s\t%s" % (i,fusion['cluster_id'],fwd_seq,rev_seq,fusion_type,fusion['gene_name1'],fusion['gene_name2'])
  inputFile.close()

  ## Process Trinity data and compare to deFuse
  matched_transcripts = dict()
  matched_orfs = dict()
  transcript_orfs = dict()
  fusions_with_transcripts = set()
  fusions_with_orfs = set()
  ## fusion['transcripts'][tx_id] { revcompl:?, bkpt:n, seq1: ,  seq2: , match1:n, match2:n}
  n = 0
  if options.transcripts: 
    with open(options.transcripts) as fp:
      for tx_full_id, seq in read_fasta(fp):
        n += 1
        for i,fusion in enumerate(fusions):
          if fusion['fwd_seq'] in seq or fusion['rev_seq'] in seq:
            fusions_with_transcripts.add(i)
            fusion['Transcript'] = 'Yes'
            tx_id = tx_full_id.lstrip('>').split()[0]
            matched_transcripts[tx_full_id] = seq
            fusion['transcripts'][tx_id] = dict()
            fusion['transcripts'][tx_id]['seq'] = seq
            fusion['transcripts'][tx_id]['full_id'] = tx_full_id
            pos = seq.find(fusion['fwd_seq'])
            if pos >= 0:
              tx_bkpt = pos + fusion['fwd_off']
              # fusion['transcripts'][tx_full_id] = tx_bkpt
              if tx_bkpt > fusion['split_max_lens'][0]:
                fusion['split_max_lens'][0] = tx_bkpt 
              len2 = len(seq) - tx_bkpt
              if len2 > fusion['split_max_lens'][1]:
                fusion['split_max_lens'][1] = len2 
              fusion['transcripts'][tx_id]['bkpt'] = tx_bkpt
              fusion['transcripts'][tx_id]['revcompl'] = False
              fusion['transcripts'][tx_id]['seq1'] = seq[:tx_bkpt]
              fusion['transcripts'][tx_id]['seq2'] = seq[tx_bkpt:]
            else: 
              pos = seq.find(fusion['rev_seq'])
              tx_bkpt = pos + fusion['rev_off']
              # fusion['transcripts'][tx_full_id] = -tx_bkpt
              if tx_bkpt > fusion['split_max_lens'][1]:
                fusion['split_max_lens'][1] = tx_bkpt
              len2 = len(seq) - tx_bkpt
              if len2 > fusion['split_max_lens'][0]:
                fusion['split_max_lens'][0] = len2 
              rseq = revcompl(seq)
              pos = rseq.find(fusion['fwd_seq'])
              tx_bkpt = pos + fusion['fwd_off']
              fusion['transcripts'][tx_id]['bkpt'] = tx_bkpt
              fusion['transcripts'][tx_id]['revcompl'] = True
              fusion['transcripts'][tx_id]['seq1'] = rseq[:tx_bkpt]
              fusion['transcripts'][tx_id]['seq2'] = rseq[tx_bkpt:]
            fseq = fusion['split_seqs'][0]
            tseq = fusion['transcripts'][tx_id]['seq1']
            mlen = min(len(fseq),len(tseq))
            fusion['transcripts'][tx_id]['match1'] = mlen
            for j in range(1,mlen+1):
              if fseq[-j] != tseq[-j]:
                fusion['transcripts'][tx_id]['match1'] = j - 1
                break
            fseq = fusion['split_seqs'][1]
            tseq = fusion['transcripts'][tx_id]['seq2']
            mlen = min(len(fseq),len(tseq))
            fusion['transcripts'][tx_id]['match2'] = mlen
            for j in range(mlen):
              if fseq[j] != tseq[j]:
                fusion['transcripts'][tx_id]['match2'] = j
                break
            # coverage = math.floor(float(fusion['transcripts'][tx_id]['match1'] + fusion['transcripts'][tx_id]['match2']) * 100. / len(fusion['split_seqs'][0]+fusion['split_seqs'][1]))
            coverage = int((fusion['transcripts'][tx_id]['match1'] + fusion['transcripts'][tx_id]['match2']) * 1000. / len(fusion['split_seqs'][0]+fusion['split_seqs'][1])) * .1
            # print >> sys.stderr, "%s\t%d\t%d\t%d\%s\t\t%d\t%d\t%d\t%d" % (tx_id,fusion['transcripts'][tx_id]['match1'],fusion['transcripts'][tx_id]['match2'],len(fusion['split_seqs'][0]+fusion['split_seqs'][1]),coverage,len( fusion['split_seqs'][0]),len(fusion['transcripts'][tx_id]['seq1']),len(fusion['split_seqs'][1]),len(fusion['transcripts'][tx_id]['seq2']))
            fusion['coverage'] = max(coverage,fusion['coverage'])
    print >> sys.stdout, "fusions_with_transcripts: %d  %s\n matched_transcripts: %d" % (len(fusions_with_transcripts),fusions_with_transcripts,len(matched_transcripts))
    ##for i,fusion in enumerate(fusions):
    ##  print >> sys.stdout, "%4d\t%6s\t%s\t%s\t%s\t%s\t%s\t%s" % (i,fusion['cluster_id'],fusion['fwd_seq'],fusion['rev_seq'],fusion['fusion_type'],fusion['gene_name1'],fusion['gene_name2'], fusion['transcripts'])
    ## Process ORFs and compare to matched deFuse and Trinity data.
    ## Proteins must be at least 100 aa long, starting at the first "M" and must end with an "*".
    if options.peptides: 
      with open(options.peptides) as fp:
        for orf_full_id, seq in read_fasta(fp):
          n += 1
          if len(seq) < options.min_pep_len:
            continue
          orf_type = re.match('^.* type:(\S+) .*$',orf_full_id).groups()[0]
          ## if not seq[-1] == '*' and not options.incomplete_orfs:
          ## if not orf_type 'complete' and not options.incomplete_orfs:
          if orf_type not in options.orf_type:
            continue
          for i,fusion in enumerate(fusions):
            if len(fusion['transcripts']) > 0:
              for tx_id in fusion['transcripts']:
                ## >m.196252 g.196252  ORF g.196252 m.196252 type:complete len:237 (+) comp100000_c5_seq2:315-1025(+)
                ## >m.134565 g.134565  ORF g.134565 m.134565 type:5prime_partial len:126 (-) comp98702_c1_seq21:52-429(-)
                if tx_id+':' not in orf_full_id:
                  continue
                m = re.match("^.*%s:(\d+)-(\d+)[(]([+-])[)].*" % re.sub('([|.{}()$?^])','[\\1]',tx_id),orf_full_id)
                if m:
                  if not m.groups() or len(m.groups()) < 3 or m.groups()[0] == None:
                    print >> sys.stderr, "Error:\n%s\n%s\n" % (tx_id,orf_full_id)
                  orf_id = orf_full_id.lstrip('>').split()[0]
                  if not tx_id in transcript_orfs:
                    transcript_orfs[tx_id] = []
                  alignments = "%s%s%s %s%s%s" % (fusion['genomic_strand1'], fusion['gene_strand1'], fusion['gene_align_strand1'], fusion['genomic_strand2'], fusion['gene_strand2'], fusion['gene_align_strand2'])
                  # print >> sys.stdout, "%d %s bkpt:%d %s rc:%s (%s)   %s" % (fusion['ordinal'], tx_id, int(fusion['transcripts'][tx_id]['bkpt']), str(m.groups()), str(fusion['transcripts'][tx_id]['revcompl']), alignments, orf_full_id) 
                  start = seq.find('M')
                  pep_len = len(seq)
                  if pep_len - start < options.min_pep_len:
                    continue
                  orf_dict = dict()
                  transcript_orfs[tx_id].append(orf_dict)
                  fusions_with_orfs.add(i)
                  matched_orfs[orf_full_id] = seq
                  fusion['Protein'] = 'Yes'
                  tx_start = int(m.groups()[0])
                  tx_end = int(m.groups()[1])
                  tx_strand = m.groups()[2]
                  tx_bkpt = fusion['transcripts'][tx_id]['bkpt']
                  orf_dict['orf_id'] = orf_id
                  orf_dict['tx_start'] = tx_start
                  orf_dict['tx_end'] = tx_end
                  orf_dict['tx_strand'] = tx_strand
                  orf_dict['tx_bkpt'] = tx_bkpt
                  orf_dict['seq'] = seq[:start].lower() + seq[start:] if start > 0 else seq
                  ## >m.208656 g.208656  ORF g.208656 m.208656 type:5prime_partial len:303 (+) comp100185_c2_seq9:2-910(+)
                  ## translate(tx34[1:910])
                  ## translate(tx34[1:2048])
                  ## comp99273_c1_seq1 len=3146 (-2772) 
                  ## >m.158338 g.158338  ORF g.158338 m.158338 type:complete len:785 (-) comp99273_c1_seq1:404-2758(-)
                  ##  translate(tx[-2758:-403])
                  ## comp100185_c2_seq9 len=2048 (904)
                  ## novel protein sequence
                  ## find first novel AA
                  ## get prior n AAs
                  ## get novel AA seq thru n stop codons 
                  ### tx_seq = matched_transcripts[tx_full_id] if tx_bkpt >= 0 else revcompl(tx_seq)
                  tx_seq = fusion['transcripts'][tx_id]['seq']
                  orf_dict['tx_seq'] = tx_seq
                  novel_tx_seq = tx_seq[tx_start - 1:] if tx_strand == '+' else revcompl(tx_seq[:tx_end])
                  read_thru_pep = translate(novel_tx_seq)
                  # fusion['transcripts'][tx_id]['revcompl'] = True
                  # tx_bkpt = fusion['transcripts'][tx_id]['bkpt'] 
                  # bkpt_aa_pos = tx_bkpt - tx_start - 1 
                  # bkpt_aa_pos = (tx_bkpt - tx_start - 1) / 3 if tx_strand == '+' else tx_end
                  # print >> sys.stdout, "%s\n%s" % (seq,read_thru_pep) 
                  stop_codons = get_stop_codons(novel_tx_seq)
                  if options.readthrough: 
                    readthrough = options.readthrough + 1
                    read_thru_pep = '*'.join(read_thru_pep.split('*')[:readthrough])
                    stop_codons = stop_codons[:readthrough]
                  orf_dict['read_thru_pep'] = read_thru_pep
                  orf_dict['stop_codons'] = ','.join(stop_codons)
      print >> sys.stdout, "fusions_with_orfs: %d  %s\n matched_orfs: %d" % (len(fusions_with_orfs),fusions_with_orfs,len(matched_orfs))
  ## Alignments 3 columns, seq columns padded out to longest seq, UPPERCASE_match  diffs lowercase
  ### defuse_id		pre_split_seq		post_split_seq
  ### trinity_id	pre_split_seq		post_split_seq
  ## Transcripts alignment output
  ## Peptide alignment output
  ## Write reports
  ## OS03_Matched_Rev.csv
  ## "count","gene1","gene2","breakpoint","fusion","Trinity_transcript_ID","Trinity_transcript","ID1","protein"
  if options.transcripts and options.matched: 
    #match_fields = ['ordinal','gene_name1','gene_name2','fwd_seq']
    outputMatchFile = open(options.matched,'w')
    #print >> outputMatchFile, '\t'.join(["#fusion_id","cluster_id","gene1","gene2","breakpoint","fusion","Trinity_transcript_ID","Trinity_transcript","Trinity_ORF_Transcript","Trinity_ORF_ID","protein","read_through","stop_codons"])
    print >> outputMatchFile, '\t'.join(["#fusion_id","cluster_id","gene1","gene2","breakpoint","fusion","Trinity_transcript_ID","Trinity_transcript","Trinity_ORF_Transcript","Trinity_ORF_ID","protein","stop_codons"])
    for i,fusion in enumerate(fusions):
      if len(fusion['transcripts']) > 0:
        for tx_id in fusion['transcripts'].keys():
          if tx_id in transcript_orfs:
            for orf_dict in transcript_orfs[tx_id]: 
              if 'tx_seq' not in orf_dict:
                print >> sys.stderr, "orf_dict %s" % orf_dict
              #fields = [str(fusion['ordinal']),str(fusion['cluster_id']),fusion['gene_name1'],fusion['gene_name2'],fusion['fwd_seq'],fusion['splitr_sequence'],tx_id, fusion['transcripts'][tx_id]['seq1']+'|'+fusion['transcripts'][tx_id]['seq2'],orf_dict['tx_seq'],orf_dict['orf_id'],orf_dict['seq'],orf_dict['read_thru_pep'],orf_dict['stop_codons']]
              fields = [str(fusion['ordinal']),str(fusion['cluster_id']),fusion['gene_name1'],fusion['gene_name2'],fusion['fwd_seq'],fusion['splitr_sequence'],tx_id, fusion['transcripts'][tx_id]['seq1']+'|'+fusion['transcripts'][tx_id]['seq2'],orf_dict['tx_seq'],orf_dict['orf_id'],orf_dict['read_thru_pep'],orf_dict['stop_codons']]
              print >> outputMatchFile, '\t'.join(fields)
    outputMatchFile.close()
  if options.transcripts and options.transcript_alignment: 
    if outputTxFile:
      id_fields = ['gene_name1','alignments1','gene_name2','alignments2','span_count','probability','gene_chromosome1','gene_location1','gene_chromosome2','gene_location2','fusion_type','Transcript','Protein','flags']
      fa_width = 80
      for i,fusion in enumerate(fusions):
        if len(fusion['transcripts']) > 0:
          alignments1 = "%s%s%s" % (fusion['genomic_strand1'], fusion['gene_strand1'], fusion['gene_align_strand1'])
          alignments2 = "%s%s%s" % (fusion['genomic_strand2'], fusion['gene_strand2'], fusion['gene_align_strand2'])
          alignments = "%s%s%s %s%s%s" % (fusion['genomic_strand1'], fusion['gene_strand1'], fusion['gene_align_strand1'], fusion['genomic_strand2'], fusion['gene_strand2'], fusion['gene_align_strand2'])
          fusion_id = "%s (%s) %s" % (i + 1,alignments,' '.join([str(fusion[x]) for x in report_fields]))
          for tx_id in fusion['transcripts'].keys():
            m1 = fusion['transcripts'][tx_id]['match1']
            f_seq1 = fusion['split_seqs'][0][:-m1].lower() +  fusion['split_seqs'][0][-m1:]
            t_seq1 = fusion['transcripts'][tx_id]['seq1'][:-m1].lower() + fusion['transcripts'][tx_id]['seq1'][-m1:]
            if len(f_seq1) > len(t_seq1):
              t_seq1 = t_seq1.rjust(len(f_seq1),'.')
            elif len(f_seq1) < len(t_seq1):
              f_seq1 = f_seq1.rjust(len(t_seq1),'.')
            m2 = fusion['transcripts'][tx_id]['match2']
            f_seq2 = fusion['split_seqs'][1][:m2] +  fusion['split_seqs'][1][m2:].lower()
            t_seq2 = fusion['transcripts'][tx_id]['seq2'][:m2] + fusion['transcripts'][tx_id]['seq2'][m2:].lower()
            if len(f_seq2) > len(t_seq2):
              t_seq2 = t_seq2.ljust(len(f_seq2),'.')
            elif len(f_seq2) < len(t_seq2):
              f_seq2 = f_seq2.ljust(len(t_seq2),'.')
            print >> outputTxFile, ">%s\n%s\n%s" % (fusion_id,'\n'.join(textwrap.wrap(f_seq1,fa_width)),'\n'.join(textwrap.wrap(f_seq2,fa_width)))
            print >> outputTxFile, "%s bkpt:%d rev_compl:%s\n%s\n%s" % (fusion['transcripts'][tx_id]['full_id'],fusion['transcripts'][tx_id]['bkpt'],str(fusion['transcripts'][tx_id]['revcompl']),'\n'.join(textwrap.wrap(t_seq1,fa_width)),'\n'.join(textwrap.wrap(t_seq2,fa_width)))
  """
  if options.peptides and options.orf_alignment: 
    pass
  """
  print >> outputFile,"%s\t%s" % ('#','\t'.join([report_colnames[x] for x in report_fields]))
  for i,fusion in enumerate(fusions): 
    print >> outputFile,"%s\t%s" % (i + 1,'\t'.join([str(fusion[x]) for x in report_fields]))

if __name__ == "__main__" : __main__()

