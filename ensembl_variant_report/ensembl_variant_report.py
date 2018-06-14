#!/usr/bin/env python
"""
Report variants to Ensembl Transcripts

FrameShift report line
# Gene	Variant Position	Reference	Variation	Prevalence	Sequencing Depth	Transcript	AA Position	AA change	AA Length	Stop Codon	Stop Region	AA Variation
TIGD6	chr5:149374879	AAG	AG	1.00	13	ENSG00000164296|ENST00000296736	345	Q344	522	A-TGA-T	KRWTSSRPST*

MissSense report line
# Gene	Variant Position	Reference	Variation	Prevalence	Sequencing Depth	Transcript	AA Position	AA change	AA Length	Stop Codon	Stop Region	AA Variation
FN1	chr2:216235089	G	A	1.00	7394	ENSG00000115414|ENST00000354785	2261	V2261I	2478	G-TAA	TGLTRGATYN_I_IVEALKDQQR
"""
import sys
import os.path
import re
import time
import optparse
from ensemblref import EnsemblRef
from Bio.Seq import reverse_complement, translate


def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    #I/O
    parser.add_option( '-i', '--input', dest='input', default=None, help='Tabular file with peptide_sequence column' )
    parser.add_option( '-s', '--format', dest='format', default='tabular', choices=['tabular','snpeff'], help='Tabular file with peptide_sequence column' )
    #Columns for tabular input
    parser.add_option( '-C', '--chrom_column', type='int', dest='chrom_column', default=1, help='column ordinal with Ensembl transctip ID' )
    parser.add_option( '-P', '--pos_column', type='int', dest='pos_column', default=2, help='column ordinal with Ensembl transctip ID' )
    parser.add_option( '-R', '--ref_column', type='int', dest='ref_column', default=3, help='column ordinal with Ensembl transctip ID' )
    parser.add_option( '-A', '--alt_column', type='int', dest='alt_column', default=4, help='column ordinal with Ensembl transctip ID' )
    parser.add_option( '-T', '--transcript_column', type='int', dest='transcript_column', default=1, help='column ordinal with Ensembl transctip ID' )
    parser.add_option( '-F', '--dpr_column', type='int', dest='dpr_column', default=1, help='column with VCF: DPR or AD' )
    parser.add_option( '-D', '--dp_column', type='int', dest='dp_column', default=1, help='column with VCF: DP' )
    parser.add_option( '-g', '--gene_model', dest='gene_model', default=None, help='GTF gene model file. Used to annotate NSJ peptide entries.')
    parser.add_option( '-2', '--twobit', dest='twobit', default=None, help='Reference genome in UCSC twobit format')
    #Output file
    parser.add_option( '-o', '--output', dest='output', default=None, help='The output report (else write to stdout)' )
    #filters
    parser.add_option( '-d', '--min_depth', type='int', dest='min_depth', default=None, help='Minimum read depth to report' )
    parser.add_option( '-f', '--min_freq', type='float', dest='min_freq', default=None, help='Minimum variant frequency to report' )
    #peptide options
    parser.add_option( '-l', '--leading_aa', type='int', dest='leading_aa', default=10, help='Number AAs before missense variant' )
    parser.add_option( '-t', '--trailing_aa', type='int', dest='trailing_aa', default=10, help='Number AAs after missense variant' )
    parser.add_option( '-r', '--readthrough', type='int', dest='readthrough', default=0, help='' )
    # 
    parser.add_option('--debug', dest='debug', action='store_true', default=False, help='Print debugging messages')
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

    def parse_tabular():
        ci = options.chrom_column - 1
        pi = options.pos_column - 1
        ri = options.ref_column - 1
        ai = options.alt_column - 1
        ti = options.transcript_column - 1
        di = options.dp_column - 1
        fi = options.dpr_column - 1
        for linenum,line in enumerate(inputFile):
            if options.debug:
                print >> sys.stderr, "%d: %s\n" % (linenum,line)
            if line.startswith('#'):
                continue
            if line.strip() == '':
                continue
            fields = line.rstrip('\r\n').split('\t')
            transcript = fields[ti]
            if not transcript:
                print >> sys.stderr, "%d: %s\n" % (linenum,line)
                continue
            chrom = fields[ci]
            pos = int(fields[pi])
            ref = fields[ri]
            alts = fields[ai]
            dp = int(fields[di])
            dpr = [int(x) for x in fields[fi].split(',')]
            for i,alt in enumerate(alts.split(',')):
                freq = float(dpr[i+1])/float(sum(dpr)) if dpr else None
                yield (transcript,pos,ref,alt,dp,freq)

    def parse_snpeff_vcf():
        for linenum,line in enumerate(inputFile):
            if line.startswith('##'): 
                if line.find('SnpEffVersion=') > 0:
                    SnpEffVersion = re.search('SnpEffVersion="?(\d+\.\d+)',line).groups()[0]
            elif line.startswith('#CHROM'): 
                pass
            else:
                fields = line.strip('\r\n').split('\t')
                if options.debug: print >> sys.stderr, "\n%s" % (fields)
                (chrom, pos, id, ref, alts, qual, filter, info) = fields[0:8]
                alt_list = alts.split(',')
                pos = int(pos)
                qual = float(qual)
                dp = None
                dpr = None
                for info_item in info.split(';'):
                    if info_item.find('=') < 0: continue
                    (key, val) = info_item.split('=', 1)
                    if key == 'DP':
                        dp = int(val)
                    if key == 'DPR':
                        dpr = [int(x) for x in val.split(',')]
                    if key in ['EFF','ANN']:
                        for effect in val.split(','):
                            if options.debug: print >> sys.stderr, "\n%s" % (effect.split('|'))
                            if key == 'ANN':
                                (alt,eff,impact,gene_name,gene_id,feature_type,transcript,biotype,exon,c_hgvs,p_hgvs,cdna,cds,aa,distance,info) = effect.split('|')
                            elif key == 'EFF':
                                (eff, effs) = effect.rstrip(')').split('(')
                                (impact, functional_class, codon_change, aa_change, aa_len, gene_name, biotype, coding, transcript, exon, alt) = effs.split('|')[0:11]
                            i = alt_list.index(alt) if alt in alt_list else 0
                            freq = float(dpr[i+1])/float(sum(dpr)) if dpr else None
                            yield (transcript,pos,ref,alt,dp,freq)


    #Process gene model
    ens_ref = None
    if options.gene_model != None:
        try:
            geneModelFile = os.path.abspath(options.gene_model)
            twoBitFile = os.path.abspath(options.twobit)
            print >> sys.stderr, "Parsing ensembl ref: %s %s" % (options.gene_model,options.twobit)
            time1 = time.time()
            ens_ref = EnsemblRef(geneModelFile,twoBitFile)
            time2 = time.time()
            print >> sys.stderr, "Parsing ensembl ref: %d seconds" % (int(time2-time1))
        except Exception, e:
            print >> sys.stderr, "Parsing gene model failed: %s" % e
            exit(2)
    try:
        parse_input = parse_tabular if options.format == 'tabular' else parse_snpeff_vcf
        for tid,pos1,ref,alt,dp,freq in parse_input():
            if not tid:
                continue
            if options.min_depth and dp is not None and dp < options.min_depth:
                continue
            if options.min_freq and freq is not None and freq < options.min_freq:
                continue
            ## transcript_id, pos, ref, alt, dp, dpr
            tx = ens_ref.get_gtf_transcript(tid)
            if not tx:
                continue
            coding = ens_ref.transcript_is_coding(tid)
            if not coding:
                continue
            frame_shift = len(ref) != len(alt)
            cds = ens_ref.get_cds(tid)
            pos0 = pos1 - 1 # zero based position
            spos = pos0 if tx.gene.strand else pos0 + len(ref) - 1
            alt_seq = alt if tx.gene.strand else reverse_complement(alt)
            ref_seq = ref if tx.gene.strand else reverse_complement(ref)
            cds_pos = ens_ref.genome_to_cds_pos(tid, spos)
            alt_cds = cds[:cds_pos] + alt_seq + cds[cds_pos+len(ref):] if cds_pos+len(ref) < len(cds) else '' 
            offset = 0
            if tx.gene.strand:
                for i in range(min(len(ref),len(alt))):
                    if ref[i] == alt[i]:
                        offset = i
                    else:
                        break
            else:
                for i in range(-1,-min(len(ref),len(alt)) -1,-1):
                    if ref[i] == alt[i]:
                        offset = i
                    else:
                        break
            refpep = translate(cds[:len(cds)/3*3])
            pep = translate(alt_cds[:len(alt_cds)/3*3])
            peplen = len(pep)
            aa_pos = (cds_pos + offset) / 3
            if aa_pos >= len(pep):
                print >> sys.stderr, "aa_pos %d >= peptide length %d : %s %d %s %s\n" % (aa_pos,len(pep),tid,pos1,ref,alt)
                continue
            if frame_shift:
                #find stop_codons
                nstops = 0
                stop_codons = []
                for i in range(aa_pos,peplen):
                    if refpep[i] != pep[i]:
                        aa_pos = i
                        break
                for i in range(aa_pos,peplen):
                    if pep[i] == '*':
                        nstops += 1
                        stop_codons.append("%s-%s%s" % (alt_cds[i*3-1],alt_cds[i*3:i*3+3],"-%s" % alt_cds[i*3+4] if len(alt_cds) > i*3 else ''))
                        if nstops > options.readthrough:
                            reported_peptide = pep[aa_pos:i+1]
                            reported_stop_codon = ','.join(stop_codons)
                            break
            else:
                reported_stop_codon = "%s-%s" % (alt_cds[peplen*3-4],alt_cds[peplen*3-3:peplen*3])
                reported_peptide = "%s_%s_%s" % (pep[max(aa_pos-options.leading_aa,0):aa_pos],
                                                 pep[aa_pos],
                                                 pep[aa_pos+1:min(aa_pos+1+options.trailing_aa,len(pep))])
            cs_pos = aa_pos * 3
            aa_pos = (cds_pos + offset) / 3
            ref_codon = cds[cs_pos:cs_pos+3]
            ref_aa = translate(ref_codon)
            alt_codon = alt_cds[cs_pos:cs_pos+3]
            alt_aa = translate(alt_codon)
            gene = tx.gene.names[0]
            report_fields = [tx.gene.names[0],
             '%s:%d %s' % (tx.gene.contig,pos1,'+' if tx.gene.strand else '-'),
             ref_seq,
             alt_seq,
             "%1.2f" % freq if freq is not None else '',
             str(dp),
             "%s|%s" % (tx.gene.gene_id,tx.cdna_id),
             "%d" % (aa_pos + 1),
             "%s%d%s" % (ref_aa,aa_pos + 1,alt_aa),
             "%d" % len(pep),
             reported_stop_codon,
             reported_peptide
             ]
            if options.debug:
                report_fields.append("%d %d %d %d  %s  %s" % (cds_pos, offset, cs_pos,aa_pos,ref_codon,alt_codon))
            outputFile.write('\t'.join(report_fields))
            if options.debug:
                print >> sys.stderr, "%s %s\n%s\n%s\n%s %s" % (
                         cds[cs_pos-6:cs_pos], cds[cs_pos:cs_pos+15],
                         translate(cds[cs_pos-6:cs_pos+15]),
                         translate(alt_cds[cs_pos-6:cs_pos+15]),
                         alt_cds[cs_pos-6:cs_pos], alt_cds[cs_pos:cs_pos+15])
            outputFile.write('\n')
    except Exception, e:
        print >> sys.stderr, "failed: %s" % e
        exit(1)


if __name__ == "__main__" : __main__()

