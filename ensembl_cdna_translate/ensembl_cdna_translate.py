#!/usr/bin/env python
"""
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#         Copyright 2017, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  James E Johnson
#
#------------------------------------------------------------------------------
"""

import argparse
import re
import sys
from time import sleep

from Bio.Seq import translate

import requests

from bedutil import BedEntry, bed_from_line
import digest
from ensembl_rest import get_toplevel, get_transcripts_bed, get_cds, get_cdna, max_region
from twobitreader import TwoBitFile


def __main__():
    parser = argparse.ArgumentParser(
        description='Retrieve Ensembl cDNAs and three frame translate')
    parser.add_argument(
        '-s', '--species', default='human',
        help='Ensembl Species to retrieve')
    parser.add_argument(
        '-R', '--regions', action='append', default=[],
        help='Restrict Ensembl retrieval to regions e.g.: X,2:20000-25000,3:100-500+')
    parser.add_argument(
        '-B', '--biotypes', action='append', default=[],
        help='Restrict Ensembl biotypes to retrieve')
    parser.add_argument(
        '-i', '--input', default=None,
        help='Use BED instead of retrieving cDNA from ensembl (-) for stdin')
    parser.add_argument(
        '-T', '--twobit', default=None,
        help='Genome reference sequence in 2bit format')
    parser.add_argument(
        '-t', '--transcripts', default=None,
        help='Path to output cDNA transcripts.bed (-) for stdout')
    parser.add_argument(
        '-r', '--raw', action='store_true',
        help='Report transcript exacty as returned from Ensembl')
    parser.add_argument(
        '-f', '--fasta', default=None,
        help='Path to output translations.fasta')
    parser.add_argument(
        '-b', '--bed', default=None,
        help='Path to output translations.bed')
    parser.add_argument(
        '-m', '--min_length', type=int, default=7,
        help='Minimum length of protein translation to report')
    parser.add_argument(
        '-e', '--enzyme', default=None,
        help='Digest translation with enzyme')
    parser.add_argument(
        '-a', '--all', action='store_true',
        help='Include reference protein translations')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    parser.add_argument('-d', '--debug', action='store_true', help='Debug')
    args = parser.parse_args()
    # print >> sys.stderr, "args: %s" % args
    species = args.species
    input_rdr = None
    if args.input is not None:
        input_rdr = open(args.input, 'r') if args.input != '-' else sys.stdin
    tx_wtr = None
    if args.transcripts is not None:
        tx_wtr = open(args.transcripts, 'w')\
            if args.transcripts != '-' else sys.stdout
    fa_wtr = open(args.fasta, 'w') if args.fasta is not None else None
    bed_wtr = open(args.bed, 'w') if args.bed is not None else None

    enzyme = digest.expasy_rules.get(args.enzyme,args.enzyme)

    # print >> sys.stderr, "args biotypes: %s" % args.biotypes
    biotypea = ['biotype=%s' % bt.strip() for biotype in args.biotypes for bt in biotype.split(',')] 
    # print >> sys.stderr, "args biotypes: %s" % biotypea
    biotypes = ';'.join(['biotype=%s' % bt.strip() for biotype in args.biotypes for bt in biotype.split(',') if  bt.strip()]) 
    # print >> sys.stderr, "biotypes: %s" % biotypes

    twobit = TwoBitFile(args.twobit) if args.twobit else None

    selected_regions = dict() # chrom:(start,end)
    region_pat = '^([^:]+)(?::(\d*)(?:-(\d+)([+-])?)?)?'
    if args.regions:
        for entry in args.regions:
           if not entry:
               continue
           regs = [x.strip() for x in entry.split(',') if x.strip()]
           for reg in regs:
               m = re.match(region_pat,reg)
               if m:
                   (chrom,start,end,strand) = m.groups()
                   if chrom:
                       if chrom not in selected_regions:
                           selected_regions[chrom] = []
                       selected_regions[chrom].append([start,end,strand])
        if args.debug: print >> sys.stderr, "selected_regions: %s" % selected_regions

    translations = dict()  # start : end : seq

    def unique_prot(tbed, seq):
        if tbed.chromStart not in translations:
            translations[tbed.chromStart] = dict()
            translations[tbed.chromStart][tbed.chromEnd] = []
            translations[tbed.chromStart][tbed.chromEnd].append(seq)
        elif tbed.chromEnd not in translations[tbed.chromStart]:
            translations[tbed.chromStart][tbed.chromEnd] = []
            translations[tbed.chromStart][tbed.chromEnd].append(seq)
        elif seq not in translations[tbed.chromStart][tbed.chromEnd]:
            translations[tbed.chromStart][tbed.chromEnd].append(seq)
        else:
            return False
        return True

    def get_sequence(chrom,start,end):
        if twobit:
            if chrom in twobit:
                return twobit[chrom][start:end]
            contig = chrom[3:] if chrom.startswith('chr') else 'chr%s' % chrom
            if contig in twobit:
                return twobit[contig][start:end]
        return None

    def translate_bed(bed):
        translate_count = 0
        if any([fa_wtr, bed_wtr]):
            transcript_id = bed.name
            refprot = None
            if twobit:
                bed.seq = get_sequence(bed.chrom,bed.chromStart,bed.chromEnd)
            else:
                bed.cdna = get_cdna(transcript_id)
            cdna = bed.get_cdna()
            cdna_len = len(cdna)
            if not args.all:
                try:
                    cds = bed.get_cds()
                    if cds is None:
                        cds = get_cds(transcript_id)
                    if len(cds) % 3 != 0:
                        cds = cds[:-(len(cds) % 3)]
                    refprot = translate(cds) if cds else None
                except:
                    refprot = None
            for offset in range(3):
                seqend = cdna_len - (cdna_len - offset) % 3
                aaseq = translate(cdna[offset:seqend])
                aa_start = 0
                while aa_start < len(aaseq):
                    aa_end = aaseq.find('*', aa_start)
                    if aa_end < 0:
                        aa_end = len(aaseq)
                    prot = aaseq[aa_start:aa_end]
                    if enzyme and refprot:
                        frags = digest._cleave(prot,enzyme)
                        for frag in reversed(frags):
                            if frag in refprot:
                                prot = prot[:prot.rfind(frag)]
                            else:
                                break
                    if len(prot) < args.min_length:
                        pass
                    elif refprot and prot in refprot:
                        pass
                    else:
                        tstart = aa_start*3+offset
                        tend = aa_end*3+offset
                        prot_acc = "%s_%d_%d" % (transcript_id, tstart, tend)
                        tbed = bed.trim(tstart, tend)
                        if args.all or unique_prot(tbed, prot):
                            translate_count += 1
                            tbed.name = prot_acc
                            bed_wtr.write("%s\t%s\n" % (str(tbed), prot))
                            bed_wtr.flush()
                            fa_id = ">%s\n" % (prot_acc)
                            fa_wtr.write(fa_id)
                            fa_wtr.write(prot)
                            fa_wtr.write("\n")
                            fa_wtr.flush()
                    aa_start = aa_end + 1
        return translate_count

    def translate_region(species,ref,start,stop,strand):
        translation_count = 0
        regions = range(start, stop, max_region)
        if not regions or regions[-1] < stop:
            regions.append(stop)
        for end in regions[1:]:
            bedlines = get_transcripts_bed(species, ref, start, end, strand=strand, params=biotypes)
            if args.verbose or args.debug:
                print >> sys.stderr,\
                    "%s\t%s\tstart: %d\tend: %d\tcDNA transcripts:%d"\
                    % (species, ref, start, end, len(bedlines))
            # start, end, seq
            for i, bedline in enumerate(bedlines):
                try:
                    bed = bed_from_line(bedline)\
                        if any([not args.raw, fa_wtr, bed_wtr])\
                        else None
                    if tx_wtr:
                        tx_wtr.write(bedline if args.raw else str(bed))
                        tx_wtr.write("\n")
                        tx_wtr.flush()
                    if bed:
                        translation_count += translate_bed(bed)
                except Exception as e:
                    print >> sys.stderr,\
                        "BED error (%s) : %s\n" % (e, bedline)
            start = end + 1
        return translation_count

    if input_rdr:
        translation_count = 0
        for i, bedline in enumerate(input_rdr):
            try:
                bed = bed_from_line(bedline)
                if bed is None:
                    continue
                if bed.biotype and biotypea and bed.biotype not in biotypea:
                    continue
                translation_count += translate_bed(bed)
            except:
                print >> sys.stderr, "BED format error: %s\n" % bedline
        if args.debug or (args.verbose and any([fa_wtr, bed_wtr])):
            print >> sys.stderr,\
                "%s\tcDNA translations:%d" % (species, translation_count)
    else:
        coord_systems = get_toplevel(species)
        if 'chromosome' in coord_systems:
            ref_lengths = dict()
            for ref in sorted(coord_systems['chromosome'].keys()):
                length = coord_systems['chromosome'][ref]
                ref_lengths[ref] = length
                if not any([tx_wtr, fa_wtr, bed_wtr]):
                    print >> sys.stderr,\
                        "%s\t%s\tlength: %d" % (species, ref, length)
            if selected_regions:
                translation_count = 0
                for ref in sorted(selected_regions.keys()):
                    if ref in ref_lengths:
                        for reg in selected_regions[ref]:
                            (_start,_stop,_strand) = reg
                            start = int(_start) if _start else 0
                            stop = int(_stop) if _stop else ref_lengths[ref]
                            strand = '' if not _strand else ':1' if _strand == '+' else  ':-1'
                            translation_count += translate_region(species,ref,start,stop,strand)
            else:
                strand = ''
                start = 0
                for ref in sorted(ref_lengths.keys()):
                    length = ref_lengths[ref]
                    translation_count = 0
                    if args.debug:
                        print >> sys.stderr,\
                            "Retrieving transcripts: %s\t%s\tlength: %d"\
                            % (species, ref, length)
                    translation_count += translate_region(species,ref,start,length,strand)
                    if args.debug or (args.verbose and any([fa_wtr, bed_wtr])):
                        print >> sys.stderr,\
                            "%s\t%s\tlength: %d\tcDNA translations:%d"\
                            % (species, ref, length, translation_count)


if __name__ == "__main__":
    __main__()
