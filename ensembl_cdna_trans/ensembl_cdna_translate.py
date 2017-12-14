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

import digest


server = "https://rest.ensembl.org"
ext = "/info/assembly/homo_sapiens?"
max_region = 4000000


def ensembl_rest(ext, headers):
    if True: print >> sys.stderr, "%s" % ext
    r = requests.get(server+ext, headers=headers)
    if r.status_code == 429:
        print >> sys.stderr, "response headers: %s\n" % r.headers
        if 'Retry-After' in r.headers:
            sleep(r.headers['Retry-After'])
            r = requests.get(server+ext, headers=headers)
    if not r.ok:
        r.raise_for_status()
    return r


def get_species():
    results = dict()
    ext = "/info/species"
    req_header = {"Content-Type": "application/json"}
    r = ensembl_rest(ext, req_header)
    for species in r.json()['species']:
        results[species['name']] = species
        print >> sys.stdout,\
            "%s\t%s\t%s\t%s\t%s"\
            % (species['name'], species['common_name'], 
               species['display_name'],
               species['strain'],
               species['taxon_id'])
    return results


def get_biotypes(species):
    biotypes = []
    ext = "/info/biotypes/%s?" % species
    req_header = {"Content-Type": "application/json"}
    r = ensembl_rest(ext, req_header)
    for entry in r.json():
        if 'biotype' in entry:
            biotypes.append(entry['biotype'])
    return biotypes
        

def get_toplevel(species):
    coord_systems = dict()
    ext = "/info/assembly/%s?" % species
    req_header = {"Content-Type": "application/json"}
    r = ensembl_rest(ext, req_header)
    toplevel = r.json()
    for seq in toplevel['top_level_region']:
        if seq['coord_system'] not in coord_systems:
            coord_systems[seq['coord_system']] = dict()
        coord_system = coord_systems[seq['coord_system']]
        coord_system[seq['name']] = int(seq['length'])
    return coord_systems


def get_transcripts_bed(species, refseq, start, length, strand='', params=None):
    bed = []
    param = params if params else ''
    req_header = {"Content-Type": "text/x-bed"}
    regions = range(start, length, max_region)
    if not regions or regions[-1] < length:
        regions.append(length)
    for end in regions[1:]:
        ext = "/overlap/region/%s/%s:%d-%d%s?feature=transcript;%s"\
            % (species, refseq, start, end, strand, param)
        start = end + 1
        r = ensembl_rest(ext, req_header)
        if r.text:
            bed += r.text.splitlines()
    return bed


def get_seq(id, seqtype,params=None):
    param = params if params else ''
    ext = "/sequence/id/%s?type=%s;%s" % (id, seqtype,param)
    req_header = {"Content-Type": "text/plain"}
    r = ensembl_rest(ext, req_header)
    return r.text


def get_cdna(id,params=None):
    return get_seq(id, 'cdna',params=params)


def get_cds(id,params=None):
    return get_seq(id, 'cds',params=params)


def get_transcript_haplotypes(species,transcript):
    ext = "/transcript_haplotypes/%s/%s?aligned_sequences=1" % (species,transcript)
    req_header = {"Content-Type" : "application/json"}
    r = ensembl_rest(ext, req_header)
    decoded = r.json()


def bed_from_line(line):
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) < 12:
        return None
    (chrom, chromStart, chromEnd, name, score, strand,
     thickStart, thickEnd, itemRgb,
     blockCount, blockSizes, blockStarts) = fields[0:12]
    bed_entry = BedEntry(chrom=chrom, chromStart=chromStart, chromEnd=chromEnd,
                         name=name, score=score, strand=strand,
                         thickStart=thickStart, thickEnd=thickEnd,
                         itemRgb=itemRgb,
                         blockCount=blockCount,
                         blockSizes=blockSizes.rstrip(','),
                         blockStarts=blockStarts.rstrip(','))
    if len(fields) == 20:
        bed_entry.second_name = fields[12]
        bed_entry.cds_start_status = fields[13]
        bed_entry.cds_end_status = fields[14]
        bed_entry.exon_frames = fields[15]
        bed_entry.biotype = fields[16]
        bed_entry.gene_name = fields[17]
        bed_entry.second_gene_name = fields[18]
        bed_entry.gene_type = fields[19]
    return bed_entry


class BedEntry(object):
    def __init__(self, chrom=None, chromStart=None, chromEnd=None,
                 name=None, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None,
                 blockCount=None, blockSizes=None, blockStarts=None):
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.name = name
        self.score = int(score) if score is not None else 0
        self.strand = '-' if str(strand).startswith('-') else '+'
        self.thickStart = int(thickStart) if thickStart else self.chromStart
        self.thickEnd = int(thickEnd) if thickEnd else self.chromEnd
        self.itemRgb = str(itemRgb) if itemRgb is not None else r'100,100,100'
        self.blockCount = int(blockCount)
        if isinstance(blockSizes, str) or isinstance(blockSizes, unicode):
            self.blockSizes = [int(x) for x in blockSizes.split(',')]
        elif isinstance(blockSizes, list):
            self.blockSizes = [int(x) for x in blockSizes]
        else:
            self.blockSizes = blockSizes
        if isinstance(blockStarts, str) or isinstance(blockSizes, unicode):
            self.blockStarts = [int(x) for x in blockStarts.split(',')]
        elif isinstance(blockStarts, list):
            self.blockStarts = [int(x) for x in blockStarts]
        else:
            self.blockStarts = blockStarts
        self.second_name = None
        self.cds_start_status = None
        self.cds_end_status = None
        self.exon_frames = None
        self.biotype = None
        self.gene_name = None
        self.second_gene_name = None
        self.gene_type = None
        self.seq = None
        self.pep = None

    def __str__(self):
        return '%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s' % (
            self.chrom, self.chromStart, self.chromEnd,
            self.name, self.score, self.strand,
            self.thickStart, self.thickEnd, str(self.itemRgb), self.blockCount,
            ','.join([str(x) for x in self.blockSizes]),
            ','.join([str(x) for x in self.blockStarts]))

    # (start, end)
    def get_subrange(self, tstart, tstop, debug=False):
        chromStart = self.chromStart
        chromEnd = self.chromEnd
        if debug:
            print >> sys.stderr, "%s" % (str(self))
        r = range(self.blockCount)
        if self.strand == '-':
            r.reverse()
        bStart = 0
        bEnd = 0
        for x in r:
            bEnd = bStart + self.blockSizes[x]
            if bStart <= tstart < bEnd:
                if self.strand == '+':
                    chromStart = self.chromStart + self.blockStarts[x] +\
                        (tstart - bStart)
                else:
                    chromEnd = self.chromStart + self.blockStarts[x] +\
                        self.blockSizes[x] - (tstart - bStart)
            if bStart <= tstop < bEnd:
                if self.strand == '+':
                    chromEnd = self.chromStart + self.blockStarts[x] +\
                        (tstop - bStart)
                else:
                    chromStart = self.chromStart + self.blockStarts[x] +\
                        self.blockSizes[x] - (tstop - bStart)
            if debug:
                print >> sys.stderr,\
                    "%3d %s\t%d\t%d\t%d\t%d\t%d\t%d"\
                    % (x, self.strand, bStart, bEnd,
                       tstart, tstop, chromStart, chromEnd)
            bStart += self.blockSizes[x]
        return(chromStart, chromEnd)

    # get the blocks for sub range
    def get_blocks(self, chromStart, chromEnd):
        tblockCount = 0
        tblockSizes = []
        tblockStarts = []
        for x in range(self.blockCount):
            bStart = self.chromStart + self.blockStarts[x]
            bEnd = bStart + self.blockSizes[x]
            if bStart > chromEnd:
                break
            if bEnd < chromStart:
                continue
            cStart = max(chromStart, bStart)
            tblockStarts.append(cStart - chromStart)
            tblockSizes.append(min(chromEnd, bEnd) - cStart)
            tblockCount += 1
        return (tblockCount, tblockSizes, tblockStarts)

    def trim(self, tstart, tstop, debug=False):
        (tchromStart, tchromEnd) =\
            self.get_subrange(tstart, tstop, debug=debug)
        (tblockCount, tblockSizes, tblockStarts) =\
            self.get_blocks(tchromStart, tchromEnd)
        tbed = BedEntry(
            chrom=self.chrom, chromStart=tchromStart, chromEnd=tchromEnd,
            name=self.name, score=self.score, strand=self.strand,
            thickStart=tchromStart, thickEnd=tchromEnd, itemRgb=self.itemRgb,
            blockCount=tblockCount,
            blockSizes=tblockSizes, blockStarts=tblockStarts)
        if self.seq:
            ts = tchromStart-self.chromStart
            te = tchromEnd - tchromStart + ts
            tbed.seq = self.seq[ts:te]
        return tbed


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

    def translate_bed(bed):
        translate_count = 0
        if any([fa_wtr, bed_wtr]):
            transcript_id = bed.name
            refprot = None
            if not args.all:
                try:
                    cds = get_cds(transcript_id)
                    if len(cds) % 3 != 0:
                        cds = cds[:-(len(cds) % 3)]
                    refprot = translate(cds) if cds else None
                except:
                    refprot = None
            cdna = get_cdna(transcript_id)
            cdna_len = len(cdna)
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
