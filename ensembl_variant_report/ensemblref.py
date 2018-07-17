"""
with gtf and twobit
get_bed(transcript_id)  
"""
from gtf_to_genes import gene, gene_utilities
from twobitreader import TwoBitFile
from Bio.Seq import reverse_complement, translate
import logging
logger = logging.getLogger("test")

class EnsemblRef(object):

    def __init__(self,gtf_file,twobitfile,read_now=True):
        self.gtf_file = gtf_file
        self.twobitfile = twobitfile
        self.twobit = TwoBitFile(self.twobitfile)
        self.gene_dict = None
        self.transcript_idx = None
        self.name_idx = None
        if read_now:
            self.get_transcript_idx()
        

    def get_gene_dict(self):
        if self.gene_dict is None:
            gene_structures = gene.t_parse_gtf('test')
            self.gene_dict = gene_structures.get_genes(self.gtf_file,logger=logger)
        return self.gene_dict


    def get_transcript_idx(self):
        if self.transcript_idx is None:
            self.transcript_idx = gene_utilities.index_transcripts(self.get_gene_dict(),by_prot_id=False)
        return self.transcript_idx


    def get_name_idx(self):
        if self.name_idx is None:
            self.name_idx = dict()
            for i,t in self.get_transcript_idx().items():
                for name in t.gene.names:
                   self.name_idx[name] = t.gene
                for name in t.names:
                   self.name_idx[name] = t
                if t.prot_id:
                   self.name_idx[t.prot_id] = t
        return self.name_idx


    def get_gtf_transcript(self,name):
        idx = self.get_transcript_idx()
        if name in idx:
            return idx[name]
        else:
            nidx = self.get_name_idx()
            if name in nidx:
                return nidx[name]
        return None


    def transcript_is_coding(self,transcript_id):
        tx = self.get_transcript_idx()[transcript_id]
        return len(tx.start_codons) > 0


    def get_transcript_start_codon(self,transcript_id):
        tx = self.get_transcript_idx()[transcript_id]
        return tx.start_codons[0] if len(tx.start_codons) > 0 else None


    def get_bed_line(self,transcript_id,score=0,itemRgb='0,0,0',coding=False):
        tx = self.get_transcript_idx()[transcript_id]
        chrom = tx.gene.contig
        chromStart = tx.coding_beg if coding else tx.beg 
        chromEnd = tx.coding_end if coding else tx.end
        name = transcript_id
        strand = '+' if tx.gene.strand else '-'
        thickStart = tx.coding_beg if tx.coding_beg else chromStart
        thickEnd = tx.coding_end if tx.coding_end else chromEnd
        exons = tx.get_coding_exons() if coding else tx.get_exons()
        blockCount = len(exons)
        if tx.gene.strand:
            strand = '+'
            blockSizes = [abs(e-s) for s,e in exons]
            blockStarts = [s - chromStart for s,e in exons]
        else:
            strand = '-'
            blockSizes = [abs(e-s) for s,e in reversed(exons)]
            blockStarts = [s - chromStart for s,e in reversed(exons)]
        blockSizes = ','.join([str(x) for x in blockSizes])
        blockStarts = ','.join([str(x) for x in blockStarts])
        return '%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s' % (chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts)


    def transcripts_in_range(self,chrom,startpos,endpos,strand=None):
        spos = min(startpos,endpos) if endpos else startpos
        epos = max(startpos,endpos) if endpos else startpos
        transcripts = []
        for i,t in self.get_transcript_idx().items():
            if t.gene.contig == chrom and t.beg <= epos and spos <= t.end:
                if strand and t.gene.strand != strand:
                    continue
                transcripts.append(t)
        return transcripts


    def genes_in_range(self,chrom,startpos,endpos,strand=None,gene_types=None):
        spos = min(startpos,endpos) if endpos else startpos
        epos = max(startpos,endpos) if endpos else startpos
        gene_dict = self.get_gene_dict()
        gtypes = set(gene_types)  & set(gene_dict.keys()) if gene_types else set(gene_dict.keys())
        genes = []
        for gt in gtypes:
            for gene in gene_dict[gt]:
                if gene.contig == chrom and gene.beg <= epos and spos <= gene.end:
                    if strand and gene.strand != strand:
                        continue
                    genes.append(gene)
        return genes


    def get_sequence(self,chrom,start,end):
        if self.twobit:
            if chrom in self.twobit:
                return self.twobit[chrom][start:end]
            contig = chrom[3:] if chrom.startswith('chr') else 'chr%s' % chrom
            if contig in self.twobit:
                return self.twobit[contig][start:end]
        return None


    def sequence_sizes(self):
       return self.twobit.sequence_sizes()


    def get_transcript_seq(self,transcript_id,coding=False):
        tx = self.get_transcript_idx()[transcript_id]
        chrom = tx.gene.contig
        exonbnds = tx.get_coding_exons() if coding else tx.get_exons()
        if tx.gene.strand:
            seqs = [self.get_sequence(chrom,s,e) for s,e in exonbnds]
        else:
            seqs = [reverse_complement(self.get_sequence(chrom,s,e)) for s,e in exonbnds]
        return ''.join(seqs)


    def get_cdna(self,transcript_id):
        return self.get_transcript_seq(transcript_id,coding=False)


    def get_cds(self,transcript_id):
        return self.get_transcript_seq(transcript_id,coding=True)


    def genome_to_transcript_pos(self,transcript_id,genome_pos,coding=False):
        tx = self.get_transcript_idx()[transcript_id]
        if not tx.beg <= genome_pos < tx.end:
            return None
        exonbnds = tx.get_coding_exons() if coding else tx.get_exons()
        cdna_pos = 0
        if tx.gene.strand:
            for s,e in exonbnds:
                if s <= genome_pos < e:
                    cdna_pos += genome_pos - s
                    break
                else:
                    cdna_pos += e - s
        else:
            for s,e in exonbnds:
                if s <= genome_pos < e:
                    cdna_pos += e - genome_pos - 1
                    break
                else:
                    cdna_pos += e - s
        return cdna_pos


    def genome_to_cdna_pos(self,transcript_id,genome_pos):
        return self.genome_to_transcript_pos(transcript_id,genome_pos,coding=False)


    def genome_to_cds_pos(self,transcript_id,genome_pos):
        return self.genome_to_transcript_pos(transcript_id,genome_pos,coding=True)


