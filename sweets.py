import os
import logging
import vcf

from Bio.Seq import Seq
from pathlib import Path
from dataclasses import dataclass, field


TRI = 3

log = logging.getLogger('sweets')

def sweet_logging(logname: str):
    FORMAT='[%(levelname)-6s][%(name)s] %(message)s'
    logging.basicConfig(format=FORMAT, level=logging.DEBUG)
    return logging.getLogger(logname)
    

def run(command, skip=False):
    log.info('> ' + command)
    if skip:
        log.info('Command skipped due the skip condition')
        return 0

    errcode = os.system(command)
    if errcode != 0:
        raise Exception(f'Command executed with status code {errcode}')
    return errcode



class BColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'


def colorize(str, start=0, end=None, *, color=BColors.WARNING):
    if not end:
        end = len(str)
    return str[:start] + color + str[start:end] + BColors.ENDC + str[end:]


def highlight(str, position, length=1, *, color=BColors.WARNING):
    return colorize(str, position, position+length, color=color)


def parse_vcf(input):
    input = Path(input)

    with input.open('r') as file:
        return list(vcf.Reader(file))


def log_compare(a, b, pos, *, env=20, color=BColors.WARNING):
    d = '...'
    spark = env
    if len(a) <= env*2:
        d = ''
        spark = pos

    h = lambda s: highlight(s, spark, color=color)
    log.info(f'from: {d}{h(a[pos-env:pos+env])}{d}')
    log.info(f'  to: {d}{h(b[pos-env:pos+env])}{d}')

    return a[pos] == b[pos]


@dataclass
class AnalyseRecord:
    dna: Seq
    codon: Seq
    prot: Seq = field(init=False)
    acid: Seq = field(init=False)

    def __post_init__(self):
        self.prot = self.dna.translate()
        self.acid = self.codon.translate()


def log_analysis(variant, ref:AnalyseRecord, var:AnalyseRecord):
    pos, ref_base, alt_base = variant.start, str(variant.REF), str(variant.ALT[0])
    aa_pos = pos//TRI
    codon_pos = aa_pos * TRI
    
    c = colorize
    sample, *_ = variant.samples
    freq = sample['FREQ']
    log.info(f'Analysis of the variant {c(ref_base)}->{c(alt_base)} freq={freq}')
    
    ncolor = BColors.OKGREEN
    log.info(colorize('Nucleotide sequence', color=ncolor))
    log_compare(ref.dna, var.dna, pos, color=ncolor)

    log_compare(ref.codon, var.codon, pos-codon_pos, color=ncolor)

    acolor = BColors.OKBLUE
    log.info(colorize('Amino acid sequence', color=acolor))
    log_compare(ref.prot, var.prot, aa_pos, color=acolor)

    log_compare(ref.acid, var.acid, 0, color=acolor)


def analyze_variant(ref_seq, variant):
    pos, ref, alt = variant.start, str(variant.REF), str(variant.ALT[0])

    aa_pos = pos//TRI
    codon_pos = aa_pos * TRI

    ref_record = AnalyseRecord(
        ref_seq,
        ref_seq[codon_pos:codon_pos+TRI]
    )

    assert ref_record.dna[pos] == ref

    new_seq = ref_seq[:pos] + alt + ref_seq[pos+1:]
    var_record = AnalyseRecord(new_seq, new_seq[codon_pos:codon_pos+TRI])
    
    assert ref_record.acid == ref_record.prot[aa_pos]
    assert var_record.acid == var_record.prot[aa_pos]

    return ref_record, var_record
