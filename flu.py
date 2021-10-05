""" This year, you prudently got the flu vaccine. So when your roommate, 
who forgot to get vaccinated, came down with the flu, you weren’t worried. 
But somehow, several days later, you started feeling feverish, weak, and 
all-around-awful. You knew it was the flu, but how could this be?!

You’ve heard of viral quasispecies, and you suspect that maybe a small portion
of the virus population mutated and evolved while replicating inside your 
roommate’s cells, which could explain how it was able to infect you. To find 
out, your have your friends set up a targeted deep sequencing experiment 
to analyze the HA genes in your roommate’s viral sample. They set up an
Illumina single-end sequencing run. When they send you the results,
you start analyzing your roommate’s sequence data right away.
"""
from fire import Fire

from sweets import BColors, sweet_logging
from sweets import run
from sweets import parse_vcf
from sweets import analyze_variant
from sweets import colorize

from pathlib import Path
from Bio import SeqIO


log = sweet_logging('flu')


def index(reference):
    run(f'bwa index {reference}')


def align(reference, sequence, output):
    run(f"bwa mem {reference} {sequence}"
        f" | samtools view -S -b -o - "
        f" | samtools sort - -o {output}"
    )
    run(f'samtools index {output}')


def varscan(reference, bam, output, *, freq=0.95, depth=16000):
    mpileup = bam.with_name(f'{bam.stem}_{depth}.mpileup')

    run(f'samtools mpileup -d {depth} {bam} -f {reference} -o {mpileup}',
        skip=mpileup.exists())

    pout = output.with_name(
        f'{output.stem}_{depth}_{freq}_{output.suffix}'
    )

    run(
        f'varscan mpileup2snp {mpileup} '
        f'--min-var-freq {freq} --variants '
        f'--output-vcf 1 > {pout}',
        skip=pout.exists()
    )

    variants = parse_vcf(pout)
    log.info(
        f'varscan found {len(variants)} variants (depth={depth} freq={freq})')

    return variants


class Main:

    REFERENCE = Path('data') / 'InfluenzaGene.fa'

    def _aling_(self, reads, out):
        return align(Main.REFERENCE, reads, out)

    def _varskan_(self, *, freq, depth):
        return varscan(
            Main.REFERENCE,
            self.alignment,
            self.variants,
            freq=freq,
            depth=depth
        )

    def _build_dirs(self, base_dir: Path):
        self.results = base_dir
        self.alignment = base_dir / 'alignment.sorted.bam'
        self.variants = base_dir / 'variants.vcf'
        return self.alignment, self.variants

    def prepare(self, patient, control_dir):
        patient = Path(patient)
        control_dir = Path(control_dir)

        respath = Path('results')
        alignment, _ = self._build_dirs(respath)
        self._aling_(patient, alignment)

        for srr_path in control_dir.iterdir():
            out_dir = respath / 'control' / srr_path.stem
            out_dir.mkdir(exist_ok=True, parents=True)
            out_alignment = out_dir / alignment.name

            self._aling_(srr_path, out_alignment)

    def experiment(self, base_dir, *, freq=0.001, depth=50000):
        base_dir = Path(base_dir)
        self._build_dirs(base_dir)
        gene = SeqIO.read(Main.REFERENCE, 'fasta')

        variants = self._varskan_(freq=freq, depth=depth)

        log = sweet_logging(base_dir.name)

        vars = []
        for variant in variants:
            ref, var = analyze_variant(gene.seq, variant)
            if ref.acid != var.acid:
                sample, *_ = variant.samples

                log.info(
                    f'{variant.start} {ref.acid}->{var.acid} {sample["FREQ"]}')
                vars.append((ref, var, variant))

        log.info(f'Found {len(vars)} variants which cause amino acid change')

        return base_dir.name, vars

    def _build_table(self, results):
        columns = []
        rows = set()
        table = {}

        for name, result in results:
            column = name
            columns.append(column)

            for ref, var, variant in result:
                row = variant.start
                rows.add(row)

                assert (row, column) not in table
                table[row, column] = f'{ref.acid}->{var.acid} '
                table[row, column] = f'{variant.samples[0]["FREQ"]}'

        return sorted(rows), columns, table

    def compare(self, *experiments, freq=0.001, depth=50000):
        results = list(map(self.experiment, map(Path, experiments)))

        rows, columns, table = self._build_table(results)

        width = 13

        print(BColors.BOLD + BColors.HEADER + '     #', end='')
        for column in columns:
            print(column.center(width), end='')
        print(BColors.ENDC)

        NODATA = '-'

        i = 0
        for row in rows:
            hcolor = [BColors.BOLD, BColors.ENDC][i % 2]
            head = hcolor + f'{row//3:6}' + BColors.ENDC
            line_color = BColors.ENDC
            
            data = []
            for column in columns:
                cell = table.get((row, column), None)
                data.append(cell)
            
            result, *control = data
            if not result:
                continue 

            line_color = BColors.WARNING

            if all(control):
                line_color = BColors.FAIL

            if all(map(lambda e: not e, control)):
                line_color = BColors.OKGREEN

            line = map(lambda e: (e or NODATA).center(width), data)
            print(head + colorize(''.join(line), color=line_color))
            i += 1


if __name__ == '__main__':
    Fire(Main)
