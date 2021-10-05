import unittest
from unittest.case import addModuleCleanup
import flu
import os

from pathlib import Path


class FluTestCase(unittest.TestCase):
    pass


DATA = Path('data')
CONTROL = DATA / 'control'

RESULTS = Path('results')
RESULTS.mkdir(parents=True, exist_ok=True)
RESULT_CONTROL = RESULTS / 'control'
RESULT_CONTROL.mkdir(exist_ok=True)


REFERENCE = DATA / 'InfluenzaGene.fa'
READS = DATA / 'SRR1705851.fastq'

ALIGNMENT = RESULTS / 'alignment.sorted.bam'
VARIANTS = RESULTS / 'variants.vcf'

class FluTest(FluTestCase):

    def test_ok(self):
        self.assertTrue(True)

    def test_align(self):
        flu.align(REFERENCE, READS, ALIGNMENT)
        self.assertTrue(ALIGNMENT.exists())

    def test_varscan(self):
        flu.varscan(REFERENCE, ALIGNMENT, VARIANTS, freq=0.1, depth=8000)

    def test_align_control(self):
        for srr_path in CONTROL.iterdir():
            out_dir = RESULT_CONTROL / srr_path.stem
            out_dir.mkdir(exist_ok=True, parents=True)
            out_alignment = out_dir / ALIGNMENT.name

            flu.align(REFERENCE, srr_path, out_alignment)

    def test_varscan_control(self):
        for control_dir in RESULT_CONTROL.iterdir():
            srr = control_dir.name
            alignment = control_dir / ALIGNMENT.name
            out_variants = control_dir / VARIANTS.name

            flu.varscan(REFERENCE, alignment, out_variants, freq=0.001, depth=50000)



if __name__ == '__main__':
    unittest.main()