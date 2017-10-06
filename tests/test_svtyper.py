
from .context import parsers as p
from .context import core
import unittest, os, subprocess

HERE = os.path.dirname(__file__)
in_vcf = os.path.join(HERE, "data/example.vcf")
in_bam = os.path.join(HERE, "data/NA12878.target_loci.sorted.bam")
lib_info_json = os.path.join(HERE, "data/NA12878.bam.json")
out_vcf = os.path.join(HERE, "data/out.vcf")
expected_out_vcf = os.path.join(HERE, "data/example.gt.vcf")

class TestCigarParsing(unittest.TestCase):
    def test_cigarstring_to_tuple(self):
        string1 = '5H3S2D1N5M3I2P2X1='
        self.assertEqual(p.SplitRead.cigarstring_to_tuple(string1),
                [(5, 5), (4, 3), (2, 2), (3, 1),
                    (0, 5), (1, 3), (6, 2), (8, 2),
                    (7, 1)])

    def test_get_query_pos_from_cigar(self):
        # forward
        cigar_string = '2S3M1D2M2I3M3S'
        cigar = p.SplitRead.cigarstring_to_tuple(cigar_string)
        query_pos = p.SplitRead.SplitPiece.get_query_pos_from_cigar(cigar, True)
        self.assertEqual(query_pos.query_start, 3)
        self.assertEqual(query_pos.query_end, 13)
        self.assertEqual(query_pos.query_length, 15)

        # get_query_pos_from_cigar currently modifies the cigar list in place.
        # that's why the code below doesn't work as intended.
        query_pos = p.SplitRead.SplitPiece.get_query_pos_from_cigar(cigar, False)
        self.assertEqual(query_pos.query_start, 2)
        self.assertEqual(query_pos.query_end, 12)
        self.assertEqual(query_pos.query_length, 15)

    def test_get_reference_end_from_cigar(self):
        cigar_string = '2S5M3D2M3S'
        cigar = p.SplitRead.cigarstring_to_tuple(cigar_string)
        self.assertEqual(p.SplitRead.get_reference_end_from_cigar(1, cigar), 11)

    def test_get_start_diagonal(self):
        cigar_string = '2S5M3D1I1M3S'
        split_piece = p.SplitRead.SplitPiece(1, 25, True, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        self.assertEqual(p.SplitRead.get_start_diagonal(split_piece), 23)
        split_piece2 = p.SplitRead.SplitPiece(1, 25, False, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        self.assertEqual(p.SplitRead.get_start_diagonal(split_piece2), 23)

    def test_get_end_diagonal(self):
        cigar_string = '2S5M3D2I1M3S'
        split_piece = p.SplitRead.SplitPiece(1, 25, True, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        split_piece.set_reference_end(34)
        self.assertEqual(p.SplitRead.get_end_diagonal(split_piece), 34 - (2 + 8))
        split_piece2 = p.SplitRead.SplitPiece(1, 25, False, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        split_piece2.set_reference_end(34)
        self.assertEqual(p.SplitRead.get_end_diagonal(split_piece2), 34 - (2 + 8))

class TestIntegration(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        if os.path.exists(out_vcf):
            os.remove(out_vcf)

    def test_integration(self):
        with open(in_vcf, "r") as inf, open(out_vcf, "w") as outf:
            core.sv_genotype(bam_string=in_bam,
                             vcf_in=inf,
                             vcf_out=outf,
                             min_aligned=20,
                             split_weight=1,
                             disc_weight=1,
                             num_samp=1000000,
                             lib_info_path=lib_info_json,
                             debug=False,
                             alignment_outpath=None,
                             ref_fasta=None,
                             sum_quals=False,
                             max_reads=None)

        fail_msg = "did not file output vcf '{}' after running sv_genotype".format(out_vcf)
        self.assertTrue(os.path.exists(out_vcf), fail_msg)

        fail_msg = ("output vcf '{}' "
                    "did not match expected "
                    "output vcf '{}'").format(out_vcf, expected_out_vcf)
        self.assertTrue(self.diff(), fail_msg)

    def diff(self):
        cmd = ['diff', "-I", "^##fileDate=", expected_out_vcf, out_vcf]

        rv = None
        with open(os.devnull, "w") as f:
            rv = subprocess.call(cmd, stdout=f)

        result = True if rv == 0 else False
        return result

if __name__ == '__main__':
    unittest.main(verbosity=2)
