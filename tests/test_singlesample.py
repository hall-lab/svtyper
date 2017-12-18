
from .context import singlesample as s
import unittest, os, subprocess

HERE = os.path.dirname(__file__)
in_vcf = os.path.join(HERE, "data/example.vcf")
in_bam = os.path.join(HERE, "data/NA12878.target_loci.sorted.bam")
lib_info_json = os.path.join(HERE, "data/NA12878.bam.json")
out_vcf = os.path.join(HERE, "data/out.vcf")
expected_out_vcf = os.path.join(HERE, "data/example.gt.vcf")

class TestIntegration(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        if os.path.exists(out_vcf):
            os.remove(out_vcf)

    def test_serial_integration(self):
        with open(in_vcf, "r") as inf, open(out_vcf, "w") as outf:
            s.sso_genotype(bam_string=in_bam,
                           vcf_in=inf,
                           vcf_out=outf,
                           min_aligned=20,
                           split_weight=1,
                           disc_weight=1,
                           num_samp=1000000,
                           lib_info_path=lib_info_json,
                           debug=False,
                           ref_fasta=None,
                           sum_quals=False,
                           max_reads=1000,
                           cores=None,
                           batch_size=1000)

        fail_msg = "did not find output vcf '{}' after running sv_genotype".format(out_vcf)
        self.assertTrue(os.path.exists(out_vcf), fail_msg)

        fail_msg = ("output vcf '{}' "
                    "did not match expected "
                    "output vcf '{}'").format(out_vcf, expected_out_vcf)
        self.assertTrue(self.diff(), fail_msg)

    def test_parallel_integration(self):
        with open(in_vcf, "r") as inf, open(out_vcf, "w") as outf:
            s.sso_genotype(bam_string=in_bam,
                           vcf_in=inf,
                           vcf_out=outf,
                           min_aligned=20,
                           split_weight=1,
                           disc_weight=1,
                           num_samp=1000000,
                           lib_info_path=lib_info_json,
                           debug=False,
                           ref_fasta=None,
                           sum_quals=False,
                           max_reads=1000,
                           cores=1,
                           batch_size=1000)

        fail_msg = "did not find output vcf '{}' after running sv_genotype".format(out_vcf)
        self.assertTrue(os.path.exists(out_vcf), fail_msg)

        fail_msg = ("output vcf '{}' "
                    "did not match expected "
                    "output vcf '{}'").format(out_vcf, expected_out_vcf)
        self.assertTrue(self.diff(), fail_msg)

    def diff(self):
        cmd = ['diff', "-I", "^##fileDate=", expected_out_vcf, out_vcf]

        rv = None
        with open(os.devnull, "w") as f:
#            rv = subprocess.call(cmd, stdout=f)
            rv = subprocess.call(cmd)

        result = rv == 0
        return result

if __name__ == '__main__':
    unittest.main(verbosity=2)
