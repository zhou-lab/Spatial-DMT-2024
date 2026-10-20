import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


class ReferenceToolsTest(unittest.TestCase):
    def test_lambda_sampling_manifest_is_copied_to_final_output(self):
        snakefile = (REPO / "main" / "Snakefile").read_text()
        self.assertIn(
            'lambda_sampling = OUTDIR + "/{sample}/{sample}_lambda_pileup_sampling.tsv"',
            snakefile,
        )

    def test_slurm_profile_exports_active_environment(self):
        profile = (REPO / "profiles" / "slurm" / "config.yaml").read_text()
        self.assertIn("--export=ALL", profile)

    def test_build_cpg_bed_uses_zero_based_dyad_intervals(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "reference.fa"
            fasta.write_text(">chrA description\nAACGTCG\n>lambda\nCGNN\n")
            result = subprocess.run(
                [sys.executable, str(REPO / "tools" / "build_cpg_bed.py"), str(fasta)],
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertEqual(
            result.stdout,
            "chrA\t2\t4\nchrA\t5\t7\nlambda\t0\t2\n",
        )


if __name__ == "__main__":
    unittest.main()
