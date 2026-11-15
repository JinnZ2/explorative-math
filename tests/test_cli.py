import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CLITests(unittest.TestCase):
    def run_cli(self, *args, entrypoint=None, expected=0):
        command = [sys.executable]
        command += [str(entrypoint)] if entrypoint else ['-m', 'explorative_math']
        result = subprocess.run(command + list(args), cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def test_invalid_bounds_refuse_without_dataset_creation(self):
        with tempfile.TemporaryDirectory() as root:
            dataset = Path(root)/'invalid'
            self.run_cli('project','--out',str(dataset),'--shape','hypercube',
                         '--dimension','20','--max-vertices','16',expected=2)
            self.assertFalse(dataset.exists())
            self.run_cli('geometry','--out',str(dataset),'--samples','0',expected=2)
            self.assertFalse(dataset.exists())

    def test_source_geometry_and_observer_overlay(self):
        with tempfile.TemporaryDirectory() as root:
            dataset = Path(root)/'geometry'
            result = self.run_cli('geometry','--out',str(dataset),'--dimensions','2',
                                 '--samples','2','--steps','1')
            self.assertTrue(json.loads(result.stdout)['completed'])
            result = self.run_cli('validate','--dataset',str(dataset))
            self.assertTrue(json.loads(result.stdout)['valid'])
            first = json.loads((dataset/'records.jsonl').read_text().splitlines()[0])
            self.run_cli('observe','--dataset',str(dataset),'--entry',first['id'],
                         '--slot','human_interest','--value','false','--observer','test observer',
                         '--reason','explicit synthetic regression input')
            result = self.run_cli('show','--dataset',str(dataset),'--entry',first['id'])
            self.assertIs(json.loads(result.stdout)['human_interest'],False)
            self.assertEqual((dataset/'divergence.jsonl').read_text(),'')
            unchanged = json.loads((dataset/'records.jsonl').read_text().splitlines()[0])
            self.assertIsNone(unchanged['human_interest'])

    def test_phone_build_runs_offline_with_packaged_catalog(self):
        with tempfile.TemporaryDirectory() as root:
            application = Path(root)/'explorative-math.pyz'
            result = subprocess.run([sys.executable,str(ROOT/'tools/build_phone.py'),
                                     '--output',str(application)],cwd=ROOT,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertTrue(application.is_file())
            self.run_cli('--help',entrypoint=application)
            dataset = Path(root)/'projection'
            result = self.run_cli('project','--out',str(dataset),'--dimension','3','--frames','2',
                                 '--composition','phi',entrypoint=application)
            self.assertTrue(json.loads(result.stdout)['completed'])
            result = self.run_cli('validate','--dataset',str(dataset),entrypoint=application)
            self.assertTrue(json.loads(result.stdout)['valid'])
            record = json.loads((dataset/'records.jsonl').read_text().splitlines()[0])
            self.assertEqual(len(record['params']['catalog']),4)
            for frame in record['outputs'][1:]:
                self.assertEqual(len(frame['catalog_comparisons']),4)
            for slot in ('human_interest','ai_interest','compressibility','divergence_note'):
                self.assertIsNone(record[slot])

    def test_secondary_integer_command(self):
        with tempfile.TemporaryDirectory() as root:
            dataset = Path(root)/'integer'
            result = self.run_cli('enumerate','--family','integer','--out',str(dataset),
                                 '--a-min','0','--a-max','0','--b-min','0','--b-max','0',
                                 '--modulus','2','--steps','2')
            manifest = json.loads(result.stdout)
            self.assertEqual(manifest['total_runs'],2)
            self.assertEqual(manifest['total_outputs'],6)
            self.run_cli('validate','--dataset',str(dataset))


if __name__ == '__main__':
    unittest.main()
