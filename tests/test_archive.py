import tempfile
import unittest
from pathlib import Path

from explorative_math.annotations import (
    add_note, get_annotation_history, observe, show_record,
)
from explorative_math.core import (
    EnumerationConfig, evolve_ca, iter_binary_words, iter_ca_trajectory,
    iter_integer_trajectory, iter_records, iter_substitution_rules,
)
from explorative_math.geometry import count_geometry, iter_geometry_runs, validate_geometry_config
from explorative_math.storage import (
    DatasetError, create_dataset, dataset_lock, enumerate_dataset, read_jsonl,
)
from explorative_math.validation import validate_dataset


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.path = self.root / 'dataset'
        config = validate_geometry_config({'experiment': 'trajectories', 'dimensions': [2], 'steps': 2})
        counts = count_geometry(config)
        self.manifest = create_dataset(self.path, config, iter_geometry_runs(config),
                                       counts['runs'], counts['outputs'], source='geometry')
        self.entry = next(read_jsonl(self.path / 'records.jsonl'))['id']

    def tearDown(self):
        self.temporary.cleanup()

    def log(self, name):
        return list(read_jsonl(self.path / name))

    def test_blank_start_and_validation(self):
        self.assertEqual(self.log('observers.jsonl'), [])
        self.assertEqual(self.log('divergence.jsonl'), [])
        self.assertTrue(self.manifest['completed'])
        for record in self.log('records.jsonl'):
            for field in ('human_interest','ai_interest','compressibility','divergence_note'):
                self.assertIsNone(record[field])
        result = validate_dataset(self.path)
        self.assertTrue(result['valid'], result['errors'])

    def test_explicit_disagreement_and_resolution(self):
        observe(self.path, self.entry, 'human_interest', True, 'human', 'explicit test report')
        self.assertEqual(self.log('divergence.jsonl'), [])
        observe(self.path, self.entry, 'ai_interest', False, 'test AI', 'explicit test report')
        self.assertEqual(self.log('divergence.jsonl')[0]['event'], 'opened')
        self.assertEqual(self.log('divergence.jsonl')[0]['status'], 'DERIVED')
        observe(self.path, self.entry, 'ai_interest', False, 'test AI', 'repeat kept')
        self.assertEqual(len(self.log('divergence.jsonl')), 1)
        observe(self.path, self.entry, 'ai_interest', None, 'test AI', 'withdrawn')
        self.assertEqual([x['event'] for x in self.log('divergence.jsonl')], ['opened','resolved'])
        self.assertEqual(len(get_annotation_history(self.path,self.entry)), 4)
        self.assertTrue(show_record(self.path,self.entry)['human_interest'])
        self.assertIsNone(show_record(self.path,self.entry)['divergence_note'])
        self.assertIsNone(self.log('records.jsonl')[0]['human_interest'])
        result = validate_dataset(self.path)
        self.assertTrue(result['valid'], result['errors'])

    def test_compressibility_requires_method(self):
        for value in (True, float('nan'), float('inf'), -1, 0.4):
            with self.assertRaises(ValueError):
                observe(self.path,self.entry,'compressibility',value,'test','reason')
        observe(self.path,self.entry,'compressibility',1.2,'test','explicit supplied value',
                'example method, ratio of compressed to raw byte count')
        self.assertEqual(show_record(self.path,self.entry)['compressibility'],1.2)
        self.assertEqual(self.log('divergence.jsonl'),[])
        observe(self.path,self.entry,'divergence_note','my note','test','explicit note')
        self.assertEqual(show_record(self.path,self.entry)['divergence_note'],'my note')

    def test_bad_interest_and_unknown_reference(self):
        for value in (0,1,'true'):
            with self.assertRaises(ValueError):
                observe(self.path,self.entry,'human_interest',value,'test','reason')
        with self.assertRaises(ValueError):
            add_note(self.path,'PROPOSED','try another projection',['missing'])
        note = add_note(self.path,'PROPOSED','try another projection',[self.entry])
        self.assertEqual(note['status'],'PROPOSED')
        result = validate_dataset(self.path)
        self.assertTrue(result['valid'], result['errors'])

    def test_refuse_overwrite_and_lock_contention(self):
        with self.assertRaises(DatasetError):
            create_dataset(self.path,{},[],0,0)
        with dataset_lock(self.path):
            with self.assertRaises(DatasetError):
                observe(self.path,self.entry,'human_interest',True,'test','reason')
        self.assertFalse((self.path/'.write.lock').exists())

    def test_failed_generation_not_complete(self):
        path = self.root/'failed'
        def fail():
            yield next(iter_geometry_runs({'experiment':'trajectories','dimensions':[2],'steps':0}))
            raise RuntimeError('deliberate test failure')
        with self.assertRaises(RuntimeError):
            create_dataset(path,{},fail(),2,2)
        import json
        manifest = json.loads((path/'manifest.json').read_text())
        self.assertFalse(manifest['completed'])
        self.assertEqual(manifest['state'],'failed')
        records = list(read_jsonl(path/'records.jsonl'))
        self.assertEqual(records[-1]['kind'],'error')
        self.assertEqual(records[-1]['status'],'OBSERVED')

    def test_tampering_detected(self):
        with (self.path/'divergence.jsonl').open('a') as stream:
            stream.write('{"id":"fake","status":"DERIVED"}\n')
        result = validate_dataset(self.path)
        self.assertFalse(result['valid'])
        self.assertTrue(result['errors'])


class SimpleRuleTests(unittest.TestCase):
    def test_rules_and_repetitions_retained(self):
        self.assertEqual(len(list(iter_binary_words(2))),6)
        self.assertEqual(len(list(iter_substitution_rules(2))),36)
        self.assertEqual(evolve_ca([0,1,0,0],0),[0,0,0,0])
        self.assertEqual(list(iter_ca_trajectory(0,[0,0,0],3)),[[0,0,0]]*4)
        self.assertEqual(list(iter_integer_trajectory(1,0,8,3,4)),[3]*5)

    def test_finite_counts_and_duplicate_seeds(self):
        config = EnumerationConfig(family='substitution',max_word_length=1,
                                   substitution_seeds=('0','0'),steps=2)
        records = list(iter_records(config))
        self.assertEqual(len(records),8)
        self.assertEqual(len({r['id'] for r in records}),8)
        self.assertEqual(sum(len(r['outputs']) for r in records),24)
        self.assertEqual(records[0]['outputs'],records[1]['outputs'])
        with tempfile.TemporaryDirectory() as root:
            path = Path(root)/'rules'
            manifest = enumerate_dataset(path,config)
            self.assertEqual(manifest['total_runs'],8)
            result = validate_dataset(path)
            self.assertTrue(result['valid'], result['errors'])


if __name__ == '__main__':
    unittest.main()
