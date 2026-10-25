import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('confirmation_wrapper', HERE / 'array_worker.py')
worker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(worker)


class BindingTests(unittest.TestCase):
    def prepare(self, root, name='source.py'):
        (root / 'source.py').write_text('unchanged science\n')
        obj = dict(schema_version=1, study=HERE.name,
                   files={name: worker.digest(root / 'source.py')})
        obj['identity_sha256'] = hashlib.sha256(worker.canonical(obj).encode()).hexdigest()
        path = root / 'freeze.json'
        path.write_text(json.dumps(obj))
        return path

    def test_valid_freeze_and_source_change_rejection(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t); path = self.prepare(root)
            self.assertEqual(worker.verify_freeze(path, root)['study'], HERE.name)
            (root / 'source.py').write_text('changed science\n')
            with self.assertRaisesRegex(ValueError, 'Frozen file changed'):
                worker.verify_freeze(path, root)

    def test_manifest_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t); path = self.prepare(root)
            obj = json.loads(path.read_text()); obj['study'] = 'other'
            path.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError, 'Invalid freeze identity'):
                worker.verify_freeze(path, root)

    def test_escape_path_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t); path = self.prepare(root, '../source.py')
            with self.assertRaisesRegex(ValueError, 'Unsafe source path'):
                worker.verify_freeze(path, root)

    def test_symlink_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t); path = self.prepare(root)
            (root / 'source.py').rename(root / 'actual.py')
            (root / 'source.py').symlink_to(root / 'actual.py')
            with self.assertRaisesRegex(ValueError, 'Frozen file changed'):
                worker.verify_freeze(path, root)


if __name__ == '__main__':
    unittest.main()
