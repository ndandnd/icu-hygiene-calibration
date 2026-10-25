"""Prospective launcher tests with fake scheduler responses; submits no jobs."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import launch as launcher
import submit_once as once


class SubmissionTests(unittest.TestCase):
    def setup_folder(self,td):
        root=Path(td).resolve();script=root/'fixture.sbatch';script.write_text('#!/bin/bash\nexit 0\n')
        return root,['--time=00:05:00',str(script)]

    def test_success_is_idempotent_and_receipt_gap_is_repaired_without_resubmit(self):
        with tempfile.TemporaryDirectory() as td:
            root,args=self.setup_folder(td)
            with patch.object(once.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='12345\n',stderr='')) as call:
                self.assertEqual(once.submit(root,'probe',args,'/scheduler/sbatch'),'12345')
                self.assertEqual(once.submit(root,'probe',args,'/scheduler/sbatch'),'12345')
                (root/'submissions/probe.jobid').unlink()
                self.assertEqual(once.submit(root,'probe',args,'/scheduler/sbatch'),'12345')
                self.assertEqual(call.call_count,1)
            self.assertEqual(json.loads((root/'submissions/probe.intent.json').read_text())['state'],'submitted')

    def test_unknown_outcome_never_retries(self):
        with tempfile.TemporaryDirectory() as td:
            root,args=self.setup_folder(td)
            with patch.object(once.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='unexpected output',stderr='')) as call:
                with self.assertRaisesRegex(ValueError,'No confirmed'):once.submit(root,'grid',args,'/scheduler/sbatch')
                with self.assertRaisesRegex(ValueError,'Unknown prior'):once.submit(root,'grid',args,'/scheduler/sbatch')
                self.assertEqual(call.call_count,1)
            self.assertFalse((root/'submissions/grid.jobid').exists())

    def test_timeout_preserves_unknown_intent_and_no_retry(self):
        with tempfile.TemporaryDirectory() as td:
            root,args=self.setup_folder(td)
            with patch.object(once.subprocess,'run',side_effect=subprocess.TimeoutExpired('fake',60)) as call:
                with self.assertRaises(subprocess.TimeoutExpired):once.submit(root,'selection',args,'/scheduler/sbatch')
                with self.assertRaisesRegex(ValueError,'Unknown prior'):once.submit(root,'selection',args,'/scheduler/sbatch')
                self.assertEqual(call.call_count,1)
            self.assertTrue((root/'submissions/selection.submission-error.json').exists())

    def test_declaration_drift_and_unapproved_names_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root,args=self.setup_folder(td)
            with patch.object(once.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='12345;cluster\n',stderr='')):
                once.submit(root,'policy',args,'/scheduler/sbatch')
                with self.assertRaisesRegex(ValueError,'mismatch'):once.submit(root,'policy',['--time=99:00:00',args[-1]],'/scheduler/sbatch')
            with self.assertRaisesRegex(ValueError,'Unapproved'):once.submit(root,'unrelated',args,'/scheduler/sbatch')


class LauncherTests(unittest.TestCase):
    def cfg(self):
        return dict(python='/venv/bin/python',sbatch='/scheduler/sbatch',scontrol='/scheduler/scontrol',squeue='/scheduler/squeue',
            account='fixture',partition='fixture',nodelist='node-02,node-03',executables=dict(squeue=dict(sha256='fixture')))

    def test_fixed_resource_ceiling_dependencies_and_priority(self):
        with patch.object(launcher,'guard'),patch.object(launcher,'submit',return_value='123') as submit:
            for name in ('probe','grid','selection','policy','analysis'):
                launcher.dispatch(Path('/private/study'),self.cfg(),name,'99')
                call=submit.call_args.args;args=call[2]
                self.assertIn('--nice=10000',args);self.assertIn('--dependency=afterok:99',args)
                self.assertIn('--cpus-per-task=1',args);self.assertIn('--mem=1G',args)
                if name in ('grid','policy'):
                    self.assertIn('--array=0-63%16',args);self.assertIn('--time=03:00:00',args)
                    self.assertEqual(args[-1],str(HERE/'array.sbatch'))
                else:
                    self.assertIn('--no-requeue',args)
                    self.assertIn('--time='+('00:05:00' if name=='probe' else '00:30:00'),args)
                export=next(x for x in args if x.startswith('--export='))
                self.assertNotIn('ALL',export);self.assertIn('ROBUST_PYTHON=',export)

    def test_priority_guard_blocks_only_active_opte_names(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td).resolve()
            with patch.object(launcher.subprocess,'run',return_value=SimpleNamespace(stdout='10|v2g-held|PENDING\n')):
                launcher.priority_guard(root,self.cfg())
            with patch.object(launcher.subprocess,'run',return_value=SimpleNamespace(stdout='10|OPTE-run|RUNNING\n')):
                with self.assertRaisesRegex(ValueError,'OPTE'):launcher.priority_guard(root,self.cfg())
            self.assertEqual(len(list((root/'priority_checks').glob('*.json'))),2)

    def test_executable_symlink_is_bound_without_changing_entry_path(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td).resolve();binary=root/'realpython';binary.write_text('fixture');binary.chmod(0o700)
            entry=root/'python';entry.symlink_to(binary)
            bound=launcher.executable(entry)
            self.assertEqual(bound['declared'],str(entry));self.assertEqual(bound['resolved'],str(binary))
            with self.assertRaisesRegex(ValueError,'Symlink'):launcher.path(entry)
            binary.write_text('changed');self.assertNotEqual(bound['sha256'],launcher.executable(entry)['sha256'])

    def test_export_rejects_ambiguous_values(self):
        with self.assertRaises(ValueError):launcher.exported({'ROBUST_ROOT':'bad,extra=value'})
        with self.assertRaises(ValueError):launcher.exported({'ROBUST_ROOT':'bad\nvalue'})

    def test_fast_auxiliary_start_waits_for_durable_matching_receipt(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td).resolve();folder=root/'submissions';folder.mkdir()
            (folder/'probe.intent.json').write_text(json.dumps(dict(state='pending_unknown_until_receipted')))
            def publish(unused):
                (folder/'probe.intent.json').write_text(json.dumps(dict(state='submitted',job_id='123')))
                (folder/'probe.jobid').write_text('123\n')
            with patch.object(launcher.time,'sleep',side_effect=publish) as sleep:
                launcher.await_stage_receipt(root,'probe','123')
                self.assertEqual(sleep.call_count,1)
            with self.assertRaisesRegex(ValueError,'differs'):launcher.await_stage_receipt(root,'probe','124')
            with self.assertRaisesRegex(ValueError,'no resubmission'):launcher.await_stage_receipt(root,'analysis','999',timeout=0)

    def test_receipt_appearing_after_stale_intent_read_is_rechecked(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td).resolve();folder=root/'submissions';folder.mkdir()
            intent=folder/'probe.intent.json';receipt=folder/'probe.jobid'
            pending=dict(state='pending_unknown_until_receipted')
            intent.write_text(json.dumps(pending));original_read=launcher.read;calls=[]
            def interleaved_read(p):
                calls.append(p)
                if len(calls)==1:
                    # The consumer has read pending; the submitter publishes both
                    # durable files before the consumer checks receipt.exists().
                    intent.write_text(json.dumps(dict(state='submitted',job_id='123')))
                    receipt.write_text('123\n')
                    return pending
                return original_read(p)
            with patch.object(launcher,'read',side_effect=interleaved_read),patch.object(launcher.time,'sleep') as sleep:
                launcher.await_stage_receipt(root,'probe','123')
                self.assertGreaterEqual(len(calls),2);sleep.assert_not_called()


if __name__=='__main__':unittest.main()
