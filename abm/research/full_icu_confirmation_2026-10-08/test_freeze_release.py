"""Synthetic release/archive fixtures only; no existing result banks accessed."""
import importlib.util
import ast
import contextlib
import io
import json
import os
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_confirmation_release_tests',HERE/'freeze_release.py')
release=importlib.util.module_from_spec(spec);spec.loader.exec_module(release)


def rehash(value):
    value.pop('identity',None)
    value['identity']=release.digest(release.canonical(value).encode())
    return value


class ArchiveTests(unittest.TestCase):
    def simple(self,base,n=2):
        root=base/'source';root.mkdir()
        names=[]
        for i in range(n):
            name=release.PREFIX+f'/fixture_{i:03d}.json'
            path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('{"synthetic":true}\n')
            names.append(name)
        output=base/'release'
        manifest=release.preserve(root,dict(test_fixture=True),names,[],output)
        return root,output,manifest

    def test_safe_canonical_paths_and_public_whitelist(self):
        for name in ('../x','/x','a/../x','a//x','.',r'a\x','C:/x','a space/x'):
            with self.subTest(name=name),self.assertRaises(ValueError):release.relative(name)
        for name in ('clinical/patients.csv',release.PREFIX+'/logs/x.json',release.PREFIX+'/secrets/x.txt',
                     release.PREFIX+'/.env',release.PREFIX+'/source.exe',release.RESULTS+'/calibration_local/completed/1.json'):
            with self.subTest(name=name),self.assertRaises(ValueError):release.public_name(name)
        self.assertEqual(release.public_name('src/agents.py'),'src/agents.py')

    def test_native_renderer_whitelist_and_byte_exact_restore(self):
        self.assertEqual(release.public_name(release.RENDERER),release.RENDERER)
        with self.assertRaises(ValueError):
            release.public_name(str(Path(release.RENDERER).parent/'unreviewed.py'))
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp).resolve();root=base/'source';source=root/release.RENDERER
            source.parent.mkdir(parents=True);source.write_bytes(b'# synthetic native renderer fixture\n')
            release.preserve(root,{},[release.RENDERER],[],base/'release')
            release.restore(base/'release',base/'restored')
            self.assertEqual((base/'restored'/release.RENDERER).read_bytes(),source.read_bytes())

    def test_aggregate_freeze_checked_dependency_closure(self):
        # Isolate packaging closure from scientific loaders; their complete-job
        # checks have independent fixture tests. Freeze hashes stay real here.
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve()
            files=sorted(release.REQUIRED_FILES)
            for name in files:
                path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('public fixture\n')
            dependency=release.PREFIX+'/bridge/AGGREGATION_REVIEW.md'
            (root/dependency).write_text('synthetic pre-outcome review\n')
            def make_freeze(name,dependencies):
                value=dict(study=release.STUDY,files={dep:release.digest((root/dep).read_bytes()) for dep in dependencies})
                value['identity_sha256']=release.digest(release.canonical(value).encode())
                (root/name).write_text(json.dumps(value))
            make_freeze(release.AGGREGATE_FREEZE,[release.PREFIX+'/bridge/aggregate.py',dependency])
            make_freeze(release.PREFIX+'/analyze_freeze.json',[release.PREFIX+'/analyze.py'])
            freeze_name=release.PREFIX+'/fixture_freeze.json';make_freeze(freeze_name,[])
            groups=[]
            for module in release.COUNTS:
                report=release.PREFIX+f'/{module}_fixture_report.json';(root/report).write_text('{}')
                groups.append(dict(module=module,design=release.PREFIX+'/fixture_design.json',banks=[],report=report,summaries=[],freeze=freeze_name))
            specification=dict(schema_version=1,groups=groups,files=files)
            with patch.object(release,'report_receipts',return_value=({},set())), \
                 patch.object(release,'validate_positive',side_effect=lambda root,group,*args:release.COUNTS[group['module']]), \
                 patch.object(release,'validate_bridge',return_value=256), \
                 patch.object(release,'ASSEMBLY_GUARDS',{}):
                selected,_=release.validate_spec(root,specification)
                self.assertIn(dependency,selected);self.assertIn(release.AGGREGATE_FREEZE,selected)
                without_preview=dict(specification,files=[name for name in files if name!=release.PREFIX+'/manuscript/preconfirmation_pdf/main.pdf'])
                with self.assertRaisesRegex(ValueError,'Explicit analyzer/release/assembly'):
                    release.validate_spec(root,without_preview)
                (root/dependency).write_text('changed review\n')
                with self.assertRaisesRegex(ValueError,'Frozen file changed'):
                    release.validate_spec(root,specification)

    def test_assembler_static_guards_and_finite_old_whitelist(self):
        def constants(path):
            tree=ast.parse(path.read_text())
            return {item.targets[0].id:item.value.value for item in tree.body if isinstance(item,ast.Assign)
                and len(item.targets)==1 and isinstance(item.targets[0],ast.Name) and isinstance(item.value,ast.Constant)}
        manuscript=constants(HERE/'manuscript/build_confirmation.py')
        guide=constants(HERE/'guide/build_confirmation.py')
        self.assertEqual(release.ASSEMBLY_GUARDS[release.MANUSCRIPT_BASE],manuscript['BASE_SHA'])
        self.assertEqual(release.ASSEMBLY_GUARDS[release.CONTACT_CLOCK],guide['CLOCK_SHA'])
        self.assertEqual(release.ASSEMBLY_GUARDS[release.PREFIX+'/manuscript/main.tex'],manuscript['MAIN_PRECONFIRMATION_SHA'])
        self.assertEqual(release.ASSEMBLY_GUARDS[release.PREFIX+'/manuscript/certificate_note.tex'],manuscript['NOTE_SHA'])
        self.assertEqual(release.ASSEMBLY_GUARDS[release.PREFIX+'/guide/preconfirmation/plain_english_guide.pdf'],guide['PRECONFIRMATION_PDF_SHA'])
        for name in release.ASSEMBLY_FILES:self.assertEqual(release.public_name(name),name)
        for allowed in (release.MANUSCRIPT_BASE,release.CONTACT_CLOCK):
            with self.assertRaises(ValueError):release.public_name(str(Path(allowed).parent/'unreviewed.json'))

    def test_assembly_sources_restore_without_private_input_paths(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp).resolve();root=base/'source';root.mkdir()
            for name in release.ASSEMBLY_FILES:
                path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
                # Only public Python source is read; the protected PDFs/TeX,
                # clock and layout are synthetic bytes in this archive fixture.
                path.write_bytes((HERE.parents[1]/name).read_bytes() if name.endswith('.py') else b'public synthetic fixture\n')
            with patch.object(release,'ASSEMBLY_GUARDS',{}):release.validate_assembly_dependencies(root)
            release.preserve(root,{},release.ASSEMBLY_FILES,[],base/'release')
            release.restore(base/'release',base/'restored')
            restored=base/'restored'
            guide=release.load('_restored_guide_source',restored/release.PREFIX/'guide/build_confirmation.py')
            manuscript=guide.assembly_module()
            editorial=manuscript.load('_restored_editorial_helpers',restored/release.PREFIX/'manuscript/build_editorial.py')
            self.assertEqual(manuscript.BASE,restored/release.MANUSCRIPT_BASE)
            self.assertEqual(guide.CLOCK,restored/release.CONTACT_CLOCK)
            self.assertEqual(editorial.span('first protected last','first','last'),'first protected ')
            with patch.object(release,'ASSEMBLY_GUARDS',{release.MANUSCRIPT_BASE:release.digest(b'public synthetic fixture\n')}):
                release.validate_assembly_dependencies(restored)
                (restored/release.MANUSCRIPT_BASE).write_text('changed')
                with self.assertRaisesRegex(ValueError,'Preserved assembler input changed'):
                    release.validate_assembly_dependencies(restored)

    def test_main_and_supporting_fragment_restore_and_regenerate_exactly(self):
        # Public frozen source plus wholly synthetic in-memory report fields.
        # No result directory or report file from the workspace is opened.
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp).resolve();root=base/'source';root.mkdir()
            source_names=[release.MANUSCRIPT_BASE]+[release.PREFIX+'/manuscript/'+name for name in
                ('build_confirmation.py','build_editorial.py','main.tex','certificate_note.tex')]
            for name in source_names:
                path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes((HERE.parents[1]/name).read_bytes())
            m=release.load('_synthetic_fragment_original',root/release.PREFIX/'manuscript/build_confirmation.py')
            folder=root/release.PREFIX/'fixture_figures';folder.mkdir()
            for name in (*m.FIGURES,'leading_rate_comparison'):
                (folder/(name+'.tex')).write_text('% synthetic figure/table fixture only\n')
            def uncertainty(value):return dict(interval95=[value-.01,value+.01],intervals={format(1-.05/3,'.12g'):[value-.02,value+.02]})
            rows={profile+'_target_candidate':dict(scenario=dict(p=.1),
                arms={arm:dict(mean_acquired=count,rate_per_1000_patient_days=rate,rate_uncertainty=uncertainty(rate))
                    for arm,count,rate in (('baseline',7.,3.),('policy',3.,2.))},
                mean_prevented=4.,mean_prevented_uncertainty=uncertainty(4.),absolute_prevented_rate=1.,
                absolute_prevented_rate_uncertainty=uncertainty(1.),rate_prevented_fraction=.123,
                rate_prevented_fraction_uncertainty=uncertainty(.123)) for profile in m.PROFILES}
            bridge=dict(profiles={m.BK[profile]:dict(K0=dict(mean=100.,normal_approximation_95_interval=[90.,110.]),
                K1=dict(mean=50.,normal_approximation_95_interval=[45.,55.]),limiting_response=.5,
                limiting_response_normal_95_interval=[.4,.6],stationary_q2_response_for_comparison_only=.75) for profile in m.PROFILES})
            environment=dict(primary_environment_effects=[dict(p=p,full_minus_disabled={name:dict(estimate=value,uncertainty=uncertainty(value))
                for name,value in (('F',-.2),('prevented_rate',.1),('mean_prevented',.3))}) for p in (.08,.11,.135)])
            context=dict(holdout=dict(baseline_matching=dict(established=False)),holdout_rows=rows,environment=environment,
                environment_rows={},bridge=bridge,figure_directory=folder,report_hashes=dict(synthetic_fixture=True),figure_manifest_sha256='f'*64)
            prose_name=release.PREFIX+'/manuscript/confirmation_prose.tex'
            (root/prose_name).write_text('\\begin{abstract}Synthetic archive mechanics only.\\end{abstract}\n'+'\n'.join(m.ANCHORS))
            placeholders=[]
            for name in ('holdout-report','environment-report','bridge-report','figure-manifest'):
                path=folder/(name+'.json');path.write_text('{"synthetic_fixture":true}\n');placeholders.extend(['--'+name,str(path)])
            original=m.HERE/'confirmation'/'fixture_original'
            with patch.object(m,'validate_inputs',return_value=context),contextlib.redirect_stdout(io.StringIO()):
                m.main([*placeholders,'--prose',str(root/prose_name),'--output',str(original)])
            outputs={str(path.relative_to(root)) for path in original.iterdir()}
            release.validate_final_manuscript_outputs(root,outputs)
            for missing in ('main.tex','antecedent_results.tex','manifest.json'):
                with self.subTest(missing=missing),self.assertRaisesRegex(ValueError,'requires main.tex'):
                    release.validate_final_manuscript_outputs(root,outputs-{str((original/missing).relative_to(root))})
            selected=set(source_names)|{prose_name}|outputs|{str(path.relative_to(root)) for path in folder.iterdir()}
            release.preserve(root,{},selected,[],base/'release');release.restore(base/'release',base/'restored')
            restored=base/'restored';release.validate_final_manuscript_outputs(restored,outputs)
            regenerated=release.load('_synthetic_fragment_restored',restored/release.PREFIX/'manuscript/build_confirmation.py')
            context=dict(context,figure_directory=restored/str(folder.relative_to(root)))
            document,receipt=regenerated.assemble((restored/prose_name).read_text(),context)
            editorial=regenerated.load('_synthetic_fragment_helpers',regenerated.HERE/'build_editorial.py')
            fragment,_=regenerated.antecedent_source(regenerated.BASE.read_text(),editorial)
            original_folder=restored/str(original.relative_to(root))
            self.assertEqual(document.encode(),(original_folder/'main.tex').read_bytes())
            self.assertEqual(fragment.encode(),(original_folder/'antecedent_results.tex').read_bytes())
            self.assertEqual(receipt['supporting_sources']['antecedent_results.tex']['sha256'],release.digest(fragment.encode()))
            (original_folder/'antecedent_results.tex').write_text('changed fixture')
            with self.assertRaisesRegex(ValueError,'bytehash binding differs'):
                release.validate_final_manuscript_outputs(restored,outputs)

    def test_exact_deterministic_roundtrip_and_manifest_identity(self):
        long_name=release.PREFIX+'/figures/'+('a'*100)+'.png'
        payload={long_name:b'public figure'}
        self.assertEqual(release.decode_chunk(release.encode_chunk(payload),
            {long_name:dict(size=13,sha256=release.digest(b'public figure'))}),payload)
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp).resolve();root,output,manifest=self.simple(base)
            second=base/'release2'
            release.preserve(root,manifest['specification'],manifest['files'],[],second)
            for path in output.iterdir():self.assertEqual(path.read_bytes(),(second/path.name).read_bytes())
            restored=base/'restored';release.restore(output,restored)
            for name in manifest['files']:self.assertEqual((root/name).read_bytes(),(restored/name).read_bytes())
            changed=dict(manifest,scope='changed');(output/'release_manifest.json').write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError,'canonical content hash'):release.restore(output,base/'bad')
            self.assertFalse((base/'bad').exists())

    def test_member_and_raw_limits_adaptive_compressed_split(self):
        with self.assertRaises(ValueError):release.encode_chunk({f'f{i}':b'x' for i in range(65)})
        with patch.object(release,'MAX_RAW_CHUNK_BYTES',1):
            with self.assertRaises(ValueError):release.encode_chunk({'data':b'xx'})
            with self.assertRaises(ValueError):release.decode_chunk(b'not opened',{'data':dict(size=2,sha256='0'*64)})
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp).resolve();root=base/'source';root.mkdir();names=[]
            for i in range(4):
                name=release.PREFIX+f'/fixture_{i}.png';path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
                path.write_bytes(os.urandom(2000));names.append(name)
            with patch.object(release,'MAX_COMPRESSED_BYTES',5000):
                manifest=release.preserve(root,{},names,[],base/'release')
                self.assertGreater(len(manifest['chunks']),1)
                self.assertTrue(all(chunk['size']<=5000 for chunk in manifest['chunks']))
                release.restore(base/'release',base/'restored')

    def test_corrupt_chunk_extra_file_symlink_and_path_conflict(self):
        for kind in ('bytes','extra','symlink','conflict','duplicate'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as temp:
                base=Path(temp).resolve();_,output,manifest=self.simple(base)
                chunk=output/manifest['chunks'][0]['file']
                if kind=='bytes':chunk.write_bytes(chunk.read_bytes()+b'changed')
                elif kind=='extra':(output/'unexpected.log').write_text('not public')
                elif kind=='symlink':
                    backup=base/'copy';chunk.rename(backup);chunk.symlink_to(backup)
                else:
                    if kind=='conflict':
                        name=next(iter(manifest['files']));manifest['files'][name+'/child']=dict(manifest['files'][name])
                    else:manifest['chunks'][0]['members']*=2
                    (output/'release_manifest.json').write_text(json.dumps(rehash(manifest)))
                with self.assertRaises(ValueError):release.restore(output,base/'bad')
                self.assertFalse((base/'bad').exists())

    def test_tar_links_duplicate_and_traversal_refused(self):
        for kind in ('link','duplicate','traversal'):
            stream=io.BytesIO();name='../escape' if kind=='traversal' else 'data.json'
            with tarfile.open(fileobj=stream,mode='w:gz') as archive:
                member=tarfile.TarInfo(name);member.size=1
                if kind=='link':member.type=tarfile.SYMTYPE;member.linkname='elsewhere';archive.addfile(member)
                else:archive.addfile(member,io.BytesIO(b'x'))
                if kind=='duplicate':archive.addfile(member,io.BytesIO(b'x'))
            with self.assertRaises(ValueError):release.decode_chunk(stream.getvalue(),{name:dict(size=1,sha256=release.digest(b'x'))})

    def test_privacy_markers_and_symlink_ancestors(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();name=release.PREFIX+'/data.json';path=root/name;path.parent.mkdir(parents=True)
            path.write_bytes(b'/Users/'+b'example/private')
            with self.assertRaisesRegex(ValueError,'Private-path'):release.read(root,name)
            path.write_text('public');(root/'alias').symlink_to(root,target_is_directory=True)
            with self.assertRaisesRegex(ValueError,'Symlink'):release.read(root/'alias',name)
            source=root/release.PREFIX/'markers.py';source.write_bytes(b'MARKERS='+repr(release.FORBIDDEN).encode())
            release.read(root,str(source.relative_to(root)),allow_marker_literals=True)

    def test_only_all_four_complete_modules_can_publish(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp).resolve();root=base/'source';root.mkdir()
            with self.assertRaisesRegex(ValueError,'all four'):
                release.freeze(root,dict(schema_version=1,groups=[],files=[]),base/'release')
            self.assertFalse((base/'release').exists())
            with patch.object(release,'validate_spec',side_effect=ValueError('Incomplete module report')):
                with self.assertRaisesRegex(ValueError,'Incomplete'):
                    release.freeze(root,{},base/'release')
            self.assertFalse((base/'release').exists())


class ReportReceiptTests(unittest.TestCase):
    def test_bridge_aggregate_binding_required_and_exact(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();report_name=release.PREFIX+'/bridge_fixture/summary.json'
            binding_name=str(Path(report_name).parent/'binding.json')
            directory=(root/report_name).parent;directory.mkdir(parents=True)
            report=dict(config_sha256='a'*64,seeds=[1,2],freeze_sha256='b'*64,freeze_identity='c'*64,
                scientific_sources={'src/agents.py':'d'*64},source_sha256={'src/agents.py':'d'*64},
                analysis_code_sha256='e'*64,config_path=release.PREFIX+'/bridge_config.json')
            (root/report_name).write_text(json.dumps(report))
            binding=dict(schema_version=1,module='bridge',status='complete',n_backgrounds=256,
                summary_sha256=release.digest((root/report_name).read_bytes()),summary_filename='summary.json',**report)
            (root/binding_name).write_text(json.dumps(binding))
            group=dict(report=report_name,summaries=[binding_name])
            release.validate_bridge_binding(root,group,report)
            for summaries in ([],[binding_name,binding_name]):
                with self.assertRaisesRegex(ValueError,'binding.json exactly once'):
                    release.validate_bridge_binding(root,dict(group,summaries=summaries),report)
            for key in ('summary_sha256','config_sha256','seeds','freeze_identity','scientific_sources','analysis_code_sha256'):
                changed=dict(binding);changed[key]='mismatch';(root/binding_name).write_text(json.dumps(changed))
                with self.subTest(key=key),self.assertRaisesRegex(ValueError,'binding differs'):
                    release.validate_bridge_binding(root,group,report)

    def fixture(self,root):
        bank=release.RESULTS+'/calibration_local';directory=root/bank;(directory/'completed').mkdir(parents=True)
        jobs=[dict(job_id=i,seed=100+i,scenario=dict(id='synthetic')) for i in range(2)]
        manifest=rehash(dict(jobs=jobs,resolved_config=dict(seeds=[100,101])))
        raw=json.dumps(manifest).encode();(directory/'manifest.json').write_bytes(raw)
        hashes={}
        for job in jobs:
            name=f"job_{job['job_id']:08d}.json";payload=json.dumps(dict(job,status='complete',identity=manifest['identity'])).encode()
            (directory/'completed'/name).write_bytes(payload);hashes[name]=release.digest(payload)
        freeze_name=release.PREFIX+'/fixture_freeze.json';path=root/freeze_name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('{}')
        (directory/'confirmation_binding.json').write_text(json.dumps(dict(freeze_sha256=release.digest(path.read_bytes()))))
        report=dict(module='calibration',status='complete',n_paired_jobs=2,banks=[dict(identity=manifest['identity'],manifest_sha256=release.digest(raw),jobs=2,completed_records_sha256=hashes)])
        report_name=release.PREFIX+'/fixture_report.json';(root/report_name).write_text(json.dumps(report))
        group=dict(module='calibration',banks=[bank],report=report_name,freeze=freeze_name,summaries=[])
        return group,report,directory

    def test_report_exact_payload_binding_and_no_alias_jobs(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();group,report,directory=self.fixture(root)
            validated,files=release.report_receipts(root,group,2)
            self.assertEqual(validated,report);self.assertEqual(len(files),5)
            job=next((directory/'completed').iterdir());job.write_bytes(job.read_bytes()+b' ')
            with self.assertRaisesRegex(ValueError,'validated report payload'):release.report_receipts(root,group,2)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();group,_,directory=self.fixture(root)
            next((directory/'completed').iterdir()).rename(directory/'completed'/'alias.json')
            with self.assertRaisesRegex(ValueError,'Complete exact job files'):release.report_receipts(root,group,2)

    def test_unvalidated_partial_wrong_count_and_wrong_freeze_reports_rejected(self):
        for kind in ('partial','count','manifest','binding'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as temp:
                root=Path(temp).resolve();group,report,directory=self.fixture(root)
                if kind=='partial':report['status']='partial'
                elif kind=='count':report['n_paired_jobs']=3
                elif kind=='manifest':report['banks'][0]['manifest_sha256']='0'*64
                else:(directory/'confirmation_binding.json').write_text('{"freeze_sha256":"wrong"}')
                (root/group['report']).write_text(json.dumps(report))
                with self.assertRaises(ValueError):release.report_receipts(root,group,2)


if __name__=='__main__':unittest.main()
