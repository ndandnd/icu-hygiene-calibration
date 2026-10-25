#!/usr/bin/env python3
"""Byte-exact release of four complete, validated synthetic confirmation modules.

Only explicit public files and validated manifests/bindings/records/reports are
preserved. No recursive directory/archive discovery and no simulator calls.
Archive primitives follow the earlier tested response release utility.
"""
import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys
import tarfile

STUDY = 'full_icu_confirmation_2026-10-08'
PREFIX = 'research/'+STUDY
RESULTS = 'results/'+STUDY
RENDERER = 'research/full_icu_calibration_2026-10-07/figures/build_baseline_response.py'
AGGREGATE_FREEZE = PREFIX+'/bridge/aggregate_freeze.json'
MANUSCRIPT_BASE = 'research/full_icu_response_2026-10-08/manuscript_draft/actual_04/manuscript.tex'
CONTACT_CLOCK = 'research/full_icu_calibration_2026-10-07/figures/nurse_contact_clock_native.png'
ASSEMBLY_FILES = {MANUSCRIPT_BASE,CONTACT_CLOCK, *{PREFIX+'/'+name for name in (
    'manuscript/build_confirmation.py','manuscript/build_editorial.py','manuscript/CONFIRMATION_INTERFACE.md','manuscript/confirmation_prose.tex',
    'manuscript/main.tex','manuscript/certificate_note.tex','manuscript/preconfirmation_pdf/main.pdf',
    'manuscript/preconfirmation_pdf/certificate_note.pdf',
    'guide/build_confirmation.py','guide/CONFIRMATION_INTERFACE.md','guide/content_confirmation.json',
    'guide/preconfirmation/plain_english_guide.pdf','figures/build_figures.py')}}
ASSEMBLY_GUARDS = {
    MANUSCRIPT_BASE:'6b066e5c6a3ff7529d66288b6749a44068978769f15492c72438fb3875272743',
    CONTACT_CLOCK:'0f5e4aae9c9ebc3485d037d582787511bed1a76757f2e8dd559eb90353380b44',
    PREFIX+'/manuscript/main.tex':'a49840a6113343c97458c3e1227d2a824b44778c772f2b0acc5e4001a18671de',
    PREFIX+'/manuscript/certificate_note.tex':'abb1588423a3e3455786c7c67dc6e3e6824a211c178fdc22dc696cba9917e2d8',
    PREFIX+'/guide/preconfirmation/plain_english_guide.pdf':'16921f6291a1ff8da0d141d796207ef784e8a6ae6747f8636420fbefb870c803'}
REQUIRED_FILES = ASSEMBLY_FILES|{PREFIX+'/analyze.py',PREFIX+'/analyze_freeze.json',PREFIX+'/freeze_release.py',
    PREFIX+'/test_freeze_release.py',PREFIX+'/REPRODUCIBILITY_DRAFT.md',PREFIX+'/bridge/aggregate.py',
    PREFIX+'/bridge/test_aggregate.py',AGGREGATE_FREEZE,RENDERER,'requirements-lock.txt'}
COUNTS = dict(calibration=384,holdout=1536,environment=3072,bridge=256)
MAX_MEMBERS = 64
MAX_COMPRESSED_BYTES = 19*1024*1024  # Strictly below 20 MiB, including tar/gzip.
MAX_RAW_CHUNK_BYTES = 32*1024*1024
FORBIDDEN = (b'/Users/',b'/home/',b'/scratch/',b'PRIVATE KEY',b'Bearer ',b'Authorization:')
BAD_PARTS = {'.git','.aws','.codex','.agents','.venv','__pycache__','logs','attempts','locks','progress','credentials','secrets','private','clinical_inputs','clinical_data'}
EXTENSIONS = {'.py','.md','.tex','.json','.csv','.png','.svg','.pdf','.sbatch','.txt','.toml','.lock'}
LEGACY_FILES = {'src/CPEmodel.py','src/agents.py','src/routes.py',
    RENDERER,MANUSCRIPT_BASE,CONTACT_CLOCK,
    'requirements.txt','requirements-lock.txt','requirements-research.txt','requirements-guide.txt',
    'research/full_icu_calibration_2026-10-07/run.py',
    *{'research/full_icu_response_2026-10-08/'+name for name in ('run.py','instrumented_model.py','array_worker.py','PROTOCOL.md')}}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def relative(value):
    if not isinstance(value,str):
        raise ValueError('Canonical relative path required')
    path = PurePosixPath(value)
    if not value or not path.parts or path.is_absolute() or str(path) != value or any(
        part in ('.','..') or not re.fullmatch(r'[A-Za-z0-9_.-]+',part) for part in path.parts):
        raise ValueError('Canonical relative path required')
    return value


def public_name(name,managed=False):
    relative(name)
    path = PurePosixPath(name)
    if set(path.parts)&BAD_PARTS or any(part.startswith('.') for part in path.parts) or path.suffix.lower() not in EXTENSIONS:
        raise ValueError('Unrecognized or excluded public file: '+name)
    if name not in LEGACY_FILES and not name.startswith(PREFIX+'/') and not (managed and name.startswith(RESULTS+'/')):
        raise ValueError('File outside explicit public scientific namespaces: '+name)
    return name


def safe_path(root,name=None):
    path = Path(root).absolute()
    for parent in reversed((path,*path.parents)):
        if parent.is_symlink():
            raise ValueError('Symlink root or parent refused')
    if name is not None:
        relative(name)
        for part in PurePosixPath(name).parts:
            path = path/part
            if path.is_symlink():
                raise ValueError('Symlink input refused')
    return path


def privacy_bytes(data,allow_marker_literals=False):
    if allow_marker_literals:
        for marker in FORBIDDEN:
            expression = rb'(?<![A-Za-z0-9_])(?:br|rb|b|r)?([\x27\x22])'+re.escape(marker)+rb'\1'
            data = re.sub(expression,b'PUBLIC_MARKER_LITERAL',data)
    return data


def read(root,name,privacycheck=True,allow_marker_literals=False):
    path = safe_path(root,relative(name))
    if not path.is_file():
        raise ValueError('Regular input file required: '+name)
    if not privacycheck and path.stat().st_size > MAX_COMPRESSED_BYTES:
        raise ValueError('Compressed chunk exceeds below-20-MiB limit')
    data = path.read_bytes()
    if privacycheck and any(marker in privacy_bytes(data,allow_marker_literals) for marker in FORBIDDEN):
        raise ValueError('Private-path or credential marker in '+name)
    return data


def checked_manifest(data):
    manifest = json.loads(data)
    identity_key = 'identity' if 'identity' in manifest else 'identity_sha256'
    if digest(canonical({key:value for key,value in manifest.items() if key != identity_key}).encode()) != manifest.get(identity_key):
        raise ValueError('Manifest identity is not its canonical content hash')
    return manifest


def encode_chunk(payloads):
    if len(payloads)>MAX_MEMBERS or sum(map(len,payloads.values()))>MAX_RAW_CHUNK_BYTES:
        raise ValueError('Chunk exceeds bounded member/raw-byte limit')
    stream = io.BytesIO()
    with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0) as compressed:
        with tarfile.open(fileobj=compressed,mode='w|',format=tarfile.GNU_FORMAT) as archive:
            for name,data in sorted(payloads.items()):
                relative(name)
                entry = tarfile.TarInfo(name);entry.size=len(data);entry.mode=0o644
                archive.addfile(entry,io.BytesIO(data))
    return stream.getvalue()


def decode_chunk(data,expected):
    if len(data)>MAX_COMPRESSED_BYTES or len(expected)>MAX_MEMBERS or sum(value['size'] for value in expected.values())>MAX_RAW_CHUNK_BYTES:
        raise ValueError('Chunk exceeds bounded compressed/member/raw-byte limit')
    decoded = {}
    with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as archive:
        for member in archive:
            name = relative(member.name)
            if not member.isfile() or member.pax_headers or name not in expected or name in decoded:
                raise ValueError('Unexpected, duplicate or non-file archive member')
            if member.size != expected[name]['size']:
                raise ValueError('Archive size mismatch')
            payload = archive.extractfile(member).read()
            if digest(payload) != expected[name]['sha256']:
                raise ValueError('Archive content mismatch')
            decoded[name] = payload
    if set(decoded) != set(expected):
        raise ValueError('Missing archive members')
    return decoded


def load(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frozen_files(root,path):
    freeze = checked_manifest(read(root,path))
    if freeze.get('study') != STUDY:
        raise ValueError('Wrong study freeze')
    for name,sha in freeze['files'].items():
        public_name(name,managed=True)
        data = read(root,relative(name),allow_marker_literals=PurePosixPath(name).suffix=='.py')
        if digest(data) != sha:
            raise ValueError('Frozen file changed: '+name)
    return freeze


def report_receipts(root,group,expected_count):
    report = json.loads(read(root,group['report']))
    if report.get('module') != group['module'] or report.get('status') != 'complete':
        raise ValueError('Validated complete module report required')
    count_key = 'n_backgrounds' if group['module']=='bridge' else 'n_paired_jobs'
    if report.get(count_key) != expected_count:
        raise ValueError('Complete declared report count required')
    receipts = report.get('banks')
    if not isinstance(receipts,list) or len(receipts) != len(group['banks']):
        raise ValueError('Report must bind every raw bank/shard')
    by_identity = {receipt['identity']:receipt for receipt in receipts}
    if len(by_identity) != len(receipts):
        raise ValueError('Duplicate report bank identities')
    files = set()
    identities = set()
    for bank in group['banks']:
        relative(bank)
        if not bank.startswith(RESULTS+'/'):
            raise ValueError('Unexpected result bank namespace')
        manifest_name = bank+'/manifest.json'
        manifest_raw = read(root,manifest_name)
        manifest = checked_manifest(manifest_raw)
        identity = manifest.get('identity',manifest.get('identity_sha256'))
        if identity in identities or identity not in by_identity:
            raise ValueError('Duplicate or unreported raw bank')
        identities.add(identity)
        receipt = by_identity[identity]
        if receipt['manifest_sha256'] != digest(manifest_raw):
            raise ValueError('Report manifest hash differs from raw bank')
        if group['module']=='bridge':
            names = {f'{seed}.json' for seed in manifest['resolved_config']['seeds']}
        else:
            jobs = manifest['jobs']
            names = {f"job_{job['job_id']:08d}.json" for job in jobs}
            if len(names) != len(jobs):
                raise ValueError('Duplicate declared job IDs')
        directory = safe_path(root,bank+'/completed')
        actual = {path.name for path in directory.iterdir()}
        if actual != names or receipt.get('jobs') != len(names) or set(receipt.get('completed_records_sha256',{})) != names:
            raise ValueError('Complete exact job files and report payload hashes required')
        for name in names:
            payload_name = bank+'/completed/'+name
            if digest(read(root,payload_name)) != receipt['completed_records_sha256'][name]:
                raise ValueError('Completed record differs from validated report payload')
            files.add(payload_name)
        files.add(manifest_name)
        binding = bank+('/bridge_binding.json' if group['module']=='bridge' else '/confirmation_binding.json')
        files.add(binding)
        receipt_data = json.loads(read(root,binding))
        if receipt_data.get('freeze_sha256') != digest(read(root,group['freeze'])):
            raise ValueError('Bank binding differs from declared module freeze')
        if group['module']=='bridge':
            files.add(bank+'/summary.json')
    if identities != set(by_identity):
        raise ValueError('Report contains extra raw bank receipts')
    if sum(receipt['jobs'] for receipt in receipts) != expected_count:
        raise ValueError('Raw report receipts do not exhaust declared module count')
    files.add(group['report'])
    for name in group['summaries']:
        files.add(public_name(name,managed=True))
    return report,files


def validate_positive(root,group,report,freeze):
    analyzer = load('_release_confirmation_analyzer',root/PREFIX/'analyze.py')
    design = json.loads(read(root,group['design']))
    book_path = root/PREFIX/'seedbook.json';proposal_path=root/PREFIX/'candidate_design.json'
    book = json.loads(read(root,str(book_path.relative_to(root))))
    proposal = json.loads(read(root,str(proposal_path.relative_to(root))))
    selection=calibration=None;selection_path=calibration_path=None
    if group['module']=='holdout':
        selection_path=root/PREFIX/'calibration_selection.json'
        calibration_path=root/PREFIX/'calibration_report/report.json'
        selection=json.loads(read(root,str(selection_path.relative_to(root))))
        calibration=json.loads(read(root,str(calibration_path.relative_to(root))))
    analyzer.validate_design(group['module'],design,book,proposal,selection,calibration,selection_path,calibration_path,proposal_path)
    af_path=root/PREFIX/'analyze_freeze.json';af=analyzer.verify_analysis_freeze(af_path)
    if report.get('analysis_freeze_identity') != af['identity_sha256'] or report.get('analysis_freeze_sha256') != digest(read(root,str(af_path.relative_to(root)))):
        raise ValueError('Report differs from pre-outcome analyzer freeze')
    if report.get('source_sha256') != digest(read(root,PREFIX+'/analyze.py',allow_marker_literals=True)) or report.get('design_sha256') != digest(read(root,group['design'])):
        raise ValueError('Report analyzer/design identity mismatch')
    if report.get('bootstrap_draws') != 200000 or report.get('bootstrap_seed') != analyzer.BOOTSTRAP_SEED:
        raise ValueError('Report differs from prespecified bootstrap')
    sources=analyzer.backend.source_hashes();sources['array_worker.py']=digest(read(root,'research/full_icu_response_2026-10-08/array_worker.py',allow_marker_literals=True))
    context=dict(root=root,freeze=freeze,freeze_path=root/group['freeze'],seedbook_path=book_path,module=group['module'],
                 output_scopes=json.loads(read(root,PREFIX+'/output_scopes.json')))
    records,receipts=analyzer.load_banks([root/bank for bank in group['banks']],design,sources,book,group['module'],context)
    if sorted(report['banks'],key=lambda item:item['identity']) != sorted(receipts,key=lambda item:item['identity']):
        raise ValueError('Report receipts disagree with fully validated raw banks')
    cases={case['id']:case for case in analyzer.normalized(design)['scenarios']}
    if len(report['scenarios']) != len(cases):
        raise ValueError('Report condition cardinality mismatch')
    def point_equal(actual,expected):
        if expected is None:return actual is None
        return type(actual) in (int,float) and math.isfinite(actual) and math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-10)
    for row in report['scenarios']:
        seeds=design['seeds'];identifier=row['scenario']['id']
        if row['scenario'] != cases.get(identifier):
            raise ValueError('Report condition settings differ from declared design')
        if row['seeds'] != seeds or row['n'] != len(seeds):
            raise ValueError('Report condition seed membership mismatch')
        rates=[];counts=[]
        for name in ('baseline','policy'):
            observed=[analyzer.arm(records[(identifier,seed)],name)['observation'] for seed in seeds]
            acquired=sum(value['acquired_colonizations'] for value in observed)
            exposure=math.fsum(value['patient_days'] for value in observed)
            actual=row['arms'][name]
            if actual['acquired_total'] != acquired or not math.isclose(actual['patient_days_total'],exposure,rel_tol=1e-12,abs_tol=1e-7):
                raise ValueError('Report point totals disagree with raw records')
            rate=None if exposure==0 else 1000*acquired/exposure
            if not point_equal(actual['rate_per_1000_patient_days'],rate):
                raise ValueError('Report pooled rate differs from raw totals')
            rates.append(rate);counts.append(acquired)
        rate_fraction=None if rates[0] in (None,0) or rates[1] is None else 1-rates[1]/rates[0]
        count_fraction=None if counts[0]==0 else 1-counts[1]/counts[0]
        checks={'rate_prevented_fraction':rate_fraction,'count_prevented_fraction':count_fraction,
                'mean_prevented':(counts[0]-counts[1])/len(seeds),
                'absolute_prevented_rate':None if None in rates else rates[0]-rates[1]}
        if any(not point_equal(row[key],value) for key,value in checks.items()):
            raise ValueError('Report response/prevented point estimate mismatch')
    if {row['scenario']['id'] for row in report['scenarios']} != {case['id'] for case in design['scenarios']}:
        raise ValueError('Report omits or adds conditions')
    return len(records)


def validate_bridge(root,group,report,freeze):
    frozen_files(root,AGGREGATE_FREEZE)
    aggregate_path=root/PREFIX/'bridge/aggregate.py'
    aggregate=load('_release_bridge_aggregate_validator',aggregate_path)
    saved_path=list(sys.path)
    try:
        rows,config,book,provenance=aggregate.collect(root,root/group['freeze'])
        expected=aggregate.summarize_complete(rows,config,book,provenance,digest(read(root,PREFIX+'/bridge/aggregate.py',allow_marker_literals=True)))
    finally:
        sys.path[:]=saved_path
    if not aggregate.equal_payload(report,expected):
        raise ValueError('Bridge aggregate report fails validated individual-row reconciliation')
    parts={part['output_scope'] for values in book['partitions'].values() for part in values}
    if set(group['banks']) != parts or len(group['banks']) != 38:
        raise ValueError('Expected all 38 frozen bridge partition scopes')
    if json.loads(read(root,group['design'])) != dict(book['scientific_config'],seeds=book['all_seeds']):
        raise ValueError('Bridge aggregate design differs from seedbook')
    validate_bridge_binding(root,group,report)
    return len(rows)


def validate_bridge_binding(root,group,report):
    binding_name=str(PurePosixPath(group['report']).parent/'binding.json')
    if PurePosixPath(group['report']).name != 'summary.json' or group['summaries'].count(binding_name) != 1:
        raise ValueError('Bridge aggregate requires its report-sibling binding.json exactly once')
    expected={'schema_version':1,'module':'bridge','status':'complete','n_backgrounds':256,
        'summary_sha256':digest(read(root,group['report'])),
        **{key:report[key] for key in ('config_sha256','seeds','freeze_sha256','freeze_identity','scientific_sources','source_sha256','analysis_code_sha256')},
        'summary_filename':'summary.json','config_path':report['config_path']}
    if json.loads(read(root,binding_name)) != expected:
        raise ValueError('Bridge aggregate binding differs from validated summary')


def validate_assembly_dependencies(root):
    for name in ASSEMBLY_FILES:
        data=read(root,public_name(name),allow_marker_literals=PurePosixPath(name).suffix=='.py')
        if name in ASSEMBLY_GUARDS and digest(data) != ASSEMBLY_GUARDS[name]:
            raise ValueError('Preserved assembler input changed: '+name)


def validate_final_manuscript_outputs(root,names):
    prefix=PREFIX+'/manuscript/confirmation/'
    folders={str(PurePosixPath(name).parent) for name in names if name.startswith(prefix)
        and PurePosixPath(name).name in ('main.tex','antecedent_results.tex')}
    for folder in folders:
        main=folder+'/main.tex';fragment=folder+'/antecedent_results.tex';receipt_name=folder+'/manifest.json'
        if not {main,fragment,receipt_name} <= names:
            raise ValueError('Final manuscript requires main.tex, antecedent_results.tex and manifest.json together')
        receipt=json.loads(read(root,receipt_name))
        expected_fragment={'status':'uncompiled_supporting_body_fragment','sha256':digest(read(root,fragment)),
            'automatically_included_in_main':False,'pdf_exported':False}
        if receipt.get('main_sha256') != digest(read(root,main)) or receipt.get('supporting_sources') != {'antecedent_results.tex':expected_fragment}:
            raise ValueError('Final manuscript main/supporting source bytehash binding differs')
        if receipt.get('base_sha256') != ASSEMBLY_GUARDS[MANUSCRIPT_BASE] or receipt.get('builder_sha256') != digest(read(root,PREFIX+'/manuscript/build_confirmation.py')):
            raise ValueError('Final manuscript base/builder binding differs')
        for name in ('antecedent_27_conditions','antecedent_18_contrasts'):
            if receipt.get('protected_spans_location',{}).get(name) != 'antecedent_results.tex':
                raise ValueError('Antecedent protected-span location differs')


def validate_spec(root,spec):
    if not isinstance(spec,dict) or set(spec) != {'schema_version','groups','files'} or spec['schema_version'] != 1:
        raise ValueError('Explicit schema-1 groups/files required')
    if not isinstance(spec['groups'],list) or not isinstance(spec['files'],list):
        raise ValueError('Explicit group/file lists required')
    names={public_name(name) for name in spec['files']}
    if len(names) != len(spec['files']):
        raise ValueError('Duplicate explicit files')
    if len(spec['groups']) != 4 or {group.get('module') for group in spec['groups']} != set(COUNTS):
        raise ValueError('Release requires all four declared complete modules')
    verified=[];banks=set()
    for group in spec['groups']:
        if set(group) != {'module','design','banks','report','summaries','freeze'} or not isinstance(group['banks'],list) or not isinstance(group['summaries'],list):
            raise ValueError('Explicit module/design/banks/report/summaries/freeze required')
        for name in (group['design'],group['report'],group['freeze']):
            public_name(name,managed=True)
        if len(group['banks']) != len(set(group['banks'])) or banks.intersection(group['banks']):
            raise ValueError('Duplicate bank/shard in release specification')
        banks.update(group['banks'])
        freeze=frozen_files(root,group['freeze'])
        # Frozen scientific declarations already form a public whitelist;
        # preserve their full dependency closure so restored freezes verify.
        for name in freeze['files']:
            names.add(public_name(name,managed=True))
        report,group_files=report_receipts(root,group,COUNTS[group['module']])
        n=validate_bridge(root,group,report,freeze) if group['module']=='bridge' else validate_positive(root,group,report,freeze)
        if n != COUNTS[group['module']]:
            raise ValueError('Raw complete module count mismatch')
        for name in group_files:
            names.add(public_name(name,managed=True))
        names.update((group['design'],group['freeze']))
        verified.append(dict(module=group['module'],records=n,report=group['report'],report_sha256=digest(read(root,group['report'])),
                             freeze_identity=freeze['identity_sha256'],banks=group['banks']))
    for freeze_name in (PREFIX+'/analyze_freeze.json',AGGREGATE_FREEZE):
        for name in frozen_files(root,freeze_name)['files']:
            names.add(public_name(name,managed=True))
    if not REQUIRED_FILES <= names:
        raise ValueError('Explicit analyzer/release/assembly/reproducibility inputs required')
    validate_assembly_dependencies(root)
    validate_final_manuscript_outputs(root,names)
    return names,verified


def preserve(root,spec,files,verified,output):
    """Internal byte/chunk stage, separately testable without scientific records."""
    metadata={}
    for name in sorted(files):
        data=read(root,name,allow_marker_literals=PurePosixPath(name).suffix=='.py')
        if len(data)>MAX_RAW_CHUNK_BYTES:
            raise ValueError('Single input exceeds bounded raw chunk limit')
        metadata[name]=dict(size=len(data),sha256=digest(data))
    # Produce all bounded chunks before publishing a manifest. A later failure
    # leaves no valid release manifest and can never masquerade as completion.
    output.mkdir(parents=True)
    manifest=dict(schema_version=1,scope='Four complete synthetic confirmation modules and explicit public scientific files.',
        study=STUDY,specification=spec,verified_groups=verified,files={},chunks=[],
        chunk_limits=dict(max_compressed_bytes=MAX_COMPRESSED_BYTES,max_raw_bytes=MAX_RAW_CHUNK_BYTES,max_members=MAX_MEMBERS))
    def publish(selected):
        payloads={name:read(root,name,allow_marker_literals=PurePosixPath(name).suffix=='.py') for name in selected}
        if any(digest(data)!=metadata[name]['sha256'] for name,data in payloads.items()):
            raise ValueError('Input changed during preservation')
        encoded=encode_chunk(payloads)
        if len(encoded)>MAX_COMPRESSED_BYTES:
            if len(selected)==1:
                raise ValueError('Single input cannot fit below-20-MiB release chunk')
            middle=len(selected)//2;publish(selected[:middle]);publish(selected[middle:]);return
        filename=f'payload_{len(manifest["chunks"]):04d}.tar.gz'
        with (output/filename).open('xb') as stream:stream.write(encoded)
        manifest['chunks'].append(dict(file=filename,sha256=digest(encoded),size=len(encoded),members=selected))
        for name in selected:manifest['files'][name]=dict(metadata[name],chunk=filename)
    selected=[];raw_size=0
    for name in sorted(metadata):
        size=metadata[name]['size']
        if selected and (len(selected)==MAX_MEMBERS or raw_size+size>MAX_RAW_CHUNK_BYTES):
            publish(selected);selected=[];raw_size=0
        selected.append(name);raw_size+=size
    if selected:publish(selected)
    manifest['identity']=digest(canonical(manifest).encode())
    encoded_manifest=json.dumps(manifest,indent=2,sort_keys=True)+'\n'
    if len(encoded_manifest.encode())>MAX_COMPRESSED_BYTES:
        raise ValueError('Release manifest exceeds below-20-MiB public-file limit')
    with (output/'release_manifest.json').open('x') as stream:stream.write(encoded_manifest)
    return manifest


def freeze(root,spec,output):
    root=safe_path(root);output=safe_path(output)
    if output.exists():raise FileExistsError(output)
    files,verified=validate_spec(root,spec)
    return preserve(root,spec,files,verified,output)


def verify_release(release):
    release=safe_path(release)
    manifest=checked_manifest(read(release,'release_manifest.json'))
    if manifest.get('schema_version') != 1 or manifest.get('study') != STUDY:
        raise ValueError('Unknown release schema/study')
    for name,metadata in manifest['files'].items():
        public_name(name,managed=True)
        if type(metadata['size']) is not int or metadata['size']<0 or not re.fullmatch('[0-9a-f]{64}',metadata['sha256']):
            raise ValueError('Invalid release file metadata')
        if any(str(parent) in manifest['files'] for parent in PurePosixPath(name).parents if str(parent)!='.'):
            raise ValueError('Conflicting file and parent paths')
    seen=set();chunks=set()
    for chunk in manifest['chunks']:
        if chunk['file'] in chunks or len(chunk['members'])!=len(set(chunk['members'])):
            raise ValueError('Duplicate chunk or member declaration')
        relative(chunk['file']);chunks.add(chunk['file'])
        data=read(release,chunk['file'],privacycheck=False)
        if digest(data)!=chunk['sha256'] or len(data)!=chunk['size']:
            raise ValueError('Compressed chunk changed')
        expected={name:manifest['files'][name] for name in chunk['members']}
        if seen.intersection(expected) or any(value['chunk']!=chunk['file'] for value in expected.values()):
            raise ValueError('Overlapping or incorrectly assigned chunk members')
        decode_chunk(data,expected);seen.update(expected)
    if seen!=set(manifest['files']):raise ValueError('Incomplete release')
    if {path.name for path in release.iterdir()} != chunks|{'release_manifest.json'}:
        raise ValueError('Unexpected release directory member')
    return manifest


def restore(release,output):
    release=safe_path(release);output=safe_path(output)
    if output.exists():raise FileExistsError(output)
    manifest=verify_release(release)
    output.mkdir(parents=True)
    for chunk in manifest['chunks']:
        data=read(release,chunk['file'],privacycheck=False)
        if digest(data)!=chunk['sha256']:raise ValueError('Release changed during restoration')
        payloads=decode_chunk(data,{name:manifest['files'][name] for name in chunk['members']})
        for name,payload in payloads.items():
            path=output/relative(name);path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('xb') as stream:stream.write(payload)
    return manifest


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='mode',required=True)
    create=sub.add_parser('freeze');create.add_argument('--root',type=Path,required=True);create.add_argument('--spec',type=Path,required=True);create.add_argument('--output',type=Path,required=True)
    unpack=sub.add_parser('restore');unpack.add_argument('--release',type=Path,required=True);unpack.add_argument('--output',type=Path,required=True)
    check=sub.add_parser('verify');check.add_argument('--release',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.mode=='freeze':result=freeze(args.root,json.loads(args.spec.read_text()),args.output)
    elif args.mode=='restore':result=restore(args.release,args.output)
    else:result=verify_release(args.release)
    print(json.dumps(dict(mode=args.mode,files=len(result['files']),chunks=len(result['chunks']),identity=result['identity'])))
    return 0


if __name__=='__main__':sys.exit(main())
