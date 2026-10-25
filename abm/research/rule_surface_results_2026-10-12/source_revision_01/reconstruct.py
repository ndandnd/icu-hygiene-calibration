#!/usr/bin/env python3
"""Reconstruct bound manuscript sources; never simulate, analyze, or compile."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys

HEX = re.compile(r'[0-9a-f]{64}')
MAX_BYTES = 8_000_000


def require(ok, message):
    if not ok: raise ValueError(message)


def sha(raw): return hashlib.sha256(raw).hexdigest()


def identity(value):
    return sha(json.dumps({k:v for k,v in value.items() if k != 'identity_sha256'},
                         sort_keys=True, separators=(',', ':'), allow_nan=False).encode())


def parse(raw):
    def pairs(rows):
        result = {}
        for k,v in rows:
            require(k not in result, 'Duplicate JSON field'); result[k] = v
        return result
    def nonfinite(value): raise ValueError('Nonfinite JSON constant')
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=nonfinite)


def safe(path):
    path = Path(os.path.abspath(path))
    require(not any(p.is_symlink() for p in (path, *path.parents)), 'Symlink path refused')
    return path


def relative(name):
    require(isinstance(name,str) and name and '\\' not in name and ':' not in name,
            'Portable relative name required')
    p = PurePosixPath(name)
    require(not p.is_absolute() and p.as_posix() == name and all(x not in ('.','..','') for x in p.parts),
            'Canonical relative name required')
    return name


def read(root, name):
    path = safe(root/relative(name))
    require(path.is_file() and stat.S_ISREG(path.stat().st_mode) and path.stat().st_size <= MAX_BYTES,
            'Bounded regular file required')
    return path.read_bytes()


def inventory(root):
    files, directories = set(), set()
    for path in root.rglob('*'):
        safe(path)
        name = path.relative_to(root).as_posix()
        if path.is_file():
            require(stat.S_ISREG(path.stat().st_mode), 'Special member refused'); files.add(name)
        else:
            require(path.is_dir(), 'Special member refused'); directories.add(name)
    return files, directories


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle', type=Path, default=Path(__file__).absolute().parent)
    p.add_argument('--expected-identity', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); bundle = safe(a.bundle); out = safe(a.output)
    require(HEX.fullmatch(a.expected_identity), 'External manifest identity required')
    require(not out.exists() and not out.is_relative_to(bundle) and not bundle.is_relative_to(out),
            'Fresh output outside bundle required')
    manifest_raw = read(bundle, 'source_binding_manifest.json'); manifest = parse(manifest_raw)
    require(manifest['kind'] == 'bound_rule_surface_document_source_bundle' and
            type(manifest['schema_version']) is int and manifest['schema_version'] == 1 and
            identity(manifest) == manifest['identity_sha256'] == a.expected_identity,
            'Externally trusted canonical manifest required')
    require(manifest['statistics_recomputed'] is False and manifest['canonical_write_authorized'] is False and
            manifest['pdf_compilation_included'] is False, 'Source-only scope required')
    descriptors = manifest['files']; require(isinstance(descriptors, list) and descriptors, 'Explicit files required')
    names, files = [], {}
    for descriptor in descriptors:
        require(set(descriptor) == {'path','sha256','bytes'} and HEX.fullmatch(descriptor['sha256']) and
                type(descriptor['bytes']) is int and 0 <= descriptor['bytes'] <= MAX_BYTES, 'Exact descriptor schema')
        name = relative(descriptor['path']); names.append(name); raw = read(bundle, name)
        require(sha(raw) == descriptor['sha256'] and len(raw) == descriptor['bytes'], 'Changed bound file: '+name)
        files[name] = raw
    require(len(names) == len(set(names)) == len({n.casefold() for n in names}), 'Unique member names required')
    actual, directories = inventory(bundle)
    expected_names = set(names)|{'source_binding_manifest.json'}
    expected_dirs = {str(parent) for n in expected_names for parent in PurePosixPath(n).parents if str(parent) != '.'}
    require(actual == expected_names and directories == expected_dirs, 'Exact bundle inventory required')
    require(files['reconstruct.py'] == Path(__file__).read_bytes(), 'Running adapter source binding')
    scientific = parse(files['input/evidence/scientific_review.json'])
    arithmetic = parse(files['input/evidence/full_analysis_arithmetic_acceptance_projection.json'])
    selector = parse(files['input/evidence/selector_arithmetic_acceptance_projection.json'])
    for r in (scientific, arithmetic, selector):
        require(r['status'] == 'pass' and identity(r) == r['identity_sha256'], 'Actual review identity and status')
    for projection, kind, origin in (
        (arithmetic, 'public_full_analysis_scientific_acceptance_projection', 'af2c825de679bbf84291864e9c61dcee13ee7a6b450493004880195a7113b128'),
        (selector, 'public_selector_scientific_acceptance_projection', '65676c46e371a1b6779026fa97030e26bbad3c96267a1cb0db1e19d371b4fb22')):
        require(projection['kind'] == kind and projection['original_receipt_sha256'] == origin,
                'Explicit scientific projection and original receipt pin required')
        retained, omitted = projection['retained_original_fields'], projection['omitted_original_fields']
        hidden = projection['omitted_original_field_name_sha256']
        hidden_pointers = projection['omitted_sensitive_pointer_sha256']
        require(isinstance(retained,list) and isinstance(omitted,list) and
                len(retained) == len(set(retained)) and len(omitted) == len(set(omitted)) and
                not set(retained).intersection(omitted) and
                type(projection['original_top_level_field_count']) is int and
                isinstance(hidden,list) and isinstance(hidden_pointers,list) and len(hidden) == len(set(hidden)) == len(hidden_pointers) and
                all(HEX.fullmatch(x) for x in hidden+hidden_pointers) and
                not set(hidden).intersection(sha(x.encode()) for x in retained+omitted) and
                len(retained)+len(omitted)+len(hidden) == projection['original_top_level_field_count'],
                'Complete declared retained/omitted original field partition')
        metadata = {'schema_version','kind','identity_sha256','original_receipt_sha256','original_receipt_identity',
                    'original_receipt_kind','original_top_level_field_count','retained_original_fields',
                    'omitted_original_fields','omitted_original_json_pointers','omission_pointer_scope',
                    'omitted_original_field_name_sha256','omitted_sensitive_pointer_sha256',
                    'complete_omitted_field_inventory_sha256','complete_omitted_pointer_inventory_sha256',
                    'omitted_entire_map_entry_counts','projection_notice'}
        require(projection['omitted_original_json_pointers'] == ['/'+x.replace('~','~0').replace('/','~1') for x in omitted],
                'Explicit whole-field JSON omission pointers')
        require(HEX.fullmatch(projection['complete_omitted_field_inventory_sha256']) and
                HEX.fullmatch(projection['complete_omitted_pointer_inventory_sha256']), 'Complete omission inventory digest')
        require(set(projection) == set(retained)|metadata and
                {'external_pins','files_sha256'}.issubset(omitted) and
                set(projection['omitted_entire_map_entry_counts']).issubset(omitted) and
                all(type(n) is int and n >= 0 for n in projection['omitted_entire_map_entry_counts'].values()),
                'No undeclared projection field or operational map value')
    require(scientific['report_sha256'] == arithmetic['actual_report_sha256'] == sha(files['input/report.json']) and
            arithmetic['actual_result_review'] is True and arithmetic['bootstrap_draws'] == 50000 and
            arithmetic['bootstrap_seed'] == 20261011 and arithmetic['primary_family_size'] == 6,
            'Exact accepted complete report binding')
    require(arithmetic['full_independent_selector_review_sha256'] == selector['original_receipt_sha256'] and
            selector['actual_selection_review'] is True and selector['records'] == 2432 and selector['conditions'] == 19 and
            selector['producer_exact_output_sha256']['research/rule_surface_followup_2026-10-11/selection.json'] == sha(files['input/selection.json']),
            'Accepted complete selection binding')
    require(read(bundle, 'source_binding_manifest.json') == manifest_raw and
            all(read(bundle,n) == raw for n,raw in files.items()), 'Input changed before reconstruction')
    completed = subprocess.run([sys.executable, '-B', str(bundle/'build_documents.py'), '--root', str(bundle),
                                '--output', str(out)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    require(completed.returncode == 0, 'Source builder failed: '+completed.stderr.decode(errors='replace'))
    expected = manifest['outputs']; require(isinstance(expected, dict) and set(expected) == {
        'manuscript.tex','supplement.tex','figures/weekly_control.png','manuscript.tex.diff','supplement.tex.diff','SOURCE_MANIFEST.json'},
        'Exactly six source outputs required')
    outputs = {}
    for name, descriptor in expected.items():
        require(set(descriptor) == {'sha256','bytes'} and HEX.fullmatch(descriptor['sha256']) and type(descriptor['bytes']) is int,
                'Exact output descriptor required')
        raw = read(out, name)
        require(sha(raw) == descriptor['sha256'] and len(raw) == descriptor['bytes'], 'Reconstructed byte mismatch: '+name)
        outputs[name] = descriptor
    require(inventory(out) == (set(expected), {'figures'}), 'Exact reconstructed inventory')
    require(read(bundle, 'source_binding_manifest.json') == manifest_raw and
            all(read(bundle,n) == raw for n,raw in files.items()) and inventory(bundle) == (actual,directories),
            'Input changed during reconstruction')
    receipt = dict(schema_version=1, kind='bound_document_source_reconstruction_receipt',
        status='pass_exact_source_bytes_only', source_binding_manifest_sha256=sha(manifest_raw),
        source_binding_manifest_identity=a.expected_identity, adapter_sha256=sha(files['reconstruct.py']),
        builder_sha256=sha(files['build_documents.py']), outputs=outputs,
        report_sha256=sha(files['input/report.json']), selection_sha256=sha(files['input/selection.json']),
        input_snapshot_sha256=sha(files['INPUTS.json']), statistics_recomputed=False, simulations_run=False,
        pdf_compiled=False, canonical_written=False, clinical_validation=False,
        scope='Exact source reconstruction; numerical, byte-replay, PDF, clinical and release acceptances remain separate.')
    receipt['identity_sha256'] = identity(receipt)
    with (out/'RECONSTRUCTION_RECEIPT.json').open('x') as f:
        f.write(json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=receipt['status'], outputs=outputs, identity_sha256=receipt['identity_sha256'])))


if __name__ == '__main__': main()
