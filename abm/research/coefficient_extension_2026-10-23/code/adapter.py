"""Load the verified workflow-coefficient adapter by path, unchanged and hash-checked.

Study: research/coefficient_extension_2026-10-23 (see ../PROTOCOL.md).

The adapter ``research/workflow_coefficient_2026-10-21/code/workflow_background.py``
(sha256 939aabd5...) is executed from its own location under the module name
``workflow_background``. It is never copied or edited. Its directory is NOT put
on ``sys.path`` (that folder also holds a ``run_batch.py`` and a
``summarize.py``), and its ``STUDY`` constant points at the earlier study, so
this study takes every output path from its own ``__file__``.

``sys.dont_write_bytecode`` is set before anything is executed, so no
``__pycache__`` is created next to the adapter or any frozen file.
"""
import sys

sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import importlib.util  # noqa: E402
from pathlib import Path  # noqa: E402

CODE = Path(__file__).resolve().parent
STUDY = CODE.parent
ABM = CODE.parents[2]
ADAPTER_PATH = ABM / 'research/workflow_coefficient_2026-10-21/code/workflow_background.py'
ADAPTER_SHA256 = '939aabd55a64ce59f1c9111052cdd36b7545fbc848578b48eef1dbb37d2e8c67'
ADAPTER_RELATIVE = 'research/workflow_coefficient_2026-10-21/code/workflow_background.py'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load():
    existing = sys.modules.get('workflow_background')
    if existing is not None:
        if Path(existing.__file__).resolve() != ADAPTER_PATH.resolve():
            raise RuntimeError('A different workflow_background module is already imported')
        return existing
    if sha256(ADAPTER_PATH) != ADAPTER_SHA256:
        raise RuntimeError('Verified adapter changed: ' + str(ADAPTER_PATH))
    spec = importlib.util.spec_from_file_location('workflow_background', ADAPTER_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules['workflow_background'] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop('workflow_background', None)
        raise
    return module


wb = _load()
if wb.adapter_sha256() != ADAPTER_SHA256:
    raise RuntimeError('Adapter hash changed after loading')
if wb.ABM.resolve() != ABM.resolve():
    raise RuntimeError('Adapter and study disagree on the ABM root')

__all__ = ['ABM', 'ADAPTER_PATH', 'ADAPTER_RELATIVE', 'ADAPTER_SHA256', 'CODE', 'STUDY', 'sha256', 'wb']
