"""Immutable source snapshots for long-running research experiments."""
import hashlib
import importlib.metadata
import platform
import subprocess
import sys
from pathlib import Path


def source_snapshot(root,entrypoint,extra_scripts=()):
    """Archive exact source bytes at process start, including dirty development code."""
    root=Path(root);paths=[Path(entrypoint),*[root/p for p in extra_scripts]]
    paths.extend(sorted((root/'src/fourier_splats').glob('*.py')))
    directory=root/'provenance/uncertainty/source-snapshots';directory.mkdir(parents=True,exist_ok=True)
    records={}
    for path in paths:
        content=path.read_bytes();digest=hashlib.sha256(content).hexdigest();dest=directory/(digest+'.txt')
        if dest.exists():
            if dest.read_bytes()!=content:raise AssertionError('Source snapshot digest collision')
        else:dest.write_bytes(content)
        records[str(path.relative_to(root))]={'sha256':digest,'snapshot':str(dest.relative_to(root))}
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
    packages={}
    for package in ['numpy','scipy','finufft','torch','cvxpy','mrcfile']:
        try:packages[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:packages[package]=None
    return {'git_head_at_snapshot':head,'interpretation':'Exact archived source bytes are authoritative; HEAD may precede uncommitted development edits.',
            'sources':records,'python':sys.version,'platform':platform.platform(),'packages':packages}
