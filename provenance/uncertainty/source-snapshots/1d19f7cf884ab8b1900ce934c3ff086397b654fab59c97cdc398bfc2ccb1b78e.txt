#!/usr/bin/env python3
"""Separate, bounded CPU continuations; preserve all parent fit outcomes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

import starfile
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
BIN = ROOT/'tmp/relion-env/bin'
PROTOCOL = ROOT/'research/uncertainty/RELION-CONTINUATION-PROTOCOL.md'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def save(path, value):
    temp = path.with_suffix('.json.tmp')
    temp.write_text(json.dumps(value, indent=2)+'\n'); temp.replace(path)


def last_checkpoint(folder, stage):
    candidates = []
    for file in folder.glob(f'{stage}_it???_optimiser.star'):
        match = re.fullmatch(rf'{stage}_it(\d{{3}})_optimiser.star', file.name)
        if match:
            iteration = int(match.group(1))
            if starfile.read(file)['rlnCurrentIteration'] != iteration:
                raise ValueError('Checkpoint name/state mismatch')
            candidates.append((iteration, file))
    if not candidates: raise RuntimeError('No saved optimizer to continue')
    return max(candidates)


def execute(ds, frozen_files):
    for name, digest in frozen_files.items():
        if sha(ROOT/name) != digest: raise RuntimeError(f'Continuation source changed: {name}')
    parent = BASE/'relion-reconstruction-v2'/ds
    parent_path = parent/'record.json'; previous = json.loads(parent_path.read_text())
    if not previous.get('complete'): raise RuntimeError('Parent is still running')
    out = BASE/'relion-reconstruction-v3'/ds
    if out.exists(): raise RuntimeError('Preserve all earlier continuation attempts')
    out.mkdir(parents=True); path = out/'record.json'
    result = {'complete':False, 'dataset':ds, 'parent_record':str(parent_path.relative_to(ROOT)),
        'parent_record_sha256':sha(parent_path), 'pixel_size_A':previous['pixel_size_A'],
        'field_A':previous['field_A'], 'box':previous['box'], 'seed':previous['seed'],
        'counts':previous['counts'], 'input_hashes':previous['input_hashes'],
        'source_snapshot':source_snapshot(ROOT,Path(__file__),[
            str(PROTOCOL.relative_to(ROOT)), 'research/uncertainty/RELION-BASELINE-PROTOCOL.md']),
        'continuation_source_hashes':frozen_files, 'runs':[],
        'refinement_converged':False, 'half_labels_verified':False,
        'scope':'Bounded continuation after a compute timeout; parent failures are not replaced.'}
    save(path,result)
    try:
        if previous.get('error') or not previous.get('runs') or previous['runs'][-1].get('status') != 'wall_limit':
            result.update(complete=True, stage='not_applicable',
                reason='Parent did not terminate at a declared wall limit')
            save(path,result); return
        for name, digest in previous['input_hashes'].items():
            if sha(ROOT/name) != digest: raise RuntimeError(f'Parent input changed: {name}')
        for name, digest in previous['executables'].items():
            if sha(BIN/name) != digest: raise RuntimeError('Parent executable changed')
        result['executables'] = previous['executables']
        stage = previous['runs'][-1]['stage']; iteration, checkpoint = last_checkpoint(parent, stage)
        family = sorted(parent.glob(checkpoint.name.removesuffix('_optimiser.star')+'*'))
        result['checkpoint'] = {'stage':stage,'iteration':iteration,
            'optimizer':str(checkpoint.relative_to(ROOT)),
            'family_hashes':{str(f.relative_to(ROOT)):sha(f) for f in family if f.is_file()}}
        save(path,result)
        data = ROOT/'data/uncertainty/relion-reconstruction-v2'/ds
        def run(which, command, limit, threads):
            log = ROOT/f'logs/uncertainty/relion-v3-{ds}-{which}.log'
            if log.exists(): raise RuntimeError('Preserve earlier process log')
            row = {'stage':which,'command':command,'complete':False,'wall_limit_seconds':limit}
            result['stage'] = which; result['runs'].append(row); save(path,result)
            begin = time.perf_counter()
            with log.open('w') as stream:
                process = subprocess.Popen(command,cwd=data,stdout=stream,stderr=subprocess.STDOUT,
                    start_new_session=True,env={**os.environ,'OMP_NUM_THREADS':str(threads),'OPENBLAS_NUM_THREADS':'1'})
                try:
                    code = process.wait(timeout=limit); status = 'finished' if code == 0 else 'failed'
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGTERM)
                    try: process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid,signal.SIGKILL); process.wait()
                    code = process.returncode; status = 'wall_limit'
            row.update(complete=True,returncode=code,status=status,seconds=time.perf_counter()-begin,
                log=str(log.relative_to(ROOT)),log_sha256=sha(log))
            save(path,result); print(ds,which,status,row['seconds'],flush=True)
            return status == 'finished'
        execution = ['--preread_images','--dont_combine_weights_via_disc','--pool','30']
        if stage == 'initial':
            ok = run('initial',[str(BIN/'relion_refine'),'--continue',str(checkpoint),
                '--o',str(out/'initial'),'--iter','100','--j','8',*execution],14400,8)
            if ok:
                reference = out/'initial_it100_class001.mrc'
                if not reference.exists(): raise RuntimeError('Completed initializer has no final map')
                common = ['--ctf','--K','1','--sym','C1','--flatten_solvent','--zero_mask',
                    *execution,'--particle_diameter',str(.8*previous['field_A']),
                    '--oversampling','1','--offset_step','2','--j','4','--random_seed',str(previous['seed'])]
                run('refine',[str(BIN/'mpirun'),'-np','3',str(BIN/'relion_refine_mpi'),
                    '--i','refinement.star','--o',str(out/'refine'),'--auto_refine','--split_random_halves',
                    '--ref',str(reference),'--firstiter_cc','--ini_high','50','--iter','50',
                    '--healpix_order','2','--auto_local_healpix_order','4','--offset_range','5',
                    '--low_resol_join_halves','40','--norm','--scale','--pad','2',*common],21600,4)
        elif stage == 'refine':
            run('refine',[str(BIN/'mpirun'),'-np','3',str(BIN/'relion_refine_mpi'),
                '--continue',str(checkpoint),'--o',str(out/'refine'),'--iter','50','--j','4',*execution],21600,4)
        else: raise RuntimeError(f'Unexpected parent stage: {stage}')
        result['optimiser_status'] = []
        for file in sorted(out.glob('*_optimiser.star')):
            meta = starfile.read(file)
            result['optimiser_status'].append({'path':str(file.relative_to(ROOT)),'sha256':sha(file),
                **{key:meta.get(key) for key in ['rlnCurrentIteration','rlnHasConverged','rlnGradHasConverged']}})
        final = out/'refine_optimiser.star'
        result['refinement_converged'] = bool(final.exists() and starfile.read(final).get('rlnHasConverged') == 1
            and result['runs'][-1]['status'] == 'finished')
        expected = starfile.read(data/'refinement.star')['particles'].set_index('rlnImageName').rlnRandomSubset.to_dict()
        result['half_label_checks'] = []
        for file in sorted(out.glob('refine*_data.star')):
            actual = starfile.read(file)['particles'].set_index('rlnImageName').rlnRandomSubset.to_dict()
            result['half_label_checks'].append({'path':str(file.relative_to(ROOT)),
                'sha256':sha(file),'same_particles_and_halves':actual == expected})
        result['half_labels_verified'] = bool(result['half_label_checks']) and all(
            r['same_particles_and_halves'] for r in result['half_label_checks'])
        result.update(complete=True,stage='fit_execution_complete',evaluation_complete=False)
    except Exception as error:
        result.update(complete=True,stage='execution_error',error=repr(error),refinement_converged=False)
        save(path,result); raise
    save(path,result)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--wait-hours',type=float,default=4)
    args = parser.parse_args()
    frozen = {str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),PROTOCOL]}
    for name in frozen:
        committed = subprocess.check_output(['git','show',f'HEAD:{name}'],cwd=ROOT)
        if hashlib.sha256(committed).hexdigest() != frozen[name]:
            raise RuntimeError('Commit continuation source/protocol before fitting')
    for ds in ['10028','10049','10076']:
        begin = time.monotonic(); parent = BASE/'relion-reconstruction-v2'/ds/'record.json'
        while not json.loads(parent.read_text()).get('complete'):
            if time.monotonic()-begin > args.wait_hours*3600: raise TimeoutError('Parent wait exceeded')
            time.sleep(30)
        try: execute(ds,frozen)
        except Exception as error: print(ds,'ERROR',repr(error),flush=True)


if __name__ == '__main__': main()
