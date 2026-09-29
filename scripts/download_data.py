#!/usr/bin/env python3
"""Fetch actual deposited particle bytes with reproducible, verified HTTP ranges.

The selection is stratified contiguous blocks over the original particle ordering,
then the published cryoDRGN particle filter (where supplied) is applied. This is
not a uniform independent sample and is explicitly recorded in the manifest.
"""
import argparse, concurrent.futures as cf, hashlib, json, pickle, struct, time
from pathlib import Path
import numpy as np
import requests
from scipy.fft import fft2, fftshift, ifftshift, ifft2

ROOT=Path(__file__).resolve().parents[1]
BASE='https://ftp.ebi.ac.uk/empiar/world_availability'

def get_range(url,start,end):
    for attempt in range(6):
        try:
            # A range-specific query avoids intermediate caches reusing a header-only request.
            response=requests.get(url+f'?range={start}-{end}',headers={'Range':f'bytes={start}-{end}'},timeout=(20,180))
            response.raise_for_status()
            if response.status_code!=206 or response.headers.get('Content-Range','').split('/')[0]!=f'bytes {start}-{end}':
                raise RuntimeError(f'Range not honored: {response.status_code}, {response.headers.get("Content-Range")}')
            if len(response.content)!=end-start+1: raise RuntimeError('Truncated download')
            return response.content,response.headers.get('ETag'),response.headers.get('Content-Range')
        except Exception:
            if attempt==5: raise
            time.sleep(2**attempt)

def run(args):
    directory=ROOT/'data'/args.dataset; directory.mkdir(parents=True,exist_ok=True)
    source=ROOT/'background/cryodrgn_empiar'/('empiar'+args.dataset)/'inputs'
    cs=np.load(next(source.glob('*.cs')),allow_pickle=False)
    # Author-provided metadata repository, never arbitrary user-supplied pickle files.
    poses=pickle.load(open(source/'poses.pkl','rb')); ctf=pickle.load(open(source/'ctf.pkl','rb'))
    n=len(cs); rng=np.random.default_rng(args.seed)
    block=64; starts=np.arange(0,n,block); rng.shuffle(starts)
    valid=np.ones(n,dtype=bool)
    if (source/'filtered.ind.pkl').exists():
        valid[:]=False; valid[np.asarray(pickle.load(open(source/'filtered.ind.pkl','rb')),dtype=int)]=True
    selected=[]
    for start in starts:
        selected.extend([i for i in range(start,min(start+block,n)) if valid[i]])
        if len(selected)>=args.count: break
    selected=np.sort(selected[:args.count]); np.save(directory/'indices.npy',selected)
    N=len(selected); D=args.box
    images=np.lib.format.open_memmap(directory/'images.npy',mode='w+',dtype='float32',shape=(N,D,D))
    np.savez(directory/'metadata.npz',rotations=poses[0][selected],translations=poses[1][selected],ctf=ctf[selected],indices=selected)
    jobs=[]
    for path in np.unique(cs['blob/path'][selected]):
        positions=np.flatnonzero(cs['blob/path'][selected]==path)
        idx=cs['blob/idx'][selected[positions]].astype(int)
        order=np.argsort(idx); idx=idx[order]; positions=positions[order]
        relative=path.decode().split('/imported/')[-1]
        if args.dataset=='10028': relative='Particles/'+relative
        url=f'{BASE}/{args.dataset}/data/{relative}'
        splits=np.flatnonzero(np.diff(idx)>16)+1
        for segment in np.split(np.arange(len(idx)),splits):
            # Bound each response to at most ~64MB for low peak memory.
            for part in np.array_split(segment,max(1,int(np.ceil(len(segment)/128)))):
                if len(part): jobs.append((url,idx[part],positions[part]))
    headers={}
    def work(job):
        url,ids,pos=job
        if url not in headers:
            raw,_,_=get_range(url,0,1023); headers[url]=raw
        header=headers[url]; nx,ny,nz,mode=struct.unpack('<4i',header[:16]); nsymbt=struct.unpack('<i',header[92:96])[0]
        if mode!=2 or nx!=ny: raise ValueError((nx,ny,nz,mode))
        if ids.max()>=nz: raise ValueError('Particle index exceeds source stack')
        stride=nx*ny*4; start=1024+nsymbt+int(ids.min())*stride; end=1024+nsymbt+(int(ids.max())+1)*stride-1
        raw,etag,content_range=get_range(url,start,end)
        x=np.frombuffer(raw,dtype='<f4').reshape(-1,ny,nx)[ids-ids.min()]
        # Crop the centered Fourier transform with symmetric Nyquist handling.
        f=fftshift(fft2(ifftshift(x,axes=(-2,-1)),workers=1),axes=(-2,-1))
        c=nx//2; h=D//2; f=f[:,c-h:c+h,c-h:c+h].copy()
        # Remove Nyquist row/column, which otherwise require folding both endpoints.
        f[:,0,:]=0; f[:,:,0]=0
        xsmall=fftshift(ifft2(ifftshift(f,axes=(-2,-1)),workers=1),axes=(-2,-1)).real*(D/nx)**2
        images[pos]=xsmall.astype('float32')
        return {'url':url,'start':start,'end':end,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'etag':etag,'content_range':content_range,'header_sha256':hashlib.sha256(header).hexdigest(),'source_shape':[nz,ny,nx],'source_indices':ids.tolist(),'output_indices':pos.tolist()}
    records=[]; t=time.time()
    with cf.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures=[executor.submit(work,j) for j in jobs]
        for future in cf.as_completed(futures):
            records.append(future.result()); images.flush()
            if len(records)%10==0 or len(records)==len(jobs): print(args.dataset,len(records),'/',len(jobs),'ranges;',round(sum(x['bytes'] for x in records)/1e9,2),'GB;',round(time.time()-t),'s',flush=True)
            (directory/'download-progress.json').write_text(json.dumps(records))
    manifest={'dataset':args.dataset,'count':N,'available_particles':n,'published_filtered_count':int(valid.sum()),'box':D,'seed':args.seed,'selection':'randomly permuted contiguous source blocks of 64, applying published filter; sorted before output','raw_pixel_size_A':float(cs['blob/psize_A'][0]),'raw_box':int(cs['blob/shape'][0,0]),'pixel_size_A':float(cs['blob/psize_A'][0])*int(cs['blob/shape'][0,0])/D,'data_sign':1 if args.dataset=='10076' else -1,'source':'experimental extracted particle images; not raw movies','metadata_repository':'https://github.com/zhonge/cryodrgn_empiar','elapsed_seconds':time.time()-t,'downloaded_bytes':sum(x['bytes'] for x in records),'ranges':records}
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print('COMPLETE',args.dataset,N,'particles',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('dataset',choices=['10028','10049','10076']); p.add_argument('--count',type=int,default=8192);p.add_argument('--box',type=int,default=64);p.add_argument('--seed',type=int,default=20260928);p.add_argument('--workers',type=int,default=4)
    run(p.parse_args())
