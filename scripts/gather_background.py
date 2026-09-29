from pathlib import Path
import requests,json,concurrent.futures,hashlib
root=Path(__file__).resolve().parents[1]
items={
'cryodrgn_nature_methods_2021':'https://www.cs.princeton.edu/courses/archive/fall22/cos597N/papers/zhong_cryodrgn_nmeth.pdf',
'cryodrgn2_iccv_2021':'https://openaccess.thecvf.com/content/ICCV2021/papers/Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_ICCV_2021_paper.pdf',
'cryospire_2025':'https://arxiv.org/pdf/2506.09063',
'cryogs_2025':'https://arxiv.org/pdf/2508.04929',
'gem_2025':'https://arxiv.org/pdf/2509.25075',
'cryosplat_iclr_2026':'https://openreview.net/pdf?id=dLaUZKBzta',
'cryodrgn_protocol_2023':'https://bpb-us-e1.wpmucdn.com/sites.mit.edu/dist/5/91/files/2025/01/Uncovering-structural-ensembles.pdf',
'cryobench_2024':'https://papers.nips.cc/paper_files/paper/2024/file/a2ef5ba272df8f168dc38037cc946be0-Paper-Datasets_and_Benchmarks_Track.pdf'}
def get(item):
 name,url=item
 try:
  r=requests.get(url,timeout=90);r.raise_for_status(); assert r.content.startswith(b'%PDF')
  (root/'background/papers'/f'{name}.pdf').write_bytes(r.content)
  return dict(name=name,url=url,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),status='downloaded locally; not redistributed')
 except Exception as e:return dict(name=name,url=url,status=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: results=list(ex.map(get,items.items()))
(root/'background/sources.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
