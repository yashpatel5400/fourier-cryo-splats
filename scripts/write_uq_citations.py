#!/usr/bin/env python3
"""Render selected bibliographic records from saved publisher/index metadata."""
import html,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'research/uncertainty/citation-metadata.json').read_text());entries=[]
def clean(text):return ' '.join(html.unescape(re.sub('<[^>]+>','',text)).split()).replace('&',r'\&').replace('_',r'\_')
for key,record in data.items():
    if key=='batlle2025':continue
    if 'metadata' in record:
        m=record['metadata'];title=clean(m['title'][0]);author=' and '.join(clean(a['family'])+', '+clean(a.get('given','')) for a in m['author'])
        journal=clean(m['container-title'][0]);year=m.get('published',m.get('issued'))['date-parts'][0][0]
    elif record.get('europepmc'):
        m=record['europepmc'][0];title=clean(m['title']).rstrip('.');author=' and '.join(clean(a['lastName'])+', '+clean(a.get('firstName','')) for a in m['authorList']['author'])
        journal=clean(m['journalInfo']['journal']['title']);year=m['pubYear']
    else:raise ValueError('Missing verified bibliographic record: '+key)
    entries.append('@article{'+key+',title={{'+title+'}},author={'+author+'},journal={'+journal+'},year={'+str(year)+'},doi={'+record['doi']+'}}')
entries.append(r'@article{batlle2025,title={{Simultaneous Frequentist Calibration of Confidence Regions for Multiple Functionals in Constrained Inverse Problems}},author={Batlle, Pau and Patil, Pratik and Stanley, Michael and Ruiz Lupon, Javier and Owhadi, Houman and Kuusela, Mikael},journal={arXiv:2510.11708v1},year={2025},url={https://arxiv.org/abs/2510.11708v1}}')
(ROOT/'paper/survey-references.bib').write_text('\n'.join(entries)+'\n')
