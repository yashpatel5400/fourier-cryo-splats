#!/usr/bin/env python3
"""Preserve final external text and audit metadata without publishing reasoning."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(data):return hashlib.sha256(data).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);directory=p.parse_args().directory
    manifest=json.loads((directory/'manifest.json').read_text())
    assert manifest['complete'] and manifest['requested_model']=='claude-fable-5-1'
    raw=(directory/'response.jsonl').read_bytes();assert digest(raw)==manifest['raw_response_sha256']
    events=[json.loads(line) for line in raw.decode().splitlines()]
    results=[e for e in events if e.get('type')=='result'];assert len(results)==1
    result=results[0];assert not result.get('is_error') and set(result.get('modelUsage',{}))=={'claude-fable-5-1'}
    targets=[directory/n for n in ['public-events.jsonl','assistant-text.md','publication-audit.json']]
    if any(p.exists() for p in targets):raise ValueError('Preserve original publication audit')
    redacted=[]
    def clean(value):
        if isinstance(value,list):return [clean(x) for x in value]
        if not isinstance(value,dict):return value
        if value.get('type') in ['thinking','redacted_thinking','thinking_delta','signature_delta','reasoning','reasoning_delta']:
            serialized=json.dumps(value,sort_keys=True,ensure_ascii=False).encode()
            record=dict(type='redacted_reasoning',original_type=value['type'],sha256=digest(serialized),bytes=len(serialized))
            redacted.append(record);return record
        out={}
        for k,v in value.items():
            if k in ['thinking','reasoning','signature'] and isinstance(v,str):
                record=dict(sha256=digest(v.encode()),bytes=len(v.encode()))
                out[k+'_redacted']=record;redacted.append(dict(field=k,**record))
            else:out[k]=clean(v)
        return out
    safe=[]
    private_events=0
    for event in events:
        if event.get('type') in ['assistant','result']:
            safe.append(clean(event))
        else:
            # Provider-internal continuations can contain unrequested context;
            # retain their identity/hash without exposing their contents.
            serialized=json.dumps(event,sort_keys=True,ensure_ascii=False).encode()
            safe.append(dict(type=event.get('type'),subtype=event.get('subtype'),
                private_event_content_omitted=True,sha256=digest(serialized),bytes=len(serialized)))
            private_events+=1
    text=[]
    for event in events:
        if event.get('type')=='assistant':
            for block in event.get('message',{}).get('content',[]):
                if block.get('type')=='text':text.append(block['text'])
    # Preserve every visible assistant text block, including continuations.
    joined='\n\n'.join(text)
    targets[0].write_text(''.join(json.dumps(e,ensure_ascii=False)+'\n' for e in safe))
    targets[1].write_text(joined+'\n')
    audit=dict(raw_response_sha256=digest(raw),raw_response_published=False,
        exact_model='claude-fable-5-1',event_count=len(events),assistant_text_blocks=len(text),
        assistant_text_matches_terminal_result=joined.strip()==result['result'].strip(),
        reasoning_blocks_or_fields_redacted=len(redacted),
        non_assistant_provider_events_redacted=private_events,
        final_result_sha256=digest(result['result'].encode()),
        public_files={p.name:dict(bytes=p.stat().st_size,sha256=digest(p.read_bytes())) for p in targets[:2]})
    targets[2].write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
