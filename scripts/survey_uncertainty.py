#!/usr/bin/env python3
"""Save a reproducible recent-literature search using the Europe PMC API.

Raw abstracts/full text remain local in the ignored third-party material folder.
The public ledger contains bibliographic metadata, queries, and inclusion flags.
Human screening and reading notes are separate; a hit is not a reviewed paper.
"""
import concurrent.futures
import datetime
import hashlib
import json
import time
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "background/papers/uncertainty"
OUT = ROOT / "research/uncertainty"
QUERIES = {
    "density_uncertainty": '("cryo-EM" OR "cryo electron microscopy") AND (uncertainty OR confidence OR Bayesian OR posterior OR bootstrap) AND FIRST_PDATE:[2023-01-01 TO 2026-09-29]',
    "validation": '("cryo-EM" OR "cryo electron microscopy") AND (validation OR benchmark OR "Fourier shell correlation") AND FIRST_PDATE:[2024-01-01 TO 2026-09-29]',
    "heterogeneity": '("cryo-EM" OR "cryo electron microscopy") AND ("heterogeneity" OR "conformational landscape") AND (method OR reconstruction OR algorithm) AND FIRST_PDATE:[2024-01-01 TO 2026-09-29]',
    "learned_priors": '("cryo-EM" OR "cryo electron microscopy") AND ("diffusion model" OR "deep learning" OR "neural network") AND (reconstruction OR validation OR uncertainty) AND FIRST_PDATE:[2024-01-01 TO 2026-09-29]',
}


def search(item):
    name, query = item
    cursor='*'; hits=[]; hashes=[]; failures=[]; total=None
    for page in range(100):
        result=None
        for attempt in range(3):
            try:
                r = requests.get(
                    "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                    params={"query":query,"format":"json","resultType":"lite",
                            "pageSize":200,"cursorMark":cursor},timeout=30)
                r.raise_for_status();result=r.json();break
            except (requests.RequestException,ValueError) as error:
                failures.append({'page':page,'attempt':attempt,'error':str(error)})
                time.sleep(attempt+1)
        if result is None:
            break
        content=r.content
        (LOCAL / f'{name}-page-{page}.json').write_bytes(content)
        hashes.append(hashlib.sha256(content).hexdigest())
        total=result['hitCount'];batch=result['resultList']['result'];hits.extend(batch)
        following=result.get('nextCursorMark')
        if not batch or not following or following==cursor or len(hits)>=total:
            break
        cursor=following
    return name,{'resultList':{'result':hits},'hitCount':total,
                 'complete':total is not None and len(hits)>=total,
                 'failures':failures},hashes


def main():
    LOCAL.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    entries, queries = {}, {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for name, result, digest in ex.map(search, QUERIES.items()):
            hits = result["resultList"]["result"]
            queries[name] = {"query": QUERIES[name], "hit_count": result["hitCount"], "returned": len(hits), "raw_page_sha256": digest,"complete":result['complete'],'request_failures':result['failures']}
            for p in hits:
                key = p.get("doi") or (p.get("source", "") + ":" + p["id"])
                if key not in entries:
                    entries[key] = {k: p.get(k) for k in ["id", "source", "title", "authorString", "pubYear", "firstPublicationDate", "doi", "pmcid"]}
                    entries[key]["queries"] = []
                    entries[key]["review_status"] = "unscreened search hit"
                entries[key]["queries"].append(name)
            print(name, queries[name]["hit_count"], "hits", flush=True)
    ledger = {"accessed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "endpoint": "https://www.ebi.ac.uk/europepmc/webservices/rest/search", "queries": queries, "unique_hits": len(entries), "entries": sorted(entries.values(), key=lambda x: (x.get("firstPublicationDate") or "", x.get("title") or ""), reverse=True)}
    (OUT / "europepmc-search-ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
    (OUT / "search-screening.tsv").write_text("date\ttitle\tdoi\tpmcid\n" + "\n".join("\t".join(str(e.get(k) or "").replace("\t", " ").replace("\n", " ") for k in ["firstPublicationDate", "title", "doi", "pmcid"]) for e in ledger["entries"]) + "\n")
    print("Unique search hits", len(entries))


if __name__ == "__main__":
    main()
