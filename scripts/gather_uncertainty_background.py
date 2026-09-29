#!/usr/bin/env python3
"""Fetch primary research material locally and publish only source/hash metadata."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path
import requests
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "background/papers/uncertainty"
OUT = ROOT / "research/uncertainty"


def fetch(entry):
    entry = dict(entry)
    if not entry["pdf_url"]:
        entry["download_status"] = "primary page listed; no direct PDF configured"
        return entry
    path = LOCAL / (entry["id"] + ".pdf")
    try:
        if not path.exists():
            r = requests.get(entry["pdf_url"], timeout=60)
            r.raise_for_status()
            if not r.content.startswith(b"%PDF"):
                raise ValueError("Response was not a PDF")
            path.write_bytes(r.content)
        content = path.read_bytes()
        reader = PdfReader(path)
        text_path = path.with_suffix(".txt")
        if not text_path.exists():
            text_path.write_text("\n\n".join("[PAGE %d]\n%s" % (i + 1, p.extract_text() or "") for i, p in enumerate(reader.pages)))
        entry.update(download_status="downloaded locally; not redistributed", bytes=len(content), sha256=hashlib.sha256(content).hexdigest(), pages=len(reader.pages))
    except Exception as e:
        entry["download_status"] = str(e)
    print(entry["id"], entry["download_status"], flush=True)
    return entry


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ids", help="Comma-separated subset; default all")
    args = p.parse_args()
    LOCAL.mkdir(parents=True, exist_ok=True)
    with (OUT / "reading-list.tsv").open() as f:
        entries = list(csv.DictReader(f, delimiter="\t"))
    if args.ids:
        selected = set(args.ids.split(","))
        entries = [e for e in entries if e["id"] in selected]
    manifest_path = OUT / "background-manifest.json"
    old = {e["id"]: e for e in json.loads(manifest_path.read_text())} if manifest_path.exists() else {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for e in ex.map(fetch, entries):
            old[e["id"]] = e
            manifest_path.write_text(json.dumps(list(old.values()), indent=2) + "\n")


if __name__ == "__main__":
    main()
