#!/usr/bin/env python3
"""Fetch real genomes from any reachable source; clean to ACGT; record provenance.

Datasets are named by what they really are (taken from the FASTA header), never assumed.
Writes data/<name>.txt (raw ACGT, one line), data/SOURCES.txt, data/fetch_log.txt.
data/genome.txt = largest real bacterial genome found, else the synthetic one (labeled in SOURCES.txt).
"""
import hashlib
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RAW = os.path.join(DATA, "raw")
os.makedirs(RAW, exist_ok=True)
GH = "https://raw.githubusercontent.com"

# name -> (kind, candidate URLs in priority order)
SOURCES = {
    "ecoli_k12_mg1655": ("bacterial", [
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta",
        "https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3",
        "https://hgdownload.soe.ucsc.edu/goldenPath/eschColi_K12/bigZips/eschColi_K12.fa.gz",
    ]),
    "ecoli536": ("bacterial", [f"{GH}/BenLangmead/bowtie/master/genomes/NC_008253.fna"]),
    "sarscov2": ("viral", [
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta",
        f"{GH}/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta",
    ]),
    "lambda": ("viral", [f"{GH}/BenLangmead/bowtie2/master/example/reference/lambda_virus.fa"]),
    "mt_human": ("organelle", [f"{GH}/lh3/minimap2/master/test/MT-human.fa"]),
    "mt_orang": ("organelle", [f"{GH}/lh3/minimap2/master/test/MT-orang.fa"]),
}
LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def download(url, out):
    r = subprocess.run(["curl", "-sS", "-m", "120", "-L", "-o", out, "-w", "%{http_code}", url],
                       capture_output=True, text=True)
    code = r.stdout.strip()
    ok = code == "200" and os.path.exists(out) and os.path.getsize(out) > 1000
    return ok, code if code else (r.stderr.strip()[:80] or "fail")


def parse_fasta(path):
    import gzip
    op = gzip.open if open(path, "rb").read(2) == b"\x1f\x8b" else open
    header, seq = None, []
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                else:
                    break  # first record only
            else:
                seq.append(line.strip())
    return header, "".join(seq)


def main():
    found = {}
    for name, (kind, urls) in SOURCES.items():
        for u in urls:
            out = os.path.join(RAW, name + ".fa")
            side = out + ".url"
            if os.path.exists(out) and os.path.getsize(out) > 1000 and os.path.exists(side):
                u = open(side).read().strip()
                log(f"[{name}] using cached {out} (from {u})")
                ok = True
            else:
                ok, code = download(u, out)
                log(f"[{name}] {'OK  ' if ok else 'FAIL'} {code} {u}")
                if ok:
                    open(side, "w").write(u + "\n")
            if ok:
                header, seq = parse_fasta(out)
                raw_len = len(seq)
                clean = "".join(c for c in seq.upper() if c in "ACGT")
                gc = (clean.count("G") + clean.count("C")) / max(1, len(clean))
                with open(os.path.join(DATA, name + ".txt"), "w") as f:
                    f.write(clean + "\n")
                found[name] = dict(kind=kind, url=u, header=header, raw_len=raw_len, len=len(clean), gc=gc,
                                   sha=hashlib.sha256(clean.encode()).hexdigest())
                break
        else:
            log(f"[{name}] ALL SOURCES FAILED")
    # K-12 claim sanity: only keep the label if header + length confirm
    if "ecoli_k12_mg1655" in found:
        e = found["ecoli_k12_mg1655"]
        if not ("K-12" in e["header"] and abs(e["len"] - 4_641_652) < 5000):
            log("ecoli_k12_mg1655 header/length do not confirm K-12 MG1655; removing label")
            os.rename(os.path.join(DATA, "ecoli_k12_mg1655.txt"), os.path.join(DATA, "unlabeled_ecoli.txt"))
            found["unlabeled_ecoli"] = found.pop("ecoli_k12_mg1655")
    # main text = largest real bacterial genome, else synthetic
    bact = [(v["len"], n) for n, v in found.items() if v["kind"] == "bacterial" or n.endswith("ecoli")]
    synth = os.path.join(DATA, "genome_synth.txt")  # created once by simulate.py
    main_name = max(bact)[1] if bact else None
    if main_name:
        shutil.copy(os.path.join(DATA, main_name + ".txt"), os.path.join(DATA, "genome.txt"))
        log(f"MAIN TEXT = {main_name} (largest real bacterial genome found)")
    else:
        if os.path.exists(synth):
            shutil.copy(synth, os.path.join(DATA, "genome.txt"))
        log("NO real bacterial genome reachable: main text = SYNTHETIC (labeled)")
    with open(os.path.join(DATA, "SOURCES.txt"), "w") as f:
        f.write(f"main text: {main_name or 'SYNTHETIC genome_synth.txt'}\n\n")
        for n, v in found.items():
            f.write(f"{n}\n  url:    {v['url']}\n  header: {v['header']}\n  length: {v['len']} (raw {v['raw_len']}), "
                    f"GC {v['gc']:.4f}\n  sha256: {v['sha']}\n")
        if os.path.exists(synth):
            f.write("genome_synth\n  SYNTHETIC random genome with injected repeats (scripts/simulate.py)\n")
    open(os.path.join(DATA, "fetch_log.txt"), "w").write("\n".join(LOG) + "\n")
    return 0 if found else 1


if __name__ == "__main__":
    sys.exit(main())
