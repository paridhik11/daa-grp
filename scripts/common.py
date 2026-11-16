"""Shared helpers for benchmarks."""
import json
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
BIN = os.path.join(ROOT, "approx")
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)


def load_cm():
    p = os.path.join(RES, "costmodel.json")
    if os.path.exists(p):
        return json.load(open(p))
    return None


def cm_env():
    cm = load_cm()
    return f"{cm['c_hit']},{cm['c_win']},{cm['c_bp']},{cm['a_scan']},{cm['c_sel']}" if cm else None


def run_bin(algo, text_path, pats, k, nofallback=False, dump=True, tag="x", times=False):
    pf = os.path.join(DATA, f"_p_{tag}_{algo}.txt")
    open(pf, "w").write("\n".join(pats) + "\n")
    env = dict(os.environ)
    env["APPROX_NOFALLBACK"] = "1" if nofallback else "0"
    c = cm_env()
    if c:
        env["APPROX_CM"] = c
    args = [BIN, algo, text_path, pf, str(k)] + (["--dump"] if dump else []) + (["--times"] if times else [])
    r = subprocess.run(args, capture_output=True, text=True, check=True, env=env)
    lines = r.stdout.splitlines()
    info = {a: float(b) for a, b in re.findall(r"(\w+)=([\d.]+)", lines[0])[1:]}
    sets, tq = {}, []
    for l in lines[1:]:
        t = l.split()
        if t[0] == "M":
            sets.setdefault(int(t[1]), set()).add(int(t[2]))
        elif t[0] == "T":  # per-pattern: ns hits verified_bp pred_filter_ns pred_scan_ns fallbacks
            tq.append(dict(ns=float(t[2]), hits=float(t[3]), bp=float(t[4]), pf=float(t[5]), ps=float(t[6]), fb=float(t[7])))
    return (info, sets, tq) if times else (info, sets)
