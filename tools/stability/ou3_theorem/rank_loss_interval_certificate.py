"""Emit the fail-closed literal rank-loss interval certificate."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .rank_loss_literal_boxes import rank_loss_factor_certificate

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);p.add_argument("--max-depth",type=int,default=10)
    a=p.parse_args();z=rank_loss_factor_certificate(max_depth=a.max_depth)
    Path(a.output).write_text(json.dumps(z,indent=2,sort_keys=True)+"\n")
    print(json.dumps(z,sort_keys=True))
if __name__=="__main__":main()
