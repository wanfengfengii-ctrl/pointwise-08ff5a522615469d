import argparse
import json
from pathlib import Path
from matrix import encode
from resample import materialize
from summary import summarize, merge
from trend import fit

def run(package, plan):
    source, windows = materialize(package)
    lookup = {w["id"]: w for w in windows}
    flat = [wid for group in plan for wid in group]
    if not plan or any(not group for group in plan) or len(flat) != len(set(flat)) or set(flat) != set(lookup):
        raise ValueError("blocks must partition all window ids")
    blocks = [summarize(source, [lookup[wid] for wid in group]) for group in plan]
    total = merge(blocks)
    return {"result": fit(total), "trace": {"windows": windows, "blocks": blocks, "merged": total}}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("package")
    parser.add_argument("plan")
    args = parser.parse_args()
    package = json.loads(Path(args.package).read_text())
    output = run(package, package["plans"][args.plan])
    print(json.dumps(encode(output), ensure_ascii=False, sort_keys=True))

if __name__ == "__main__":
    main()
