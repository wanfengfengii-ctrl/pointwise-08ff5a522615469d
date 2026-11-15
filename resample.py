from matrix import unpack

def materialize(package):
    source = {}
    for row in package["observations"]:
        ident = row["id"]
        if ident in source:
            raise ValueError("duplicate source id")
        source[ident] = {k: unpack(row[k]) for k in ("time", "value", "variance")}
        if source[ident]["variance"] <= 0:
            raise ValueError("nonpositive variance")
    windows = []
    seen = set()
    for spec in package["windows"]:
        if spec["id"] in seen:
            raise ValueError("duplicate window id")
        seen.add(spec["id"])
        weights = {k: unpack(v) for k, v in spec["weights"].items()}
        if not weights or any(k not in source or v < 0 for k, v in weights.items()):
            raise ValueError("invalid weight")
        if sum(weights.values()) != 1:
            raise ValueError("weights must sum to one")
        weights = {k: v for k, v in weights.items() if v}
        windows.append({
            "id": spec["id"], "weights": weights,
            "time": sum(v * source[k]["time"] for k, v in weights.items()),
            "value": sum(v * source[k]["value"] for k, v in weights.items())
        })
    return source, windows
