#!/usr/bin/env python3
"""Fetch gnomAD v4 variant annotations from the official GraphQL API."""
import csv
import datetime as dt
import json
import pathlib
import sys
import urllib.error
import urllib.request

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
VCF = ROOT / "student_31.vcf"
RAW = SCRIPT_DIR / "gnomad_raw"
OUT = SCRIPT_DIR / "gnomad_student31.tsv"
ENDPOINT = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"


def post(query, variables=None):
    payload = json.dumps({"query": query, "variables": variables or {}}).encode()
    req = urllib.request.Request(ENDPOINT, data=payload, headers={"Content-Type": "application/json", "User-Agent": "bioinformatics-homework-annotation/1.0"})
    with urllib.request.urlopen(req, timeout=45) as response:
        return response.status, json.loads(response.read())


def read_variants():
    rows = []
    with VCF.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, input_id, ref, alt = line.rstrip("\n").split("\t")[:5]
            rows.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, "input_id": input_id})
    return rows


def main():
    RAW.mkdir(exist_ok=True)
    schema_query = '{ __type(name: "VariantDetails") { fields { name type { name kind ofType { name kind ofType { name kind } } } } } }'
    schema_status, schema = post(schema_query)
    detail = schema.get("data", {}).get("__type", {})
    detail_fields = {f["name"] for f in detail.get("fields", [])}
    nested_types = {}
    for field in detail.get("fields", []):
        t = field["type"]
        while t and not t.get("name"):
            t = t.get("ofType")
        if t and field["name"] in {"exome", "genome", "joint"}:
            nested_types[field["name"]] = t["name"]
    aliases = " ".join(f'{name}: __type(name: "{typename}") {{ fields {{ name }} }}' for name, typename in nested_types.items())
    nested_query = "{ " + aliases + " }" if aliases else None
    nested_data = post(nested_query)[1] if nested_query else {"data": {}}
    type_fields = {name: {x["name"] for x in nested_data.get("data", {}).get(name, {}).get("fields", [])} for name in nested_types}
    (RAW / "schema_introspection.json").write_text(json.dumps({"variant_details": {"http_status": schema_status, "query": schema_query, "body": schema}, "nested_types_query": nested_query, "nested_types": nested_data, "fields_by_dataset_slice": {k: sorted(v) for k, v in type_fields.items()}}, indent=2) + "\n")
    has_joint = "joint" in nested_types
    joint_fields = type_fields.get("joint", set())
    joint_query_field = " joint { ac an" + (" filters" if "filters" in joint_fields else "") + " }" if has_joint else ""
    filter_fields = {name: (" filters" if "filters" in type_fields.get(name, set()) else "") for name in ("exome", "genome")}
    q = "query VariantLookup($variantId: String!, $dataset: DatasetId!) { variant(variantId: $variantId, dataset: $dataset) { variant_id reference_genome exome { ac an af" + filter_fields["exome"] + " } genome { ac an af" + filter_fields["genome"] + " }" + joint_query_field + " } }"
    now = dt.datetime.now(dt.timezone.utc).isoformat()
    settings = {"retrieved_at_utc": now, "endpoint": ENDPOINT, "dataset": DATASET, "reference": "GRCh38", "source_vcf": str(VCF.relative_to(ROOT)), "query_template": q, "joint_available_in_variant_details_schema": has_joint, "available_slice_fields": {k: sorted(v) for k, v in type_fields.items()}}
    (RAW / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    results = []
    for i, row in enumerate(read_variants(), 1):
        vid = f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}"
        query_id = f"{i:02d}_{vid}"
        request = {"query": q, "variables": {"variantId": vid, "dataset": DATASET}}
        (RAW / f"{query_id}.request.json").write_text(json.dumps(request, indent=2) + "\n")
        status, data = "ok", {}
        try:
            http_status, data = post(q, request["variables"])
            (RAW / f"{query_id}.response.json").write_text(json.dumps({"http_status": http_status, "body": data}, indent=2) + "\n")
            if not data.get("data", {}).get("variant"):
                errors = data.get("errors", [])
                status = "not_found" if errors and all(e.get("message") == "Variant not found" for e in errors) else "api_error"
            elif data.get("errors"):
                status = "api_error"
            elif data["data"]["variant"].get("variant_id") != vid or data["data"]["variant"].get("reference_genome") != "GRCh38":
                status = "identity_mismatch"
            else:
                status = "found"
        except urllib.error.HTTPError as e:
            status = "api_error"
            try:
                data = json.loads(e.read())
            except Exception:
                data = {"error_type": type(e).__name__, "error": str(e)}
            (RAW / f"{query_id}.response.json").write_text(json.dumps({"http_status": e.code, "body": data}, indent=2) + "\n")
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            status = "request_error"
            data = {"error_type": type(e).__name__, "error": str(e)}
            (RAW / f"{query_id}.response.json").write_text(json.dumps(data, indent=2) + "\n")
        results.append((row, status, data.get("data", {}).get("variant") if status == "found" else None))
        print(f"{i}/50 {vid}: {status}", flush=True)
    fields = ["chrom", "pos", "ref", "alt", "input_id", "gnomad_status", "gnomad_dataset", "gnomad_exome_af", "gnomad_exome_ac", "gnomad_exome_an", "gnomad_genome_af", "gnomad_genome_ac", "gnomad_genome_an", "gnomad_exome_filters", "gnomad_genome_filters"] + (["gnomad_joint_af", "gnomad_joint_ac", "gnomad_joint_an", "gnomad_joint_filters"] if has_joint else []) + ["gnomad_url"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row, status, variant in results:
            exome = variant.get("exome") or {} if variant else {}
            genome = variant.get("genome") or {} if variant else {}
            joint = variant.get("joint") or {} if variant else {}
            joint_af = (joint.get("ac") / joint.get("an")) if joint.get("ac") is not None and joint.get("an") else None
            record = {**row, "gnomad_status": status, "gnomad_dataset": DATASET,
                      "gnomad_exome_af": exome.get("af"), "gnomad_exome_ac": exome.get("ac"), "gnomad_exome_an": exome.get("an"),
                      "gnomad_genome_af": genome.get("af"), "gnomad_genome_ac": genome.get("ac"), "gnomad_genome_an": genome.get("an"),
                      "gnomad_exome_filters": json.dumps(exome.get("filters")) if exome.get("filters") is not None else None,
                      "gnomad_genome_filters": json.dumps(genome.get("filters")) if genome.get("filters") is not None else None,
                      **({"gnomad_joint_af": joint_af, "gnomad_joint_ac": joint.get("ac"), "gnomad_joint_an": joint.get("an"), "gnomad_joint_filters": json.dumps(joint.get("filters")) if joint.get("filters") is not None else None} if has_joint else {}),
                      "gnomad_url": f"https://gnomad.broadinstitute.org/variant/{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}?dataset={DATASET}"}
            w.writerow(record)
    if any(status not in {"found", "not_found"} for _, status, _ in results):
        print("Annotation incomplete: inspect error statuses and saved responses before using the TSV.", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
