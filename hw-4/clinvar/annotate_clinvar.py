import argparse
import csv
import datetime as dt
import json
import pathlib
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
INPUT = ROOT / "student_31.vcf"
RAW = SCRIPT_DIR / "clinvar_raw"
OUT = SCRIPT_DIR / "clinvar_student31.tsv"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "bioinformatics-homework/1.0 (NCBI E-utilities)"})
    with urllib.request.urlopen(req, timeout=45) as response:
        data = response.read()
    time.sleep(0.36)  # NCBI unauthenticated rate guidance: <= 3 requests/sec.
    return data


def parse_vcf(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line:
            continue
        chrom, pos, ident, ref, alt, *_ = line.split("\t")
        if "," in alt:
            raise ValueError(f"Multiallelic input not supported: {line}")
        rows.append({"chrom": chrom.removeprefix("chr"), "pos": int(pos), "ref": ref.upper(), "alt": alt.upper(), "input_id": ident})
    return rows


def saved_search(term, stem):
    path = RAW / f"{stem}_esearch.json"
    if path.exists():
        parsed = json.loads(path.read_text())
    else:
        url = BASE + "esearch.fcgi?" + urllib.parse.urlencode({"db": "clinvar", "term": term, "retmode": "json", "retmax": "1000"})
        data = fetch(url)
        path.write_bytes(data)
        parsed = json.loads(data)
    result = parsed.get("esearchresult", {})
    err = result.get("errorlist", {})
    if "count" not in result or "idlist" not in result:
        raise RuntimeError(f"Incomplete ESearch response for {term}; saved in {path}")
    if err.get("fieldsnotfound"):
        raise RuntimeError(f"NCBI did not recognize search field for {term}: {err['fieldsnotfound']}; saved in {path}")
    if int(result["count"]) != len(result["idlist"]):
        raise RuntimeError(f"ESearch returned {len(result['idlist'])} of {result['count']} results for {term}; increase retmax before annotation")
    return parsed


def fetch_record(uid):
    path = RAW / f"uid_{uid}.xml"
    parsed = None
    if path.exists():
        try:
            parsed = ET.parse(path).getroot()
        except ET.ParseError:
            parsed = None
    # ESearch's numeric IDs must be passed to EFetch as Variation IDs. Without
    # is_variationid=1, EFetch may return an unrelated VCV. Validate that the
    # returned VariationID equals the requested ID before reusing a cache.
    archives = ([parsed] if parsed is not None and parsed.tag == "VariationArchive" else []) + (parsed.findall(".//VariationArchive") if parsed is not None else [])
    if parsed is None or not any(a.get("VariationID") == str(uid) for a in archives):
        url = BASE + "efetch.fcgi?" + urllib.parse.urlencode({"db": "clinvar", "id": uid, "rettype": "vcv", "retmode": "xml", "is_variationid": "1"})
        path.write_bytes(fetch(url))
    try:
        root = ET.parse(path).getroot()
        archives = ([root] if root.tag == "VariationArchive" else []) + root.findall(".//VariationArchive")
        if not any(a.get("VariationID") == str(uid) for a in archives):
            raise RuntimeError(f"NCBI returned no matching VariationArchive (VariationID={uid}) for ESearch UID {uid}; saved in {path}")
        return root
    except ET.ParseError as exc:
        raise RuntimeError(f"Invalid saved NCBI XML {path}: {exc}") from exc


def text_values(elem, tag):
    return list(dict.fromkeys((x.text or "").strip() for x in elem.iter(tag) if (x.text or "").strip()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=pathlib.Path, default=INPUT)
    ap.add_argument("--output", type=pathlib.Path, default=OUT)
    args = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    rows = parse_vcf(args.input)
    results = []
    for row in rows:
        rsid = row["input_id"]
        if not rsid.startswith("rs"):
            raise ValueError(f"Expected dbSNP rsID: {rsid}")
        # Search both the supplied dbSNP alias and the exact chromosomal
        # position. Either route can surface a VCV candidate; neither is
        # considered a match until GRCh38 position and REF/ALT are confirmed.
        search = saved_search(rsid, rsid)
        coord_term = f"{row['chrom']}[chr] AND {row['pos']}[chrpos38]"
        coord_search = saved_search(coord_term, f"chrpos38_{row['chrom']}_{row['pos']}_{row['ref']}_{row['alt']}")
        sr = search.get("esearchresult", {})
        cr = coord_search.get("esearchresult", {})
        uid_list = list(dict.fromkeys(sr["idlist"] + cr["idlist"]))
        candidates = []
        for uid in uid_list:
            root = fetch_record(uid)
            archives = ([root] if root.tag == "VariationArchive" else []) + root.findall(".//VariationArchive")
            for va in archives:
                varid = va.get("VariationID", "")
                vcv = va.get("Accession", "")
                name = va.get("VariationName", "")
                # Match only the primary called allele, or an IncludedRecord
                # allele for evidence-only VCVs. Do not traverse SCV assertions
                # or nested haplotypes and accidentally borrow their locations.
                locs = va.findall("./ClassifiedRecord/SimpleAllele/Location/SequenceLocation")
                locs.extend(va.findall(".//IncludedRecord/SimpleAllele/Location/SequenceLocation"))
                for loc in locs:
                    if loc.get("Assembly") != "GRCh38":
                        continue
                    # ClinVar's VCF-normalized allele fields are explicit for SNVs.
                    pos = loc.get("positionVCF") or loc.get("start")
                    ref = loc.get("referenceAlleleVCF") or loc.get("referenceAllele")
                    alt = loc.get("alternateAlleleVCF") or loc.get("alternateAllele")
                    chrom = (loc.get("Chr") or "").removeprefix("chr")
                    if (chrom, pos, ref, alt) == (row["chrom"], str(row["pos"]), row["ref"], row["alt"]):
                        # Use the overall classified significance where present; include all distinct statements.
                        # Aggregate VCV classification (e.g. conflicting classifications)
                        # rather than conflating it with all individual SCV assertions.
                        classes = va.findall(".//ClassifiedRecord/Classifications/*")
                        preferred = next((x for x in classes if x.tag in ("GermlineClassification", "OncogenicityClassification", "SomaticClinicalImpact")), None)
                        if preferred is None and classes:
                            preferred = classes[0]
                        # Use the aggregate classification's own direct children;
                        # descendant Description nodes include per-condition results.
                        descriptions = [preferred.findtext("Description", default="").strip()] if preferred is not None and preferred.findtext("Description", default="").strip() else []
                        review = [preferred.findtext("ReviewStatus", default="").strip()] if preferred is not None and preferred.findtext("ReviewStatus", default="").strip() else []
                        conditions = []
                        for rcv in va.findall(".//RCVList/RCVAccession"):
                            title = rcv.get("Title", "")
                            if " AND " in title:
                                conditions.append(title.rsplit(" AND ", 1)[1])
                        candidates.append({
                            "clinvar_variation_id": varid,
                            "clinvar_significance": "; ".join(descriptions),
                            "clinvar_review_status": "; ".join(review),
                            "clinvar_conditions": "; ".join(dict.fromkeys(conditions)),
                            "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{varid}/",
                            "clinvar_accession": vcv,
                            "clinvar_name": name,
                            "ncbi_uid": uid,
                        })
        if candidates:
            # Preserve all distinct VCVs matching the exact normalized allele.
            dedup = {c["clinvar_variation_id"]: c for c in candidates}
            matches = sorted(dedup.values(), key=lambda c: c["clinvar_variation_id"])
            combined = {}
            for key in matches[0]:
                vals = list(dict.fromkeys(c[key] for c in matches if c.get(key)))
                combined[key] = "; ".join(vals)
            results.append({**row, "clinvar_status": "found", **combined})
        else:
            search_url = "https://www.ncbi.nlm.nih.gov/clinvar/?" + urllib.parse.urlencode({"term": coord_term})
            results.append({**row, "clinvar_status": "not_found", "clinvar_variation_id": "", "clinvar_significance": "", "clinvar_review_status": "", "clinvar_conditions": "", "clinvar_url": search_url, "clinvar_search_url": search_url, "clinvar_accession": "", "clinvar_name": "", "ncbi_uid": ""})
    for result in results:
        result.setdefault("clinvar_search_url", "https://www.ncbi.nlm.nih.gov/clinvar/?" + urllib.parse.urlencode({"term": f"{result['chrom']}[chr] AND {result['pos']}[chrpos38]"}))
    fields = ["chrom", "pos", "ref", "alt", "input_id", "clinvar_status", "clinvar_variation_id", "clinvar_significance", "clinvar_review_status", "clinvar_conditions", "clinvar_url", "clinvar_search_url", "clinvar_accession", "clinvar_name", "ncbi_uid"]
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    settings = {
        "source": "NCBI ClinVar via Entrez E-utilities (esearch.fcgi by rsID and chromosome/position [chrpos38], then efetch.fcgi rettype=vcv retmode=xml is_variationid=1)",
        "input": str(args.input.relative_to(ROOT) if args.input.is_relative_to(ROOT) else args.input),
        "assembly": "GRCh38",
        "coordinate_search_field": "NCBI ClinVar ESearch [chrpos38], confirmed in ESearch querytranslation; EInfo also lists C37 as the separate GRCh37 field.",
        "matching": "Exact chromosome, VCF position, REF and ALT compared to ClinVar SequenceLocation Assembly=GRCh38 VCF-normalized fields; rsID and explicit chrpos38 searches only retrieve candidates.",
        "empty_search": "count=0 is treated as not_found; HTTP/XML/JSON errors abort the run and do not become not_found.",
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "output": str(args.output.relative_to(ROOT) if args.output.is_relative_to(ROOT) else args.output),
        "rows": len(results),
        "found": sum(r["clinvar_status"] == "found" for r in results),
        "not_found": sum(r["clinvar_status"] == "not_found" for r in results),
    }
    (RAW / "settings.json").write_text(json.dumps(settings, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(settings, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
