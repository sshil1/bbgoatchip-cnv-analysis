#!/usr/bin/env python3
"""
07_rebuild_string_network.py

PURPOSE
  No saved script for the manuscript's original STRING network figure
  (CD244-CD48 edge, C3/CD48 hubs) could be found anywhere on the server --
  it was almost certainly built by pasting a gene list into the string-db.org
  web interface interactively, and which gene list (254-gene fixed set vs.
  2,856-gene genome-wide set) is not recorded or remembered.

  Rather than guess at the original figure's provenance, this rebuilds the
  network directly and reproducibly via STRING's REST API, using the
  CORRECTED, reciprocal-overlap-confirmed 114-gene set from
  05_rebuild_fixed_enrichment.py. Whatever this produces is defensible on
  its own terms and replaces reliance on the untraceable original.

  Also queries the genome-wide gene list (if you provide one) as a
  comparison, so you can see whether CD244-CD48 survive there even if they
  don't in the corrected fixed set (we already know CD244 and CD48 are both
  LOST from the corrected fixed-114 set -- see 05's output -- so this run's
  main job against that list is to see what DOES remain, not to expect the
  original edge to reappear).

INPUTS
  gene_list_file   plain text, one gene symbol per line (e.g.
                    corrected_fixed_genes_246cnvr.txt from script 05)

OUTPUT
  string_network_edges_<label>.tsv   -- all interactions above the score
                                         threshold, with combined_score
  string_network_summary_<label>.txt -- node degree ranking (hub genes),
                                         and explicit CD244/CD48/C3 status
  A PNG network image is also downloaded if --save-image is passed.

USAGE:
  python3 07_rebuild_string_network.py corrected_fixed_genes_246cnvr.txt \
      --label corrected_114 --species 9925 --score 400 --save-image

  # for comparison, if you have/rebuild the genome-wide gene list too:
  python3 07_rebuild_string_network.py genome_wide_genes_2856.txt \
      --label genome_wide_2856 --species 9925 --score 400

  --species 9925 = Capra hircus NCBI taxonomy ID (confirmed correct for
  this project's reference genome, GCF_001704415.2 ARS1.2).
  --score 400 = STRING's standard "medium confidence" threshold (0-1000
  scale); this is the same default STRING's own web interface uses, so a
  network built this way is comparable to what the original web-interface
  session most likely produced.
"""
import argparse
import io
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

STRING_API = "https://string-db.org/api"


def api_post(endpoint, params, fmt="tsv"):
    url = f"{STRING_API}/{fmt}/{endpoint}"
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode()


def get_string_ids(genes, species):
    """Map gene symbols to STRING identifiers. Returns dict: input_gene -> string_id"""
    params = {
        "identifiers": "\r".join(genes),
        "species": species,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": "bbg_cnv_manuscript_rebuild",
    }
    raw = api_post("get_string_ids", params)
    mapping = {}
    unmapped = []
    lines = [l for l in raw.strip().split("\n") if l]
    if not lines:
        return mapping, genes
    header = lines[0].split("\t")
    for line in lines[1:]:
        f = line.split("\t")
        row = dict(zip(header, f))
        query = row.get("queryItem") or row.get("queryIndex")
        sid = row.get("stringId")
        if query and sid:
            mapping[query] = sid
    unmapped = [g for g in genes if g not in mapping]
    return mapping, unmapped


def get_network(string_ids, species, score):
    params = {
        "identifiers": "\r".join(string_ids),
        "species": species,
        "required_score": score,
        "caller_identity": "bbg_cnv_manuscript_rebuild",
    }
    raw = api_post("network", params)
    lines = [l for l in raw.strip().split("\n") if l]
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        f = line.split("\t")
        rows.append(dict(zip(header, f)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gene_list_file")
    ap.add_argument("--label", default="run")
    ap.add_argument("--species", default="9925", help="NCBI taxid, 9925=Capra hircus")
    ap.add_argument("--score", default="400", help="STRING required_score, 0-1000")
    ap.add_argument("--save-image", action="store_true")
    args = ap.parse_args()

    with open(args.gene_list_file) as f:
        genes = [l.strip() for l in f if l.strip()]
    print(f"Loaded {len(genes)} genes from {args.gene_list_file}")

    print("Mapping gene symbols to STRING IDs...")
    mapping, unmapped = get_string_ids(genes, args.species)
    print(f"  Mapped: {len(mapping)}/{len(genes)}  Unmapped: {len(unmapped)}")
    if unmapped:
        print(f"  Unmapped genes (first 20): {unmapped[:20]}")

    if not mapping:
        sys.exit("No genes mapped to STRING IDs -- check --species and gene symbols.")

    time.sleep(1)  # be polite to the API between calls
    print(f"Querying network at required_score={args.score}...")
    edges = get_network(list(mapping.values()), args.species, args.score)
    print(f"  {len(edges)} edges returned")

    edge_path = f"string_network_edges_{args.label}.tsv"
    with open(edge_path, "w") as out:
        if edges:
            cols = list(edges[0].keys())
            out.write("\t".join(cols) + "\n")
            for e in edges:
                out.write("\t".join(e.get(c, "") for c in cols) + "\n")
    print(f"Edge list written: {edge_path}")

    # degree / hub ranking
    degree = defaultdict(int)
    named_edges = []
    for e in edges:
        a = e.get("preferredName_A", e.get("stringId_A", "?"))
        b = e.get("preferredName_B", e.get("stringId_B", "?"))
        degree[a] += 1
        degree[b] += 1
        named_edges.append((a, b, e.get("score", "")))

    summary_path = f"string_network_summary_{args.label}.txt"
    with open(summary_path, "w") as out:
        out.write(f"STRING network rebuild: {args.label}\n")
        out.write(f"Gene list: {args.gene_list_file} ({len(genes)} genes, "
                   f"{len(mapping)} mapped)\n")
        out.write(f"Species: {args.species}  Score threshold: {args.score}\n")
        out.write(f"Edges: {len(edges)}\n\n")
        out.write("Node degree ranking (potential hubs, highest first):\n")
        for name, d in sorted(degree.items(), key=lambda x: -x[1])[:20]:
            out.write(f"  {name}: {d}\n")
        out.write("\nKey marker gene status:\n")
        for marker in ["CD244", "CD48", "C3"]:
            in_list = marker in genes
            in_net = marker in degree
            out.write(f"  {marker}: in gene list={in_list}, "
                      f"appears in network={in_net}"
                      + (f", degree={degree[marker]}" if in_net else "") + "\n")
        out.write("\nDirect CD244-CD48 edge check:\n")
        cd_edge = [e for e in named_edges
                   if {e[0], e[1]} == {"CD244", "CD48"}]
        out.write(f"  CD244-CD48 edge present: {bool(cd_edge)}"
                   + (f" (score={cd_edge[0][2]})" if cd_edge else "") + "\n")

    print(f"Summary written: {summary_path}")
    with open(summary_path) as f:
        print("\n" + f.read())

    if args.save_image:
        img_url = (f"{STRING_API}/image/network?identifiers="
                    + urllib.parse.quote("\r".join(mapping.values()))
                    + f"&species={args.species}&required_score={args.score}"
                    + "&caller_identity=bbg_cnv_manuscript_rebuild")
        img_path = f"string_network_{args.label}.png"
        urllib.request.urlretrieve(img_url, img_path)
        print(f"Network image saved: {img_path}")


if __name__ == "__main__":
    main()
