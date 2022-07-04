
# build LCA ANI table from prefetch csv(s)
import sys
import argparse
import sqlite3
import csv
import numpy as np
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify
from sourmash.distance_utils import containment_to_distance


def get_lineage(ident, tax_assign):
    try:
        lineage = tax_assign[ident]
    except KeyError:
        if "GCF" in ident:
            ident = ident.replace("GCF", "GCA")
        elif "GCA" in ident:
            ident = ident.replace("GCA", "GCF")
        #print("new ident: ", ident)
        # ARGH, having 202 --> 207 lineage matching issues. just let it be None and notify
        lineage = tax_assign.get(ident, None)
    return lineage

def lineages_lca(linA, linB):
    for rank in tax_utils.ascending_taxlist():
        if lca_utils.is_lineage_match(linA, linB, rank):
            lca = lca_utils.pop_to_rank(linA, rank)
            return lca

def get_lca(id1, id2, lineages):
    lin1 = get_lineage(id1, lineages)
    lin2 = get_lineage(id2, lineages)
    if lin1 is None or lin2 is None:
        if lin1 is None:
            notify(f"{id1} is not in taxonomy files")
        if lin2 is None:
            notify(f"{id2} is not in taxonomy files")
        return None
    return lineages_lca(lin1,lin2)

#def get_avg_contanment_ani(query_bp, ):

def main(args):
    # load in taxonomy
    tax_assign = tax_utils.MultiLineageDB.load(args.taxonomy_csvs, keep_identifier_versions=False)

    # handle file input
    prefetch_csvs = args.prefetch_csvs
    if args.from_file:
        ff_csvs = [x.strip() for x in open(args.from_file, 'r')]
        prefetch_csvs += ff_csvs

    # read in each file and load into table
    with open(args.output_csv, 'w') as outF:
        fields = ["comparison", "query_name", "match_name", "lca_rank",
                  "lca_lineage", "query_ani", "match_ani", "avg_ani",
                  "query_containment", "match_containment", "avg_containment"]

        writer = csv.writer(outF)
        writer.writerow(fields)

        comparisons = set()
        for inF in prefetch_csvs:
            with open(inF, 'r') as pf:
                prefetch_r = csv.DictReader(pf)
                for n, row in enumerate(prefetch_r):
                    if n % 10000 == 0:
                        notify(f"row {n}")
                    query_name = row['query_name']
                    query_id = tax_utils.get_ident(query_name)
                    match_name = row['match_name']
                    match_id = tax_utils.get_ident(match_name)

                    comparison_name = f"{query_id}_x_{match_id}"
                    reverse_name = f"{match_id}_x_{query_id}"

                    if any([comparison_name in comparisons, reverse_name in comparisons]):
                        # avoid dupes
                        continue
                    comparisons.add(comparison_name)

                    lca_lin = get_lca(query_id, match_id, tax_assign)
                    if lca_lin is None:
                        # if missing lineage, can't get LCA. Skip.
                        continue
                    q_containment = float(row['f_match_query'])
                    m_containment = float(row['f_query_match'])
                    # depending on version of prefetch, might not have avg contain -- recalc here.
                    avg_containment = np.mean([q_containment, m_containment])

                    if args.recalculate_ani:
                        ksize = int(row['ksize'])
                        scaled = int(row['scaled'])
                        n_unique_kmers = int(row['query_bp'])
                        query_dist = containment_to_distance(q_containment, ksize, scaled, n_unique_kmers=n_unique_kmers).dist
                        match_dist = containment_to_distance(m_containment, ksize, scaled, n_unique_kmers=n_unique_kmers).dist
                        # don't let any ANI values get zeroed out --> estimate independtly
                        query_ani = 1-query_dist
                        match_ani = 1-match_dist
                    else:

                        query_ani = row['query_ani']
                        match_ani = row['match_ani']

                    avg_ani = np.mean([float(query_ani), float(match_ani)])

                    # write csv
                    writer.writerow([comparison_name, query_id, match_id, lca_lin[-1].rank, lca_lin[-1].name, query_ani, match_ani, avg_ani, q_containment, m_containment, avg_containment])



if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('prefetch_csvs', nargs='*')
    p.add_argument('--from-file', '--prefetch-from-file', help="file containing paths to prefetch csvs")
    p.add_argument('-t', '--taxonomy-csvs', nargs='+', help='taxonomy information', required=True)
    p.add_argument('-o', '--output-csv', required=True, help='output csv')
    p.add_argument('-r', '--recalculate-ani', action='store_true')
    args = p.parse_args()
    sys.exit(main(args))
