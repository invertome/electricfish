import os
import argparse
import pandas as pd
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from gseapy import enrichr, barplot
from Bio import KEGG

def read_ncbi_gene2go(gene2go, taxids, **kws):
    return Gene2GoReader(gene2go, taxids=taxids, **kws).get_id2gos(namespace='BP')

def run_enrichr(gene_list, description, gene_sets, organism, output_dir):
    return enrichr(gene_list=gene_list, description=description, gene_sets=gene_sets, organism=organism, outdir=output_dir)

def plot_pathway(kegg_id, output_dir):
    pathway = KEGG.get(kegg_id).read()
    KEGG.draw(pathway, filename=os.path.join(output_dir, f'{kegg_id}_pathway.png'))

def extract_blast_hit_ids(blast_results_file):
    # Load the BLAST results into a DataFrame
    blast_results = pd.read_csv(blast_results_file, header=None)

    # Extract the hit IDs
    hit_ids = blast_results[1].unique()

    # Return hit IDs as a list
    return hit_ids.tolist()

if __name__ == '__main__':
    # Handle argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, nargs='+', required=True, help='Input files (CSV format)')
    parser.add_argument('-o', type=str, required=True, help='Output directory')
    parser.add_argument('--alpha', type=float, default=0.05, help='Significance level')
    parser.add_argument('--method', type=str, default='fdr_bh', help='Method for multiple test correction')
    parser.add_argument('--species', type=str, required=True, help='Species for analysis')
    parser.add_argument('--taxid', type=int, required=True, help='Taxonomy ID for species')
    parser.add_argument('--background', type=str, help='Path to background genes file')
    parser.add_argument('--blast_results', type=str, help='Path to BLAST results file')
    parser.add_argument('-kegg', action='store_true', help='Conduct KEGG pathway enrichment analysis')
    parser.add_argument('--plot_pathway', type=str, help='KEGG pathway ID to plot')
    args = parser.parse_args()

    # Prepare the output directory
    os.makedirs(args.o, exist_ok=True)

    # Download necessary files if they don't exist
    go_obo = os.path.join(args.o, "go-basic.obo")
    gene2go = os.path.join(args.o, "gene2go.gz")
    if not os.path.exists(go_obo):
        os.system(f'wget http://purl.obolibrary.org/obo/go/go-basic.obo -O {go_obo}')
    if not os.path.exists(gene2go):
        os.system(f'wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz -O {gene2go}')

    obodag = GODag(go_obo)
    geneid2gos_species = read_ncbi_gene2go(gene2go, taxids=[args.taxid])

    # Extract hit IDs from BLAST results file if provided and no other background file is given
    if args.blast_results and not args.background:
        args.background = extract_blast_hit_ids(args.blast_results)

    # For each input file
    for input_file in args.i:
        # Read the CSV file
        df = pd.read_csv(input_file)
        df = df.drop_duplicates('ProteinID')  # remove duplicate entries
        df = df.dropna(subset=['ProteinID'])  # remove entries with missing ProteinIDs

        gene_ncbi = df['ProteinID'].tolist()

        # Running the GO enrichment analysis
        goeaobj = GOEnrichmentStudyNS(
            gene_ncbi,
            geneid2gos_species,
            obodag,
            propagate_counts=False,
            alpha=args.alpha,
            methods=[args.method]
        )

        goea_results_all = goeaobj.run_study(gene_ncbi)
        goea_results_sig = [r for r in goea_results_all if getattr(r, 'p_' + args.method) < args.alpha]

        # Save significant results as a CSV
        goeaobj.wr_tsv(os.path.join(args.o, f'{os.path.basename(input_file).split(".")[0]}_goea_significant.tsv'), goea_results_sig)

        # Check if KEGG analysis should be conducted
        if args.kegg:
            print("Running KEGG pathway enrichment analysis...")
            enr_res = run_enrichr(gene_list=gene_ncbi, description='pathway', gene_sets='KEGG_Pathways', organism=args.species, outdir=args.o)
            enr_res.res2d.to_csv(os.path.join(args.o, 'KEGG_enrichment.csv'))
            barplot(enr_res.res2d, title='KEGG Pathway Enrichment', ofname=os.path.join(args.o, 'KEGG_enrichment.png'))

        # Check if a KEGG pathway should be plotted
        if args.plot_pathway:
            plot_pathway(args.plot_pathway, args.o)
