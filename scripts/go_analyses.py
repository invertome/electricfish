import os
import argparse
import pandas as pd
import logging
import gzip
import shutil
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from gseapy import enrichr, barplot
from gseapy.plot import gseaplot
from Bio import KEGG
from Bio.KEGG.KGML import KGML_parser

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Function to read gene to GO term associations
def read_ncbi_gene2go(gene2go, taxids, **kws):
    return Gene2GoReader(gene2go, taxids=taxids, **kws).get_id2gos(namespace='BP')

# Function to run Enrichr for gene list
def run_enrichr(gene_list, description, gene_sets, taxid, output_dir):
    return enrichr(gene_list=gene_list, description=description, gene_sets=gene_sets, organism=str(taxid), outdir=output_dir)

# Function to plot KEGG pathway
def plot_pathway(kegg_id, output_dir):
    pathway = KEGG.get(kegg_id).read()
    KEGG.draw(pathway, filename=os.path.join(output_dir, f'{kegg_id}_pathway.png'))

# Function to extract hit IDs from BLAST results
def extract_blast_hit_ids(blast_results_file):
    # Load the BLAST results into a DataFrame
    blast_results = pd.read_csv(blast_results_file, sep="\t", header=None)

    # Extract the hit IDs from the second column
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
    parser.add_argument('--taxid', type=int, required=True, help='Taxonomy ID for species')
    parser.add_argument('--background', type=str, help='Path to background genes file')
    parser.add_argument('--blast_results', type=str, help='Path to BLAST results file')
    parser.add_argument('-kegg', action='store_true', help='Conduct KEGG pathway enrichment analysis')
    parser.add_argument('--plot_pathway', type=str, help='KEGG pathway ID to plot')
    args = parser.parse_args()

    # Create the output directory if it doesn't exist
    os.makedirs(args.o, exist_ok=True)

    # Set up logging to file
    logging.basicConfig(filename=os.path.join(args.o, 'run_log.log'), level=logging.INFO)

    # Start logging
    logger.info('Script started.')

    # Download necessary files if they don't exist
    go_obo = os.path.join(args.o, "go-basic.obo")
    gene2go = os.path.join(args.o, "gene2go.gz")
    if not os.path.exists(go_obo):
        os.system(f'wget http://purl.obolibrary.org/obo/go/go-basic.obo -O {go_obo}')
    if not os.path.exists(gene2go):
        os.system(f'wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz -O {gene2go}')

    # Load the GO DAG
    obodag = GODag(go_obo)
    logger.info('GO DAG loaded.')

    # Read gene to GO term associations
    # Unzip the gzipped file
    with gzip.open(gene2go, 'rb') as f_in:
        with open(gene2go[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    geneid2gos_species = read_ncbi_gene2go(gene2go[:-3], taxids=[args.taxid])
    logger.info('Gene to GO term associations read.')

    # Extract hit IDs from BLAST results file if provided and no other background file is given
    if args.blast_results and not args.background:
        args.background = extract_blast_hit_ids(args.blast_results)

    # For each input file
    for input_file in args.i:
        # Read the CSV file
        df = pd.read_csv(input_file)
        df = df.drop_duplicates('ProteinID')  # remove duplicate entries
        df = df.dropna(subset=['ProteinID'])  # remove entries with missing ProteinIDs

        # Create a list of gene NCBI IDs
        gene_ncbi = df['ProteinID'].tolist()

        # Conduct GO enrichment analysis
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
        logger.info(f'Significant results saved for file {input_file}.')

        # Conduct KEGG analysis if specified
        if args.kegg:
            logger.info('Running KEGG pathway enrichment analysis...')
            enr_res = run_enrichr(gene_list=gene_ncbi, description='pathway', gene_sets='KEGG_Pathways', taxid=args.taxid, output_dir=args.o)
            enr_res.res2d.to_csv(os.path.join(args.o, 'KEGG_enrichment.csv'))

            # Plot KEGG pathway enrichment
            gseaplot(enr_res.ranking, term='KEGG Pathway Enrichment', ofname=os.path.join(args.o, 'KEGG_enrichment.png'))
            logger.info('KEGG pathway enrichment analysis completed.')

        # Plot specific KEGG pathway if specified
        if args.plot_pathway:
            plot_pathway(args.plot_pathway, args.o)
            logger.info(f'KEGG pathway {args.plot_pathway} plotted.')

    logger.info('Script finished.')
