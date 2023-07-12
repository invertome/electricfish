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
from Bio import KEGG, Entrez
from Bio.KEGG.KGML import KGML_parser

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Function to read gene to GO term associations
def read_ncbi_gene2go(gene2go, taxids, gene_info_df, **kws):
    id2gos = Gene2GoReader(gene2go, taxids=taxids, **kws).get_id2gos(namespace='BP')
    id2gos = {ncbi_id: go_terms for ncbi_id, go_terms in id2gos.items() if ncbi_id in gene_info_df.index}
    return id2gos

def run_enrichr(gene_list, gene_sets, organism, output_dir):
    # Run Enrichr
    return enrichr(gene_list=gene_list, gene_sets=gene_sets, organism=organism, outdir=output_dir)

# Function to download and process gene2refseq file
def get_gene2refseq_file(refseq_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)  
    if not os.path.exists(refseq_file):
        logger.info('Downloading gene2refseq.gz...')
        os.system(f'wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz -O {refseq_file}')
    if not os.path.exists(os.path.join(output_dir, 'gene2refseq')):
        logger.info('Extracting gene2refseq.gz...')
        with gzip.open(refseq_file, 'rb') as f_in:
            with open(os.path.join(output_dir, 'gene2refseq'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    logger.info('gene2refseq file is ready.')

# Function to convert ProteinID to GeneID using gene2refseq file
def protein_to_gene(protein_id, gene2refseq_df):
    try:
        gene_id = gene2refseq_df.loc[protein_id, 'GeneID']
    except KeyError:
        gene_id = None
    return gene_id



# Function to plot KEGG pathway
def plot_pathway(kegg_id, output_dir):
    pathway = KEGG.get(kegg_id).read()
    KEGG.draw(pathway, filename=os.path.join(output_dir, f'{kegg_id}_pathway.png'))

# Function to extract hit IDs from BLAST results
def extract_blast_hit_ids(blast_results_file, gene_info_df):
    # Load the BLAST results into a DataFrame
    blast_results = pd.read_csv(blast_results_file, sep="\t", header=None)

    # Extract the hit IDs from the second column
    hit_ids = blast_results[1].unique()

    # Convert the hit IDs from NCBI IDs to Entrez IDs
    hit_ids = [ncbi_id for ncbi_id in hit_ids if ncbi_id in gene_info_df.index]

    # Return hit IDs as a list
    return hit_ids
    
# Function to download and extract gene_info file
# Define the mapping outside the function
taxid_to_organism = {
    9606: 'Homo_sapiens',
    10090: 'Mus_musculus',
    7955: 'Danio_rerio',
    7227: 'Drosophila_melanogaster',
    # Add more if needed
}

taxid_to_path = {
    9606: 'Mammalia/Homo_sapiens.gene_info.gz',
    10090: 'Mammalia/Mus_musculus.gene_info.gz',
    7955: 'Non-mammalian_vertebrates/Danio_rerio.gene_info.gz',
    7227: 'Invertebrates/Drosophila_melanogaster.gene_info.gz',
    # Add more if needed
}

def get_gene_info_file(taxid, gene_info_gz, output_dir):
    os.makedirs(output_dir, exist_ok=True)  
    if not os.path.exists(gene_info_gz):
        logger.info('Downloading gene_info.gz...')
        path = taxid_to_path.get(taxid, 'Invertebrates/All_Invertebrates.gene_info.gz')  # use a default path if taxid is not in the dictionary
        os.system(f'wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/{path} -O {gene_info_gz}')
    if not os.path.exists(os.path.join(output_dir, f'{organism}.gene_info')):
        logger.info('Extracting gene_info.gz...')
        with gzip.open(gene_info_gz, 'rb') as f_in:
            with open(os.path.join(output_dir, f'{organism}.gene_info'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    logger.info('Gene_info file is ready.')





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
    parser.add_argument('--ref_dir', type=str, required=True, help='Reference files directory')
    args = parser.parse_args()


    # Create the output directory if it doesn't exist
    os.makedirs(args.o, exist_ok=True)

    # Set up logging to file
    logging.basicConfig(filename=os.path.join(args.o, 'run_log.log'), level=logging.INFO)

    # Start logging
    logger.info('Script started.')


    # Define the reference files paths
    ref_dir = args.ref_dir
    go_obo = os.path.join(ref_dir, "go-basic.obo")
    gene2go = os.path.join(ref_dir, "gene2go.gz")
    if not os.path.exists(go_obo):
        os.system(f'wget http://purl.obolibrary.org/obo/go/go-basic.obo -O {go_obo}')
    if not os.path.exists(gene2go):
        os.system(f'wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz -O {gene2go}')

    organism = taxid_to_organism.get(args.taxid, 'Other')
    gene_info_gz = os.path.join(ref_dir, f"{organism}.gene_info.gz")

    # Get the gene_info file
    get_gene_info_file(args.taxid, gene_info_gz, ref_dir)

    # Load gene_info into a DataFrame
    gene_info_df = pd.read_csv(os.path.join(ref_dir, f'{taxid_to_organism.get(args.taxid, "Other")}.gene_info'), sep='\t', index_col='GeneID')

    # After reading the gene_info_df DataFrame
    print(gene_info_df.head())
    
    # Get the gene2refseq file
    gene2refseq_gz = os.path.join(ref_dir, "gene2accession.gz")
    get_gene2refseq_file(gene2refseq_gz, ref_dir)
    # Load gene2refseq into a DataFrame
    gene2refseq_df = pd.read_csv(os.path.join(ref_dir, 'gene2refseq'), sep='\t', index_col='protein_accession.version')

    # Convert a test ProteinID to a GeneID
    test_protein_id = 'XP_001920640.3'  
    print(protein_to_gene(test_protein_id, gene2refseq_df))


    # Load the GO DAG
    obodag = GODag(go_obo)
    logger.info('GO DAG loaded.')


    # Mapping NCBI Taxonomy IDs to Enrichr organism names
    taxonomy_id_to_organism = {
        9606: 'Human',
        10090: 'Mouse',
        10116: 'Rat',
        7227: 'Fly',
        6239: 'Worm',
        7955: 'Fish',
        4932: 'Yeast',
        9986: 'Rabbit',
        # Add more if needed
    }
    
    # Convert the taxid from an integer to a string for Enrichr
    enrichr_organism = taxonomy_id_to_organism.get(args.taxid, 'Other')

    
    # Read gene to GO term associations
    # Unzip the gzipped file
    with gzip.open(gene2go, 'rb') as f_in:
        with open(gene2go[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    geneid2gos_species = read_ncbi_gene2go(gene2go[:-3], taxids=[args.taxid], gene_info_df=gene_info_df)
    logger.info('Gene to GO term associations read.')

    # Extract hit IDs from BLAST results file if provided and no other background file is given
    if args.blast_results and not args.background:
        hit_ids = extract_blast_hit_ids(args.blast_results, gene_info_df)
        geneid2gos_species = {gene_id: go_terms for gene_id, go_terms in geneid2gos_species.items() if gene_id in hit_ids}

    # For each input file
    for input_file in args.i:
        # Read the CSV file
        df = pd.read_csv(input_file)
        df = df.drop_duplicates('ProteinID')  # remove duplicate entries
        df = df.dropna(subset=['ProteinID'])  # remove entries with missing ProteinIDs

        # Create a list of gene NCBI IDs
        gene_ncbi = df['ProteinID'].tolist()

        # Convert NCBI IDs to Entrez IDs
        gene_entrez = [protein_to_gene(id, gene2refseq_df) for id in gene_ncbi if protein_to_gene(id, gene2refseq_df) is not None]

        # Conduct GO enrichment analysis
        goeaobj = GOEnrichmentStudyNS(
            gene_entrez,
            geneid2gos_species,
            obodag,
            propagate_counts=False,
            alpha=args.alpha,
            methods=[args.method]
        )

        # Print for debugging
        print("gene_entrez: ", gene_entrez)
        print("geneid2gos_species: ", geneid2gos_species)
        print("Number of entries in geneid2gos_species: ", len(geneid2gos_species))
        #print("obodag: ", obodag)

        goea_results_all = goeaobj.run_study(gene_ncbi)

        # Print the total number of enrichment results
        print("Total number of enrichment results: ", len(goea_results_all))

        goea_results_sig = [r for r in goea_results_all if getattr(r, 'p_' + args.method) < args.alpha]

        # Print the number of significant results
        print("Number of significant results: ", len(goea_results_sig))

        # Save significant results as a CSV
        if len(goea_results_sig) > 0:
            goeaobj.wr_tsv(os.path.join(args.o, f'{os.path.basename(input_file).split(".")[0]}_goea_significant.tsv'), goea_results_sig)
            logger.info(f'Significant results saved for file {input_file}.')
        else:
            logger.warning(f'No significant results for file {input_file}.')


        # Conduct KEGG analysis if specified
        if args.kegg:
            logger.info('Running KEGG pathway enrichment analysis...')
            enr_res = run_enrichr(gene_list=gene_ncbi, gene_sets='KEGG_2016', organism=enrichr_organism, output_dir=args.o)
            enr_res.res2d.to_csv(os.path.join(args.o, 'KEGG_enrichment.csv'))

            # Plot KEGG pathway enrichment
            gseaplot(enr_res.ranking, term='KEGG Pathway Enrichment', ofname=os.path.join(args.o, 'KEGG_enrichment.png'))
            logger.info('KEGG pathway enrichment analysis completed.')

        # Plot specific KEGG pathway if specified
        if args.plot_pathway:
            plot_pathway(args.plot_pathway, args.o)
            logger.info(f'KEGG pathway {args.plot_pathway} plotted.')

    logger.info('Script finished.')

