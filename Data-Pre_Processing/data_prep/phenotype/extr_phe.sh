#SBATCH --output=phenotype_generation.out
#SBATCH --error=phenotype_generation.err
#SBATCH --time=02:00:00      # Increased time limit to 2 hours
#SBATCH --mem=8G            # Increased memory allocation to 8 GB
#SBATCH --cpus-per-task=1

echo "Starting script..."

# Create and activate a new Conda environment
echo "Creating and activating new Conda environment..."
conda create -n new_env python=3.9 -y
source ~/miniconda3/bin/activate
conda activate new_env || { echo "Failed to activate Conda environment"; exit 1; }
echo "Conda environment activated."

# Add Bioconda and Conda-forge channels
echo "Adding Bioconda and Conda-forge channels..."
conda config --add channels bioconda
conda config --add channels conda-forge

# Install required packages
echo "Installing required packages..."
conda install -c bioconda -c conda-forge numpy pandas bcftools gsl -y || { echo "Failed to install required packages"; exit 1; }
echo "Required packages installed."

# Set library path (if needed)
export LD_LIBRARY_PATH=~/miniconda3/lib:$LD_LIBRARY_PATH
echo "Library path set."

# Create output directory if it doesn't exist
OUTPUT_DIR=""
echo "Creating output directory if it doesn't exist..."
mkdir -p $OUTPUT_DIR
echo "Output directory is ready."

# Set paths
RNA_TAR_FILE=""
RNA_OUTPUT_DIR=""
BED_FILE="$OUTPUT_DIR/phenotypes.bed"
NORMALIZED_COUNTS_FILE="$RNA_OUTPUT_DIR/normalizedCountMatrix_all.tsv"
GTF_FILE=""
ANNOTATION_FILE="$OUTPUT_DIR/gene_annotation.tsv"

# Extract RNA-seq data
echo "Extracting RNA-seq data..."
tar -xzvf $RNA_TAR_FILE -C $RNA_OUTPUT_DIR
echo "RNA-seq data extracted."

# Extract gene locations from GTF file
echo "Extracting gene locations from GTF file..."
python3 -c "
import pandas as pd
import gzip

gtf_file = '$GTF_FILE'
output_file = '$ANNOTATION_FILE'

annotations = []
with gzip.open(gtf_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] == 'gene':
            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            info = fields[8]
            gene_id = info.split('gene_id \"')[1].split('\"')[0]
            annotations.append([gene_id, chrom, start, end])

annotations_df = pd.DataFrame(annotations, columns=['gene_id', 'chr', 'start', 'end'])
annotations_df.to_csv(output_file, sep='\t', index=False)
"
echo "Gene locations extracted."

# Convert normalized RNA-seq data to BED format
echo "Converting normalized RNA-seq data to BED format..."
python3 -c "
import pandas as pd

normalized_counts = pd.read_csv('$NORMALIZED_COUNTS_FILE', sep='\t', index_col=0)
annotation = pd.read_csv('$ANNOTATION_FILE', sep='\t', index_col=0)

bed_data = []
header = ['#chr', 'start', 'end', 'phenotype_id'] + list(normalized_counts.columns)
bed_data.append(header)
processed_genes = 0
skipped_genes = 0

for gene in normalized_counts.index:
    if gene in annotation.index:
        chr = annotation.loc[gene, 'chr']
        start = annotation.loc[gene, 'start']
        end = annotation.loc[gene, 'end']
        row = [chr, start, end, gene] + list(normalized_counts.loc[gene])
        bed_data.append(row)
        processed_genes += 1
    else:
        skipped_genes += 1
        print(f'Skipped gene: {gene}')

bed_df = pd.DataFrame(bed_data[1:], columns=bed_data[0])
bed_df.to_csv('$BED_FILE', sep='\t', index=False)
print(f'Processed genes: {processed_genes}')
print(f'Skipped genes: {skipped_genes}')
print(f'Output file path: $BED_FILE')
"
echo "RNA-seq phenotype generation completed successfully."
