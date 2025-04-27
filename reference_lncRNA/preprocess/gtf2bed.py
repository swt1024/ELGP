import os
import csv
import sys

def gtf_to_bed(gtf_file, bed_file):
    """Convert a GTF file to BED format and save the output."""
    with open(gtf_file, 'r') as infile, open(bed_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = "chr" + fields[0]
            feature = fields[2]
            start = int(fields[3]) - 1  # Convert to 0-based indexing for BED format
            end = fields[4]
            strand = fields[6]
            attributes = fields[8]
            
            if feature == 'gene':
                gene_id = None
                name = None

                for attr in attributes.split(';'):
                    attr = attr.strip()
                    
                    # extract gene_id
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                        gene_id = gene_id.split('.')[0]

                    # extract gene_name
                    elif attr.startswith('gene_name'):
                        name = attr.split('"')[1]

                # Write to BED format: chrom, start, end, name, gene_id, strand
                writer.writerow([chrom, start, end, name, gene_id, strand])

def process_all_gtf_files(input_folder, output_folder):
    """Find all GTF files in the current folder and convert them to BED format."""
    
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # List all GTF files in the current directory
    gtf_files = [f for f in os.listdir(input_folder) if f.endswith('.gtf')]

    if not gtf_files:
        print("No GTF files found in the current directory.")
        return

    for gtf_file in gtf_files:
        input_gtf_path = os.path.join(input_folder, gtf_file)
        output_bed_path = os.path.join(output_folder, gtf_file.replace('.gtf', '.bed'))
        
        print(f"Processing: {gtf_file} → {output_bed_path}")
        gtf_to_bed(input_gtf_path, output_bed_path)

    print(f"Conversion complete. BED files saved to {output_folder}")

if __name__ == "__main__":

    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python gtf2bed.py <input_directory> <output_directory>")
        sys.exit(1)

    # Retrieve the input and output directory paths from command line arguments
    input_directory = sys.argv[1]
    output_directory = sys.argv[2]

    # Process all GTF files in the input directory
    process_all_gtf_files(input_directory, output_directory)

