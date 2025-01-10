This script processes fastq files directly from the output of the MinKNOW software without requiring preprocessing on other software. Modifications to the code may be necessary to adapt to specific sequence requirements.
1.	Prepare input and output paths:
Specify the paths for the input fastq files and the output directory through command-line arguments:
--fastq (-f): Path to the folder containing fastq files (either directly or within subfolders).
--output (-o): Path to the output folder.
--index: Index IDs (e.g., R1,R2,R3) to match during the sequence analysis.

2.	Main options:
--match: Logical variable indicating whether to perform sequence matching (default: TRUE).
--stat: Logical variable indicating whether to generate statistical outputs (default: TRUE).
--min_mis: Minimum number of mismatched bases allowed for pattern matching (default: 0).
--max_mis: Maximum number of mismatched bases allowed for pattern matching (default: 0).
--cores (-t): Number of CPU cores to use for parallel processing (default: 1).

3.	Major steps of the analysis:
i. File detection and reading
The script detects all fastq files within the input directory (including subdirectories) and reads them into memory.
The resulting data table contains:
Barcode ID (barcode)
Sequence length (width)
Raw sequence data (seq)

ii. Sequence matching
If --match = TRUE, the script matches the raw sequences to the specified index sequences provided via the --index argument (e.g., R1,R2,R3).
Matches are performed based on user-defined mismatch criteria (--min_mis, --max_mis).
Outputs include:
Match counts for each sequence pattern.
Metadata for matched sequences.

iii. Statistical analysis
If --stat = TRUE, the script generates the following:
Total number of matched sequences for each barcode.
Percentage of matches for each index sequence relative to a reference (S1).

4.	Output files:
Matched sequence information grouped by barcode and saved as .csv files.
Summary statistics:
summary_counts.csv: Total match counts for each barcode.
summary_percent.csv: Percentage match statistics for each barcode.

5.	How to run the script:
Example command to execute the script:
bash
Rscript nanopore.r --fastq /path/to/fastq_folder --output /path/to/output_folder --index R1,R2,R3 --cores 4

6.	Dependencies:
This script requires the following R packages:
optparse
ShortRead
data.table
Biostrings
parallel

7.	Key outputs:
Match results: Detailed information on sequence matching for each barcode.
Statistical summaries: CSV files for match counts and percentages.
