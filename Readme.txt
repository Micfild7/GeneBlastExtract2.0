This script is meant to search a group of genome assemblies for a single reference gene match, using BLAST and extract the gene, even if fragmented, if the matches are in the same contig. It will also give a report of allelic profiles based on complete nucleotide identity.

Dependencies:
python 3.11.2 or higher
BLASTn 2.7.1+

USAGE
main.py [assemblies_directory] [reference_gene_file] [reference_gene_extended_file] [Output_folder] [Start_offset] [End_offset]

All arguments have to be given, none can be skipped or left blank

- Assemblies_directory
[Path] - Directory path where the assembly fasta files are located

- Reference_gene_file 
[Path] - File path to the reference gene fasta file.

- Reference_gene_extended_file
[Path] - File path to the reference gene extended fasta file. If you do not have an extended gene file you want to search, give the same path as the reference_gene_file and the program will skip this step

[Output_folder] 
[Path] - Directory where the output will be placed

[Start_offset] 
[Positive Integer] - Quantity of nucleotides before the gene to be included in the extracted gene output. Setting it to 0 will include none.

[End_offset]
[Positive Integer] - Quantity of nucleotides after the gene to be included in the extracted gene output. Setting it to 0 will include none.
