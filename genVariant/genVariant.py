# Query is the reference sequence
# Subject is the bacterial whole genome

import os
import sys
import subprocess
import csv
from pathlib import Path
from operator import itemgetter
from collections import defaultdict

def varFinder(Directory,Outdir,Query):

  genes_path = Path(Directory)
  output_path = Path(Outdir)
  
  try:                                 # If you don't offer a query path or an invalid one, then it won't add a query gene to the mix
    query_path = Path(Query)    
    if not query_path.is_file():
       query_path = ''       
  except TypeError:
    query_path = ''
 
  output_path.mkdir(parents=True, exist_ok=True)

  #--------------For the Update Bar---------
  total_files = len(os.listdir(genes_path))
  current_files = 1
  #-----------------------------------------

  genes = {}
  print ("Reading gene files...")
  for file_path in genes_path.rglob('*.fasta'):
    if file_path.is_file():
       name, gene = genGetter(file_path)
       genes[name] = gene
       UpdateBar(current_files,total_files)
       current_files +=1

  grouped = defaultdict(list)            # its a factory function designed to create a dictionary whose values are type list
  if query_path:                                          # If we have a queried gene to add then
    print ("Adding Query to report!")
    query_name, query_gene = genGetter(query_path)   # We start by initiating with the Queried Gene
    grouped[query_gene].append(query_name)           # Then we add the queried gene to the final dictionary
  
  
  for key, value in genes.items():       # dict.items() returns tuples of the key-value pairs of a dictionary
    grouped[value].append(key)           # it will create a dictionary with the gene as a key and a list of filenames as values


  print ("Creating Alleles.tsv...")
  with open(output_path / "Alleles.tsv", "a", newline='', encoding='utf-8') as allele_file:
      allele_writer = csv.writer(allele_file, delimiter='\t')  # Creater a writer for the report file
      allele_writer.writerow(["Allele","Representative","Quantity","All_genomes","Gene"])
      counter = 1
      total_alleles = len(grouped)
      for key in grouped:
        allele_writer.writerow([counter,grouped[key][0],len(grouped[key]),grouped[key],key])
        UpdateBar(counter,total_alleles)
        counter += 1

  
  print ("Creating Alleles.fasta...")
  with open(output_path / "Alleles.fasta", "a", newline='', encoding='utf-8') as fasta_file:
     total_alleles = len(grouped)
     counter = 1  
     for key in grouped:
        header = [">",grouped[key][0],"\n"]
        gene = [key,"\n"]
        fasta_file.write(''.join(header))
        fasta_file.write(''.join(gene))
        fasta_file.write('\n')
        UpdateBar(counter,total_alleles)
        counter += 1

  print ("All done!")       
         
'''
AUXILIARY FUNCTIONS
'''

def genGetter(file):
  #File format is first line header then content
  with open(file, "r", encoding='utf-8') as fastafile: # Open the file
    header = next(fastafile)                           # skips the header
    filename = header.split()[0]                       # Gets the genome name from the header
    filename = filename.replace(">","")                # Replaces any ">" from fasta files
    data = []
    for line in fastafile:                             # Reads file line by line 
      if line:
        data.append(line.strip())                      # Appends lines to data list
  gene = ''.join(data)                                 # Joins all the lines into one string
  return filename, gene                                # Outputs the filename and the gene



def UpdateBar(Current,Total):
  percent = Total//20
  try:
      print(f"{'Working... '}{'#' * (Current//percent)}{'-' * (20-(Current//percent))} file: {Current}'/'{Total}", end="\r", flush=True)
  except ZeroDivisionError:
      print(f'Working...', end="\r", flush=True)
  if Current == Total:
    print('\n')

# this is so that you can call the program from the command line and feed it the arguments there
if __name__== "__main__":
    #Query,Subject,Outdir
    varFinder(sys.argv[1], sys.argv[2], sys.argv[3])

'''
    for file in Path(sys.argv[2]).iterdir():
      varFinder(sys.argv[1], file, sys.argv[3])
      UpdateBar(current_files,total_files)
      current_files += 1
'''
#Open file, gather info, output report, send to extractor
#Extract sequence and output it to fasta

