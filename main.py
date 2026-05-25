#pipeline
#Run BLAST on the genomes assemblies for each reference sequence (from python use subprocess.run to run a cmd command)
#output in a folder named after the reference sequence
#Output should several files containing each the location of the gene in question on each genome

#verify if gene is present or not

#Run the extractor algorithm on each genome

# genblast requires the following arguments (Query,Subject,Outdir)
from genBlast.genBlast import genBlast
from genXtractor.genXtractor import SeqEx
from genVariant.genVariant import varFinder
from pathlib import Path
from operator import itemgetter
import csv
import sys
import os

def genFindXtract(Genome,Query,NewQuery,Output,Before_offset,After_offset):
  
  #create path objects
  genome_path = Path(Genome)       # directory path of the Genome assemblies
  query_path = Path(Query)         # File path of the reference gene you want to look for  
  output_path = Path(Output) / "Results"      # Directory path where the results of are to be put in
  output_path.mkdir(parents=True, exist_ok=True)

  report_path = output_path/ str(query_path.stem) / "Report.tsv"
  if NewQuery:
    new_query_path = Path(NewQuery)  # File path of the new reference gene in case the whole gene is present
    ext_report_path = output_path/ str(new_query_path.stem) / "Report.tsv"
  else:
    new_query_path = query_path
    ext_report_path = report_path

  #Start Counter  
  total_files = len(os.listdir(genome_path))
  current_files = 1

# RUN GENBLAST AND PRODUCE OUTPUT FILES
  for file in genome_path.iterdir():
    if genBlast(query_path,file,output_path) and NewQuery:
       genBlast(new_query_path,file,output_path)

    UpdateBar(current_files,total_files)
    current_files += 1
  
# RUN ANALYSIS
  with open(report_path, "r") as report_file: # Open the report file
    report_reader = csv.reader(report_file,delimiter='\t')                         # Initiate the reader variable
    header = next(report_reader)                                                   # Skip the first line (header)
    current_files = 1                          # Loop ended, reset the variable
    # Subject	Query	Presence	Single_contig	Alignment_length	Extended_ref
    for line in report_reader:
      UpdateBar(current_files,total_files)                                         # Creates an Progress Bar
      current_files += 1                                                           # Update Progress Bar Counter      
      line = dict(zip(header,line))
      if line["Presence"] == "False":                                # Case 1 - No matches        
        continue  # in case 1 (empty) we skip to the next one
     
      elif line["Single_contig"] == "False":                         # Case 2 - Multiple contigs
        #print ("extract multiple contigs")
        True

      elif int(line["Alignment_length"]) < 200:                           # Case 3 - Single contig with less than 200 bp
        # print ("extract single contig with less than 200 bp")        
        data = gather_data(output_path / str(query_path.stem) / 'BLAST' / (line['Subject'] + ".tsv"))  # Gives a list of each line in the blast tsv file sorted by match position on query (qstart)
        # Skip if file already exists (avoids resource waste)
        test_genome_path = Path(genome_path/ line["Subject"])
        test_output_path = Path(output_path/str(query_path.stem)/"Extracted" / "Short" / (str(test_genome_path.stem) + "_extracted" + ".fasta"))
        if test_output_path.is_file(): # If the file already exists, skip it
          continue
        # since its on the same contig, its either all sense or antisense and that can be detemrined by start and send
        # Element keywords = query scontig pident length mismatch gapopen qstart qend sstart send
        if data[0]['sstart'] > data[0]['send']:  # This check if the sequence is antisense
          start_point = data[-1]['send']         # Sequence is antisense, so start and end points are reversed
          end_point = data[0]['sstart']
          sense = '-'
        else:
          end_point = data[-1]['send']           # Sequence is sense, so start and end points are unchanged
          start_point = data[0]['sstart']
          sense = '+'

        SeqEx(genome_path/ line["Subject"],                                  #  Genome file path
                output_path/str(query_path.stem)/"Extracted" / "Short",      #  Output directory path
                new_query_path,                                              #  Query file path
                data[0]['scontig'],                                          #  Contig header name (in the file)
                int(start_point),                                            #  Sequence starting position
                int(end_point),                                              #  Sequence ending position                
                Before_offset,                                               #  How many base pairs before the sequence start
                After_offset,                                                #  How many base pairs after sequence ends
                sense)                                                       #  Sequence sense (for the purposes or contig header)
      
  with open(ext_report_path, "r") as report_file: # Open the report file
    report_reader = csv.reader(report_file,delimiter='\t')                         # Initiate the reader variable
    header = next(report_reader)                                                   # Skip the first line (header)
    ext_total_files = len(os.listdir(ext_report_path.parent / "BLAST"))
    current_files = 1                                                              # Loop ended, reset the variable
    # Subject	Query	Presence	Single_contig	Alignment_length	Extended_ref
    for line in report_reader:
      UpdateBar(current_files,ext_total_files)                                         # Creates an Progress Bar
      current_files += 1                                                           # Update Progress Bar Counter      
      line = dict(zip(header,line))      

      if line["Single_contig"] == "True":                           # Case 4 - Full extraction with new ref
        # print ("full extraction of the new ref")
        data = gather_data(output_path / str(new_query_path.stem) / 'BLAST' / (line['Subject'] + ".tsv"))  # Gives a list of each line in
        # Skip if file already exists (avoids resource waste)
        test_genome_path = Path(genome_path/ line["Subject"])
        test_output_path = Path(output_path/str(new_query_path.stem)/"Extracted" / "Full" / (str(test_genome_path.stem) + "_extracted" + ".fasta"))
        if test_output_path.is_file(): # If the file already exists, skip it
          continue
        # the blast tsv file sorted by match position on query (qstart)
        # since its on the same contig, its either all sense or antisense and that can be detemrined by start and send
        # Element keywords = query scontig pident length mismatch gapopen qstart qend sstart send
        if data[0]['sstart'] > data[0]['send']:  # This check if the sequence is antisense
          start_point = data[-1]['send']         # Sequence is antisense, so start and end points are reversed
          end_point = data[0]['sstart']
          sense = '-'
        else:
          end_point = data[-1]['send']           # Sequence is sense, so start and end points are unchanged
          start_point = data[0]['sstart']
          sense = '+'

        SeqEx(genome_path/ line["Subject"],                                  #  Genome file path
                output_path/str(new_query_path.stem)/"Extracted",            #  Output directory path
                new_query_path,                                              #  Query file path
                data[0]['scontig'],                                          #  Contig header name (in the file)
                int(start_point),                                            #  Sequence starting position
                int(end_point),                                              #  Sequence ending position                
                Before_offset,                                               #  How many base pairs before the sequence start
                After_offset,                                                #  How many base pairs after sequence ends
                sense)                                                       #  Sequence sense (for the purposes or contig header)
      else:
        continue  
  # Case 4 Full extraction with new ref 
  final_report(report_path)
  final_report(ext_report_path)
  varFinder(output_path/str(new_query_path.stem)/"Extracted",ext_report_path.parent,new_query_path)

'''
HERE STARTS THE AUXILIARY FUNCTIONS
'''

def UpdateBar(Current,Total):
  percent = Total//20
  try:
      print(f"{'Working... '}{'#' * (Current//percent)}{'-' * (20-(Current//percent))} file: {Current}'/'{Total}", end="\r", flush=True)
  except ZeroDivisionError:
      print(f'Working...', end="\r", flush=True)
  if Current == Total:
    print('\n')

def gather_data(file):
  with open(file, "r") as tsvfile: # Open the file
    tsvreader = csv.reader(tsvfile,delimiter='\t') # Creates a reader object to work with the tsv file
    data = []    
    for row in tsvreader:                          # For each row in the file
      entry = {"query": row[0], "scontig":row[1], "pident":row[2], "length":int(row[3]), "mismatch":int(row[4]), "gapopen":int(row[5]), "qstart":int(row[6]), "qend":int(row[7]), "sstart":int(row[8]), "send":int(row[9])}
      data.append(entry)                             # stores information as a list in a the data list variable
    # sort data by query match position (from lowest to highest query start), index = 6
    if not data:
      return data
    else:
      data_sorted = sorted(data, key=itemgetter("qstart"))      
      return data_sorted

def gather_element(data,name):   # creates a list with all the values of a column (an element)
  elements = {"query": 0, "scontig":1, "pident":2, "length":3, "mismatch":4, "gapopen":5, "qstart":6, "qend":7, "sstart":8, "send":9}
  if name in elements:
    if data:      
      sublist = [] 
      for item in data:              
       sublist.append(item[elements[name]])
      return sublist                             # returns a list with all the desired elements
    else:                                        # if the data is empty, exit with an error.
      print ("Empty data list")
      sys.exit()
  else:                                          # if the name isn't in the list of elements, exit with an error
      print ("Element name does not exist")
      sys.exit()

def reverse_complement_fast(seq):
    mapping = str.maketrans('ATCG', 'TAGC') # this creates the mapping table to translate the nucleotides into their complements
    return seq.translate(mapping)[::-1] # Will execute the translation on the sequence given and output it backwards (reverse)

def final_report(file):
  print('Writing final report...')
  empty_counter = 0
  multiple_contigs_counter = 0
  single_contig_counter = 0
  bp200_counter = 0
  extended_ref_counter = 0
  total_files_counter = 0
  query = ""  

  with open(file, "r") as report_file: # Open the report file
    report_reader = csv.reader(report_file,delimiter='\t')                         # Initiate the reader variable
    header = next(report_reader)                                                   # Skip the first line (header)
    for line in report_reader:      
      # Subject	Query	Presence	Single_contig	Alignment_length	Extended_ref
      line = dict(zip(header,line))
      query = line['Query']
      total_files_counter += 1
      if line['Presence'] == 'False':
        empty_counter += 1 
      elif line['Single_contig'] == "True":
        single_contig_counter += 1
        if line["Extended_ref"] == 'True':
          extended_ref_counter += 1
        else:
          bp200_counter += 1
      else:
        multiple_contigs_counter += 1        

  with open( file.parent / "Final_Report.tsv", "w" , newline='', encoding="utf-8") as final:
    final_writer = csv.writer(final,delimiter='\t')                         # Initiate the reader variable
    newheader = ["Query", "No_matches", "Multiple_Contigs", "Single_Contig", "Less_than_200bp", "Full_sequence", "Total" ]
    data = [query, empty_counter, multiple_contigs_counter, single_contig_counter, bp200_counter, extended_ref_counter, total_files_counter ]
    final_writer.writerow(newheader)  
    final_writer.writerow(data)       
    

# this is so that you can call the program from the command line and feed it the arguments there
if __name__== "__main__":
    #Genome, Query, Output, Before_offset, After_offset
    genFindXtract(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]))