# Query is the reference sequence
# Subject is the bacterial whole genome

import os
import sys
import subprocess
import csv
from pathlib import Path
from operator import itemgetter

def genBlast(Query,Subject,Outdir):
  #Turns the paths into path objects with windows or Posix path format depending on OS
  query_path = Path(Query) #File path
  subject_path = Path(Subject) #File path
  base_outdir_path = Path(Outdir) #Directory path

  #Verify if paths lead to files and/or directories 
  if not Path.is_dir(base_outdir_path):
    sys.exit("An error occurred: Outdir is not a directory path.")
  if not Path.is_file(query_path):
    sys.exit("An error occurred: Query is not a file path.")
  if not Path.is_file(subject_path):
    sys.exit("An error occurred: Subject is not a file path.")
  
  # Output directory file Structure
  #./Outdir
  #./Outdir/Query_name/Results/SubjectName_Result.tsv
  #./Outdir/Query_name/Report.tsv

  # establish the final_outdir_path
  outdir_path = base_outdir_path / query_path.stem / 'BLAST' # query_path.stem returns the filename without extensions
  # create final_outdir_path
  outdir_path.mkdir(parents=True, exist_ok=True)
  # parents = true will create all intermediary directories also
  # exist_ok = true will not raise an error if the directory already exists
  
  #Create output_filename_path
  output_filename = subject_path.name + ".tsv" # to put the .tsv on the filename
  output_filename_path = outdir_path / output_filename # Creates the final output path to go on the blast command
  if output_filename_path.is_file(): # If the file already exists, no need to run the analysis again.
    return

  #BLAST the query onto the Subject
  #There will be two queries. The first one is to identify if the gene is present or absent (base ref gene)
  #Second query will be done with the Gene + 500bp upstream region as ref. (Actual ref gene for extraction)
  #Skip Second query is first returns empty
  
  subprocess.run(["cmd", "/c", f'blastn -query {query_path} -out {output_filename_path} -subject {subject_path} -perc_identity 80 -outfmt 6 -word_size 100 -dust {"no"}'])  # i know this runs on windows, but i don't think it can run on linux in this current form
  #subprocess.run(["cmd", "/c", f'blastn -query {query_path} -out {outdir_path / (subject_path.name + ".xml")} -subject {subject_path} -perc_identity 80 -outfmt 5 -word_size 100 -dust {"no"}'])  # i know this runs on windows, but i don't think it can run on linux in this current form -outfmt 5 = XML

  #Generate Report
  #query_acc.ver	subject_acc.ver	%_identity	alignment_length	mismatches	gap_opens	q.start	q.end	s.start	s.end s.strand
  if not (outdir_path.parent / "Report.tsv").is_file():        # Checks if the file exists. If not, create file and write header
    with open(outdir_path.parent / "Report.tsv", "a", newline='', encoding='utf-8') as report_file:
      output_writer = csv.writer(report_file, delimiter='\t')  # Creater a writer for the report file
      output_writer.writerow(["Subject","Query","Presence","Single_contig","Alignment_length","Extended_ref"])           
      
  #Gather data and sort it based on query match
  filedata = gather_data(output_filename_path) # The data in the file is read and stored in a list of lists
                                               # The data then is sorted by qstart in ascending order
                                               # If the file is empty, filedata will be an empty list
  # Start decision tree
  blast_extended = False                       # Flag if another blast will be needed after this one (case 4)
  # Case 1: There are no matches (filadata is empty)
  if not filedata:
    # ["Subject","Query","Presence","Single_contig"," alignment_length","extended_ref"]
    # [string,string,Bool,Bool,int,Bool]
    data = [output_filename_path.stem,output_filename_path.parent.parent.stem,False,False,"",False]    
    #print ("The file is empty")

  # Case 2: There is a match, but its in multiple contigs
  elif not is_single_contig(filedata):
    # ["Subject","Query","Presence","Single_contig"," alignment_length","extended_ref"]
    # [string,string,Bool,Bool,int,Bool]
    data = [output_filename_path.stem,output_filename_path.parent.parent.stem,True,False,alignment_length(filedata),False]    
    #print ('is not a single contig')      

  # Case 3: There is a match in a single contig, but the size is inferior to 200 bp      
  elif alignment_length(filedata) < 200:
    data = [output_filename_path.stem,output_filename_path.parent.parent.stem,True,True,alignment_length(filedata),False]
    #print ('alingment is smaller than 200 bp')
    
  # Case 4: There is a match in a single contig bigger than 200 bp    
  else:
    data = [output_filename_path.stem,output_filename_path.parent.parent.stem,True,True,alignment_length(filedata),True]
    blast_extended = True
    #print ("Base Sequence to extract")
     
  with open(outdir_path.parent / "Report.tsv", "a", newline='', encoding='utf-8') as report_file:
    output_writer = csv.writer(report_file, delimiter='\t')  # Creater a writer for the report file
    output_writer.writerow(data)                             # Writes the data to the report file 

  return blast_extended    
  # All done

#Auxiliary functions
def gather_data(file):
  with open(file, "r") as tsvfile: # Open the file
    tsvreader = csv.reader(tsvfile,delimiter='\t') # Creates a reader object to work with the tsv file
    data = []
    for row in tsvreader:                          # For each row in the file      
      data.append(row)                             # stores information as a list in a the data list variable
    # sort data by query match position (from lowest to highest query start), index = 6
    if not data:
      return data
    else:
      data_sorted = sorted(data, key=itemgetter(6))      
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
      
def is_single_contig(data): #Returns True if the file is empty, False otherwise    
  if len(set(gather_element(data,"scontig"))) == 1:
    return True
  else:
    return False 
         
def alignment_length(data):  # Return the sum of the length of all matched fragments
  total = 0
  for item in gather_element(data,"length"):
    total += int(item)  
  return total

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
    total_files = len(os.listdir(sys.argv[2]))
    current_files = 1
    #Query,Subject,Outdir

    for file in Path(sys.argv[2]).iterdir():
      genBlast(sys.argv[1], file, sys.argv[3])
      UpdateBar(current_files,total_files)
      current_files += 1

#Open file, gather info, output report, send to extractor
#Extract sequence and output it to fasta

