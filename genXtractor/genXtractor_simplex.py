# genome_path -> (Assembled Fasta path) The genome of the organism you want to search things in
# output_path -> (path) Where the output fasta file will be placed
# Contig -> (string) name of the header of contig before any spaces
# SeqStart -> (positive int) Position that the sequence sequence starts on the contig (Example: 15405)
# SeqEnd -> (positive int) Position that the sequence sequence ends on the contig (Example: 32500)
# Before_offset -> (positive int) How many base pairs before the sequence would you like to extract
# After_offset -> (positive int) How many base pairs after the sequence would you like to extract
import os
import sys
from pathlib import Path

def SeqEx(Genome,Output,Query,Contig,SeqStart,SeqEnd,Before_offset,After_offset,Sense):
  #Turn paths into Path objects with windows or Posix path format depending on OS
  genome_path = Path(Genome) # Genome file path
  output_path = Path(Output) # Output dir path
  query_path = Path(Query)   # Query file path
  
  #Counters 
  line_cnt = 0             # lines parsed
  
  #Flags
  is_contig = False        # Marks if we reached the contig of interest
  
  # If output path doesn't exist. Make it
  output_path.mkdir(parents=True, exist_ok=True)

  #  
  #MAIN LOOP
  #Read each line
  #Find the right contig by counting ">" occurrences
  #Run algorith when reaching the right contig
  #If contig end or reached the end of the sequence, break loop end program
  
  #==============START======================
  output_file_name = output_path / (str(genome_path.stem) + "_extracted" + ".fasta.temp")
  test =  output_path / (str(genome_path.stem) + "_extracted" + ".fasta")
  if test.is_file(): # If the file already exists, no need to run the analysis again.    
    return
      
  with open(genome_path, "r") as genome:  # Genome to parse
    with open (output_file_name, "a") as output:   # Output file to write
      output.write(">" + str(genome_path.stem) + " Query=" + str(query_path.stem) +" Contig=" + str(Contig) + " SeqStart=" + str(SeqStart) + " SeqEnd=" + str(SeqEnd) + " BeforeOffset=" + str(Before_offset) + " AfterOffset=" + str(After_offset) + " SeqSize=" + str(SeqEnd-SeqStart) + " Sense=" + str(Sense) + "\n")

      char = 0             # Initiate the char counting variable
      start_flag = False   # Flag (bool) so we know when we reach the start of the sequence to be extracted
      SeqStart = SeqStart - Before_offset
      SeqEnd = SeqEnd + After_offset

      for line in genome:   # Parse the genome line by line
         
        if ">" in line:             # if we reached ">" header) need to check if its the contig of interest          
          if str(Contig) in line:   # if it is
            is_contig = True        # change is_contig to true
            continue                # go to next line
          else:                     # if it isn't, then we are either before or after it
            if not is_contig:       # if is_contig is false, we are before and thus can go to the next line
              continue              
            else:                   # if is_contig is true, we are after our contig of interest and can break the loop
              break
        elif not is_contig:         # if we haven't reached our contig and its not a header line, skip
          continue

        # Will only get here if we are after the header of our contig of interest, but before the next one
        linesize = len(line)-1
        char = char + linesize                     # Start counting the characters so we know when we reach the start
                
        # A bit less intuitive but more computationaly efficient. We'll start by checking if we need to write full lines to the output (middle of the sequence)
        if start_flag:
          if char < SeqEnd:
            output.write(line)                  # Write the whole line
          else:
            output.write(line[:linesize-(char-SeqEnd)])
            #output.write(str(char-(char-SeqEnd)))
            break
        else:
          if char >= SeqStart:
            output.write(line[linesize-(char-SeqStart)-1:])
            #output.write(str(char-(char-SeqStart)-1))
            start_flag = True
          else:
            continue


    # At this point there is nothing else to do. The sequence has been extracted. The files will close automaticaly and the program is done.
  #=================END======================
  # One last step remains, if the file is anti-sense, we need to output a extracted sequence in the correct direction
  if Sense == '-':
    reverse_complement_file(output_file_name,output_file_name.with_suffix(''))
    os.remove(output_file_name)
  else:
    os.rename(output_file_name,output_file_name.with_suffix(''))

def reverse_complement_file(inputfile, outputfile):
    inputpath = Path(inputfile)
    outputpath = Path(outputfile)

    with open(inputpath, 'r') as read_file:  # open the file as binary (can't work straight with UTF encoded files)
        positions = []
        while True:
            pos = read_file.tell()               # saves the position of the cursor at the start of the line
            line = read_file.readline()          # moves the cursor to the next line
            if not line:                         # if end of file, exit the loop
                break
            positions.append(pos)                # Save Start_of_line position to the the list of positions

        with open(outputpath, 'a', encoding='utf-8') as output_file:
            read_file.seek(positions[0])                              # positions the cursor at the start of the header
            output_file.write(read_file.readline())   # writes the header to the new file          
            
            for index in range(len(positions)-1, -1, -1):        # Iterate through indexes backwards
                if index == 0:                                   # if index is 0, we hit the header and should stop
                  break
                else:                                            # else we still have lines to go through
                    read_file.seek(positions[index])             # find the Start_line_pos and place the cursor there
                    line_to_write = reverse_complement(read_file.readline())                                      
                    output_file.write(line_to_write)  # write that line to a new file, since it goes by reverse index its going to write the file in reverse

def reverse_complement(seq):
    mapping = str.maketrans('ATCG', 'TAGC') # this creates the mapping table to translate the nucleotides into their complements
    newseq = seq.translate(mapping)
    return newseq[::-1]

# this is so that you can call the program from the command line and feed it the arguments there
if __name__== "__main__":
    #genome_path,output_path,Contig,SeqStart,SeqEnd,SeqSize,Before_offset,After_offset,Sense
    SeqEx(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), sys.argv[9])
