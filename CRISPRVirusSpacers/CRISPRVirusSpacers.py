#!/usr/bin/python3

import os.path
from os import path
import sys
import zipfile
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from subprocess import Popen, PIPE

cd_path = os.path.abspath("")

path = ""
outpath = ""
kingdom = "AB"
# dbpath = cd_path + "/DB/SpacersDB.fa"
dbpath = cd_path + "/NewDB/Spacers_ShortnameDB.fa"
user_dbpath = ""
fasta_path = ""
e_cutoff_s = 1e-7
max_target_seqs = 500
dbsize = 1e4
b_cutoff = 50
pr_cutoff = 0.7
min_nocshftsbg = 1
non_crisprbank = 0
mask_arrays = 0
ind = 0




    
####################FERTIG########################        
def for_loop(ARGV):
    for i in range(len(ARGV)-1):
        if ARGV[i] == '-in':
            global path
            path = ARGV[i+1]
            # print("in : ", path)
        if ARGV[i] == '-out':
            global outpath
            outpath = ARGV[i+1]
        if ARGV[i] == 'kingdom':
            global kingdom
            kingdom = ARGV[i+1]
            # print("kingdom : ", kingdom)
        if ARGV[i] == '-db':
            global dbpath
            dbpath= ARGV[i+1]
            user_dbpath = ARGV[i+1]
        if ARGV[i] == '-dbsize':
            global dbsize
            dbsize = int(ARGV[i+1])
        if ARGV[i] == '-e_cutoff_s':
            global e_cutoff_s
            e_cutoff_s = ARGV[i+1]
        if ARGV[i] == '-b_cutoff':
            global b_cutoff
            b_cutoff = ARGV[i+1]
        if ARGV[i] == '-bp_cutoff':
            global pr_cutoff
            pr_cutoff = ARGV[i+1]
        if ARGV[i] == '-min_nosch':
            global min_nocshftsbg
            min_nocshftsbg = ARGV[i+1]
        if ARGV[i] == '-fasta_path':
            global fasta_path
            fasta_path = ARGV[i+1]
        if ARGV[i] == '-non_crisprbank':
            global non_crisprbank
            non_crisprbank = 1
        if ARGV[i] == '-max_target_seqs':
            global max_target_seqs
            max_target_seqs = ARGV[i+1]
        if ARGV[i] == '-mask_arrays':
            global mask_arrays
            mask_arrays = 1
        if (ARGV[i] == '-h') or (ARGV[i] == '-help'):
            print("\n")
            os.system("cat help.txt")
            print("\n")
            sys.exit()
####################FERTIG########################

def read_map(file_name):
    # Read the Map
    Map = open(cd_path + "/NewDB/" + file_name, "r")
    contents = list(Map.readlines())
    Map_Sort_Full_Name = {}
    for i in range(0,len(contents),2):
        key = contents[i].strip("\n")
        value = contents[i+1].strip("\n")
        Map_Sort_Full_Name[key] = value
    Map.close
    return Map_Sort_Full_Name

def map_function(file_name, Map_Sort_Full_Name):
    # read_output.txt
    read_output = open(outpath + "/" + file_name, "r")
    contents = list(read_output.readlines())
    read_output.close

    # write to new file.
    WRITE = open(outpath + "/new_output.txt", "w")
    for j in range(0, len(contents)):
        head = (contents[j].rstrip("\n").split("\t"))
        process = Map_Sort_Full_Name[">" + head[1]]
        process = process.split("|",1)[1]
        process = process.replace(" ", "-")
        process = process.replace("|#","#")
        # print(process)
        l = []
        l.append(head[0]) 
        l.append(process) 
        l += head[2:]

        line = ""
        for i in l:
            line += i
            line += "\t"

        # l = "\t".join(str(x) for x in l)
        WRITE.write(line + "\n")
    WRITE.close()
    return "new_output.txt"

if __name__ == "__main__":
   
    if len(sys.argv) < 3:

        print("python CRISPRVirusSpacers.py -in <fasta_path> -out <dir_path> <options>")
        sys.exit()
    else:
        ARGV = sys.argv
        for_loop(ARGV)
        # unzip DB.zip
        dir_list = os.listdir()

        if "DB" not in dir_list and "DB.zip" in dir_list:
            os.system("unzip " + cd_path + "/DB.zip")

        # Delete outpath Folder.
        # TODO: not sure delete outpath or not ? 
        if outpath in dir_list:
            os.system("rm -rf " + outpath)
        os.system("mkdir " + outpath)

        # Check if Fasta path exist.
        # TODO: CHECK IF IT WORKS...
        if (len(fasta_path) > 0 and fasta_path in dir_list):
            os.system("cp " + fasta_path + " " + outpath + "/input_spacers.fa")
            cmd_fadb = "makeblastdb -in " + outpath + "/input_spacers.fa -out " + outpath + "/input_spacers.fa -dbtype nucl -parse_seqids"
            flag = os.system(cmd_fadb)
            dbpath = outpath + "/input_spacers.fa"

        path_2 = outpath + "/input_genomes.fa"
        path_3 = outpath + "/input_genomes_masked.fa"
        path_mask = outpath+ "/c_arrays.gff"
        
        # Bio SeqIO
        # TODO: ERROR getting alot of lines.
        f = open(path_2, "w")
        for record in SeqIO.parse(ARGV[2], "fasta"):
            # this works but a lot of lines
            SeqIO.write(record, path_2, "fasta")
        f.close()

        # TODO: CHECK if it works
        # Check if minced file exist.
        if "minced" not in dir_list: 
            if mask_arrays > 0:
                print("CRISPR-array predictor \'minced\' not found, skip the array-masking step.\n")
            mask_arrays = 0
        

        flag = os.system("makeblastdb -in " + path_2 + " -out " + path_2 + " -dbtype nucl -parse_seqids")

        qfile = path_2

        if ((mask_arrays > 0) and (os.path.isfile(path_3))):
            qfile = path_3

        flag = os.system("blastn -query " + qfile +" -db " + dbpath + " -word_size 7 -out " + outpath + "/Output.txt -gapopen 10 -max_target_seqs " + str(max_target_seqs) + " -gapextend 2 -penalty -1 -reward 1 -task blastn-short -lcase_masking -outfmt \'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq nident slen qlen sstrand\' -dbsize " + str(int(dbsize)) + " -evalue " + str(e_cutoff_s) + " -num_threads 24 >/dev/null 2>&1")
        ff= "blastn -query " + qfile +" -db " + dbpath + " -word_size 7 -out " + outpath + "/Output.txt -gapopen 10 -max_target_seqs " + str(max_target_seqs) + " -gapextend 2 -penalty -1 -reward 1 -task blastn-short -lcase_masking -outfmt \'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq nident slen qlen sstrand\' -dbsize " + str(int(dbsize)) + " -evalue " + str(e_cutoff_s) + " -num_threads 24 >/dev/null 2>&1"
        # print(ff)
        
        # TODO : output = map_function("Output.txt")
        Map_Sortname_FullnameDB = read_map("Map_Sortname_FullnameDB.fa")
        output = map_function("Output.txt", Map_Sortname_FullnameDB)
        # exit()
        # output = "Output.txt"
        x = "python3 " + cd_path + "/ProcessOutput.py " + outpath + "/" + output + " " + outpath + "/Full_HostVirst_Interaction_Results.csv "+ path_2  + " " + dbpath + " " + str(b_cutoff) + " " + str(pr_cutoff) + " " + str(min_nocshftsbg) + " " + str(non_crisprbank) + " " + outpath + "/Spacer_matched.fasta"
        # print(x)
        flag = os.system("python3 " + cd_path + "/ProcessOutput.py " + outpath + "/" + output + " "  + outpath + "/Full_HostVirst_Interaction_Results.csv "+ path_2  + " " + dbpath + " " + str(b_cutoff) + " " + str(pr_cutoff) + " " + str(min_nocshftsbg) + " " + str(non_crisprbank) + " " + outpath + "/Spacer_matched.fasta")

        # RAWDAT
        RAWDAT = open(outpath + "/Full_HostVirst_Interaction_Results.csv", "r")
        contents = list(RAWDAT.readlines())
        RAWDAT.close
        head = contents[0].rstrip("\n")

        # TEMP
        TEMP = open(outpath + "/Output.temp", "w")
        for j in range(1, len(contents)):
            line = contents[j].strip("\n")
            TEMP.write(line + "\n") 
        TEMP.close()

        # SDAT
        SDAT = open(outpath + "/Full_HostVirst_Interaction_Results.csv.temp2", "w")
        SDAT.write(head + "\n")
        SDAT.close()

        flag = os.system("cat " + outpath + "/Output.temp | sort -t \',\' -k17nr,17 -k10nr,10 >> " + outpath + "/Full_HostVirst_Interaction_Results.csv.temp2")

        unlink = os.system("rm -f " + outpath + "/Output.temp")
        unlink = os.system("rm -f " + outpath + "/Full_HostVirst_Interaction_Results.csv")
        
        flag = os.system("mv " + outpath + "/Full_HostVirst_Interaction_Results.csv.temp2 " + outpath + "/Full_HostVirst_Interaction_Results.csv")

    

        # SUMM
        # print(outpath + "/summ.csv")
        SUMM = open(outpath + "/summ.csv", "r")
        contents = list(SUMM.readlines())
        SUMM.close()

        # SUMMOUT     
        SUMMOUT = open(outpath + "/Summary_VirusResults.csv", "w")
        for k in range(len(contents)):
            line = contents[k].strip("\n")
            line = line.replace("\"", "")
            SUMMOUT.write(line + "\n")
        SUMMOUT.close()
        unlink = os.system("rm -f " + outpath + "/summ.csv")


        # FULL
        FULL = open(outpath + "/full.csv", "r")
        contents = list(FULL.readlines())
        FULL.close()

        # FULLOUT     
        FULLOUT = open(outpath + "/Fasta_Full_Results.csv", "w")
        for k in range(len(contents)):
            line = contents[k].strip("\n")
            line = line.replace("\"", "")
            FULLOUT.write(line + "\n")
        FULLOUT.close()
        unlink = os.system("rm -f " + outpath + "/full.csv")
        unlink = os.system("rm -f " + outpath + "/input_*.*")
