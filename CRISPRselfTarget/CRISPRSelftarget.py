import shutil
import os
import re
import shutil
import sys
import argparse

try:
	from urllib.request import HTTPError
	from urllib.error import URLError
except ImportError:
	from urllib2 import HTTPError, URLError

try:
	from http.client import IncompleteRead
except ImportError:
	from httplib  import IncompleteRead


def is_complete(file_path):
	with open (file_path,"r") as file:
		first_line = file.readline()
		if "complete" in first_line:
			return True

def is_empty(file_path):
	if os.stat(file_path).st_size == 0:
		return True

def phaster(fasta_path, fasta_file, Results):
	# Generates the file with phage location
	# Also generate list of organism caused index error.
	# Index error problem arose after update of phaster site.
	identified_phages = Results + "/phages_identified_by_phaster.txt"
	error_list = Results + "/phaster_error_list.txt"
	phaster_output_file = Results + "/phaster_output_file.txt"
	phast_output_file = Results + "/phast_output_file.txt"
	phast_output_file2 = Results + "/phast_output_file2.txt"

	open(identified_phages, "w").close()
	open(error_list, "w").close()
	open(phaster_output_file, "w").close()
	open(phast_output_file, "w").close()
	open(phast_output_file2, "w").close()
	acc_no = fasta_file.split(".")[0]
	print("============================================================ ")
	print("	Phage Search using Phaster webservice for " + acc_no)
	print("============================================================ ")

	os.system('wget "http://phaster.ca/phaster_api?acc=' + acc_no + '"' + ' -O ' + phaster_output_file)
	if not is_complete(phaster_output_file):
		with open(fasta_path, "r") as file:
			first_line = file.readline()
			if first_line.startswith(">"):
				acc_no = first_line.split()[0].split("|")[0].split(".")[0].strip(">")
				print("acc_no: " + acc_no)
				os.system('wget "http://phaster.ca/phaster_api?acc=' + acc_no + '"' + ' -O ' + phaster_output_file)
	if not is_complete(phaster_output_file):
		os.system('wget --post-file="%s"  "http://phaster.ca/phaster_api"' %fasta_path + ' -O ' + phaster_output_file)

	if is_complete(phaster_output_file) and is_empty(phaster_output_file):
		print ("-----------There are results from Phaster-----------\n")
	try:
		data_file = open(phaster_output_file, "r")
		text = data_file.readline()
		loc = [m.start() for m in re.finditer('--', text)]
		# upto the index where -- exist and loc[-1]+2: ] will give us text
		# loc[-1] gives the last location
		text = text[loc[-1]+2:].strip().split()
		# text = text.strip().split()
		# \\n postions found in text list bcz next position gives us elements
		# which contains phages information.
		pos = [pos for pos, char in enumerate(text) if char == '\\n']

		with open(identified_phages, 'a') as f:
			for i in range(len(pos)):
				k = pos[i]
				wr = ("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (acc_no, text[k+1],
				text[k+2], text[k+3],
				text[k+7], text[k+5],
				text[k+17]))
				f.write(wr + '\n')
		data_file.close()
	except:
		"""
		os.system(
		'wget "http://phast.wishartlab.com/cgi-bin/phage_command_line.cgi?acc='
		+ acc_no + '"' + ' -O ' + phast_output_file)
		"""
		try:
			line = open(phast_output_file, "r").readlines()[3]
			line = line.replace(' ', '')[:-5]
			os.system('wget ' + line + ' -O ' + phast_output_file2)
			with open(phast_output_file2, "r") as fr:
				for content in fr:
					if "--" in content:
						loc = content.index("--")
			data_file = open(phast_output_file2, "r")
			text = data_file.readlines()
			l = len(text)
			for i in range(l):
				if i > loc:
					text = text[i].strip().split()
					with open(identified_phages, 'a') as f:
						k = 0
						wr = ("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
						acc_no, text[k + 0], text[k + 1], text[k + 2], text[k + 6], text[k + 4],
						text[k + 16]))
						f.write(wr + '\n')
			data_file.close()
		except (ValueError, IndexError):
			with open(error_list, 'a') as wr:
				wr.write(acc_no + '\n')
	return identified_phages

def count_spacer_lenght(genom_folder, Results, spacer_num):
	counter = 0
	spacer_lenght = open(Results + "/temp-Spacers%s-len" % spacer_num, "w")
	with open(genom_folder + "/Spacers_%s.fa" % spacer_num, "r") as file:
		for line in file:
			counter += 1
			if line.startswith(">"):
				spacer_name = line.strip(">").strip()
			else:
				lenght = len(line)
			if counter%2 == 0 :
				spacer_lenght.write(spacer_name + "\t" + str(lenght) + "\n")
	spacer_lenght.close()

def find_simi_spacers(genome_path, genom_folder, Results, spacer_num):
	count_spacer_lenght(genom_folder,Results, spacer_num)
	collected_spacers = []
	# calculate similarity between spacers and the whole genome
	cmd = "Bin/fasta36 " + genom_folder + "/Spacers_%s.fa " % str(spacer_num) + "%s -m 8 > "% genome_path + Results + "/Fasta-simi_Spacers_%s_genome.fastab" % str(spacer_num)
	print(cmd)
	res = os.system(cmd)
	cmd = "cat " + Results + "/Fasta-simi_Spacers_%s_genome.fastab | Bin/hashcol " % str(spacer_num)  + Results + "/temp-Spacers%s-len | Bin/tabcol2 'abs($13-$4)+$5 $13 $1 $7 $8 $9 $10' | sort -n > " % str(spacer_num) + Results + "/Spacers_%s_result2.txt" %spacer_num
	res = os.system(cmd)
	print(cmd)
	# ignore the regions which belong to any CRISPR-arrays
	with open(Results + "/Spacers_%s_result3.txt" % spacer_num, "w") as f_out:
		with open(Results + "/Spacers_%s_result2.txt" % spacer_num, "r") as f_in:
			for i, line in enumerate(f_in):
				sim_start = int(line.strip().split()[5])
				sim_end = int(line.strip().split()[6])
				f_out.write(line)

	# filter 100% similarity
	# 0	29	CRT_391564_391797_0	29	1	1266914	1266942
	cmd = "awk '{if ($1 == 0) print }' " + Results + "/Spacers_%s_result3.txt > "%spacer_num + Results + "/Spacers_%s_100.tab"% spacer_num
	res = os.system(cmd)

	# filter 90-99% similarity
	# 1	29	CRT_391875_394476_32	29	1	1258486	1258514
	cmd = "awk '{if ((($2-$1)/$2 >= 0.9) && (($2-$1)/$2 < 1.0)) print }' " + Results + "/Spacers_%s_result3.txt > "%spacer_num + Results + "/Spacers_%s_90.tab"% spacer_num
	res = os.system(cmd)

	# filter 80-89% similarity
	# 1	29	CRT_391875_394476_32	29	1	1258486	1258514
	cmd = "awk '{if ((($2-$1)/$2 >= 0.8) && (($2-$1)/$2 < 0.9)) print }' " + Results + "/Spacers_%s_result3.txt > " %spacer_num + Results + "/Spacers_%s_80.tab"% spacer_num
	res = os.system(cmd)

def main():
	cmdline_parser = argparse.ArgumentParser('Self-target')
	cmdline_parser.add_argument('-i', '--input_files',
	nargs='+',
	default= ['FP565176.fasta'],
	help='Filenames of the input data.',
	type=str)
	cmdline_parser.add_argument('-a', '--annotation',
	nargs='+',
	default= ['FP565176/annotations_details.tab'],
	help='Filenames of the annotations',type=str)
	cmdline_parser.add_argument('-tdir', '--fold',
	nargs='+',
	default= ['FP565176'],
	help='folder of files',
	type=str)
	cmdline_parser.add_argument('-out', '--result_folder',
	nargs='+',
	default= ['Results'],
	help='Results',
	type=str)
	cmdline_parser.add_argument('-top', '--topspacers',
	nargs='+',
	default= ['0'],
	help='top k spacers',
	type=int)
	args, unknowns = cmdline_parser.parse_known_args()

	fasta_path = args.input_files[0]
	fasta_file = fasta_path.split("/")[-1]
	genome_name = fasta_file.split(".")[0]
	Results = args.result_folder[0]
	annotation_file = args.annotation[0]
	genom_folder = args.fold[0]
	top_num = int(args.topspacers[0])
	virsorter_cpu = 4
	virsorter_data_dir = "virsorter_data"
	phage_finder_dir = "Phage_finder_cluster"
	phage_finder_sh = "Phage_finder_cluster/phage_finder_customised.sh"
	pwd_path = ""

	if not os.path.isfile(fasta_path):
		print("The fasta file: %s isn't exist" % fasta_file)
		return
	if not os.path.isdir(genom_folder):
		print("The Genome folder: %s isn't exist" %genom_folder)
		return
	if not os.path.isdir(genom_folder):
		print("The annotation file: %s isn't exist" %annotation_file)
		return

	if os.path.isdir(Results):
		shutil.rmtree(Results)
	os.mkdir(Results)
	# phaster
	phaster_identified_phages = phaster(fasta_path,fasta_file, Results)
	print("phaster_identified_phages", phaster_identified_phages)
	#find simi spacers
	Spacer_files = []
	path = genom_folder +"/"
	for i in os.listdir(path):
		if os.path.isfile(os.path.join(path,i)) and 'Spacers' in i:
			Spacer_files.append(i)
	Spacer_files.sort()
	spacer_number = len(Spacer_files)
	for j in range(0,spacer_number):
		spacer_num = int(Spacer_files[j].split("_")[1].split(".")[0])
		print("find_simi....spacer_number:", spacer_num)
		find_simi_spacers(fasta_path, genom_folder, Results, spacer_num)
	#the end result
	End_Result_path = Results + "/End_result.tab"
	End_Result = open (End_Result_path,"w")
	End_Result.write("")
	counter_crispr = 0
	spacers = 0
	casgene_number = 0
	prophadge_number = 0
	genomic_number = 0
	#iterate on spaceres Files
	for i in range(0,spacer_number):
		spacer_num = int(Spacer_files[i].split("_")[1].split(".")[0])
		if not is_empty(Results + "/Spacers_%s_result2.txt" % spacer_num):
		    with open (Results + "/Spacers_%s_result2.txt" % spacer_num, "r") as file:
		    	#iterate on spacer in spaceres file
		    	for spacer_line in file:
		    		spacers += 1
		    		founded = False
		    		a_CRISPR = False
		    		spacer_list = spacer_line.split()
		    		spacer_id = spacer_list[0]
		    		spacer_first = spacer_list[1]
		    		spacer_second = spacer_list[2]
		    		spacer_third = spacer_list[3]
		    		spacer_forth = spacer_list[4]
		    		spacer_start = int(spacer_list[5].strip())
		    		spacer_end = int(spacer_list[6].strip())
		    		#search in annotation_file
		    		if not founded:
		    			with open(annotation_file, "r") as anno_file:
		    				for anno_line in anno_file:
		    					if not anno_line.startswith(("#","%","!")):
		    						anno_start = int(anno_line.split()[2])
		    						anno_end = int(anno_line.split()[3])
		    						anno_type = anno_line.split()[4]
		    						anno_label =  anno_line.split()[5]
		    						if ((spacer_start <= anno_end) and (spacer_start >= anno_start)) or ((spacer_end <= anno_end) and (spacer_end >= anno_start)):
		    							if anno_type != "CRISPR":
		    								result_line = spacer_id +"\t" + spacer_first +"\t" + spacer_second +"\t" +"Array_%s " %spacer_num + "\t" + spacer_third + "\t" + spacer_forth + "\t" + str(spacer_start) + "\t" + str(spacer_end) + "\tCasgeneRegion\t" + anno_label + "\n"
		    								founded = True
		    								casgene_number += 1
		    								break
		    							else:
		    								founded = True
		    								a_CRISPR = True
		    								counter_crispr +=1
		    								break
		    		#search in phages_identified_by_phaster.txt
		    		if not founded:
		    			with open(Results + "/phages_identified_by_phaster.txt", "r") as phaster_file:
		    				for phaster_line in phaster_file:
		    					phaster_start = int(phaster_line.split()[5].split("-")[0])
		    					phaster_end = int(phaster_line.split()[5].split("-")[1])
		    					if phaster_start > phaster_end:
		    						holder = phaster_start
		    						phaster_start = phaster_end
		    						phaster_end = holder
		    					if ((spacer_start <= phaster_end) and (spacer_start >= phaster_start)) or ((spacer_end <= phaster_end) and (spacer_end >= phaster_start)):
		    						result_line = spacer_id +"\t" + spacer_first + "\t" + spacer_second +"\t" +"Array%s " %spacer_num + "\t" + spacer_third + "\t" + spacer_forth + "\t" + str(spacer_start) + "\t" + str(spacer_end) + "\tProphageRegion\tProphageRegion\n"
		    						founded = True
		    						prophadge_number += 1
		    						break
		    		if not founded:
		    			result_line = spacer_id +"\t" + spacer_first +"\t" + spacer_second +"\t" +"Array_%s " %spacer_num + "\t" + spacer_third + "\t" + spacer_forth +"\t" + str(spacer_start) +"\t" + str(spacer_end) +"\tGenomicRegion\tGenomicRegion\n"
		    			genomic_number += 1
		    		if not a_CRISPR:
		    			if top_num == 0:
		    				End_Result.write(result_line)
		    			elif top_num >= int(result_line.split()[0]):
		    				End_Result.write(result_line)
	print("\n\n-------------------------------")
	print(genome_name)
	print("The number of all spacers: ", str(spacers))
	print("-------------------------------")
	print("The number of spacers from type CRISPR: ", str(counter_crispr))
	print("The number of spacers in CasgeneRegion: " , str(casgene_number))
	print("The number of spacers in ProphageRegion: ", str(prophadge_number))
	print("The number of spacers in GenomicRegion: ", str(genomic_number))
	print("The total number of spacers written in the end result: ",str(casgene_number + prophadge_number + genomic_number))
	End_Result.close()
	os.system("sort -k1n " + End_Result_path + " | uniq > "+ Results + "/Summary_Selftarget_Results.tab")
	os.system("rm " + End_Result_path)


if __name__ == "__main__":
	main()
