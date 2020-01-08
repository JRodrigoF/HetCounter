#!/usr/bin/python
# het_counter.py
# Feb 18

# ./het_counter.py /mmg/jflores/Het_P/data/1K_phase3_20130502/SNPs.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz <window size> 
# ./het_counter.py /mmg/jflores/Het_P/data/1K_phase3_20130502/SNPs.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz <50000>
# ./het_counter.py /mmg/jflores/Het_P/data/1K_phase3_20130502/SNPs.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 0
 

import sys
import re
import gzip
import os

out_dir_path = os.getcwd()

if sys.stdin.isatty():

	try:
		input_file = str(sys.argv[1])
		# file_handle = str(sys.argv[1])
	except IndexError:
		sys.exit("Input VCF is missing\n")

	try:
		W_length = sys.argv[2]
	except:
		sys.exit("Window size value is missing as the second argument\n")
		pass

	if input_file.endswith('.gz'):

		file_handle = gzip.open(input_file, 'r')
		gzflag = True

	elif input_file.endswith('.vcf'):
		file_handle = open(input_file, 'r')

else:

    file_handle = sys.stdin
    W_length = sys.argv[1]
    input_file = "stdin"
    gzflag = False


chrom_dict = {

	'chr1':249250621,
	'chr2':243199373,
	'chr3':198022430,
	'chr4':191154276,
	'chr5':180915260,
	'chr6':171115067,
	'chr7':159138663,
	'chr8':146364022,
	'chr9':141213431,
	'chr10':135534747,
	'chr11':135006516,
	'chr12':133851895,
	'chr13':115169878,
	'chr14':107349540,
	'chr15':102531392,
	'chr16':90354753,
	'chr17':81195210,
	'chr18':78077248,
	'chr19':59128983,
	'chr20':63025520,
	'chr21':48129895,
	'chr22':51304566,
	'chrX':155270560,
	'chrY':59034049
}

flag_FIRST_LINE = True
flag_HEAD_WINDOWS = True
flag_haploids = False
flag_missingness = False
gzflag = gzflag
total_VCF_positions = 0
biallelic_SNPs_counter = 0
SNPs_per_window_counter = 0
window_counter = 0

# Tests
# W_right_position = int(50900000)
# W_left_position = int(W_right_position) - int(W_length)
# W_step = int(W_length)
#


for rawline in file_handle:

	if gzflag:
		line = rawline.decode('utf8').rstrip('\n')
	else:
		line = rawline.rstrip('\n')

	header_INFO = re.search(r'^##', line)			# skips header
	header = re.search(r'^#CHROM', line)		

	if header_INFO:
		pass

	elif header:
		
		temp = re.split(r'\t', line)			
		new_header = "Window\tTotal_SNPs\tVCF_line\t" + "\t".join(temp[9:])

	else:											 

		field_line = re.split(r'\t', line)			
														


		#############################################
		if flag_FIRST_LINE:						  
			
			chr_in_file = field_line[0]

			try:
				chromosome_size = int(chrom_dict[chr_in_file])

			except:
				
				chr_in_file = "chr" + chr_in_file
				chromosome_size = int(chrom_dict[chr_in_file])

			if W_length == "all":
				W_length = int(3000000000)
				chr_in_file = "all"

			if W_length == "0":
				W_length = chromosome_size

			W_right_position = W_length
			W_left_position = W_length - W_length
			W_step = W_length

			if W_length == 3000000000:
				W_length = "all"

			file_name = out_dir_path + "/" + str(chr_in_file) + "_window_" + str(W_length) + "_het_stats.tab"
			with open("list_window_" + str(W_length) + "_het_stats_files.tab", 'a') as out_list:
			    out_list.write(file_name + "\t" + str(chr_in_file) + "\n")
			out_list.close()

			out = open(file_name, 'w')
			out.write(new_header + "\n")

			num_samples = len(field_line)-9

			list_counters = [float(0.0)] * num_samples
			counter_string = "\t".join(str(value) for value in list_counters)
			flag_FIRST_LINE = False

		#############################################

		index = 9

		# for now it will ignore variation other than SNPs,
		# an update should ideally consider whatever and rely rather on previous filtering of VCF input
		# if len(field_line[3]) > 1 or len(field_line[4]) > 1:
		alt_alleles = re.split(',', field_line[4])	
		MNP = False
		for allele in alt_alleles:
		    if len(allele)>1:
		        MNP = True

		if len(field_line[3]) > 1 or MNP: # 
			total_VCF_positions += 1
			continue
			
		conditon = True
		while conditon:

			if int(field_line[1]) > W_right_position and flag_HEAD_WINDOWS == True:

				out.write(str(W_right_position) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n")
				W_right_position += W_step
				W_left_position += W_step
				window_counter += 1

			elif int(field_line[1]) >= W_left_position and int(field_line[1]) < W_right_position:

				for sample in range(0, num_samples):			
			
				# 1K phase 3 -> 0|0 1|0 0|1 1|1 .. 2|0 and so on if multi-allelic are taken
				# Pagani et al -> 0|0 1|0 0|1 1|1 ..
				# neanderthal -> ./. 0/0 0/1 1/1 
				# chimp -> . 0/0 0/1 1/1

					genotype_field = re.split(r':', field_line[sample + index])
					genotype = re.split(r'[\|/]', genotype_field[0])

					# if re.search(r'[^01]', genotype[0]) or re.search(r'[^01]', genotype[1]):
						# pass # only passes one sample, continue you can be used here and re-start the loop
					
					# elif genotype[0] != genotype[1]:		
					if len(genotype) != 2:
						flag_haploids = True
						pass
						# print(field_line)
						# sys.exit()

					elif genotype[0] == '.' or genotype[1] == '.':
						flag_missingness = True
						pass
						# print(field_line)
						# sys.exit()

					elif genotype[0] != genotype[1]:		

						list_counters[sample] += 1.0
				
				SNPs_per_window_counter += 1
				total_VCF_positions += 1				
				flag_HEAD_WINDOWS = False
				conditon = False


			elif int(field_line[1]) >= W_right_position and int(field_line[1]) < chromosome_size: # saving equal to chrom size or "incomplete" final windows for last 

				# print W_right_position
				counter_string = "\t".join(str(value) for value in list_counters)
				out.write(str(W_right_position) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n")
				# print str(W_right_position) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n"
				list_counters = [float(0.0)] * num_samples
				SNPs_per_window_counter = 0
				W_right_position += W_step
				W_left_position += W_step
				window_counter += 1


		biallelic_SNPs_counter += 1

try:
	file_handle.close()
except:
	pass


counter_string = "\t".join(str(value) for value in list_counters)
out.write(str(W_right_position) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n")
# print str(chromosome_size) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n"
# print "Already OUT\n\n"
# print str(W_right_position) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n"

list_counters = [float(0.0)] * num_samples
counter_string = "\t".join(str(value) for value in list_counters)
SNPs_per_window_counter = 0
# total_VCF_positions = 0

if W_length != "all":

	while int(W_right_position) < chromosome_size:

		W_right_position += W_step
		window_counter += 1
		if int(W_right_position) > chromosome_size:
			# print "Aqui\t" + str(chromosome_size)
			out.write(str(chromosome_size) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n")
			# print str(chromosome_size) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n"
		else:
			out.write(str(W_right_position) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n")
			# print str(W_right_position) + "\t" + str(SNPs_per_window_counter) + "\t" + str(total_VCF_positions) + "\t" + counter_string + "\n"
			# print "Aca\t" + str(W_right_position)

out.close()

# print "\n\n\n"
print("file: " + str(input_file))
print("chromosome: " + str(chr_in_file))
print("chromosome size: " + str(chromosome_size))
print("window size: " + str(W_step))
print("windows total number: " + str(window_counter))
print("VCF positions total: " + str(total_VCF_positions))
print("SNPs total: " + str(biallelic_SNPs_counter))

if flag_haploids:
	print("\n\nWarning mesage: There are haploid genotypes in the dataset\n")

if flag_missingness:
	print("\n\nWarning mesage: There are missing genotypes in the dataset\n")



# Voila
