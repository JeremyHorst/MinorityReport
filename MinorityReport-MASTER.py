#!/usr/bin/python
# MinorityReport-MASTER.py
# 20151221 JAHorst 
# 20160115 last update

from __future__ import print_function
import sys
import os
from sys import argv
from subprocess import call, Popen
import tempfile

tempdir = 'variant_temp_dir/'
clear_temp_files = False

################################################################################
# master script to split GFF, FASTA, and parent SAM & variant SAM into chromosomes, run in parallel.

try:
	fasta = argv[1]
	gff = argv[2]
	parent_sam = argv[3]
	mutant_sam = argv[4]
	MinorityReport = argv[5]
	print('\nBEWARE: this script will start 1 job per chromosome in your FASTA & GFF files in ~1 minute!', file=sys.stderr)
except:
	print('Usage:  MinorityReport-MASTER.py <genome.fasta> <genome.gff> <parent_mapped_reads.sam> <variant_mapped_reads.sam> <path_to_MinorityReport.py> [options for MinorityReport.py]', file=sys.stderr)
	print('e.g.:   ~/scripts/MinorityReport-MASTER.py Pfalciparum3D7-20160101.fasta Pfalciparum3D7-20160101.gff SensitiveStrain.sam ResistantStrain.sam ~/scripts/MinorityReport.py -vp 0.95 -wp 0.05 -vc 20 -wc 1 -wtc 0 -pe', file=sys.stderr)
	exit()

pass_arguments = ""
if len(argv)>5:
	for arg in argv[6:]:
		if arg[0]=='-':
			pass_arguments += arg +' '
			value_index = argv.index(arg)+1
			if len(argv) >= value_index:
				try:
					value = argv[value_index]
					if value[0] != '-':
						pass_arguments += value +' '
				except: single_flag_argument = True

################################################################################

try:  os.mkdir(tempdir)
except:  already_there = True

########
#chromosome_files = {}

# Methods
def get_sequence_from_fasta(fasta, name):										# bug fix 20151217: old solution runs into problems if one chromosome name is subset of another, e.g. chr1 and chr10
	#return ''.join('>'.join(open(fasta).read().split('>')[1:]).split(name)[1].split('\n>')[0].split('\n')[1:])
	real_entry = '\n' # gives blank result if real sequence not found
	for entry in open(fasta).read()[1:].split('\n>'):
		if entry.split('\n')[0].split('|')[0].strip() == name:
			real_entry=entry
			break
	return ''.join(real_entry.split('\n')[1:])

#def get_sequence_from_fasta(fasta, name):
#	return ''.join(open(fasta).read().split('>'+name+'\n')[1].split('\n>')[0].split('\n'))

def call_bash(command):
	call(command, shell=True, executable='/bin/bash')


################################################################################
# FASTA
print('processing FASTA file:',fasta, file=sys.stderr)
# get chromosomes from FASTA
chromosomes = []
for entry in open(fasta).read()[1:].split('\n>'):
	chromosome_name = entry.split()[0]
	sequence = ''.join(entry.split('\n')[1:]).strip()
	chromosomes += [[ len(sequence), chromosome_name ]]

# sort by length
chromosomes.sort(reverse=True)
for i in range(len(chromosomes)):
	chromosomes[i] = chromosomes[i][1]

# write chromosome-specific FASTA
for chromosome in chromosomes:
	chromosome_fasta = tempdir + fasta.split('/')[-1].split('.fa')[0]+'-'+chromosome+'.fasta'
	w=open(chromosome_fasta,'w')
	w.write('>'+chromosome+'\n')
	w.write(get_sequence_from_fasta(fasta, chromosome)+'\n')
	w.close()



################################################################################
# GFF
print('processing GFF file:',gff, file=sys.stderr)
chromosome_gff_data={}
for chromosome in chromosomes:
	chromosome_gff_data[chromosome] = [[],[]]  # [IDs, GFF lines]

for line in open(gff).readlines():
	chromosome = line.split()[0]
	if chromosome in chromosomes:
		ID=None
		# find lines of interest in GFF
		chromosome_gff_data[chromosome][1] += [line.strip()]
		ID = line.split('ID=')[1].split(';')[0]
		if ID:
			# store IDs to retrieve their sequences later, which may be at the end of the GFF file
			if ID not in chromosome_gff_data[chromosome][0]:  chromosome_gff_data[chromosome][0] += [ID]

# grabbing chromosome-specific actual GFF lines is fast & easy. 
# PROBLEM: getting the end-of-file sequences into the chromosome-specific GFFs is SLOW.

# make grep-able 1-line FASTA from GFF
grepable_GFF_FASTA = tempdir + gff.split('/')[-1]+'.fasta'
w=open(grepable_GFF_FASTA,'w')
for entry in open(gff).read().split('\n>')[1:]:
	w.write( '>'+ entry.split('\n')[0] +'\n'+ ''.join(entry.split('\n')[1:]) +'\n')
w.close()


# write chromosome-specific GFF files with all entries from end of GFF file
for chromosome in chromosomes:
	
	chromosome_gff = tempdir + gff.split('/')[-1].split('.gff')[0]+'-'+chromosome+'.gff'
#	print 'gff retrieval for chromosome:',chromosome, 'into:', chromosome_gff
	
	# write chromosome-specific GFF entries to file
	w=open(chromosome_gff,'w')
	for gff_entry in chromosome_gff_data[chromosome][1]:
		w.write(gff_entry +'\n')
	w.close()
	
	# prepare file of chromosome gene/protein sequence entries from GFF lines to search in sequence
	genelist_handle = tempdir+chromosome+'.genelist.temp'
	w=open(genelist_handle,'w')
	w.write('\n'.join(chromosome_gff_data[chromosome][0])+'\n')
	w.close()
	
	# append chromosome-specific sequence entries
	call_bash("""grep -A 1 -f %s %s >> %s""" % (genelist_handle, grepable_GFF_FASTA, chromosome_gff))
	


################################################################################
# SAM parent & mutant
print('processing SAM files:',parent_sam,'&',mutant_sam, file=sys.stderr)
if len(chromosomes) > 1:
	for chromosome in chromosomes[2:]:
	#	call_bash("""echo grep -P \"%s\\t\" %s \">\" %s &""" % (chromosome, parent_sam, tempdir + parent_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam') )
		call_bash("""grep -P \"%s\\t\" %s > %s &""" % (chromosome, parent_sam, tempdir + parent_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam') )
	#	call_bash("""echo grep -P \"%s\\t\" %s \">\" %s &""" % (chromosome, mutant_sam, tempdir + mutant_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam'), shell=True)
		call_bash("""grep -P \"%s\\t\" %s > %s &""" % (chromosome, mutant_sam, tempdir + mutant_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam') )

	# run the last pair sequentially, __not in background__, to allow extra time for the above to finish
	# the first chromosome tends not to be the shortest, and lets do the last one as well, just in case
	for chromosome in chromosomes[:2]:
		call_bash("""grep -P \"%s\\t\" %s > %s""" % (chromosome, parent_sam, tempdir + parent_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam') )
		call_bash("""grep -P \"%s\\t\" %s > %s""" % (chromosome, mutant_sam, tempdir + mutant_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam') )
else:
	chromosome = chromosomes[0]
	call_bash("""grep -P \"%s\\t\" %s > %s""" % (chromosome, parent_sam, tempdir + parent_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam') )
	call_bash("""grep -P \"%s\\t\" %s > %s""" % (chromosome, mutant_sam, tempdir + mutant_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam') )
	

################################################################################
#print tempdir,'files complete:', os.listdir(tempdir)
print('\nLaunching MinorityReport.py for each chromosome:\n', file=sys.stderr)
subprocesses = []
stdouts = []
for chromosome in chromosomes:
	command = MinorityReport +' '+ tempdir+fasta.split('/')[-1].split('.fa')[0]+'-'+chromosome+'.fasta ' + tempdir+gff.split('/')[-1].split('.gff')[0]+'-'+chromosome+'.gff '+ tempdir+parent_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam '+ tempdir+mutant_sam.split('/')[-1].split('.sam')[0]+'-'+chromosome+'.sam '
	command += pass_arguments +' '
	print(command, file=sys.stderr)
	
	# if the output (stdout) is more than 64KB, the subprocess will hang forever.
	# usually this is handled by communicate(). e.g. Popen(command...).communicate()
	# however, communicate() does not let us launch multiple jobs in parallel.
	# so, this is a work around.
	stdout_tempfile = tempfile.TemporaryFile()
	stderr_tempfile = tempfile.TemporaryFile()
	subprocesses += [Popen(command, stdout = stdout_tempfile.fileno(), stderr = stderr_tempfile.fileno(), close_fds = True, shell=True, executable='/bin/bash')]
	stdouts += [stdout_tempfile]

# wait for all processes to finish
print('waiting...', file=sys.stdout)
wait_until_done = [p.wait() for p in subprocesses]
print('Done.', file=sys.stdout)
# report output from all the separate runs
# handle the first job separately, remove the explanatory first output line from all others
first = True
for stdout_tempfile in stdouts:
	# go to the beginning of the file (NECESSARY)
	stdout_tempfile.seek(0)
	if first:
		printme = stdout_tempfile.read().strip() +'\n'
		first = False
	else:
		add_to_printme = '\n'.join(stdout_tempfile.read().strip().split('\n')[1:])
		if add_to_printme:  printme += add_to_printme+'\n'
if printme.strip():  print(printme, file=sys.stdout)

# clear up temp files
if clear_temp_files:
	# remove sub- SAM, GFF, & FASTA files
	call_bash('rm -rf ./'+tempdir)
	

