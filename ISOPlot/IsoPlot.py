#this script organizes the exons coordinates for sashimi plots

import re,os,sys

#get the list of exons
#EM.out.organized.FDR5
ifile=open(sys.argv[1])
ifile.readline();
list=ifile.readlines();
ASM_list=[];
for i in list:
	elements=re.findall('[^\t\n ]+',i);
	ASM_list.append(elements[0]);

#get the coordinates of all events from annotations
#PC3E-1Aligned.sort.bam.IsoExon
ifile=open(sys.argv[2]);
list=ifile.readlines();
ASM2coord={};ASM2name={};
iscoord = 0;
for i in list:
	elements=re.findall('[^\t\n]+',i);
	if iscoord == 1:
		start = re.findall('[^,]+',elements[0]); start = start[0];
		end = re.findall('[^,]+',elements[-1]); end = end[-1];
		query = chr+':'+strand+':'+start+':'+end;
		ASM2coord[ASM] = query;
		ASM2name[ASM] = gene;
	if 'ASM#' in elements[0]:
		ASM = elements[0];
		strand = elements[3];
		chr = elements[4];
		gene = elements[5]+'_'+elements[6];
		iscoord = 1;
	else:
		iscoord = 0;

#write the coordinates and 	gene names
ofile=open(sys.argv[3],'w');
for i in ASM_list:
	ofile.write(ASM2coord[i]+'\n');
ofile.close();

#write the coordinates and 	gene names
ofile=open(sys.argv[4],'w');
for i in ASM_list:
	ofile.write('ASM'+i[4:]+'_'+ASM2name[i]+'\n');
ofile.close();
