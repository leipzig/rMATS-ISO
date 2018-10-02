#this script classifies the ASM modules
#Shihao Shen. 05/18/18

import re,os,sys

#read the IsoExon format
ifile = open(sys.argv[1]);
ilines = ifile.readlines();
asm_class = {};
this_pattern = {};
asm_list=[];class_list=[];
for i in ilines:
	elements=re.findall('[^\t\n]+',i);
	if 'ASM#' in elements[0]:
		asm_list+=[elements[0]];
		this_pattern=list(this_pattern.keys());
		this_pattern.sort();
		write_pattern = '';
		for k in this_pattern:
			write_pattern = write_pattern + k + ':';
		if not(write_pattern in asm_class):
			asm_class[write_pattern]=1;
		else:
			asm_class[write_pattern]=asm_class[write_pattern]+1;
		class_list+=[write_pattern];
		this_pattern={};
		skip = 1;
		exon_count=int(elements[1]);
	else:
		if skip == 1:
			skip = 0;
			exon_newindex=0;exon_add='A';
			exon_newindex_dic={};
			k = 0;left=0;
			while k < exon_count-1:
				first=re.findall('[^,]+',elements[k])
				second=re.findall('[^,]+',elements[k+1])
				left=max([left,int(first[1])]);				
				#overlapped exons
				if left > int(second[0]):
					if not(k in exon_newindex_dic):
						exon_newindex_dic[k]=str(exon_newindex)+exon_add;
					exon_add+='A'
					if not(k+1 in exon_newindex_dic):
						exon_newindex_dic[k+1]=str(exon_newindex)+exon_add;
					left=max([left,first[1],second[1]]);
				else:
					left=0;
					if not(k in exon_newindex_dic):
						exon_newindex_dic[k]=str(exon_newindex)
					exon_add='A';
					exon_newindex+=1;
				k+=1;
			if not(k in exon_newindex_dic):
				exon_newindex_dic[k]=str(exon_newindex)
		else:
			#print(exon_newindex_dic);
			for juc in range(len(elements[1:-1])):
				this_pattern[exon_newindex_dic[int(elements[juc+1])]+'-'+exon_newindex_dic[int(elements[juc+2])]]=0;
ifile.close();

#write the pattern summary
ofile=open(sys.argv[2],'w'); skip = 0;
for i in asm_class.keys():
	if skip == 0:
		skip = 1; continue;
	ofile.write(i+'\t'+str(asm_class[i])+'\n');
ofile.close();

#write the pattern for each asm
ofile=open(sys.argv[3],'w');
for i in range(len(asm_list)-1):
	ofile.write(asm_list[i]+'\t'+class_list[i+1]+'\n');
ofile.close();

