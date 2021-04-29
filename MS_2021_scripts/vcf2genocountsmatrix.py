#!/usr/bin/env python
# TL-051120
#last maj 091120

import sys

# usage: 
# python vcf2genocountsmatrix.py Bromeliad1776_mappedAnanas_allsamples_primitives_withmono_filter4_DPindelgap_maxmissing20_10klines.vcf pineapple.20150427.annotat.gff3.LG01 Bromeliad1776.genes.in.set.pineapple.20150427.genes.gff3 

vcf = open(sys.argv[1])
gff = open(sys.argv[2])
targetedgenes = open(sys.argv[3])
flanking= int(500) # define the size of the flanking regions
outputfile=sys.argv[1]+".matrixcounts.summaryvcf"
output_structure = open(outputfile,"w")


### 1st: create different dictionnaries corresponding to different features of the total gff (to speed up the script, but requires a lot of memory, subset the gff per chromosomes in cases of memory issues)
gff_CDS_dict = {} 
gff_exon_dict = {} 
gff_gene_dict = {} 
gff_flanking_dict = {} 
print("Step1:reading the whole gff")
# read the vcf
for line1 in gff.readlines():
	if line1[0] != "#": 
		line1 = line1.replace('\n','')
		line1 = line1.replace(';','')
		splitted_line1 = line1.split('\t')
		scaffID = splitted_line1[0] 
		feature = splitted_line1[2] 
		posstart = int(splitted_line1[3])
		posend = int(splitted_line1[4])
		geneID=  splitted_line1[8] 
		genepos=1 
		if feature == "CDS":
			for positiongenegenome in range(int(posstart), int(posend), 1):
				genepos+=1
				key = (scaffID + '-' + str(positiongenegenome)) 
				gff_CDS_dict[key] = str(geneID)+"\t"+str(genepos)
				#print(key)
		if feature == "exon":
			for positiongenegenome in range(int(posstart), int(posend), 1):
				genepos+=1
				key = (scaffID + '-' + str(positiongenegenome)) 
				gff_exon_dict[key] = str(geneID)+"\t"+str(genepos)
				#print(key)
		if feature == "gene":
			# position gene (if not already considered in the previous CDS dict) will be considered as intronic later.
			for positiongenegenome in range(int(posstart), int(posend), 1):
				genepos+=1
				key = (scaffID + '-' + str(positiongenegenome)) 
				gff_gene_dict[key] = str(geneID)+"\t"+str(genepos)
				#print(key)
			# flanking regions of genes
			for positiongenegenome in range(int(posstart-flanking),int(posstart-1)):
				genepos+=1
				key = (scaffID + '-' + str(positiongenegenome)) 
				gff_flanking_dict[key] = str(geneID)+"\t"+str(genepos)
				#print(key)
			for positiongenegenome in range(int(posend+1),int(posend+flanking)):
				genepos+=1
				key = (scaffID + '-' + str(positiongenegenome)) 
				gff_flanking_dict[key] = str(geneID)+"\t"+str(genepos)
				#print(key)

### 2nd: create dictionnaries for the targeted genes (to speed up the script later)
gff_ontarget_dict = {} 
gff_flankingtarget_dict = {} 
# read the vcf
print("Step2:reading the gff with the targeted genes")
for line2 in targetedgenes.readlines():
	if line2[0] != "#": 
		line2 = line2.replace('\n','')
		line2 = line2.replace(';','')
		splitted_line2 = line2.split('\t')
		scaffID = splitted_line2[0] 
		feature = splitted_line2[2] 
		posstart = int(splitted_line2[3])
		posend = int(splitted_line2[4])
		geneID=  splitted_line2[8] 
		genepos=1 
		for positiongenegenome in range(int(posstart), int(posend), 1):
			genepos+=1
			key = (scaffID + '-' + str(positiongenegenome)) 
			gff_ontarget_dict[key] = str(geneID)+"\t"+str(genepos)
			#print(key)
			# flanking regions of genes
		for positiongenegenome in range(int(posstart-flanking),int(posstart-1)):
			genepos+=1
			key = (scaffID + '-' + str(positiongenegenome)) 
			gff_flankingtarget_dict[key] = str(geneID)+"\t"+str(genepos)
			#print(key)
		for positiongenegenome in range(int(posend+1),int(posend+flanking)):
			genepos+=1
			key = (scaffID + '-' + str(positiongenegenome)) 
			gff_flankingtarget_dict[key] = str(geneID)+"\t"+str(genepos)
			#print(key)


### 3rd: then read the vcf line by line and report summary statistics for each line of the vcf
nameInd = []
IDs="locus\tallref\tallalt"
print("Step 3: reading the vcf")
print("Generating the '"+ outputfile + "' file")
for line in vcf:
	#print(line)
	line = line.rstrip()
	
	# count number of individuals and store names in array
	if line[0:6] == "#CHROM":
		allalleles=0
		arrline = line.split()
		for i in arrline[9:]:
			nameInd.append(i)
			IDs=IDs+"\t"+i+"_all1\t"+i+"_all2"
			allalleles+=2
		#print(IDs)
		header="Scaffold"+"\t"+"Position"+"\t"+"TotalNbAlleles"+"\t"+"GeneticVariationFreeBayes"+"\t"+"Ref"+"\t"+"Alt"+"\t"+"MonoOrPolymorphicBasedOnCounts"+"\t"+"GenicOrIntergenic"+"\t"
		header=header+"CDSOrUTROrIntronsOrIntergenic"+"\t"+"OnOrOffTargets"+"\t"+"CountsMissing"+"\t"+"CountsRef"+"\t"+"CountsAlt1"+"\t"+"CountsAlt2"+"\t"+"CountsAlt3"+"\n"
		output_structure.write(header)
	if line[0] != "#":
		NAalleles=0
		allele1=0
		allele2=0
		allele3=0
		allele4=0
		multiallelicstate=0
		allGT=""
		status0=""
		status1=""
		status2=""
		status3=""
		status4=""
		arrline = line.split()
		chromosome = arrline[0]
		pos = int(arrline[1])
		ref = str(arrline[3])
		alt = str(arrline[4])
		# nb of alleles & status0 (=nature of the genetic variation)
		if alt == ".":
			nballeles=1
			status0="monomorphic"
		else:
			nballeles = len(alt.split(","))+1
			currentallelelength=0
			maxallelelength=0
			for allalleles in range(0, nballeles-1):
				currentlength=len(alt.split(",")[allalleles])
				if currentlength > maxallelelength:
					maxallelelength=currentlength
			#status
			if maxallelelength == 1:
				status0="SNP"
			elif maxallelelength > 1:
				status0="indel"	
			else:
				status0="Unexpected situation"	
		# loop over the individuals calls and keep only the GT info
		GT = []
                GTfieldExists = 0
		GTPos = 0
		info = arrline[8].split(":")
		for i in info:
			if i == "GT":
				GTfieldExists = 1
                                break
			GTPos= GTPos+1
		for i in range(0, len(nameInd)):
			ID = nameInd[i]
			ind = arrline[i+9]
			arrInd = ind.split(":")
			PASStagExists = 0 # test number of info
			if  GTfieldExists == 1 and len(arrInd) >= GTPos: # QC for the VCF info (looks like a vcf)
				myGT = arrInd[GTPos]
				#print(myGT)
				if len(arrInd[GTPos].split("/"))!=2 and len(arrInd[GTPos].split("|"))!=2:
					if arrInd[GTPos].split("/")[0] == "." or arrInd[GTPos].split("|")[0] == "." :
						NAalleles+=2
					else:
						print("unexpected situation")
				else:
					# First allele
					myfirstall=arrInd[GTPos].split("/")
					if len(myfirstall)==1:
						myfirstall=arrInd[GTPos].split("|")
						if len(myfirstall)==1:
							print("something weird here"+str(arrInd[GTPos]))
						else:
							myfirstall=arrInd[GTPos].split("|")[0]
					else:
						myfirstall=arrInd[GTPos].split("/")[0]
					if myfirstall == ".":
						NAalleles+=1
					elif myfirstall == "0":
						allele1+=1
					elif myfirstall == "1":
						allele2+=1
					elif myfirstall == "2":
						allele3+=1
					elif myfirstall == "3":
						allele4+=1
					else:
						multiallelicstate=1
					# Second allele
					if len(arrInd[GTPos].split("/"))==1:
						if len(arrInd[GTPos].split("|"))==1:
							print("something weird here"+str(arrInd[GTPos]))
						else:
							mysecondall=arrInd[GTPos].split("|")[1]
					else:
						mysecondall=arrInd[GTPos].split("/")[1]
					if mysecondall == ".":
						NAalleles+=1
					elif mysecondall == "0":
						allele1+=1
					elif mysecondall == "1":
						allele2+=1
					elif mysecondall == "2":
						allele3+=1
					elif mysecondall == "3":
						allele4+=1
					else:
						multiallelicstate="weird"
					#allGT  = str(allGT) + "\t"+ str(myfirstall) + "\t"+ str(mysecondall)
			else:
				continue
		# CHECK STATUS
		#status1 = monomorphic or polymorphic
		totalalleles=NAalleles+allele1+allele2+allele3+allele4
		# if calls are for only one category (all missing, all 0, all 1, etc...) or if it is the case after taken into account the missing data, consider this SNP as monomorphic
		if NAalleles == totalalleles or allele1 == totalalleles or allele2 == totalalleles or allele3 == totalalleles or allele4 == totalalleles or allele1 == ( totalalleles - NAalleles ) or allele2 == ( totalalleles - NAalleles ) or allele3 == ( totalalleles - NAalleles ) or allele4 == ( totalalleles - NAalleles ) :
			status1 = "monomorphic"
		else:
			status1 = "polymorphic"

		#status2=genic or intergenic (status 3: if exonic: CDS or UTR?)
		key_vcf = (str(chromosome) + '-' + str(pos))
		if (gff_CDS_dict.has_key(key_vcf)): 
			status2="exonic"
			status3="CDS"
		elif (gff_exon_dict.has_key(key_vcf)): 
			status2="exonic"
			status3="UTR"
		elif (gff_gene_dict.has_key(key_vcf)): 
			status2="intronic"
			status3="intron"
		elif (gff_flanking_dict.has_key(key_vcf)): 
			status2="flanking"
			status3="flanking"
		else:
			status2="intergenic"
			status3="intergenic"
		#status4=on-target or off-target?
		if (gff_ontarget_dict.has_key(key_vcf)): 
			status4="on-target"
		elif (gff_flankingtarget_dict.has_key(key_vcf)): 
			status4="flanking-targets"
		else:
			status4="off-target"
		#print(status2)
		result_analysis=str(chromosome)+"\t"+str(pos)+"\t"+str(nballeles)+"\t"+str(status0)+"\t"+str(ref)+"\t"+str(alt)+"\t"+str(status1)+"\t"+str(status2)+"\t"+str(status3)+"\t"+str(status4)+"\t"+str(NAalleles)+"\t"+str(allele1)+"\t"+str(allele2)+"\t"+str(allele3)+"\t"+str(allele4)+"\n"
		output_structure.write(result_analysis)
output_structure.close()
