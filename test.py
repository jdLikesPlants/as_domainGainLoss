#!/usr/bin/env python2.7

def main():
	args = getArgs()
	tmap, dat, dom, gtf = args.tmap, args.dat, args.dom, args.gtf
	transMap, geneDict  = storeTmap(tmap = tmap)
	exonDict, chromoDict = parseGtf(gtf = gtf)  
	asDict, asCoordDict = parseSuppa(dat = dat, geneDict = geneDict, exonDict = exonDict, transMap = transMap)
	domainsDict, geneDomains = parseDomains(dom = dom, geneDict = geneDict)
	findChange(transMap = transMap, asDict = asDict, domainsDict = domainsDict, asCoordDict = asCoordDict, geneDomains = geneDomains, exonDict = exonDict)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = 'identify gain/loss of protein domains as a result of alternative splicing')
	parser.add_argument('-tmap', action = 'store', type = str, required = True, help = 'transmap file')
	parser.add_argument('-dat', action = 'store', type = str, required = True, help = 'pasa .dat file')
	parser.add_argument('-dom', action = 'store', type = str, required = True, help = 'hmmerscan genomic domains file')
	parser.add_argument('-gtf', action = 'store', type = str, required = True, help = 'filtered gtf file')
	args = parser.parse_args()
	return(args)

def storeTmap(tmap):
	import re,io
	from collections import defaultdict
	transToGene = defaultdict(dict)
	geneToTrans = defaultdict(dict)
	geneList = []
	fh = open(tmap, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		gene, tid = arr[0], arr[1]
		if gene in transToGene:
			transToGene[gene].append(tid)
		else:
			transToGene[gene] = [tid]
		if gene not in geneList:
			geneList.append(gene)
		geneToTrans[tid] = gene
	return(transToGene, geneToTrans)

def parseGtf(gtf):
	import re,io
	from collections import defaultdict
	exonDict = defaultdict(dict)
	chromoDict = defaultdict(dict)
	fh = open(gtf, "r")
	for line in fh:
		line = line.strip()
		if(re.search('exon', line)):
			arr = line.split('\t')
			chromo = arr[0]
			tid = arr[8].split(';')[0].split('\"')[1]
			if not tid in chromoDict:
				chromoDict[tid] = chromo
			left = arr[3]
			right = arr[4]
			exons = left + ":" + right
			if tid not in exonDict:
				exonDict[tid] = [exons]
			elif tid in exonDict:
				exonDict[tid].append(exons)
	return(exonDict, chromoDict)




def parseSuppa(dat, geneDict, exonDict, transMap):
	import io,re
	from collections import defaultdict
	from as_event_metadata import parseDat
	asDict = defaultdict(dict)
	asCoordDict = defaultdict(dict)
	asDict, asCoordDict = parseDat(dat = dat, geneDict = geneDict, exonDict = exonDict, transMap = transMap)
	return(asDict, asCoordDict)

def parseDomains(dom, geneDict):
	import io,re
	from collections import defaultdict
	domainsDict, geneDomains = defaultdict(dict), defaultdict(dict)
	fh = open(dom, "r")
	for line in fh:
		line = line.strip()
		arr = line.split('\t')
		#print(line)
		domain = arr[5]
		#tid = arr[1].split(':')[2]
		tid = arr[1]
		gene = geneDict[tid]
		loc = arr[2] + ":" + arr[3]
		
		if domain in domainsDict:
			if loc in domainsDict[domain]:
				domainsDict[domain][loc].append(tid)		
			elif loc not in domainsDict[domain]:
				domainsDict[domain][loc] = [tid]
		if domain not in domainsDict:
			domainsDict[domain][loc] = [tid]
		#print(tid)	
		#print(gene)
		if gene in geneDomains:
			if domain not in geneDomains[gene]:
				geneDomains[gene].append(domain)
		elif gene not in geneDomains:
			geneDomains[gene] = [domain]
	return(domainsDict, geneDomains)

def findChange(transMap, asDict, domainsDict, asCoordDict, geneDomains, exonDict):
	from collections import defaultdict
	import re
	#print(str(geneDomains["Potri.011G058400"]))
	#print(str(domainsDict["PF02365.14"]))  gain tids in the dict
	for gene in transMap:
		#if gene == "Potri.011G058400":
			#print(gene)
			for dom in geneDomains[gene]:       ##### looping domains for the gene here
				#if gene == "Potri.011G058400":
				#	print(dom)
				hasDom = defaultdict(dict)
				lacksDom = []
				for tid in transMap[gene]:				
					for loc in domainsDict[dom]:
						if tid in domainsDict[dom][loc]: ##and tid not in hasDom:
							hasDom[tid][loc] = dom
						
						#if gene == "Potri.011G058400":
						#	print(loc)
						
							
							
							
							#if gene == "Potri.011G058400":
							#	print(tid)
							
				
				
							#if gene == "Potri.011G058400":
							#	print(str(hasDom))
				for tid in transMap[gene]:
					if tid not in hasDom:
						lacksDom.append(tid)			

				asTypes = defaultdict(dict)
				lacksAS = []
				
				for tid in hasDom:

					### exists multiple isos or iso(s) lacking domain
					if(len(transMap[gene]) < 2 or len(lacksDom) < 1):
						continue
					for domLoci in hasDom[tid]:
						### redundant conditional?
						if dom in hasDom[tid][domLoci]:
							locList = domLoci.split(':')
							domLeft, domRight = int(locList[0]), int(locList[1])
						
						
							for lossTid in lacksDom:
								if(lossTid in asDict):
									for asType in asDict[tid]:
										### maybe dont nest
										for lossType in asDict[lossTid]:   


											hasLocList = asDict[tid][asType]
											lossLocList = asDict[lossTid][lossType]
											for hasLoc in hasLocList:
												coordsArr = hasLoc.split(':')
												if(re.search('-', hasLoc)):
													hasLeft, hasRight = int(coordsArr[2]), int(coordsArr[1])
												else:						
													hasLeft, hasRight = int(coordsArr[1]), int(coordsArr[2])

												for lossLoc in lossLocList:
													coordsArr = lossLoc.split(':')
													if(re.search('-', lossLoc)):
														lossLeft, lossRight = int(coordsArr[2]), int(coordsArr[1])
													else:						
														lossLeft, lossRight = int(coordsArr[1]), int(coordsArr[2])
													hasCheck, lossCheck = domCheck(asType = asType, lossType = lossType, hasLeft = hasLeft, hasRight = hasRight, \
													lossLeft = lossLeft, lossRight = lossRight, domLeft = domLeft, domRight = domRight, \
																				   loci = hasLoc, asDict = asDict, lossTid = lossTid)															
													if(hasCheck == 1 and lossCheck == 1):
														### store in dict
														print("loss" + '\t' + gene + '\t' + dom + '\t' + str(domLeft) + '-'+ str(domRight) + '\t' + lossTid + '\t' + lossType +  '\t' + str(lossLeft) \
															  + '-' + str(lossRight))
														print("gain" + '\t' + gene + '\t' + dom + '\t' + str(domLeft) + '-'+ str(domRight) + '\t' + tid + '\t' + asType +  '\t' + str(hasLeft) \
															  + '-' + str(hasRight))
														#gainDict[asType][tid] = hasLoc
														#lossDict[lossType][lossTid] = lossLoc


												
def domCheck(asType, lossType, hasLeft, hasRight, lossLeft, lossRight, domLeft, domRight, loci, asDict, lossTid):
	hasCheck, lossCheck = 0, 0
	
####### AE ##########
	if(asType == "AE" and hasCheck == 0 and lossType == "AE"):
		#print(str(domLeft) + '-'+ str(domRight) + '\t'  + lossType +  '\t' + str(hasLeft) +'\t' + str(hasRight) + '\t' + str(lossLeft)  + '-' + str(lossRight))
		if(hasLeft >= domRight and hasRight <= domLeft):								
			hasCheck = 1
		if(hasLeft <= domLeft and hasRight >= domRight):					## 5' or 3' AE, whole domain in exon
			hasCheck = 1		
		elif(hasLeft >= domLeft and hasRight >= domRight and hasLeft <= domRight):	## 3' AE
			hasCheck = 1
		elif(hasLeft <= domLeft and hasRight <= domRight and hasRight >= domLeft):		## 5' AE
			hasCheck = 1
		if(hasCheck == 1):
			lossCheck = 1
			
			
			
####### AE ##########

#### A3 A5 ##########
	elif(asType == "A5" and lossType == "A5"):
		if(hasLeft <= domRight and hasLeft >= domLeft):								
			hasCheck = 1
		elif(hasRight >= domLeft and hasRight <= domRight):
			hasCheck = 1
	
		if(lossLeft == hasLeft or lossRight == hasRight):
			if(hasCheck == 0):
				if(lossLeft >= domLeft and lossLeft <= domRight):
					lossCheck, hasCheck = 1, 1
				elif(lossRight <= domRight and lossRight >= domLeft):
					lossCheck, hasCheck = 1, 1
			else:
				lossCheck = 1
				
	elif(asType == "A3" and lossType == "A3"):
		if(hasLeft <= domRight and hasLeft >= domLeft):								
			hasCheck = 1
		elif(hasRight >= domLeft and hasRight <= domRight):
			hasCheck = 1
	
		if(lossLeft == hasLeft or lossRight == hasRight):
			if(hasCheck == 0):
				if(lossLeft >= domLeft and lossLeft <= domRight):
					lossCheck, hasCheck = 1, 1
				elif(lossRight <= domRight and lossRight >= domLeft):
					lossCheck, hasCheck = 1, 1
			else:
				lossCheck = 1		
				
	
#### A3 A5 ##########			
			
#### SI RI ##########

### needs to confirm if in domain
	elif(asType == "SI" and lossType == "RI"):
		hasCheck, lossCheck = checkSI(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, loci = loci, riLeft = lossLeft, \
									  riRight = lossRight, domLeft = domLeft, domRight = domRight)
	elif(asType == "RI" and lossType == "SI"):
		hasCheck, lossCheck = checkRI(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, loci = loci, siLeft = lossLeft, \
									  siRight = lossRight, domLeft = domLeft, domRight = domRight)			
#### SI RI ##########

#### SE RE ##########
	elif(asType == "RE" and lossType == "SE"):
		hasCheck, lossCheck = checkRE(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, domLeft = domLeft, \
									  domRight = domRight, seLeft = lossLeft, seRight = lossRight)
	elif(asType == "SE" and lossType == "RE"):
		hasCheck, lossCheck = checkSE(asDict = asDict, leftCoord = hasLeft, rightCoord = hasRight, domLeft = domLeft, \
									  domRight = domRight, reLeft = lossLeft, reRight = lossRight)
#### SE RE ##########
	return(hasCheck, lossCheck)
################################################			
			
					
										
			

def getA3Coords(left, right, exonDict, tid, loci):
	import re
	exons = exonDict[tid]
	#print(str(exons))
	if(re.search('-', loci)):
		right -= 1
		left -= 1
		for i in range(0, len(exonDict[tid])):
			exonArr = exons[i].split(':')
			#print(tid + '\t' + str(left) + '\t' + exonArr[1])
			if(right == int(exonArr[1])):
				exRight = exons[i + 1].split(':')[0]
				right = int(exRight)
				#print(tid + '\t' + str(left) + '\t' + str(right))
						
	else:
		right += 1 
		for i in range(0, len(exonDict[tid])):
						
			exonArr = exons[i].split(':')
			if(right == int(exonArr[0])):
				exLeft = exons[i -1].split(':')[1]
				left = int(exLeft)

	return(left, right)

def getA5Coords(left, right, exonDict, tid, loci):
	import re
	exons = exonDict[tid]
	if(re.search('-', loci)):
		left += 1 
		for i in range(0, len(exonDict[tid])):
			exonArr = exons[i].split(':')
			if(left == int(exonArr[0])):
				exLeft = exons[i -1].split(':')[1]
				left = int(exLeft)
	else:
		left -= 1
		for i in range(0, len(exonDict[tid])):
			exonArr = exons[i].split(':')
			#print(tid + '\t' + str(left))
			#print(str(exonArr))
			if(left == int(exonArr[1])):
				exRight = exons[i + 1].split(':')[0]
				right = int(exRight)
	return(left, right)

def checkSI(asDict, leftCoord, rightCoord, riLeft, riRight, domLeft, domRight, loci):	
	import re
	domCheck, lossCheck = 0, 0
	#if("RI" in asTypes):
	#for lossTid in lacksDom:
	if(re.search('-', loci)):

		if(leftCoord <= domRight and leftCoord >= domLeft):                                                             
				domCheck = 1
		elif(rightCoord >= domLeft and rightCoord <= domRight):
				domCheck = 1
	else:
		if(leftCoord >= domLeft and leftCoord <= domRight):                                                             
				domCheck = 1
		elif(rightCoord <= domRight and rightCoord >= domLeft):
				domCheck = 1

	if(riLeft == leftCoord and riRight == rightCoord):
		lossCheck = 1
		#if gene == "Potri.001G008400":
		#	print(gene + '\t' + loc + '\t' + str(asTypes) + '\t' + str(lacksAS))
		
	return(domCheck, lossCheck)

def checkRI(asDict, leftCoord, rightCoord, siLeft, siRight, domLeft, domRight, loci):
	import re
	domCheck, lossCheck = 0, 0
	#for lossTid in lacksDom:			
	if(re.search('-', loci)):
		if(leftCoord <= domRight and leftCoord >= domLeft):                                                             
				domCheck = 1
		elif(rightCoord >= domLeft and rightCoord <= domRight):
				domCheck = 1
	else:
		if(leftCoord >= domLeft and leftCoord <= domRight):                                                             
				domCheck = 1
		elif(rightCoord <= domRight and rightCoord >= domLeft):
				domCheck = 1

	if(siLeft == leftCoord and siRight == rightCoord):
		lossCheck = 1
		#if gene == "Potri.001G008400":
		#	print(gene + '\t' + loc + '\t' + str(asTypes) + '\t' + str(lacksAS))
	return(domCheck, lossCheck)

def checkSE(asDict, leftCoord, rightCoord, reLeft, reRight, domLeft, domRight):
	domCheck, lossCheck = 0, 0
	if(leftCoord > domLeft and rightCoord < domRight and domCheck == 0):		
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	elif(leftCoord < domLeft and rightCoord > domRight and domCheck == 0):  ### dom contined in the RE
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	elif(leftCoord > domLeft and rightCoord > domRight and leftCoord < domRight and domCheck == 0):
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	elif(leftCoord < domLeft and rightCoord < domRight and rightCoord > domLeft and domCheck == 0):
		domCheck = reCheck(seLeft = leftCoord, seRight = rightCoord, reLeft = reLeft, reRight = reRight)	
	if(domCheck == 1):
		lossCheck = 1
	return(domCheck, lossCheck)
	
def checkRE(asDict, leftCoord, rightCoord, seLeft, seRight, domLeft, domRight):
	domCheck, lossCheck = 0, 0
	if(leftCoord > domLeft and rightCoord < domRight and domCheck == 0):		
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	elif(leftCoord < domLeft and rightCoord > domRight and domCheck == 0):  ### dom contined in the RE
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	elif(leftCoord > domLeft and rightCoord > domRight and leftCoord < domRight and domCheck == 0):
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	elif(leftCoord < domLeft and rightCoord < domRight and rightCoord > domLeft and domCheck == 0):
		domCheck = reCheck(reLeft = leftCoord, reRight = rightCoord, seLeft = seLeft, seRight = seRight)	
	if(domCheck == 1):
		lossCheck = 1
	return(domCheck, lossCheck)	

def reCheck(reLeft, reRight, seLeft, seRight):
	domCheck = 0
	if(reLeft > seLeft and reRight < seRight):		
		domCheck = 1
	return(domCheck)



	

if __name__ == "__main__":
	main()
