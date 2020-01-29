import copy
import collections
file=raw_input("Enter the name of the pdb file without any Extension: ")
filename = file + ".pdb"
data = open(filename,"r")
output = open("Side_ChainOutput.pdb","w")
#output.write("                                  SIDE-CHAIN ATOMS                              \n \n")
data.readline()
Array_Atom=[]
array=[]
AminoCount=0
count =1
NullCount=1
dict={}
def getatoms():
	global AminoCount
	global a,b
	for line in data:
		a=''
		b=''
		list1=line.split()
		if 'ATOM' in line:
			if (list1[0]=='ATOM'):
				if (list1[1]=='1'):
					AminoCount=list1[5]
				if(len(list1[4])>1):
					for u in range(0,len(list1[4])):
						if (u<1):
							a=a+(list1[4][u])
						else:
							b=b+(list1[4][u])
					list1[4]=a
					list1.insert(5,b)
					a=''
					b=''
				if(len(list1[7])>6):
					for u in range(0,len(list1[7])):
						if (u<6):
							a=a+(list1[7][u])
						else:
							b=b+(list1[7][u])
					list1[7]=a
					list1.insert(8,b)
					a=''
					b=''
				if (len(list1[9]) > 4):
					for u in range(0,len(list1[9])):
						if (u<4):
							a=a+(list1[9][u])
						else:
							b=b+(list1[9][u])
					list1[9]=a
					list1.insert(10,b)
					a=''
					b=''
				if (len(list1[2])>4):
					for u in range(0,len(list1[2])):
						if (u<4):
							a=a+(list1[2][u])
						else:
							b=b+(list1[2][u])
					list1[2]=a
					list1.insert(3,b)
				Array_Atom.append(list1)
	getatombyamino(Array_Atom,AminoCount)
	return
def gethelix():
	data = open(filename,"r")
	data.readline()
	helixarray=[]
	global count5	
	global AminoCount
	global a,b
	templist={}
	for line in data:
		a=''
		b=''
		list1=line.split()
		if 'HELIX' in line:
			if (list1[0]=='HELIX'):
				AminoCount=int(list1[5])
				lastAminoCount= int(list1[8])
				templist[AminoCount]=lastAminoCount

	templist=collections.OrderedDict(sorted(templist.items()))
	gethelixatoms(templist)

def gethelixatoms(templist):
	global count5 
	count5=1
	dictofca={}
	data = open(filename,"r")
	data.readline()
	helixarray=[]
	for k in templist:
		for y in range(k,templist[k]+1):
			for line in data:
				a=''
				b=''
				list1=line.split()
				if 'ATOM' in line:
					if (list1[0]=='ATOM'):
						if(len(list1[4])>1):
							for u in range(0,len(list1[4])):
								if (u<1):
									a=a+(list1[4][u])
								else:
									b=b+(list1[4][u])
							list1[4]=a
							list1.insert(5,b)
							a=''
							b=''
						if(len(list1[7])>6):
							for u in range(0,len(list1[7])):
								if (u<6):
									a=a+(list1[7][u])
								else:
									b=b+(list1[7][u])
							list1[7]=a
							list1.insert(8,b)
							a=''
							b=''
						if (len(list1[9]) > 4):
							for u in range(0,len(list1[9])):
								if (u<4):
									a=a+(list1[9][u])
								else:
									b=b+(list1[9][u])
							list1[9]=a
							list1.insert(10,b)
							a=''
							b=''
						if (len(list1[2])>4):
							# print"hi"
							for u in range(0,len(list1[2])):
								if (u<4):
									a=a+(list1[2][u])
								else:
									b=b+(list1[2][u])
							list1[2]=a
							list1.insert(3,b)
						if (int(list1[5])==y):
								helixarray.append(list1)
						elif (int(list1[5])>y):
							break
	for line in helixarray :
		if (line[2]=='CA'):
			dictofca[count5]=line
			count5=count5+1								
	helixfinal(dictofca)
def helixfinal(dictofca):
	helixoutput = open("HelixOutput.pdb","w")
	#helixoutput.write("                                  HELICES                              \n \n")
	totalx=0
	totaly=0
	totalz=0
	othertotal=0
	finallist=copy.deepcopy(dictofca)
	t=1
	while(t in dictofca):
		totalx=float(dictofca[t][6])+float(dictofca[t+1][6])+float(dictofca[t+2][6])+float(dictofca[t+3][6])
		avgx=round(totalx/4,3)
		totaly=float(dictofca[t][7])+float(dictofca[t+1][7])+float(dictofca[t+2][7])+float(dictofca[t+3][7])
		avgy=round(totaly/4,3)
		totalz=float(dictofca[t][8])+float(dictofca[t+1][8])+float(dictofca[t+2][8])+float(dictofca[t+3][8])
		avgz=round(totalz/4,3)
		othertotal=float(dictofca[t][10])+float(dictofca[t+1][10])+float(dictofca[t+2][10])+float(dictofca[t+3][10])
		avg=round(othertotal/4,3)
		finallist[t][1]=t
		finallist[t][2]='S'
		finallist[t][6]=float(avgx)
		finallist[t][5]=int(dictofca[t][5])
		finallist[t][7]=float(avgy)
		finallist[t][8]=float(avgz)
		finallist[t][9]=float(dictofca[t][9])
		finallist[t][10]=float(avg)
		t=t+1
		if (len(dictofca)-t < 3):
			break
	value=len(finallist)-1
	value2=len(finallist)-2
	del finallist[len(finallist)]
	del finallist[value]
	del finallist[value2]
	for r in finallist:
		a= "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(finallist[r][0],finallist[r][1],finallist[r][2],"",finallist[r][3],finallist[r][4],finallist[r][5],"",finallist[r][6],finallist[r][7],finallist[r][8],finallist[r][9],finallist[r][10],finallist[r][11],"")
		data.close()
		helixoutput.write(a+'\n')
def getatombyamino(Array_Atom,AminoCount):
	global array
	global count
	global NullCount
	for line in Array_Atom :
		list2=line
		if (list2[5]==AminoCount):
			if (list2[2] not in ('N','C','CA','O','H')):
				array.append(list2)
		else:
			dict[int(AminoCount)]=array
			AminoCount=str(int(AminoCount)+1)
			array=[]
			if (list2[2] not in ('N','C','CA','O','H')):
				array.append(list2)
	dict[int(AminoCount)]=array
	copydict=copy.deepcopy(dict)
	for j in dict:
		length = len(dict[j])
		if (length==0):
			NullCount=NullCount+1
			continue
		xcor = float(dict[j][0][6])
		ycor = float(dict[j][0][7])
		zcor = float(dict[j][0][8])
		last = float(dict[j][0][10])
		for k in range(1,length):
			xcor = float(dict[j][k][6]) + xcor 
			ycor = float(dict[j][k][7]) + ycor 
			zcor = float(dict[j][k][8]) + zcor 
			last = float(dict[j][k][10]) + last 
		averagex=round(xcor/length,3)
		averagey=round(ycor/length,3)
		averagez=round(zcor/length,3)
		averagelast=round(last/length,3)
		copydict[j][0][1]=count
		copydict[j][0][2]='S'
		copydict[j][0][6]=averagex
		copydict[j][0][7]=averagey
		copydict[j][0][8]=averagez
		copydict[j][0][5]=int(copydict[j][0][5])
		copydict[j][0][9]=float(copydict[j][0][9])
		copydict[j][0][10]=averagelast
		count=count+1
		r=copydict[j][0]
		a= "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(r[0],r[1],r[2],"",r[3],r[4],r[5],"",r[6],r[7],r[8],r[9],r[10],r[11],"")
		output.write(a+'\n')
		data.close()
	return

def getbeta():
	data = open(filename,"r")
	data.readline()
	betaarray=[]
	global count8	
	global AminoCount
	global a,b
	betalist={}
	for line in data:
		a=''
		b=''
		list1=line.split()
		if 'SHEET' in line:
			if (list1[0]=='SHEET'):
				AminoCount=int(list1[6])
				lastAminoCount= int(list1[9])
				betalist[AminoCount]=lastAminoCount
	if(len(betalist)==0):
		print "This particular file does not have any beta strands in it."
		exit()
	betalist=collections.OrderedDict(sorted(betalist.items()))

	getbetaatoms(betalist)

def getbetaatoms(betalist):
	global count8
	count8=1
	dictofbeta={}
	data = open(filename,"r")
	data.readline()
	betaarray=[]
	for line in data:
		a=''
		b=''
		list1=line.split()
		if 'ATOM' in line:
			if (list1[0]=='ATOM'):
				if(len(list1[4])>1):
					for u in range(0,len(list1[4])):
						if (u<1):
							a=a+(list1[4][u])
						else:
							b=b+(list1[4][u])
					list1[4]=a
					list1.insert(5,b)
					a=''
					b=''
				if(len(list1[7])>6):
					for u in range(0,len(list1[7])):
						if (u<6):
							a=a+(list1[7][u])
						else:
							b=b+(list1[7][u])
					list1[7]=a
					list1.insert(8,b)
					a=''
					b=''
				if (len(list1[9]) > 4):
					for u in range(0,len(list1[9])):
						if (u<4):
							a=a+(list1[9][u])
						else:
							b=b+(list1[9][u])
					list1[9]=a
					list1.insert(10,b)
					a=''
					b=''
				if (len(list1[2])>4):
					for u in range(0,len(list1[2])):
						if (u<4):
							a=a+(list1[2][u])
						else:
							b=b+(list1[2][u])
					list1[2]=a
					list1.insert(3,b)
				if (list1[2] in ('N','CA','C')):
					betaarray.append(list1)
	finalbetaprocess(betaarray,betalist)

def finalbetaprocess(betaarray,betalist):
	count=1
	finaldictofbeta={}
	for k in betalist:
		for e in range(k,betalist[k]+1):
			for some in betaarray:
				if (int(some[5])==int(e)):
					finaldictofbeta[count]=some
					count=count+1
	finalprocessofbeta(finaldictofbeta)

def finalprocessofbeta(finaldictofbeta):
	betaoutput = open("BetaSheetsOutput.pdb","w")
	#betaoutput.write("                                  SHEETS                               \n")
	totalx=0
	totaly=0
	totalz=0
	othertotal=0
	finalbetalist=copy.deepcopy(finaldictofbeta)
	t=1
	while(t in finaldictofbeta):
		totalx=float(finaldictofbeta[t][6])+float(finaldictofbeta[t+1][6])+float(finaldictofbeta[t+2][6])+float(finaldictofbeta[t+3][6])
		avgx=round(totalx/4,3)
		totaly=float(finaldictofbeta[t][7])+float(finaldictofbeta[t+1][7])+float(finaldictofbeta[t+2][7])+float(finaldictofbeta[t+3][7])
		avgy=round(totaly/4,3)
		totalz=float(finaldictofbeta[t][8])+float(finaldictofbeta[t+1][8])+float(finaldictofbeta[t+2][8])+float(finaldictofbeta[t+3][8])
		avgz=round(totalz/4,3)
		othertotal=float(finaldictofbeta[t][10])+float(finaldictofbeta[t+1][10])+float(finaldictofbeta[t+2][10])+float(finaldictofbeta[t+3][10])
		avg=round(othertotal/4,3)
		finalbetalist[t][1]=int(t)
		finalbetalist[t][2]='S'
		finalbetalist[t][6]=float(avgx)
		finalbetalist[t][5]=int(finaldictofbeta[t][5])
		finalbetalist[t][7]=float(avgy)
		finalbetalist[t][8]=float(avgz)
		finalbetalist[t][9]=float(finaldictofbeta[t][9])
		finalbetalist[t][10]=float(avg)

		value = t+1
		value1=t+2
		del finalbetalist[value]
		del finalbetalist[value1]	
		t=t+3	
		if (len(finaldictofbeta)-t < 5):
			break
	value2 = t
	value3=t+1
	value4=t+2
	del finalbetalist[value2]
	del finalbetalist[value3]	
	del finalbetalist[value4]
	for r in finalbetalist:
		a= "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(finalbetalist[r][0],finalbetalist[r][1],finalbetalist[r][2],"",finalbetalist[r][3],finalbetalist[r][4],finalbetalist[r][5],"",finalbetalist[r][6],finalbetalist[r][7],finalbetalist[r][8],finalbetalist[r][9],finalbetalist[r][10],finalbetalist[r][11],"")
		data.close()
		betaoutput.write(a+'\n')
print("\nThe output Files \n* BetaSheetsOutput.pdb\n* HelixOutput.pdb\n* Alpha_CarbonOutput.pdb\nare created")
getatoms()
gethelix()
getbeta()