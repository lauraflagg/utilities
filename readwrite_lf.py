import shutil, csv, os, time
import numpy as np

#November 3, 2017: added time to string of backupup file name





def read2cols(file_to_read,skip=0):
	#skip is integer of lines 2 skip
	f = open(file_to_read, 'r')
	lines=f.readlines()
	f.close()
	count=len(lines)

	xl=[]

	yl=[]

	c=skip
	while c<count:
		if lines[c][0]!='#':
	#		print y
	#	skipping any commented line:
			nums=lines[c].split()
	#		print nums
			xl.append(float(nums[0]))
			yl.append(float(nums[1]))
		c=c+1
	wl1=np.array(xl)
	fl1=np.array(yl)
	
	return wl1, fl1
	

def read3cols(file_to_read):

	f = open(file_to_read, 'r')
	lines=f.readlines()
	f.close()
	count=len(lines)

	xl=[]

	yl=[]
	zl=[]

	c=0
	while c<count:
		if lines[c][0]!='#':
	#		print y
	#	skipping any commented line:
			nums=lines[c].split()
	#		print nums
			xl.append(float(nums[0]))
			yl.append(float(nums[1]))
			zl.append(float(nums[2]))
		c=c+1
	wl1=np.array(xl)
	fl1=np.array(yl)
	er1=np.array(zl)
	
	return wl1, fl1, er1


def backupandwrite(filename,heads,data):
	#heads is a list
	#data is an array
	if os.path.exists(filename):
		dt=os.path.getmtime(filename)
		dts=time.strftime("%d%b%Y_%Hh%Mm",time.gmtime(dt))
		nfn=filename[0:-4]+'_'+dts+'.csv'
		shutil.copyfile(filename, nfn)
	
		
	with open(filename, 'w') as fp: 
		writer = csv.writer(fp, delimiter=',', lineterminator = '\n') # 
		writer.writerow(heads) # write header 
		writer.writerows(data)

def write(filename,heads,data):
	with open(filename, 'wb') as fp: 
		writer = csv.writer(fp, delimiter=',') # 
		writer.writerow(heads) # write header 
		writer.writerows(data)