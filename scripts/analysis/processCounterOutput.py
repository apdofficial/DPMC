from outputFileProcessors import processDPMCFile, processSSTDFile
import json
from sys import argv
from os.path import exists

tdType = int(argv[1])

def processOutput(resDir, exps, fPrefix, processor):
	'''
	Process Flowcutter output files
	'''
	oDir = resDir + '/output/'
	tDir = resDir + '/timing/'
	expectedName = None
	allStats = {}
	for j in exps:
		stats = []
		for i in range(400):
			fp = open(tDir+'/'+str(j)+'/'+fPrefix+'_array_'+str(i)+'.time')
			line = fp.readline()
			fp.close()
			maxMem = int(line.split()[-1].split('m')[0])
			fp = open(oDir+'/'+str(j)+'/'+fPrefix+'_array_'+str(i)+'.out')
			procOut = processor(fp)
			stats.append([maxMem, procOut])
			fp.close()
		allStats[j] = list(stats)
	return allStats

if tdType == 1:	
	print("Process DPMC...")
	fcDir = '/home/adi/Downloads/prob_inf/results/dpmc/arjun_no-proj/'
	allStats = processOutput(fcDir,[1,2,3,4,5,6,7,8],'dpmc',processDPMCFile)
	with open("dpmcData.json", 'w') as f:
		json.dump(allStats, f, indent=2) 
elif tdType == 2:
	print("Process sstd...")
	fcDir = '/home/adi/Downloads/prob_inf/results/sstd/arjun_no-proj/'
	allStats = processOutput(fcDir,[1],'sstd',processSSTDFile)
	with open("sstdData.json", 'w') as f:
		json.dump(allStats, f, indent=2)  