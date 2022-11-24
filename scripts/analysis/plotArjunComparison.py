import matplotlib.pyplot as plt
import json

colors=['green','blue','red','brown','orange']
mkrs= ['o','^','+','x','D']

def loadData():
	with open("arjun-mis-proj_Data.json", 'r') as f:
		arjunTWs = json.load(f)
	with open("htdData.json", 'r') as f:
		htdTWs = json.load(f)
	with open("arjun_no-proj_Data.json", 'r') as f:
		anpTWs = json.load(f)

	return arjunTWs, htdTWs, anpTWs

def plotBars(plList, thresh, plotPMC = False):
	print(plList)
	fig, ax = plt.subplots()
	ax.bar(range(0,200,2), [0 if len(el[1])==0 else el[1][-1] if el[1][-1] < thresh else thresh for el in plList], label='Primal')
	ax.bar(range(0,200,2), [0 if len(el[2])==0 else el[2][-1] if el[2][-1] < thresh else thresh for el in plList], label='Incidence',alpha=0.8)
	if plotPMC == True:
		ax.bar(range(1,200,2), [0 if len(el[3])==0 else el[3][-1] if el[3][-1] < thresh else thresh for el in plList], label='PMC Primal',alpha=0.8)
		ax.bar(range(1,200,2), [0 if len(el[4])==0 else el[4][-1] if el[4][-1] < thresh else thresh for el in plList], label='PMC Incidence',alpha=0.8)
	ax.set_ylabel('TWs')
	ax.set_title('Benchmarks')
	ax.legend()
	#plt.savefig('graphs/primal-incidence_comparison/mcc21_track1-2.png',bbox_inches='tight')
	plt.show()

def plotHTDComparison(al, hl, thresh):
	#for each of 4 sets
		#for each of hl tl fl
			#create list
			#for each of 100 benchmarks
				#add time to list
		#plot bars
	j = 0
	for r in [list(range(100)),list(range(100,200)),list(range(200,300)),list(range(300,400))]:
		fig, ax = plt.subplots()
		i = 0
		for l in [al, hl]:
			outL = []
			for ind in r:
				outL.append(thresh if len(l[ind][1])==0 else l[ind][1][-1] if l[ind][1][-1] < thresh else thresh)
			#ax.bar(range(i,300,3), outL, label=['Arjun','Vanilla'][i], alpha=0.8)
			ax.bar(range(100), outL, label=['Arjun','Vanilla'][i], alpha=0.8)
			i += 1
		ax.set_ylabel('TWs')
		ax.set_title(['MCC21_t1','MCC21_t2','MCC22_t1','MCC22_t2'][j])
		ax.legend()
		#plt.savefig('graphs/primal-incidence_comparison/mcc21_track1-2.png',bbox_inches='tight')
		plt.show()		
		j += 1

def plotArjunComparison(hl, anpl, thresh):
	j = 0
	for r in [list(range(100)),list(range(100,200)),list(range(200,300)),list(range(300,400))]:
		fig, ax = plt.subplots()
		i = 0
		for l in [al, anpl]:
			outL = []
			for ind in r:
				outL.append(thresh if len(l[ind][1])==0 else l[ind][1][-1] if l[ind][1][-1] < thresh else thresh)
			#ax.bar(range(i,300,3), outL, label=['Arjun','Vanilla'][i], alpha=0.8)
			ax.bar(range(100), outL, label=['Vanilla','Arjun_no-proj'][i], alpha=0.8)
			i += 1
		ax.set_ylabel('TWs')
		ax.set_title(['MCC21_t1','MCC21_t2','MCC22_t1','MCC22_t2'][j])
		ax.legend()
		#plt.savefig('graphs/primal-incidence_comparison/mcc21_track1-2.png',bbox_inches='tight')
		plt.show()		
		j += 1
		
al, hl, anpl = loadData()
#plotFC(fl)
#plotHTDComparison(al, hl, 300)
plotArjunComparison(hl, anpl, 300)
