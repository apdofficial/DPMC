import sys,os
from cnfParser import parseMCC21CNF
from heapq import nlargest

k = int(sys.argv[1])
inF = sys.argv[2]
tmp = os.environ['TMP']

nVars, nCls, clauses, probType, wEncountered, litWts, projVars, inds = parseMCC21CNF(inF)

counts = {}

for cl in clauses:
	for l in cl:
		if l in counts.keys():
			counts[l] += 1
		else:
			counts[l] = 1

minLCounts = []
for v in range(1,nVars+1):
	cp = 0
	cn = 0
	if v in counts.keys():
		cp = counts[v]
	if -v in counts.keys():
		cn = counts[-v]
	minLCounts += [min(cp,cn)]

vList = nlargest(k,range(len(minLCounts)),key = lambda idx: minLCounts[idx]) 

print(vList)
# arjun weighted command:
#/usr/bin/time -o /projects/vardi/prob_inf/results/arjun/timing/mcc/arjun_array_${SLURM_ARRAY_TASK_ID}.time timeout -k 1s 1000s ~/arjun/arjun-1190e8a --renumber=0 --input $FNAME /projects/vardi/prob_inf/benchmarks/cnf/arjun/arjun_array_${SLURM_ARRAY_TASK_ID}.cnf
for i in range(2**k):
	print('Doing '+str(i)+' out of 2^'+str(k))
	'''
	ls=[]
	for j,v in enumerate(vList,start=0):
		print((i>>j)%2,end='')
		l = v if (i>>j)%2 == 1 else -v
		ls.append(l)
	print(' '+str(ls))
	
	'''
	cnfF = 'tmp_'+str(i)+'.cnf'
	if os.path.exists(tmp+'/'+cnfF):
		os.remove(tmp+'/'+cnfF)
	f = open(tmp+'/'+cnfF,'w')
	f.write('p cnf '+str(nVars)+' '+str(nCls+k)+'\n')
	for cl in clauses:
		for l in cl:
			f.write(str(l)+' ')
		f.write('0\n')
	for j,v in enumerate(vList,start=0):
		l = v if (i>>j)%2 == 1 else -v
		f.write(str(l)+' 0\n')
	f.close()
	if os.path.exists(tmp+'/arjun_'+cnfF):
		os.remove(tmp+'/arjun_'+cnfF)
	arjun_cmd = 'timeout -k 1s 120s /home/adi/Downloads/arjun_mis/arjun-1190e8a --input '+tmp+'/'+cnfF+' '+tmp+'/arjun_'+cnfF+' >/dev/null'
	os.system(arjun_cmd)
	dpmc_cmd = '/home/adi/Downloads/dpmc//lg/build/lg "/home/adi/Downloads/dpmc/lg/solvers/htd-master/bin/bin/htd_main -s 1234567 --print-progress --strategy challenge  --opt width --iterations 0  --preprocessing full --patience 20" < '+tmp+'/arjun_'+cnfF+' | /home/adi/Downloads/dpmc/dmc/dmc --pw=10 --cf '+tmp+'/arjun_'+cnfF+' --dp=s --lc=1 --wc=0 --vs=2 --dv=5 --dy=0 --mm=6500 --jp=s --sa=2 --aa=1 --er=0 --tc=7 --tr=3 --ir=5'
	os.system(dpmc_cmd)

