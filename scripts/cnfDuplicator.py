from cnfParser import parseMCC21CNF
from sys import argv

inF = argv[1]
outF = argv[2]
numRep = int(argv[3]) #number of times to replicate cnf

nVars, nCls, clauses, probType, wEncountered, litWts, projVars, inds = parseMCC21CNF(inF)

nnVars = nVars * numRep
nnCls = nCls * numRep

of = open(outF,'w')
of.write('p cnf '+str(nnVars)+' '+str(nnCls)+'\n')
for i in range(numRep):
	for cl in clauses:
		for v in cl:
			sign = v>0
			outVar = abs(v)+nVars*i
			outLit = outVar if sign else -outVar
			of.write(str(outLit)+' ')
		of.write('0\n')
of.close()
