def processTDFile(fp, expectedFName=None):
	fp.readline()
	fp.readline()
	fp.readline()
	fp.readline()
	fName = fp.readline()
	if expectedFName!=None:
		assert(fName == expectedFName)
	fp.readline()
	twList = []
	nvars = 0
	for line in fp:
		if line.startswith("s td"):
			tw = int(line.split()[3])
			nvars = int(line.split()[4])
			twList.append(tw)
	return nvars, twList

