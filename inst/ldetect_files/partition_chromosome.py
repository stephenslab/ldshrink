#!/usr/bin/env python3

import sys, os, gzip, math

infile = gzip.open(sys.argv[1])
nind = float(sys.argv[2])
outfile = open(sys.argv[3], "w")
N = 5000


pos2gpos = dict()
poss = list()
line = infile.readline()
while line:
	line = line.strip().split()
	pos = int(line[1])
	gpos = float(line[2])
	pos2gpos[pos] = gpos
	poss.append(pos)
	line = infile.readline()

print(len(poss))
nsnp = len(poss)

chunk = float(nsnp)/float(N)
print(chunk)
chunk = int(math.floor(chunk))
print(chunk)

for i in range(chunk):
	start = i*N
	end = i*N + N
	if i == chunk-1:
		end = len(poss)-1
		startpos = poss[start]
		endpos = poss[end]
		# print >> outfile, startpos, endpos
		print(str(startpos)+' '+str(endpos), file=outfile)
		continue
	startpos = poss[start]
	endpos = poss[end-1]
	endgpos = pos2gpos[endpos]
	test =end+1
	testpos = poss[test]
	stop = False
	while stop == False:
		if test == len(poss):
			stop = True
			continue
		testpos = poss[test]
		testgpos = pos2gpos[testpos]
		df = testgpos-endgpos
		tmp = math.exp(	-4.0*11418.0*df / (2.0*nind))
		if tmp < 1.5e-8:
			stop = True
		else:
			test = test+1
	# print >> outfile, startpos, testpos
	print(str(startpos)+' '+str(testpos), file=outfile)



