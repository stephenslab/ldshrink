#!/usr/bin/env python3

import sys, os, gzip, math
import numpy as np

# calculate Wen/Stephens shrinkage LD estimate
panelf = open(sys.argv[1])
gmapfile = gzip.open(sys.argv[2]) # genetic map 
indfile = open(sys.argv[3]) #list of individuals
# NE = 11418.0
NE = float(sys.argv[4])
# CUTOFF = 1e-7
CUTOFF = float(sys.argv[5])
outfile = gzip.open(sys.argv[6], "wt") # outfile file

inds = list()
# for line in indfile.xreadlines():
for line in indfile:
        line = line.strip().split()
        inds.append(line[0])

haps= list()
theta = 0
nind = len(inds)
s = 0
for i in range(1, 2*nind):
        print(i)
        s = s+ 1.0/float(i)

s = 1/s
#print "s", s
theta = s/(2.0*float(nind)+s)
print(theta)

pos2gpos = dict()
line = gmapfile.readline()
while line:
        line = line.strip().split()
        pos = int(line[1])
        gpos = float(line[2])
        pos2gpos[pos] = gpos
        line = gmapfile.readline()

line = panelf.readline()

allpos = list()
allrs = list()
ind2index = dict()
while line:
        line = line.strip().split()
        if line[0][0] == "#" and line[0][1] == "#":
                line = panelf.readline()
                continue
        if line[0] == "#CHROM":
                for i in range(9, len(line)):
                        tmp = line[i]
                        if tmp in inds:
                                ind2index[tmp] = i
                line = panelf.readline()
                continue
        h = list()
        allpos.append(int(line[1]))
        allrs.append(line[2])
        for ind in inds:
        	index = ind2index[ind]
        	tmpg = line[index]
        	tmpg = tmpg.split(":")[0]
        	g = tmpg.split("|")
        	h.append(g[0])
        	h.append(g[1])
        haps.append(h)
        line = panelf.readline()

for i in range(len(allpos)):
        pos1 = allpos[i]
        gpos1 = pos2gpos[pos1]
        toofar = False
        j = i
        while j < len(allpos) and toofar == False:
                pos2 = allpos[j]
                gpos2 = pos2gpos[pos2]
                df = gpos2-gpos1
                ee = math.exp( - df*NE*4.0 / (2.0*float(len(inds))))
                if ee < CUTOFF:
                        toofar = True
                        j = j+1
                        continue
                g1 = haps[i]
                g2 = haps[j]
                n11 = 0
                n10 = 0
                n01 = 0
                for k in range(len(g1)):
                        if g1[k] == "1" and g2[k] == "1":
                                n11 = n11+1
                        elif g1[k] == "0" and g2[k] == "1":
                                n01 = n01 +1
                        elif g1[k] == "1" and g2[k] == "0":
                                n10 = n10 +1
                f11 = float(n11)/float(len(g1))
                f1 = (float(n11)+float(n10))/float(len(g1))
                f2 = (float(n11)+float(n01))/float(len(g1))
                D = f11 - f1*f2
                Ds = D*ee
                Ds2 = (1-theta)*(1-theta)*Ds
               	if math.fabs(Ds2) < CUTOFF:
                        j = j+1
                        continue
                if i == j:
                        Ds2 = Ds2 + (theta/2.0)*(1-theta/2.0)
                # print >> outfile, allrs[i], allrs[j], pos1, pos2, gpos1, gpos2, D, Ds2
                print(str(allrs[i])+' '+str(allrs[j])+' '+str(pos1)+' '+str(pos2)+' '+str(gpos1)+' '+str(gpos2)+' '+str(D)+' '+str(Ds2), file=outfile)
                j = j+1
