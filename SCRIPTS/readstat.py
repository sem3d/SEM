#!/usr/bin/env python

# Read and synthetize stats.log

import sys
from math import sqrt

stats = []

def parse_stats(fn):
    f=file(fn)
    res = {}
    for line in f:
        l = line.split()
        if len(l)!=3:
            continue
        k = l[0]
        v = float(l[2])
        r = res.setdefault(k, [])
        r.append(v)
    stat = {}
    for key, values in res.items():
        N = len(values)
        avg = sum(values)/N
        ect = [ (v-avg)**2 for v in values ]
        std = sqrt(sum(ect)/(N))
        stat[key] = (N, avg, std)
    return stat

for fn in sorted(sys.argv[1:]):
    name = fn
    if name.endswith("/stat.log"):
       name = name[:-9]
    stats.append( (name, parse_stats(fn)) )


keys = stats[0][1].keys()

cols = [ max( len(k) for k in keys ) ]
for fn, st in stats:
    cols.append(max(len(fn),13))

for n in cols:
    print "="*n,
print
    
print " "*cols[0],
for i, (fn, st) in enumerate(stats):
    print fn.ljust(cols[i+1]),
print
for n in cols:
    print "-"*n,
print

for k in keys:
    print k.ljust(cols[0]),
    for n,st in enumerate(stats):
        #print st, k
        print (("%.1f" % st[1][k][1])+" ~"+("%5.1f" % st[1][k][2])).rjust(cols[n+1]),
    print

for n in cols:
    print "="*n,

