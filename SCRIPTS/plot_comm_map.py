import sys
from pylab import *

# Fichier genere par read_comm_map
fmap = sys.argv[1]
# Nombre de process / noeud
npc = int(sys.argv[2])


procs = {}
nodes = {}
f =open(fmap)
for line in f:
    line = line.strip()
    if not line:
        continue
    SRC, DST = line.split(":")
    SRC = int(SRC)
    DST = array([int(d) for d in DST.split(",")])
    procs[SRC] = DST
    SRC = SRC//npc
    DST = DST//npc
    node = nodes.setdefault(SRC, set())
    for nn in DST:
        node.add(nn)

def nodeplot(nodes):
    NP = len(nodes)
    for s in nodes.keys():
        for d in nodes[s]:
            if s>d:
                continue
            thetaS = 2*s*pi/NP
            thetaD = 2*d*pi/NP
            plot([cos(thetaS),cos(thetaD)],[sin(thetaS),sin(thetaD)],"o-")
    show()

def matplot(nodes):
    NP = len(nodes)
    M = zeros( (NP,NP), uint8 )
    for s in nodes.keys():
        for d in nodes[s]:
            c = 100
            if s//npc == d//npc:
                c = 50
            M[s,d] = c
    imshow(M)
    show()

#nodeplot(nodes)
matplot(procs)
    
