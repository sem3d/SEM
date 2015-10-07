# Exemple de transformation de coordonnees

import pyproj
import sys

defp1 = sys.argv[1]
defp2 = sys.argv[2]
input = sys.argv[3]
output = sys.argv[4]

p1 = pyproj.Proj(defp1)
p2 = pyproj.Proj(defp2)

f1 = open(input,"r")
f2 = open(output,"w")

for line in f1:
    (x,y,z) = [ float(u) for u in line.split() ]
    u,v,w = pyproj.transform(p1,p2,x,y,z)
    f2.write("%f %f %f\n" % (u,v,w) )

