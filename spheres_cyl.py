# Calculate sphere centers on surface of cylinder.
#
# Usage:
#    python spheres_cyl.py <cyl_radius> <cyl_height>  <sphere_radius>  <name-of-celltype>
#    e.g., confirm working
#    python spheres_cyl.py 20 40 8.0 epi
#
#    then, save to a file in PhysiCell cells.csv format:
#    python spheres_cyl.py 400.0 100 epi > fib3D.csv

import sys,string
# import math
import numpy as np

argc=len(sys.argv)
# print('argc=',argc)
# print('argv=',sys.argv)
# print('argv[0]=',sys.argv[0])
if argc < 5:
    print("Usage: <cyl_R> <cyl_height>  <sphere_R>  <name-of-celltype>")
    sys.exit(-1)

#R = string.atof(sys.argv[1])
cyl_R = float(sys.argv[1])
cyl_height = float(sys.argv[2])
sphere_R = float(sys.argv[3])
celltype_name = sys.argv[4]

x0 = 0.0
y0 = 0.0

delta_theta = np.arcsin(sphere_R / cyl_R)

o1_value = 0.0
o2_value = 360.0
rmod = 1
start_radians = o1_value * np.pi/180.
end_radians = o2_value * np.pi/180.

z0 = 0.0
zdel = 1.5*sphere_R
print(f"x,y,z,type,volume,cycle entry,custom:GFP,custom:sample")
kstack = 1
start_radians_del = 0.3
for zval in np.arange(z0,cyl_height,zdel):
    for theta in np.arange(start_radians, end_radians, rmod*2*delta_theta):
    # print("theta= ",theta)
        xval = x0 + cyl_R * np.cos(theta)
        yval = y0 + cyl_R * np.sin(theta)
        print(f"{xval},{yval},{zval},{celltype_name}")
    start_radians += kstack * start_radians_del
    kstack *= -1

