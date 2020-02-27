import numpy as np

from Functions_Interpolation import find_interpolant,compute_interpolant,integrate_polynomial,integrate_area

ref_file = open("aerodynamicloadf100.dat", "r")
lines = ref_file.readlines()
ref_file.close()

aerodata = np.mat(np.zeros((81, 41)))

for line in lines:
    idx = lines.index(line)
    line = line.replace("\n", "")
    values = line.split(',')
    values = np.matrix(values)

    aerodata[idx, :] = values

print(aerodata[0])