import numpy as np

#Opening and reading the file text of the Fokker 100
file_f100 = open("aerodynamicloadf100.dat","r")
lines = file_f100.readlines()

#Inserting the values of the text file to a matrix 81 times 41
matrix = np.zeros((81,41))
for line in lines:
    row = lines.index(line)
    line = line.replace("\n","")
    magnitude_list = line.split(",")
    matrix[row,:] = magnitude_list

print(matrix)