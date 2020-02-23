import numpy as np
from scipy.integrate import quad

'''
Define fun with lamda constructor to define variables
'''
def integral(f, a, b, dx=10**(-5), steps=None):

    if steps is None and dx is not None:
        steps = (b-a)/dx

    elif dx is None and steps is not None:
        dx = (b-a)/steps

    else: raise ValueError

    I = 0 
    for i in range(int(steps)):
        I += f(a+i*dx)*dx

    return sum(I)


f1 = lambda x: np.sin(x)
x = np.arange(0,1,0.1)
x = np.transpose(x)
f1 = lambda x: np.arange(0, 1, 0.1)*x
print(integral(f1, 0, 1, 0.1))