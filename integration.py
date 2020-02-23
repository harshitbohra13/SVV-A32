import numpy as np
from scipy.integrate import quad

    
def integral(f, a, b, dx=10**(-5), steps=None):

    if steps is None and dx is not None:
        steps = (b-a)/dx

    elif dx is None and steps is not None:
        dx = (b-a)/steps

    else: raise ValueError

    I = 0

    for i in range(int(steps)):
        I += f(a+i*dx)*dx

    return I

f1 = lambda x: np.sin(x)
print(integral(f1, 0, 3.14))