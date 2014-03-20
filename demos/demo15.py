"""
Generate the Sobol quasirandom sequence.

Author:
    Ilias Bilionis

Date:
    3/19/2014

"""


import design
import numpy as np
import matplotlib.pyplot as plt

# The number of input dimensions
num_dim = 2
# The number of points
num_points = 100
# How many elements to skip
skip = 10
# Draw the design
X = design.sobol(num_points, num_dim, skip=10)
print X
# And plot it
plt.plot(X[:, 0], X[:, 1], '.', markersize=10)
plt.xlabel('$x_1$', fontsize=16)
plt.ylabel('$x_2$', fontsize=16)
plt.title('Sobol Quasirandom Sequence', fontsize=16)
plt.show()
