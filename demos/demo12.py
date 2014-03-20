"""
Generate the Faure quasirandom sequence.

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
# Draw the design
X = design.faure(num_points, num_dim)
print X
# And plot it
plt.plot(X[:, 0], X[:, 1], '.', markersize=10)
plt.xlabel('$x_1$', fontsize=16)
plt.ylabel('$x_2$', fontsize=16)
plt.title('Faure Quasirandom Sequence', fontsize=16)
plt.show()
