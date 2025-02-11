"""
Adjust an D dimensional dataset of N points so that it forms a Latin hypercube.

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
# The number of points you want
num_points = 50
# Draw some random points
X0 = np.random.rand(num_points, num_dim)
# Latinize it
X = design.latinize(X0)
# Look at it
print('Before:')
print(X0)
print('After:')
print(X)
# And plot it
plt.plot(X0[:, 0], X[:, 1], 'b.', markersize=10)
plt.plot(X[:, 0], X[:, 1], 'r.', markersize=10)
plt.legend(['Original', 'Latinized'], loc='best')
plt.xlabel('$x_1$', fontsize=16)
plt.ylabel('$x_2$', fontsize=16)
plt.title('Latinize Demonstration', fontsize=16)
plt.show()
