"""
Construct a centered Latin Square design.

Author:
    Ilias Bilionis

Date:
    3/19/2014

"""


import design
import matplotlib.pyplot as plt

# The number of input dimensions
num_dim = 2
# The number of points you want
num_points = 10
# Create the design
X = design.latin_center(num_points, num_dim)
# Look at it
print X
# And plot it
plt.plot(X[:, 0], X[:, 1], '.', markersize=10)
plt.xlabel('$x_1$', fontsize=16)
plt.ylabel('$x_2$', fontsize=16)
plt.title('Centered Latin Square Design', fontsize=16)
plt.show()
