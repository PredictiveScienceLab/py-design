Design of Random Experiments in Python
======================================

Description
-----------

The ``py-design`` package defines the Python module ``design`` which implements
several routines for the design of random experiments. Basically, it serves as
a wrapper the Fortran 90 codes for experimental design written by 
[John Burkardt](http://people.sc.fsu.edu/~jburkardt/). I have collected, 
probably, all of them here.


Related Packages
----------------

To the best of my knowledge, there is also another Python package implementing
several designs called [PyDOE](http://pythonhosted.org/pyDOE/index.html). I
concentrate more on what is known as **randomized designs** used in sampling
models in order to create surrogate surfaces as well as performing Monte Carlo
tasks.


Demos
-----

Here are some demos demonstrating how to use the package:
+ [``demos/demo1.py``](demos/demo1.py): Centered Latin Square Design.
+ [``demos/demo2.py``](demos/demo2.py): Latin Edge Square Design.
