Design of Experiments in Python
===============================

Description
-----------

The ``py-design`` package defines the Python module ``design`` which implements
several routines for the design of experiments. Basically, it serves as
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
+ [``demos/demo3.py``](demos/demo3.py): Latin Random Square Design.
+ [``demos/demo4.py``](demos/demo4.py): Adjust a ``D`` dimensional dataset of ``N`` points so that it forms a Latin hypercube.
+ [``demos/demo5.py``](demos/demo5.py): Sparse Grid: Clenshaw Curtis Closed Fully Nested rule.
+ [``demos/demo6.py``](demos/demo6.py): Sparse Grid: Fejer 1 Open Fully Nested rule.
+ [``demos/demo7.py``](demos/demo7.py): Sparse Grid: Fejer 2 Open Fully Nested rule.
+ [``demos/demo8.py``](demos/demo8.py): Sparse Grid: Gauss Patterson Open Fully Nested rule.
+ [``demos/demo9.py``](demos/demo9.py): Sparse Grid: Gauss Legendre Open Weakly Nested rule.
+ [``demos/demo10.py``](demos/demo10.py): Sparse Grid: Gauss Hermite Open Weakly Nested rule.
+ [``demos/demo11.py``](demos/demo11.py): Sparse Grid: Gauss Laguerre Open Non Nested rule.
+ [``demos/demo12.py``](demos/demo12.py): Generate the Faure quasirandom sequence.
+ [``demos/demo13.py``](demos/demo13.py): Generate the Halton quasirandom sequence.
+ [``demos/demo14.py``](demos/demo14.py): Generate the Hammersley quasirandom sequence.
+ [``demos/demo15.py``](demos/demo15.py): Generate the Sobol quasirandom sequence.
+ [``demos/demo16.py``](demos/demo16.py): Generate the Lambert quasirandom sequence.
+ [``demos/demo17.py``](demos/demo17.py): Generate the Improved Distributed Hypercube Sequence.


TODO
----
+ Add references to each algorithm.
