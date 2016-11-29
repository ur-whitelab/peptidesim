
Goal
====
We want to be able to post-process explore/analyze the data. We need to know the following then:

   1. The peptides/amounts, density, etc (traits)
   2. The mdp parameters
   3. the trajectory locations
   4. metadata system of simulations


Developer Environment
===

First, prepare the docker image in the test-docker folder by running
the build script. Then, run the test script in the root directory. It
is only necessary to rebuild the docker script when newer gromacs,
gromacswrapper packages are available.
