# potential-inference
Infers the gravitational potential/mass distribution of a galaxy from observations of stars' positions and velocities.

Needs boost, [Eigen](http://eigen.tuxfamily.org), and [nanoflann](https://github.com/jlblancoc/nanoflann/blob/master/include/nanoflann.hpp).

To generate 300 stars with parameters 1.5 and 0.8, use

    ./generate 300 1.5 0.8 > stars.txt
     
To infer parameters, use

    ./infer < stars.txt
    
You can provide MCMC with reasonable starting parameters using

    ./infer 1.4 0.9 < stars.txt
    
By default, MCMC starts with parameters all equal to 1.
