# ITensor excited states

Finding excited states using the ITensor library <br />
Examples: 1d Ising and 1d Heisenberg model

ITensor - Software Library for Tensor Network Calculations <br />
by Matthew Fishman and Steven R. White and E. Miles Stoudenmire <br />
arxiv.org/abs/2007.14822 <br />
http://itensor.org/

In collaboration with Heitor Peres Casagrande <br />
https://github.com/heitorc7


## Starting guide
http://itensor.org/docs.cgi?vers=cppv3&page=install <br />
$ git clone https://github.com/ITensor/ITensor itensor <br />
$ cd itensor <br />
$ cp options.mk.sample options.mk <br />
Edit options.mk <br />
$ make <br />
Makefile -> Change path to ITensor in line 3 <br />
$ cd eigenvalue-dmrg/Examples <br />
$ make <br />
$ ./ExcitedDMRG <br />

Optional: Run plot.py

