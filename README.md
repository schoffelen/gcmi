# gcmi : Gaussian-Copula Mutual Information

Functions for calculating mutual information and other information theoretic quantities using a parametric Gaussian copula. 

## Installation

#### Matlab

Add the contents of the `matlab` directory to your Matlab path.

#### Python

...

## Usage

This is an overview of the functions available. Fuller descriptions of function arguments are provided in the in-line documentaion (Python docstrings, Matlab help strings).

The function names follow a format where the first part of the name is the quantity which the function calculates (eg mi, cmi, gcmi) followed by an underscore, followed by the types of variables that function acts on, in order corresponding to the arguments - g: Gaussian, d: Discrete, c: Continuous (any distribution). 

Please note a crucial difference in the interface between the Python and Matlab versions for multidimensional variables. In keeping with the convention of each platform in Matlab samples should be the first axis of the array (data stored in Fortran order); in Python samples should be the last axis of the array (data stored in C order by default).

| Platform | Multivariate data layout    |
| :------: | :----------------------:    |
| Matlab   | samples FIRST axis (rows)   |
| Python   | samples LAST axis (columns) |

Discrete inputs are passed as two arguments, the vector of values over samples y, and an integer parameter Ym specifying the size of the discrete space. In Python discrete variables should be stored in an integer data type array; in Matlab a double is used but should contain only integer values. y takes values between `0` and `Ym-1` inclusive. Empty classes are not supported.

Care should be taken with continuous variables that contain many repeated values. The copula transform which depends on a rank ordering will not be well defined in such cases. Possible approaches include repeated calculations while jittering the data with low amplitude noise to avoid the numerically equivalent values, or using binning and discrete methods. 

For functions with a `biascorrect` option, this is an optional true or false switch (default true) which indicates whether analytic bias correction for the entropy of Gaussian variables is applied. The bias correction increases computation time and is not needed when combined with permutation testing.

#### GCMI functions

`gcmi` functions estimate the Gaussian Copula Mutual Information, including input data checking and the copula transform step. We suggest new who are trying out the measure start with these functions.

*  `I = gcmi_cc(x,y)` 

    Calculate GCMI between two (possibly multidimensional) continuous variables x and y. x and y can have any marginal distribution but should not contain repeated values. 

*  `I = gcmi_cd(x,y,Ym)` 

    Calculate GCMI between a (possibly multidimensional) continuous variables x and a discrete y (with values between 0 and Ym-1 inclusive).

*  `I = gccmi_ccc(x,y,z)` 

    Calculate GCMI between two (possibly multidimensional) continuous variables x and y conditioned on a third, z. x / y / z can have any marginal distribution but should not contain repeated values.

*  `I = gccmi_ccd(x,y,z,Zm)` 

    Calculate GCMI between two (possibly multidimensional) continuous variables x and y conditioned on a third discrete variable, z (with values between 0 and Zm-1 inclusive). x and y can have any marginal distribution but should not contain repeated values.

#### Low level functions

These functions implement the different steps for the GCMI calculation. They are provided separately for computational efficiency (e.g. copula transform only needs to be performed once prior to permutation testing).

##### Copula transformation functions

*  `c = ctransform(x)`

    Compute empirical CDF value (copula transform) for each value. If x is >2D transformation is performed on each dimension separately. 

*  `cx = copnorm(x)` 

    Perform copula normalisation (equivalent to `norminv(ctransform(x))`). Returns standard normal samples with rank ordering preserved. If x is >2D normalization is performed on each dimension separately.  

##### Information theoretic quantities for Gaussian variables

These functions calculate information theoretic quantities (mi: mutual information, cmi: conditional mutual information) for Gaussian variables. Together with copula normalization above they implement the GCMI estimator: `gcmi(x,y) = mi_gg(copnorm(x),copnorm(y))`.

*  `I = mi_gg(x,y,biascorrect)` 

    Calculate MI between two (possibly multidimensional) Gaussian variables x and y. 

*  `I = mi_gd(x,y,Ym,biascorrect)` 

    Calculate MI between a (possibly multidimensional) Gaussian variable x and a discrete y (with values between 0 and Ym-1 inclusive). 

*  `I = cmi_ggg(x,y,biascorrect)` 

    Calculate CMI between two (possibly multidimensional) Gaussian variables x and y conditioned on a third, z. 

##### Miscellaneous functions

*  `H = ent_g(x, biascorrect)`

    Return analytic entropy of a (possibly multidimensional) Gaussian variable. 

