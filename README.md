# Root-Finding-Using-Matrices

Implimentation of root-finding algorithm for non-linear problems using companion matrices. Project completed as a part of final project for
an undergraduate class - Computational Physics 1.
 
# Motivation and Method: 
  Any root finding problem requires iterative techniques that converge to a single root; and is not particularly reliable to provide all
  roots of a function (even with different Initial conditions). In this project, we take the given root finding problem, and approximate
  it by its taylor series upto a provided tolerance. 
  Following the taylor approximation, we create the so called companion matrix assosciated to the polynomial. 
  The companion matrix M of a polynomial f has the property that the characteristic polynomial of M is f. 
  Upon obtaining the companion matrix, we run M through an eigenvalue solver (a QR scheme) to obtain the eigenvalues of M -
  Through this process, by definition, we then obain the roots of $f$ to degree 2 less than the degree of the taylor approximation. 

