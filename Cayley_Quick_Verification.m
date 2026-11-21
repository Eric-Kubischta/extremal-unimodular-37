////////////////////////////////////////////////////////////////////////
// Cayley lattice: quick verification from the explicit Gram matrix
// arXiv:2611.07138
// Requires gram_matrix_Cayley.m in the same directory.
////////////////////////////////////////////////////////////////////////

load "gram_matrix_Cayley.m";

print "\n=== Cayley lattice: quick verification ===";
print "Gram matrix is 37 x 37:", Nrows(G) eq 37 and Ncols(G) eq 37;
print "Gram matrix is symmetric:", G eq Transpose(G);

// Here Lambda is defined solely by the stored Gram matrix.
Lambda := LatticeWithGram(G);

print "Rank(Lambda):", Rank(Lambda);
print "Lambda is integral:", IsIntegral(Lambda);
print "Lambda is odd:", not IsEven(Lambda);
print "det(Lambda):", Determinant(Lambda);
print "Lambda is unimodular:",
      IsIntegral(Lambda) and Determinant(Lambda) eq 1;

print "Computing minimum norm...";
time m := Minimum(Lambda);
print "Minimum(Lambda):", m;
print "Extremal (minimum 4):", m eq 4;

print "Computing kissing number...";
time k := KissingNumber(Lambda);
print "KissingNumber(Lambda):", k;
print "Kissing number equals 35594:", k eq 35594;

assert Nrows(G) eq 37 and Ncols(G) eq 37;
assert G eq Transpose(G);
assert Rank(Lambda) eq 37;
assert IsIntegral(Lambda) and not IsEven(Lambda);
assert Determinant(Lambda) eq 1;
assert m eq 4 and k eq 35594;

print "\n=== Quick verification passed ===";
