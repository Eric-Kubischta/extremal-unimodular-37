////////////////////////////////////////////////////////////////////////
// Cayley lattice: full construction and certification
// arXiv:2611.07138
//
// Constructs Lambda in Q^40 from the incidence vectors v_ij,
// kappa, theta, and tau (Sections 2-3 of the paper).
// Then verifies its invariants, shadow, Z/4Z-code, and automorphisms.
// No external data files are required.
////////////////////////////////////////////////////////////////////////

Q := Rationals();

// Standard orthonormal basis of Q^40:
// e, e_i^R, e_i^C, e_i^S, with i in Z/13Z.

E := [Vector(Q, [j eq i select 1 else 0 : j in [1..40]])
      : i in [1..40]];

e  := E[1];
eR := E[2..14];
eC := E[15..27];
eS := E[28..40];

// ---------------------------------------------------------------
// 1. The incidence lattice L (Section 2)
// ---------------------------------------------------------------

// v_ij = e + e_i^R + e_j^C + e_{i+j}^S.

vij := [];
for i in [0..12] do
    for j in [0..12] do
        Append(~vij, e + eR[i+1] + eC[j+1]
                           + eS[((i+j) mod 13)+1]);
    end for;
end for;

X := Matrix(Q, #vij, 40, &cat[Eltseq(x) : x in vij]);
L := Lattice(X);

print "\n(1) Incidence lattice L";
print "Number of cell vectors:", #vij;
print "All v_ij have norm 4:",
      &and[Norm(x) eq 4 : x in vij];
print "Rank(L):", Rank(L);
print "L is integral:", IsIntegral(L);
print "L is even:", IsEven(L);
print "det(L):", Determinant(L);

// ---------------------------------------------------------------
// 2. The intermediate lattice Gamma (Section 2)
// ---------------------------------------------------------------

// vbar = (1/169) sum_{i,j} v_ij.
// kappa = (13/4) vbar.

vbar := (&+vij)/169;
kappa := (13/4)*vbar;

XGamma := Matrix(Q, #vij+1, 40,
                 &cat[Eltseq(x) : x in (vij cat [kappa])]);
Gamma := Lattice(XGamma);

print "\n(2) Intermediate lattice Gamma";
print "Norm(kappa):", Norm(kappa);
print "4*kappa = sum_i v_ii:",
      4*kappa eq &+[vij[14*i+1] : i in [0..12]];
print "Rank(Gamma):", Rank(Gamma);
print "Gamma is integral:", IsIntegral(Gamma);
print "Gamma is odd:", not IsEven(Gamma);
print "det(Gamma):", Determinant(Gamma);
print "Index [Gamma:L]:", Index(Gamma,L);

// Check that kappa is characteristic in Gamma.
kg := Gamma ! Eltseq(kappa);
print "kappa is characteristic in Gamma:",
      &and[
        (Integers()!(InnerProduct(kg,b)-Norm(b))) mod 2 eq 0
        : b in Basis(Gamma)
      ];

// ---------------------------------------------------------------
// 3. The glue vectors theta and tau (Section 3)
// ---------------------------------------------------------------

// Fundamental weights omega_X of the three A_12 blocks.

omegaR := eR[1] - (&+eR)/13;
omegaC := eC[1] - (&+eC)/13;
omegaS := eS[1] - (&+eS)/13;

// Weyl vector rho, with the sign reversed in the symbol block.

rho := &+[(i-6)*(eR[i+1] + eC[i+1] - eS[i+1])
         : i in [0..12]];

// The two glue vectors of Lemmas 3.3 and 3.6.

s1 := eS[2] - eS[1];

theta := omegaR + 3*omegaC - 4*omegaS;

tau := -4*rho/13 - 2*(omegaR+omegaC-omegaS) - s1;

print "\n(3) Glue vectors";
print "Norm(theta):", Norm(theta);
print "Norm(tau):", Norm(tau);
print "Inner product (theta,tau):",
      InnerProduct(theta,tau);

// Check theta,tau belong to the dual of Gamma.
gammaGens := vij cat [kappa];

print "theta belongs to Gamma^*:",
      &and[Denominator(InnerProduct(theta,x)) eq 1
           : x in gammaGens];

print "tau belongs to Gamma^*:",
      &and[Denominator(InnerProduct(tau,x)) eq 1
           : x in gammaGens];

// ---------------------------------------------------------------
// 4. The completed Cayley lattice Lambda (Section 3)
// ---------------------------------------------------------------

// Lambda = Gamma + Z*theta + Z*tau.

gens := vij cat [kappa,theta,tau];

X := Matrix(Q, #gens, 40, &cat[Eltseq(x) : x in gens]);
Lambda := Lattice(X);

print "\n(4) Cayley lattice Lambda";
print "Rank(Lambda):", Rank(Lambda);
print "Lambda is integral:", IsIntegral(Lambda);
print "Lambda is odd:", not IsEven(Lambda);
print "det(Lambda):", Determinant(Lambda);
print "Lambda is unimodular:",
      IsIntegral(Lambda) and Determinant(Lambda) eq 1;

print "Index [Lambda:Gamma]:", Index(Lambda,Gamma);
print "Index [Lambda:L]:", Index(Lambda,L);

// Check that kappa remains characteristic.
kl := Lambda ! Eltseq(kappa);
print "kappa is characteristic in Lambda:",
      &and[
        (Integers()!(InnerProduct(kl,b)-Norm(b))) mod 2 eq 0
        : b in Basis(Lambda)
      ];

// ---------------------------------------------------------------
// 5. Minimum and kissing number (Proposition 3.7)
// ---------------------------------------------------------------

print "\n(5) Extremality and kissing number";

print "Computing minimum of Lambda...";
time m := Minimum(Lambda);
print "Minimum(Lambda):", m;
print "Extremal (minimum 4):", m eq 4;

print "Computing kissing number...";
time k := KissingNumber(Lambda);
print "KissingNumber(Lambda):", k;
print "Kissing number equals 35594:", k eq 35594;


// ---------------------------------------------------------------
// 6. Shadow and s-extremality (Proposition 3.7)
// ---------------------------------------------------------------

// Since kappa is characteristic,
// S(Lambda) = kappa/2 + Lambda.
// Compute its minimum by the closest-vector algorithm.

print "\n(6) Shadow and s-extremality";
print "Computing minimum norm in the shadow...";

time shadowVectors, shadowMin := ClosestVectors(
    Lambda, -kappa/2 : Max := 1
);

print "Minimum norm of the shadow:", shadowMin;
print "Shadow minimum equals 13/4:", shadowMin eq 13/4;

print "8*min(Lambda) + 4*min(shadow):", 8*m+4*shadowMin;
print "rank(Lambda) + 8:", Rank(Lambda)+8;

print "Lambda is s-extremal:",
      8*m+4*shadowMin eq Rank(Lambda)+8;

// ---------------------------------------------------------------
// 7. A 4-frame and its Z/4Z-code (Section 5, Appendix A)
// ---------------------------------------------------------------

// Nine representatives from Table 1, multiplied by 26.
// Coordinates are (a; r, c, s).

z := [0 : i in [1..13]];
u := [26] cat [0 : i in [1..12]];

reps := [
<26, u, u, u>,
<13, [-23] cat [3 : i in [1..12]],
 [9,-17,9,-17,9,9,9,9,9,9,-17,9,-17],
 [1 : i in [1..13]]>,
<0, [0,0,0,0,0,26,-26,-26,26,0,0,0,0], z, z>,
<0,
 [0,2,4,6,8,10,-14,14,-10,-8,-6,-4,-2],
 [0,2,4,6,8,10,-14,14,-10,-8,-6,-4,-2],
 [0,-2,-4,-6,-8,16,-12,12,-16,8,6,4,2]>,
<0,
 [0,4,8,12,-10,-6,-2,2,6,10,-12,-8,-4],
 [0,4,8,-14,16,-6,-2,2,6,-16,14,-8,-4],
 [0,-4,-8,14,10,6,2,-2,-6,-10,-14,8,4]>,
<0,
 [0,26,0,-26,0,0,0,0,0,0,0,0,0],
 [0,0,-26,0,26,0,0,0,0,0,0,0,0], z>,
<0,
 [-14,-10,20,-2,2,6,10,-12,-8,-4,0,4,8],
 [10,14,-8,-4,0,4,8,-14,-10,-6,-2,2,6],
 [4,0,-4,-8,-12,10,6,2,-2,-6,-10,12,8]>,
<0,
 [6,-10,0,10,-6,4,14,-2,8,-8,2,-14,-4],
 [-8,2,-14,-4,6,-10,0,10,20,4,-12,-2,8],
 [2,-8,8,-2,-12,4,-6,10,0,-10,6,-4,12]>,
<13,
 [-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13],
 [-7,-5,-3,-1,1,3,5,7,9,-15,13,-11,17],
 [5,3,1,-1,-3,-5,-7,-9,-11,13,11,9,7]>
];

// Index reversal i -> -i (mod 13).
function ReverseBlock(x)
    return [x[1]] cat [x[14-i] : i in [1..12]];
end function;

// Generate the orbits under <C, M^6>, up to sign.
frame := [];
orbitSizes := [];

for q in reps do
    a := q[1]; r := q[2]; c := q[3]; s := q[4];
    oldSize := #frame;

    for j in [0..2] do
        for eps in [0..1] do
            rr := r; cc := c; ss := s;

            if eps eq 1 then
                rr := ReverseBlock(rr);
                cc := ReverseBlock(cc);
                ss := ReverseBlock(ss);
            end if;

            f := Vector(Q, [a] cat rr cat cc cat ss)/26;

            if not (f in frame or -f in frame) then
                Append(~frame, f);
            end if;
        end for;

        // C(a;r,c,s) = (a;c,s^vee,r^vee).
        oldr := r; oldc := c; olds := s;
        r := oldc;
        c := ReverseBlock(olds);
        s := ReverseBlock(oldr);
    end for;

    Append(~orbitSizes, #frame-oldSize);
end for;

F := Matrix(Q, #frame, 40,
            &cat[Eltseq(f) : f in frame]);

// Pairings of lattice generators with frame vectors.
P := X*Transpose(F);

print "\n(7) 4-frame and self-dual Z/4Z-code";
print "Orbit sizes:", orbitSizes;
print "Orbit sizes agree with Appendix A:",
      orbitSizes eq [1,3,3,3,3,6,6,6,6];
print "Number of frame vectors:", #frame;
print "First frame vector is v_00:", frame[1] eq vij[1];

print "Frame has Gram matrix 4*I_37:",
      #frame eq 37 and
      F*Transpose(F) eq 4*IdentityMatrix(Q,37);

// Integral pairings imply frame membership in Lambda,
// since Lambda is unimodular.
frameIntegral := &and[
    Denominator(t) eq 1 : t in Eltseq(P)
];

print "Frame vectors pair integrally with lattice generators:",
      frameIntegral;
print "Frame vectors lie in Lambda (by unimodularity):",
      frameIntegral and IsIntegral(Lambda)
                    and Determinant(Lambda) eq 1;

// C = { ((x,f_i) mod 4)_i : x in Lambda }.
Z4 := IntegerRing(4);

PZ4 := ChangeRing(ChangeRing(P, Integers()), Z4);
C4 := LinearCode(PZ4);

print "Length of Z/4Z-code:", Length(C4);
print "Number of codewords:", #C4;
print "Code has size 2^37:", #C4 eq 2^37;
print "Code is self-dual:", IsSelfDual(C4);
print "Code is Type I (because Lambda is odd):",
      IsSelfDual(C4) and not IsEven(Lambda);

// v_01 supplies a word of Euclidean weight 16.
// The lattice minimum supplies the matching lower bound.
witnessWeight := EuclideanWeight(PZ4[2]);

print "Euclidean weight of codeword from v_01:",
      witnessWeight;
print "Minimum Euclidean weight is 16 (bound + witness):",
      m eq 4 and witnessWeight eq 16;

// ---------------------------------------------------------------
// 8. Automorphism group (Section 4)
// ---------------------------------------------------------------

// The Cayley-table symmetries:
// T: (i,j) -> (i+1,j+3)
// M: (i,j) -> (2i,2j)
// C: (i,j) -> (j,-i-j)
// Indices are taken modulo 13.

S40 := Sym(40);

T := S40 ! ([1]
    cat [2+((i+1) mod 13) : i in [0..12]]
    cat [15+((i+3) mod 13) : i in [0..12]]
    cat [28+((i+4) mod 13) : i in [0..12]]);

M := S40 ! ([1]
    cat [2+((2*i) mod 13) : i in [0..12]]
    cat [15+((2*i) mod 13) : i in [0..12]]
    cat [28+((2*i) mod 13) : i in [0..12]]);

C := S40 ! ([1]
    cat [28+((-i) mod 13) : i in [0..12]]
    cat [2+i : i in [0..12]]
    cat [15+((-i) mod 13) : i in [0..12]]);

// Apply a coordinate permutation to an ambient vector.
function ActAmbient(x, g)
    y := Vector(Q, [0 : i in [1..40]]);
    for i in [1..40] do
        y[i^g] := x[i];
    end for;
    return y;
end function;

// Check the incidence generators, kappa, and glue vectors.
function PreservesLambda(g)
    fixesKappa := ActAmbient(kappa,g) eq kappa;
    permutesCells := &and[ActAmbient(x,g) in vij : x in vij];
    preservesGlue := &and[
        &and[Denominator(InnerProduct(ActAmbient(x,g),y)) eq 1
             : y in gens] : x in [theta,tau]
    ];
    return fixesKappa and permutesCells and preservesGlue;
end function;

print "\n(8) Automorphism group";
print "T preserves Lambda:", PreservesLambda(T);
print "M preserves Lambda:", PreservesLambda(M);
print "C preserves Lambda:", PreservesLambda(C);

Haff := sub<S40 | T,M>;
H := sub<S40 | T,M,C>;
Z := M^8*C;

print "Orders of T, M, C:", [Order(g) : g in [T,M,C]];
print "Order <T,M> (AGL(1,13)):", #Haff;
print "Order <T,M,C>:", #H;

// Magma uses a right-action convention.
print "M^-1*T*M = T^2:", M^-1*T*M eq T^2;
print "C^-1*T*C = T^3:", C^-1*T*C eq T^3;
print "M*C = C*M:", M*C eq C*M;
print "M^8*C has order 3:", Order(Z) eq 3;
print "M^8*C is central:",
      Z*T eq T*Z and Z*M eq M*Z and Z*C eq C*Z;

print "Adjoining -I gives a subgroup of order:", 2*#H;

// Independent full automorphism-group computation.
print "Computing the full automorphism group of Lambda...";
time AutLambda := AutomorphismGroup(Lambda);

print "Order Aut(Lambda):", #AutLambda;
print "Aut(Lambda) = <T,M,C,-I>:",
      PreservesLambda(T) and PreservesLambda(M) and
      PreservesLambda(C) and #AutLambda eq 2*#H;

print "Aut(Lambda) = AGL(1,13) x C6:",
      #Haff eq 156 and #H eq 468 and #AutLambda eq 936 and
      M^-1*T*M eq T^2 and C^-1*T*C eq T^3 and
      M*C eq C*M and Order(Z) eq 3 and
      Z*T eq T*Z and Z*M eq M*Z and Z*C eq C*Z;


// ---------------------------------------------------------------
// 9. Gram matrix
// ---------------------------------------------------------------

print "\n(9) Gram matrix of Lambda";
print "37 x 37 integral Gram matrix:";
print GramMatrix(Lambda);

// Stop with an error if any of the principal invariants fails.
assert Rank(L) eq 37 and IsEven(L) and Determinant(L) eq 16*13^4;
assert Norm(kappa) eq 13 and Index(Gamma,L) eq 4;
assert IsIntegral(Lambda) and Determinant(Lambda) eq 1;
assert Rank(Lambda) eq 37 and not IsEven(Lambda);
assert m eq 4 and k eq 35594 and shadowMin eq 13/4;
assert #frame eq 37 and frameIntegral and IsSelfDual(C4);
assert #C4 eq 2^37 and witnessWeight eq 16;
assert #AutLambda eq 936 and 2*#H eq #AutLambda;

print "\n=== Full certification passed ===";
