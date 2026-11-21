# The Cayley lattice: an extremal unimodular lattice in rank 37

Magma code accompanying [arXiv:2611.07138](https://arxiv.org/abs/2611.07138).

These files describe the **Cayley lattice** $\Lambda$ from the paper: an odd, extremal unimodular lattice of rank 37.

## Files

| File | Description |
| --- | --- |
| [`gram_matrix_Cayley.m`](gram_matrix_Cayley.m) | The **exact $37\times37$ Gram matrix** printed by the full construction in the author's Magma V2.30-1 run, with all 37 rows written explicitly. The file defines `G`. |
| [`Cayley_Quick_Verification.m`](Cayley_Quick_Verification.m) | Constructs $\Lambda$ directly from `G` and checks rank, integrality, oddness, determinant, minimum, and kissing number; no Cayley construction is needed. |
| [`Cayley_Construction_Full_Certification.m`](Cayley_Construction_Full_Certification.m) | Reconstructs $\Lambda$ **from scratch in $\mathbb{Q}^{40}$**, using the vector coordinates $e,e_i^R,e_j^C,e_k^S$ and $v_{ij}$ from the paper, followed by $\kappa,\theta,\tau$. It verifies the intermediate lattices $L$ and $\Gamma$, minimum and kissing number, shadow and $s$-extremality, a 4-frame and self-dual $\mathbb{Z}_4$-code, and the full automorphism group. It also prints `GramMatrix(Lambda)`. |

The explicit matrix uses **the same row and column order** as the matrix printed by the author's successful full Magma run; it is not merely an isometric replacement obtained by a basis change. A different Magma version may select a different lattice basis during the $\mathbb{Q}^{40}$ construction, even though it produces the same lattice.

## Results

| Invariant | Value |
| --- | --- |
| Rank | 37 |
| Determinant | 1 |
| Parity | Odd |
| Minimum squared norm | 4 |
| Kissing number | 35,594 |
| Shadow minimum squared norm | $13/4$ |
| Order of full automorphism group | 936 |
| Code | Type I self-dual $\mathbb{Z}_4$-code of length 37, minimum Euclidean weight 16 |

## Running with Magma

**Local installation:** Put all three `.m` files in the same directory, change the Magma working directory to that directory, and run:

```magma
load "Cayley_Quick_Verification.m";
```

For the full construction and certification, run:

```magma
load "Cayley_Construction_Full_Certification.m";
```

**[Free Magma Online Calculator](https://magma-maths.org/calc/):** Each run starts a fresh session and cannot load your local files by path.

- **Quick check:** Paste the entire contents of `gram_matrix_Cayley.m` into the input, followed by `Cayley_Quick_Verification.m` **without** its `load "gram_matrix_Cayley.m";` line. Run both together in one submission.
- **Full check:** Paste `Cayley_Construction_Full_Certification.m` directly and run it. It needs no other files.

The author's full run with **Magma V2.30-1** took approximately **151 seconds**, including about **146 seconds** for the independent full automorphism-group computation. The online calculator currently limits each run to five minutes and 1 GiB of memory, so runtime and success may vary.

### Optional exact matrix comparison

To verify that the stored matrix is *entry-for-entry* the matrix produced by the construction, run the following in a **single** local Magma session:

```magma
load "Cayley_Construction_Full_Certification.m";
load "gram_matrix_Cayley.m";
assert GramMatrix(Lambda) eq G;
print "Stored Gram matrix matches the construction exactly.";
```

On the online calculator, concatenate the full construction code and the matrix data, then append the last two lines (omitting the `load` commands). The exact equality test is basis-dependent; in different Magma versions, the Gram matrices may differ by an integral change of basis.

## Scope

These scripts reproduce the **Cayley lattice construction** and verify the listed computational invariants. Proofs of the other results in the paper are contained in the manuscript.
