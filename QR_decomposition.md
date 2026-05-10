# Computing R in QR Decomposition

In QR decomposition, **R** is computed via **Gram-Schmidt orthogonalization** applied to the columns of A.

---

## Setup

Given an m×n matrix A with columns **a₁, a₂, ..., aₙ**, we want:
- **Q** — orthonormal columns **q₁, q₂, ..., qₙ**
- **R** — upper triangular matrix

---

## Step-by-step (Classical Gram-Schmidt)

For each column j = 1, 2, ..., n:

1. Start with **uⱼ = aⱼ**
2. Subtract its projections onto all previous orthonormal vectors:
   ```
   uⱼ = aⱼ - Σᵢ﹤ⱼ (qᵢᵀ aⱼ) qᵢ
   ```
3. Normalize: **qⱼ = uⱼ / ‖uⱼ‖**

---

## How R is Built

The entries of R come directly from those projections and norms:

| Entry | Formula | Meaning |
|-------|---------|---------|
| R[i,j] for i < j | qᵢᵀ aⱼ | dot products (off-diagonal) |
| R[j,j] | ‖uⱼ‖ | norms (diagonal) |
| R[i,j] for i > j | 0 | upper triangular by construction |

---

## Intuition

R encodes how much of each original column **aⱼ** lies along each orthonormal basis vector **qᵢ**. Because we process columns left to right and only subtract earlier basis vectors, the result is necessarily upper triangular — column j has no component in directions qⱼ₊₁ onwards.

---

## Worked Example

**A = matrix(1:6, nrow = 2)** fills column-by-column in R:

```
A = [1  3  5]
    [2  4  6]
```

Columns: **a₁ = [1, 2]**, **a₂ = [3, 4]**, **a₃ = [5, 6]**

---

**Column 1**

u₁ = a₁ = [1, 2]

R[1,1] = ‖u₁‖ = √(1² + 2²) = √5 ≈ 2.2361

q₁ = u₁ / √5 = [1/√5, 2/√5]

---

**Column 2**

R[1,2] = q₁ᵀ a₂ = (1/√5)(3) + (2/√5)(4) = 11/√5 ≈ 4.9193

u₂ = a₂ − R[1,2] q₁ = [3, 4] − (11/5)[1, 2] = [4/5, −2/5]

R[2,2] = ‖u₂‖ = √(16/25 + 4/25) = 2/√5 ≈ 0.8944

q₂ = u₂ / (2/√5) = [2/√5, −1/√5]

---

**Column 3**

R[1,3] = q₁ᵀ a₃ = (1/√5)(5) + (2/√5)(6) = 17/√5 ≈ 7.6026

R[2,3] = q₂ᵀ a₃ = (2/√5)(5) + (−1/√5)(6) = 4/√5 ≈ 1.7889

u₃ = a₃ − R[1,3] q₁ − R[2,3] q₂ = [0, 0]  *(a₃ is linearly dependent)*

---

**Result**

```
Q = [1/√5   2/√5 ]     R = [√5    11/√5   17/√5]
    [2/√5  −1/√5 ]         [ 0     2/√5    4/√5]

  ≈ [0.4472  0.8944]      ≈ [2.2361  4.9193  7.6026]
    [0.8944 −0.4472]         [0       0.8944  1.7889]
```

Verify: Q %*% R recovers A exactly. The third column of A being linearly dependent on the first two is reflected by u₃ = 0 — R's third column still exists but carries no new diagonal entry.

---

## Worked Example 2

**B = matrix(1:9, nrow=3, byrow=FALSE)** fills column-by-column (R's default):

```
B = [1  4  7]
    [2  5  8]
    [3  6  9]
```

Columns: **b₁ = [1, 2, 3]**, **b₂ = [4, 5, 6]**, **b₃ = [7, 8, 9]**

Note: b₃ = 2b₂ − b₁, so the columns are linearly dependent (rank 2).

---

**Column 1**

u₁ = b₁ = [1, 2, 3]

R[1,1] = ‖u₁‖ = √(1²+2²+3²) = √14 ≈ 3.7417

q₁ = [1/√14, 2/√14, 3/√14]

---

**Column 2**

R[1,2] = q₁ᵀ b₂ = (4+10+18)/√14 = 32/√14 ≈ 8.5522

u₂ = b₂ − R[1,2] q₁ = [4,5,6] − (32/14)[1,2,3] = [12/7, 3/7, −6/7]

R[2,2] = ‖u₂‖ = 3√21/7 ≈ 1.9640

q₂ = [4/√21, 1/√21, −2/√21]

---

**Column 3**

R[1,3] = q₁ᵀ b₃ = (7+16+27)/√14 = 50/√14 ≈ 13.3631

R[2,3] = q₂ᵀ b₃ = (28+8−18)/√21 = 18/√21 ≈ 3.9279

u₃ = b₃ − R[1,3] q₁ − R[2,3] q₂ = [0, 0, 0]  *(linear dependence)*

---

**Result**

```
Q = [1/√14   4/√21]     R = [√14    32/√14   50/√14]
    [2/√14   1/√21]         [ 0    3√21/7    18/√21]
    [3/√14  -2/√21]         [ 0       0          0 ]

  ≈ [0.2673   0.8729]      ≈ [3.7417   8.5522  13.3631]
    [0.5345   0.2182]         [0        1.9640   3.9279]
    [0.8018  -0.4364]         [0        0        0     ]
```

Q is 3×2 (thin QR — only two independent columns), R is 2×3. The zero row in R signals rank deficiency: B has no third independent direction, so Gram-Schmidt produces u₃ = 0 and no third diagonal entry.

---

## R's `qr()` Implementation

R uses **Householder reflections** rather than Gram-Schmidt — numerically more stable for large or ill-conditioned matrices. Each reflection zeroes out the sub-diagonal entries of a column using a carefully chosen orthogonal transformation, building R column by column from left to right.
