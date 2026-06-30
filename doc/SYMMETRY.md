# Symmetries and Global Frame Determination

The position calibration problem estimates antenna coordinates
$\{\mathbf{p}_i = (x_i, y_i)\}$ from three kinds of data:

1. **Radial distances** $r_i = |\mathbf{p}_i|$ — distance from each antenna
   to the phase centre.
2. **Pairwise distances** $\delta_{ij} = |\mathbf{p}_i - \mathbf{p}_j|$ —
   distance between every pair of antennas.
3. **Reference bearing** $\theta_k$ — geographic bearing (degrees east of
   north) of one reference antenna $k$.

The phase centre is defined as the origin $(0,0)$.

## Symmetries of the distance data

Both radial and pairwise distances are invariant under rigid
transformations that fix the origin:

| Symmetry              | Type                | DOF     |
|-----------------------|---------------------|---------|
| Translation           | Continuous (2)      | $x,y$   |
| Rotation about origin | Continuous (1)      | $\phi$  |
| Reflection            | Discrete ($\mathbb{Z}_2$) | —   |

- **Translation** is resolved by fixing the phase centre at $(0,0)$. All
  radial distances are measured relative to this point, so the origin is
  not free to move.
- **Rotation** is resolved by the reference bearing $\theta_k$, which
  pins one antenna to a known geographic direction.
- **Reflection** is *not* resolved by the distance data or the single
  bearing constraint.

## The residual chirality ambiguity

Consider a reflection $R$ across the line through the origin at the
reference bearing $\theta_k$. Under $R$:

- The reference antenna $k$ lies on (or near) the axis of reflection,
  so its bearing is unchanged: $\theta_k' = \theta_k$.
- Every other antenna is reflected to the opposite side of the axis.
- All radial distances $r_i$ are preserved.
- All pairwise distances $\delta_{ij}$ are preserved.

The true configuration and its mirror image across the $\theta_k$ ray
produce **identical** values of the least-squares cost $L(x,y)$. They
are degenerate global minima — the measurement data alone cannot
distinguish them.

### Why this matters

If the initial position estimate in the design file were somehow closer
(in parameter space) to the reflected configuration than to the true
one, the local optimiser (L-BFGS-B or TRF) would converge to the mirror
image. All residuals would look equally good, and the error would be
undetectable from the calibration output alone.

In practice this is unlikely because:

1. The design positions are close to the as-built configuration and
   share its chirality.
2. The coordinate bounds ($\pm 4.2$ m default) typically exclude the
   reflected solution for antennas far from the symmetry axis.

But neither of these is a *guarantee*.

## The chirality constraint

To eliminate the reflection degeneracy definitively, the code supports a
**chirality constraint** via the `--chirality-index` CLI argument (or
the `chirality_index` parameter in the Python API).

### How it works

1. Pick a second antenna $c$ (the *chirality antenna*) that is **not**
   collinear with the reference antenna $k$. Any antenna on a different
   spiral arm works well.
2. Compute the cross product of the two position vectors from the
   *initial guess*:

   \[
   C = \mathbf{p}_k \times \mathbf{p}_c = x_k y_c - y_k x_c
   \]

3. Record the sign $s = \operatorname{sign}(C)$. This encodes the
   chirality of the initial configuration.
4. During optimisation, add a one-sided penalty that is zero when the
   sign matches, and quadratic when it flips:

   \[
   P(\mathbf{p}_k, \mathbf{p}_c) =
   \min\!\bigl(0,\; s \cdot (\mathbf{p}_k \times \mathbf{p}_c) / (r_k r_c)\bigr)
   \]

   The cross product is normalised by the measured radial distances so
   the penalty is dimensionless ($\approx \sin$ of the angle between the
   vectors).

5. Add $w \cdot P^2$ to the objective, where $w = 10^6$ by default.
   This acts as an effective hard constraint: the cost of flipping
   chirality is $\sim 10^{12}$ vs. a converged objective of $\sim 10^2$.

### Choosing the chirality antenna

The chirality antenna must not be collinear with the reference antenna
and the origin. If $\mathbf{p}_k \times \mathbf{p}_c \approx 0$, the
code raises an error because the sign is indeterminate.

A good choice is an antenna on a different spiral arm from the reference
antenna, ideally roughly perpendicular. For the UNAM array with
`--rot-index 3`, antenna 10 or 15 would be suitable (they lie on
different arms and are far from collinear with antenna 3).

### Usage

```bash
tart-position-cal \
    --measurements antenna_measurements.ods \
    --positions antenna_positions.json \
    --output calibrated_antenna_positions.json \
    --rot-index 3 \
    --rot-degrees 18.1125 \
    --chirality-index 10
```

If the initial position file has the correct chirality (which is always
the case for TART design files), this constraint is satisfied
automatically and adds zero cost. Its only effect is to prevent the
optimiser from crossing into the mirror solution.

### In IRLS mode

In IRLS mode (`--irls`), the chirality residual is excluded from the
per-iteration reweighting — it uses a fixed weight of 1.0 throughout.
Only the radial and pairwise distance residuals are subject to
Tukey/Huber downweighting. This ensures the chirality constraint is
never relaxed by the robust weighting scheme.

## Summary

| Constraint              | CLI flag            | Symmetry resolved |
|-------------------------|---------------------|-------------------|
| Phase centre at origin  | (implicit)          | Translation       |
| Reference bearing       | `--rot-index`, `--rot-degrees` | Rotation |
| Chirality (cross product sign) | `--chirality-index` | Reflection |

With all three constraints active, the global frame is fully determined.
