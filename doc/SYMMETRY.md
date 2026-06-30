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

If the initial position estimate in the design file is closer (in
parameter space) to the reflected configuration than to the true one,
the local optimiser (TRF) will converge to the mirror image. All
residuals look equally good, and the error is undetectable from the
calibration output alone — the cost function cannot distinguish the
true array from its reflection.

**This is not a hypothetical edge case.** Antenna arrays are sometimes
physically assembled as a mirror image of the nominal design — e.g. a
spiral arm wound the opposite way round during construction. The UNAM
deployment is a real example: the as-built array turned out to have the
*opposite* handedness to the nominal design file used as the initial
guess. Nothing in the radial or pairwise distance measurements, nor the
bearing constraint, can catch this, because all three are invariant
under reflection by construction (see above). The coordinate bounds
($\pm 4.2$ m default) only help when the reflected antenna position
happens to fall outside the box — for antennas near the symmetry axis
it usually doesn't, so it is not a guarantee either.

The consequence: **you cannot assume the design/initial-guess file has
the correct chirality.** It must be checked against the true, as-built
array (see "Determining the true chirality" below) before it is trusted.

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

3. Specify the sign $s$ of the as-built array via the
   `--chirality-sign` CLI argument (`positive` = $+1$, `negative` =
   $-1$).  The sign must come from an independent observation of the
   physical array; **it is not inferred from the `initial_guess`
   file.**  See "Determining the true chirality" below for how to
   determine the correct sign on site.
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

### Determining the true chirality

Choosing *which* antenna to use as $c$ is a separate question from
knowing *what sign it should produce*. The sign that matters is the
sign of the **actual, as-built** array — not whatever the nominal
design/initial-guess file says.

The distance measurements cannot tell you this: $r_i$ and $\delta_{ij}$
are exactly the same for the true array and its mirror image, and the
bearing constraint only pins down rotation, not reflection (see
"Symmetries of the distance data" above). So the true handedness has to
come from an **independent observation of the physical array**, for
example:

- Standing at the phase centre, facing the reference antenna $k$, and
  noting by eye whether antenna $c$ is physically to your left or your
  right.
- A site photo or sketch taken during installation that records which
  way the spiral arms actually wind, compared against the winding
  direction implied by the design file.

Once you know the true sign, pass it directly to the software via
`--chirality-sign`:

- If antenna $c$ is to the **left** of the reference bearing line (as seen
  standing at the phase centre facing $k$), use `--chirality-sign positive`.
- If antenna $c$ is to the **right**, use `--chirality-sign negative`.

The software will enforce exactly that sign regardless of what the
`initial_guess` file says.  You do **not** need to mirror or otherwise
manipulate the input file — the sign you specify is the sign the
optimiser respects.

If the design file agrees with reality, the chirality constraint adds
zero cost; if the design file disagrees, the constraint actively
steers the optimiser away from the (wrong) initial-guess chirality
toward the (correct) as-built chirality.

### Intuition: picking it without doing any maths

You don't need to compute a cross product by hand to choose a good
chirality antenna — a couple of mental pictures are enough.

1. **Draw the reference ray.** Imagine standing at the phase centre and
   facing antenna $k$ (the reference bearing direction). Extend that
   line all the way through the origin and out the other side. This
   line is the mirror axis: it's the one place a reflection could pivot
   around without anyone noticing from the distances alone.

2. **Pick an antenna clearly off to one side.** Any antenna $c$ that
   sits clearly to the *left* or *right* of that line — not hugging it
   — works as your $c$. But picking $c$ only chooses *which* antenna to
   use; it does not choose the *correct sign*. That sign must come from
   the actual, physically built array (see "Determining the true
   chirality" above) — the code will happily lock onto whatever
   handedness your initial-guess file has, even if it's the wrong one,
   so don't assume the design file is automatically correct.

3. **The further around the array, the better.** An antenna roughly
   90° away around the array from the reference antenna (e.g. on a
   spiral arm pointing in a very different direction) gives the
   strongest signal, because the penalty is normalised by $\approx
   \sin$ of the angle between $\mathbf{p}_k$ and $\mathbf{p}_c$ — that's
   largest near 90° and vanishes near 0°/180°. So as a rule of thumb:
   **pick an antenna on a different arm, as close to a right angle from
   the reference antenna as you can get.**

4. **Avoid the danger zone.** Antennas close to the reference ray
   itself, or close to its 180° opposite, are the ones to avoid — not
   because they're forbidden (only *exact* collinearity raises an
   error), but because a near-collinear choice gives a fragile,
   easily-flipped signal early in the optimisation, before the
   positions have settled.

In short: stand at the centre, look along the reference bearing, and
pick anything confidently off to one side and roughly perpendicular —
ideally on a different spiral arm. That gives you a robust *choice* of
$c$. It does not, by itself, tell you which sign is correct — for that
you still need to check against the actual built array, not the
design file.

### Usage

```bash
tart-position-cal \
    --measurements antenna_measurements.ods \
    --positions antenna_positions.json \
    --output calibrated_antenna_positions.json \
    --rot-index 3 \
    --rot-degrees 18.1125 \
    --chirality-index 10 \
    --chirality-sign positive
```

If the observed chirality sign matches this, the constraint adds zero
cost during optimisation; its only effect is to prevent the optimiser
from crossing into the mirror solution.  If the as-built array has the
opposite handedness (as happened at UNAM), pass `--chirality-sign
negative` instead — the software will enforce the correct, as-built
sign regardless of what the `initial\_guess` file says.

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
| Chirality (cross product sign) | `--chirality-index`, `--chirality-sign` | Reflection |

With all three constraints active, the global frame is fully determined.
