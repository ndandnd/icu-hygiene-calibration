# Exact full-ICU bridge coefficient estimator

This directory implements `K_h` from the clean-reservoir `p=0` background of the
frozen ICU engine. It does not estimate the coefficient by finite differencing
positive-p outcomes. The frozen model, instrumentation, study outputs, and
canonical mathematical notes are read-only. The short local validation output
is an implementation check, not an estimate for the 90+28-day study horizon.

## Algorithm and estimand

One seeded background retains admissions (with unique patient IDs), imported
carriage, carrier screening, isolation and relocation, staffing, contact
episodes, environment contacts, zone boundaries, role resets and all warmup
ticks. All hands/sites start clean and remain clean at p=0; imported patients
remain carriers. The source/recipient contact weight is the actual isolation
weight immediately before screening at that contact. Susceptible-recipient
reward is counted only after warmup. A source contact can precede observation.

For each worker and each requested hygiene schedule, keep a row-vector `ν` for
the unselected completion-state law and `b` for all surviving source-weight
masses. States are `absent`, `miss`, and `complete`. At every opportunity,
including opportunities when actual hands are clean, apply `ν ← νQ` and
`b ← bQD`, where `D=diag(1,1,1-e)`. The absent-state row of Q is `(0,1-h,h)`;
the other rows implement the frozen Markov completion law. At an imported
carrier contact inject `b ← b + c_i ν`; at an observation susceptible contact
add `d_j sum(b)` to K. An implicit role reset clears b while preserving ν.
The policy changes h after the final warmup tick and retains both ν and b.

Linearity sums every same-worker source/recipient bridge in O(events × arms)
time, with constant memory per worker/arm. It does not use an automatically
stationary q₂. Clean initial sites produce no environmental injection into
K; intervening site contacts still trigger their exact physical opportunities.
There is no direct patient-to-site shedding in the frozen engine.

Optional `histogram_cap=2` yields exact K contributions for 0, 1 and 2
opportunities plus an exact aggregate tail >=3. A fixed larger cap remains
linear time; `null` disables bins and is the fast default. Tail bins advance
through the same Markov survival transition as other bins, preserving their
total mass. Raw pair histograms count pairs without an intervening implicit
reset, regardless of whether their isolation weight or survival is zero;
they are not finite-p acquisition-route frequencies. Warmup-source K is
tracked separately. Completion-state laws at the endpoint are also recorded.

The default profiles `(rho,e)=(0,1),(.5,1),(0,.8)` and both h=.9/.95 arms are
evaluated together on one background. The actual p=0 engine runs with h0 and
the first profile, solely to generate its demographic/contact trajectory.
Under ideal fresh independent draws, that background law is common across
completion profiles; hygiene has no movement or patient feedback at p=0.
Legacy iid hygiene consumes the shared epidemiological RNG. Changing hygiene
implementations can therefore alter literal same-seed paths even at p=0.
This estimator uses distributional equality, not alternative-profile PRNG
replay or a positive-p eventwise coupling. Conditional survival is integrated
exactly; uncertainty is across independently initialized background seeds.

## Portable use and freezing

Run from the ABM root with the existing Python/Mesa runtime:

```sh
.venv/bin/python -m unittest discover -s research/full_icu_confirmation_2026-10-08/bridge -p 'test_bridge.py' -v
.venv/bin/python research/full_icu_confirmation_2026-10-08/bridge/run.py --config research/full_icu_confirmation_2026-10-08/bridge/tiny_validation_config.json --output research/full_icu_confirmation_2026-10-08/bridge/tiny_validation --freeze-only
.venv/bin/python research/full_icu_confirmation_2026-10-08/bridge/run.py --config research/full_icu_confirmation_2026-10-08/bridge/tiny_validation_config.json --output research/full_icu_confirmation_2026-10-08/bridge/tiny_validation --wall-seconds 60
```

The CLI accepts `--abm-root` for portable copies. Scientific model files are
hash guarded. Config, source hashes and runtime versions bind the manifest;
an existing manifest must match exactly. Thus a reviewed local scientific
config can be reused on a cluster, but a local runtime-bound output directory
cannot be resumed under a different runtime/platform. Freeze on the target
runtime into a fresh output directory when needed. Completed seed files and
the final summary are immutable. Incomplete seeds restart from initialization.
The run lock permits one writer per output directory. Independent process
shards need separate seed configs/output directories; combine their seed rows
only after verifying common scientific config/source identities. No cluster
job or full-horizon study run was launched by this implementation subtask.

Outputs contain K0/K1, paired differences, the ratio of seed means and limiting
response `1-K1/K0`, seed-level Monte Carlo SE, paired delta-method SE of the
response, observation exposure, event counts, ordered-event SHA256 digests,
warmup-source contributions, optional opportunity mixtures, and a stationary
q₂ comparison explicitly labelled as comparison only. Reported normal/delta
intervals are asymptotic; four smoke-test seeds do not establish coverage.

## Interpretation of the fixed-horizon remainder

The theorem is `A_h(p)=p²K_h+O(p³)` at fixed warmup-plus-observation horizon.
The O constant depends on the horizon/contact count. With a deterministic
padded-slot bound L, clean-reservoir acquisition reward bounded by M, and
weights in [0,1], the proof supplies the conservative inequality, for L>=3,

```
|A_h(p)-p²K_h|
  <= p²K_h [1-(1-p)^(L-2)] + M Pr{Binomial(L,p)>=3}
  <= p³ [(L-2)K_h + M choose(L,3)].
```

This distinguishes the fixed-p asymptotic remainder from a uniform accuracy
claim over long horizons. `pL` controls the expected number of dominating
gates, and small `pL` is a conservative sufficient regime for few such gates.
It is not a necessary threshold for observable accuracy: many gates are
ineligible, or have no acquisition reward, and their effects can cancel in a
response. Nor does `pL` alone certify a small relative error if K is tiny;
the displayed error divided by `p²K` also depends on M/K. A large `pL` makes
this crude gate-count certificate uninformative; it does not prove failure of
the coefficient approximation. Empirical small-p convergence needs a separate
positive-p confirmation with uncertainty and a declared horizon.

From clean initial sites, an environment-last-pickup route requires patient
pickup, site deposition, site pickup and recipient transfer (four successes).
Repeated transfers by an original dirty hand and acquisition feedback can
already enter with three successes. A large finite-p environment-last-pickup
share warns that routes absent from the leading patient-only coefficient are
active. It is a descriptive proximate-origin share, not a quantitative lower
bound on count/rate/response error, not a causal decomposition, and not proof
that a stationary q₂ response is inaccurate by that share. Response effects
can agree or cancel even when route shares differ. Dirty initial sites would
add environmental two-success terms and are explicitly rejected here.

## Verification

Exact tiny checks enumerate completion and efficacy histories independently
of the vector recurrence, including initialized memory, policy transients,
warmup pickups, isolation weights, resets and multiple workers. Fixed cap-tail
and disabled-bin calculations conserve K. Two small transmission-gate checks
illustrate repeated hand transfer and the clean-environment four-success
path. A direct reset test covers nurse shift, doctor end-of-duty, Xray restart
and Xray final-route resets without resetting ν. Bounded 2700-tick same-seed
comparisons reproduce the unmodified p=0 patient's contact/state/exposure and
both RNG states across separate/shared/zone conventions. These tests validate
implementation, not practical asymptotic accuracy at positive p.
