# Prospective controlled-assay uncertainty benchmark

17 October 2026. This protocol is saved before drawing any assay outcomes. The
benchmark is local and separate from every ICU simulation and cluster study.
Analytical predictions are known when choosing this design; this is a method
demonstration with known truths, not an empirical test of a clinical hypothesis.

## Question and scope

How much can an independently measured two-opportunity survival probability
reduce uncertainty in hygiene response beyond one-opportunity survival, when
average completion is known? Evaluate finite-sample coverage and contraction of
a projected confidence set in a declared finite stationary Markov model family.
Do not convert arbitrary parameter ranges into priors or clinical probabilities.
Do not claim this quantifies full-ICU structural or interpolation uncertainty.

## Frozen design

- Baseline completion h0=0.90, known without error.
- Candidate family: rho=i/40, i=0,...,40; e=j/100, j=0,...,100, giving 4,141
  combinations. Bounds are sharp over this grid, not the continuous parameter
  rectangle. Retain rho=1 with its explicitly stationary, reducible-chain law.
- Four known truths, all on grid: independent (rho=0,e=1), clustered (.5,1),
  imperfect (0,.8), combined (.5,.8). These are chosen illustrative cases.
- Separate controlled assays for k=1 and k=2. Each replicate initializes
  contamination independently of a newly stationary completion process; clearing
  trials are independent conditional on completion. Record survival as Bernoulli.
  The two assay banks are independent. Draw their success counts directly as
  independent Binomial(n,q1) and Binomial(n,q2), which is distributionally identical
  to explicitly generating these independently reset Bernoulli experiments.
- n=500, 2,000 and 8,000 replicates **per assay**, with 1,000 independent repeated
  experiments per truth/sample-size combination. Total 12,000 repeated experiments;
  no data-dependent enlargement. These are simulation diagnostics, not observations.
- NumPy PCG64 generator, SeedSequence([2026101703, truth_index, n]), with truth
  indices in the order above. Draw the 1,000 q1 counts and then 1,000 q2 counts.
- Prediction targets held out from assay fitting: h1=.925,.95,.975, each with
  k=2 and k=6. Predict F=1-qk(h1)/qk(h0), holding rho/e fixed under intervention.
  No assay at a policy h is used. There is no estimate of p or a full-ICU rate.

The finite family includes exact e=0, for which all qk=1 and F=0. All baseline
qk are strictly positive for the declared finite k and h0<1.

## Confidence-set construction and estimands

For each assay use a two-sided exact Clopper–Pearson interval with marginal
coverage at least 97.5%, hence tail probability .0125. The joint two-assay set has
coverage at least 95% by the union bound. This statement does not require the
assays to be independent, although this design makes them so.

The first set retains candidate parameters whose q1 lies in its interval. The
second retains candidates satisfying both q1 and q2 intervals, using the same q1
interval so that sets are nested in every experiment. For each target, enumerate
the retained grid and report its minimum and maximum F. These are exact extrema
of the finite candidate set, to floating-point tolerance. Do not assume extrema
occur at rectangle corners. Add a numerical tolerance of 2e-12 when testing set
membership; this can only widen the set at the comparison boundary. Record this
tolerance and independently check every reported minimum/maximum against retained
candidates. If a set is empty, retain the event and report undefined bounds;
count its coverage as failure. Never substitute a nearest model.

Because each target is a deterministic function of the same parameters, inclusion
of the true grid point guarantees **simultaneous** coverage of all six targets.
No extra target multiplicity correction is needed for this projection statement.
The guarantee is conditional on the true model being on the grid and all model
and measurement assumptions. It is not a continuous-family confidence claim.

For every truth and sample size report:

1. Empty-set counts; parameter inclusion and simultaneous six-target coverage.
2. Per-target coverage, median projected interval width and 10th/90th width
   quantiles for each information set. Undefined widths remain separately counted.
3. Median width reduction and fraction with strict contraction, among experiments
   where both intervals exist. Check nestedness for every valid experiment.
4. Exact 95% binomial intervals for empirical repeated-experiment coverage. These
   describe the benchmark's Monte Carlo precision, not clinical uncertainty;
   there is no pass/fail rule based on a realized coverage estimate.

Save the raw assay counts, per-experiment confidence limits, candidate counts,
prediction extrema and coverage indicators in a compressed numerical archive.
Bind protocol, script, archive and result hashes and actual runtime versions in
a receipt. Results are never overwritten without first preserving an old run.

## Prespecified model-class counterexample

Report a deterministic non-Markov alternative with the same h, q1 and q2 as the
combined truth at each h. Let rho=.5 and e=.8. At each arm draw an episode-level
completion propensity H equal to (1-rho)h with probability
w=(1-h)/(1-h+rho*h), and 1 otherwise; conditional on H, completions are independent.
The episode retains H across all its opportunities. Then E[H]=h and
Var(H)=rho*h*(1-h), so its q1 and q2 equal the Markov model's exactly. Its longer
survival is w*(1-e*(1-rho)*h)^k+(1-w)*(1-e)^k, generally different for k>=3.
This is a constructed counterexample, known analytically before execution.

Reuse the combined truth's assay counts (their distribution is identical) to
report how often the Markov-projected intervals include these alternative k=6
targets. This has **no nominal coverage guarantee**; it illustrates dependence on
the model class. Do not call it an estimated clinical misspecification rate.
Its intervention law is explicitly arm-specific through h and w; it is not a
claim about how a real intervention changes worker heterogeneity.

## Interpretation boundaries

- Natural contamination episodes may select completion states. The controlled
  independently initialized assay assumption is necessary for these q formulas.
- h is treated as known; error in compliance measurement is not included.
- rho/e are constant across arms and clearing is homogeneous. Heterogeneity,
  higher-order dependence and a change in rho/e under policy are excluded from
  the primary family; the counterexample demonstrates one consequence.
- A smaller projected interval is useful information within the family, not
  biological validation. Discretization can omit the true parameter; these
  intervals cannot be promoted to continuous-family coverage by calling the
  grid dense. No off-grid truth is used to claim such coverage.
- This benchmark estimates hand-survival response, not the full ICU's finite-p
  policy response. Existing Monte Carlo, calibration and structural uncertainties
  must remain separate in any manuscript use.
- For the separate robustness study, bootstrap uncertainty in baseline-selected
  p is reportable. Policy intervals at a single simulated p must not pretend to
  propagate that uncertainty across unsimulated p values. A response model or
  new policy runs would be needed, with separately declared assumptions/design.

## Allowed execution repairs

Fix coding/runtime defects preserving all cases, counts, seed construction,
confidence level, projection rule and targets. Record any such correction before
rerunning. Do not tune the family, samples or reported targets to contraction or
coverage outcomes. No cluster access, simulator or manuscript changes are part
of this protocol.
