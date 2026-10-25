# Saved-data worker heterogeneity robustness

Post hoc conditional finite-horizon leading-coefficient worker reweighting only; no Monte Carlo uncertainty, new backgrounds, resampling, or clinical validation.

Exact fractional optimization treats the saved decimal worker coefficients as exact inputs. Numerical precision of the original coefficient engine is not bounded by this calculation. All seed worker sums agree with recorded K to the tolerance reported in results.json.

| Workflow | Profile | Nominal F (%) | t in [.5,2] F (%) | t in [.25,4] F (%) |
|---|---|---:|---:|---:|
| original | independent_perfect | 75.1319 | [75.0440, 75.1699] | [74.8586, 75.1813] |
| original | clustered_perfect | 53.6472 | [53.6040, 53.7520] | [53.5873, 53.8830] |
| original | independent_imperfect | 27.9939 | [27.9732, 28.0312] | [27.9621, 28.0818] |
| original | combined | 37.0118 | [36.8864, 37.3137] | [36.8372, 37.6887] |
| alternative | independent_perfect | 75.1831 | [75.1820, 75.1847] | [75.1795, 75.1861] |
| alternative | clustered_perfect | 53.6515 | [53.6052, 53.7604] | [53.5848, 53.8872] |
| alternative | independent_imperfect | 27.9861 | [27.9654, 28.0272] | [27.9538, 28.0785] |
| alternative | combined | 37.0026 | [36.8700, 37.3114] | [36.8064, 37.6728] |

## Conservative shared-multiplier contrasts

P−Q denotes prevention difference in percentage points. Each mechanism is optimized separately, so the interval contains every shared-t contrast but need not be attainable as a joint interval. A zero-containing interval is unresolved, not evidence of reversal.

| Workflow | P | Q | Bounds | Contrast interval (pp) | Guaranteed sign |
|---|---|---|---|---:|---|
| original | independent_perfect | clustered_perfect | [0.5, 2.0] | [21.2920, 21.5659] | positive |
| original | independent_perfect | clustered_perfect | [0.25, 4.0] | [20.9756, 21.5940] | positive |
| original | independent_perfect | independent_imperfect | [0.5, 2.0] | [47.0128, 47.1966] | positive |
| original | independent_perfect | independent_imperfect | [0.25, 4.0] | [46.7768, 47.2192] | positive |
| original | independent_perfect | combined | [0.5, 2.0] | [37.7303, 38.2834] | positive |
| original | independent_perfect | combined | [0.25, 4.0] | [37.1699, 38.3441] | positive |
| original | clustered_perfect | independent_imperfect | [0.5, 2.0] | [25.5728, 25.7787] | positive |
| original | clustered_perfect | independent_imperfect | [0.25, 4.0] | [25.5055, 25.9209] | positive |
| original | clustered_perfect | combined | [0.5, 2.0] | [16.2903, 16.8655] | positive |
| original | clustered_perfect | combined | [0.25, 4.0] | [15.8986, 17.0458] | positive |
| original | independent_imperfect | combined | [0.5, 2.0] | [-9.3404, -8.8552] | negative |
| original | independent_imperfect | combined | [0.25, 4.0] | [-9.7266, -8.7554] | negative |
| alternative | independent_perfect | clustered_perfect | [0.5, 2.0] | [21.4216, 21.5795] | positive |
| alternative | independent_perfect | clustered_perfect | [0.25, 4.0] | [21.2922, 21.6013] | positive |
| alternative | independent_perfect | independent_imperfect | [0.5, 2.0] | [47.1548, 47.2193] | positive |
| alternative | independent_perfect | independent_imperfect | [0.25, 4.0] | [47.1010, 47.2322] | positive |
| alternative | independent_perfect | combined | [0.5, 2.0] | [37.8706, 38.3147] | positive |
| alternative | independent_perfect | combined | [0.25, 4.0] | [37.5067, 38.3797] | positive |
| alternative | clustered_perfect | independent_imperfect | [0.5, 2.0] | [25.5780, 25.7951] | positive |
| alternative | clustered_perfect | independent_imperfect | [0.25, 4.0] | [25.5063, 25.9334] | positive |
| alternative | clustered_perfect | combined | [0.5, 2.0] | [16.2938, 16.8904] | positive |
| alternative | clustered_perfect | combined | [0.25, 4.0] | [15.9120, 17.0809] | positive |
| alternative | independent_imperfect | combined | [0.5, 2.0] | [-9.3461, -8.8428] | negative |
| alternative | independent_imperfect | combined | [0.25, 4.0] | [-9.7189, -8.7279] | negative |

## Mathematics and limits

For worker u, let A_u and B_u be pooled baseline and policy coefficients. F(t)=1−Σt_u B_u/Σt_u A_u. Sorting B_u/A_u and evaluating all threshold endpoint assignments yields the exact linear-fractional box extrema; the script verifies the root certificates with rational arithmetic. Multipliers apply to the same worker ID across all 256 seeds and both arms; their overall scale cancels. This is a realizable subclass of directional heterogeneity: assign a worker-specific pickup factor a_u and transfer factor b_u with a_u b_u=t_u. It does not vary a_u/b_u separately, because that directional ratio cancels at this order.

Workers with zero pooled baseline and policy coefficient contribute nothing and are omitted. A positive policy coefficient with zero baseline is rejected. The saved worker labels are used literally; no role interpretation is required. The boxes are mathematical stress ranges with no empirical probability interpretation. Their extrema do not incorporate background Monte Carlo uncertainty, parameter estimation, finite-transmission feedback, environmental routes, endpoint-specific heterogeneity, or arbitrary variation within a worker.

This is a post hoc analysis under the saved PLAN.md. It establishes conditional robustness only over the specified worker-multiplier class on these saved backgrounds. The leading-limit invariance to global direction multipliers remains an algebraic fact; these nontrivial bounds concern heterogeneous products. Positive or negative contrast certificates use conservative per-profile interval differences. No exact shared-t contrast optimization is claimed.

Input hashes include both manifests, all 512 analysis records, all 256 original archived reference records, source/config/seed files, this script and the saved plan. The original extension is linked per seed to archived source data using archived file hash, sealed payload, event hash and exact selected worker-dictionary equality. The alternative workflow has different event histories by construction; it is linked to the same seedbook and frozen scientific/source/configuration identities, not falsely asserted to have matching histories.
