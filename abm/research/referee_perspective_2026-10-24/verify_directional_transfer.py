"""Independent exact-arithmetic toy checks; no ABM imports or simulation."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path


def ratio(w, r, t):
    return sum(x*y*z for x, y, z in zip(w, r, t)) / sum(x*z for x, z in zip(w, t))


A = [F(1, 100), F(1)]
s0 = [F(1), F(1, 100)]
s1 = [F(1), F(1, 400)]
def coefficient(t, survival):
    return sum(a*x*s for a, x, s in zip(A, t, survival))

results = {}
for name, t, expected in [
    ('nominal', [F(1), F(1)], F(3, 8)),
    ('zero_opportunity_heavy', [F(9), F(1)], F(3, 40)),
    ('two_opportunity_heavy', [F(1), F(9)], F(27, 40)),
]:
    k0, k1 = coefficient(t, s0), coefficient(t, s1)
    prevention = 1-k1/k0
    assert prevention == expected
    # Choosing eps^2=b/K0 gives exact equal baselines in these disjoint bridges.
    target = F(1, 10**6)
    epsilon_squared = target/k0
    assert epsilon_squared*k0 == target
    assert epsilon_squared*max(t)**2 <= 1  # even put t entirely on one direction
    results[name] = {'K0': str(k0), 'K1': str(k1), 'F': str(prevention),
                     'calibrated_epsilon_squared': str(epsilon_squared)}

# Exact Bernoulli products are monomials: no common-p equality is used.
a, b, c, d = F(2), F(3), F(5), F(7)
eps = F(1, 100)
hand = (eps*a)*(eps*b)
environment = (eps*a)*(eps*c)*(eps*d)*(eps*b)
assert hand/eps**2 == a*b
assert environment/eps**4 == a*b*c*d
results['isolated_path_products'] = {'direct_order': 2, 'environment_order': 4,
                                     'direct_coefficient': str(a*b),
                                     'environment_coefficient': str(a*b*c*d)}

# Finite box: linear-fractional extrema attained at vertices; verify root signs.
w = [F(1, 5), F(3, 10), F(1, 2)]
r = [F(1, 10), F(1, 2), F(9, 10)]
m, M = F(1, 2), F(2)
vertices = list(product([m, M], repeat=3))
minimum = min(ratio(w, r, t) for t in vertices)
maximum = max(ratio(w, r, t) for t in vertices)
def root_function(x, minimize):
    return sum(wi*(ri-x)*((m if ri >= x else M) if minimize
                           else (M if ri >= x else m)) for wi, ri in zip(w, r))
assert root_function(minimum, True) == 0
assert root_function(maximum, False) == 0
assert all(minimum <= ratio(w, r, t) <= maximum for t in vertices)
results['box_vertex_check'] = {'vertices': len(vertices), 'minimum_ratio': str(minimum),
                               'maximum_ratio': str(maximum)}

# Bound sharpness: m=1, M=4, mu=2, w(low)=2/3 and w(high)=1/3.
w_sharp, t_sharp, r_sharp = [F(2, 3), F(1, 3)], [F(1), F(4)], [F(0), F(1)]
mean = sum(x*y for x, y in zip(w_sharp, t_sharp))
tv = sum(abs(x*y/mean-x) for x, y in zip(w_sharp, t_sharp))/2
shift = ratio(w_sharp, r_sharp, t_sharp)-ratio(w_sharp, r_sharp, [F(1), F(1)])
assert tv == shift == F(1, 3) == (F(2)-1)/(F(2)+1)
results['TV_sharpness'] = {'m': 1, 'M': 4, 'TV': str(tv), 'ratio_shift': str(shift)}

receipt = {'status': 'PASS', 'method': 'Python standard-library exact Fraction arithmetic',
           'scope': 'Independent toy algebra only; no ABM simulation or clinical validation',
           'checks': results}
output = Path(__file__).with_name('directional_transfer_verification.json')
output.write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
