"""CC0. Co-authored. EXPLORE diagnostic for Q9 (not pre-registered); needs runs/q9_half_mass_seed0.json first: is the +3 offset the linear-fit window?"""
import sys, math, json, statistics
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / 'experiments'):
    sys.path.insert(0, str(_p))
import q5_q9_sweeps as R
from toolkit import loghook
B = R.BAND
def exact(n, m=20000):
    f = lambda t: (1-t*t)**((n-3)/2)
    def simp(a,b):
        h=(b-a)/m; tot=f(a)+f(b)
        for i in range(1,m): tot += (4 if i%2 else 2)*f(a+i*h)
        return tot*h/3
    return simp(0,B)/simp(0,1-1e-12)
dims = R.DIMS
ex = [exact(d) for d in dims]
cross_all, slope = R.fit_cross(dims, ex)
local = [d for d in dims if 45 <= d <= 80]
cross_local, _ = R.fit_cross(local, [exact(d) for d in local])
# true crossing by bisection on integers then interpolation
dd = list(range(40,90)); pe = [exact(d) for d in dd]
true = next(d-1+(0.5-pe[i-1])/(pe[i]-pe[i-1]) for i,d in enumerate(dd) if i and pe[i-1]<0.5<=pe[i])
res = json.load(open(ROOT / 'runs' / 'q9_half_mass_seed0.json'))['outputs']['per_seed']
# empirical local fit per seed
loc = []
for r in res:
    ds=[d for d in dims if 45<=d<=80]; ps=[p for d,p in zip(dims,r['shares']) if 45<=d<=80]
    loc.append(R.fit_cross(ds,ps)[0])
out = {"exact_shares_on_grid": dict(zip(map(str,dims),ex)), "true_crossing_interpolated": true,
       "exact_curve_linear_fit_crossing_dims_30_100": cross_all, "exact_curve_linear_fit_crossing_dims_45_80": cross_local,
       "empirical_pooled_mean_fit_dims_45_80": statistics.fmean(loc), "empirical_sd_fit_dims_45_80": statistics.stdev(loc),
       "empirical_pooled_mean_fit_dims_30_100": statistics.fmean(r['fit_crossing'] for r in res),
       "mode": "EXPLORE", "status": "PROPOSED"}
print(loghook.log_run("q9_fit_diagnostic", 0, {"experiment":"q9_fit_diagnostic","mode":"EXPLORE","repo_commit":"3162c9950136a4a7f14e1e064f8d335db33968c8"}, out, out_dir=None))
for k,v in out.items():
    if k!="exact_shares_on_grid": print(k, v)
