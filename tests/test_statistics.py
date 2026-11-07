import itertools
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'))
import numpy as np
from run_experiments import exact_scan, planted_matrix, recovery_spectral, empirical_cutoff, wilson

def test_exact():
    rng=np.random.default_rng(42)
    for n,k in ((5,2),(6,3)):
        for _ in range(8):
            X=rng.normal(size=(n,n))
            ours,rows,cols=exact_scan(X,k,return_indices=True)
            best=-float('inf')
            for a in itertools.combinations(range(n),k):
                for b in itertools.combinations(range(n),k):
                    val=X[np.ix_(a,b)].sum()/k
                    if val>best: best=val
            assert abs(ours-best)<1e-12
            assert abs(X[np.ix_(rows,cols)].sum()/k-best)<1e-12
    print('Exact scan matches brute-force enumeration and returns maximizing index sets.')

def test_spectral_easy():
    rng=np.random.default_rng(188)
    for n,k in ((16,4),(64,16)):
        X,S,T=planted_matrix(n,k,20.,rng)
        a,b=recovery_spectral(X,k)
        assert np.array_equal(S,a) and np.array_equal(T,b)
    print('Spectral recovers easy planted supports exactly.')

def test_calibration():
    t=list(range(99))
    assert empirical_cutoff(t,0.05)==94
    lo,hi=wilson(0,400)
    assert lo==0 and 0<hi<.02
    print('Calibration and Wilson endpoint checks passed.')

if __name__=='__main__':
    test_exact();test_spectral_easy();test_calibration()