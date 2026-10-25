import unittest
import numpy as np
import analyze_saved as a


class SavedAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.config = dict(bootstrap_draws=1000, bootstrap_seed=20261010, chunk_size=73,
                           stream_ordinals=dict(weekly=1,slopes=2,calibration=3),
                           warmup_days=90,observation_days=56,week_length_days=7)

    def test_synchronized_seed_weights_and_pooled_exposure(self):
        values=np.array([[[1,1],[2,1]],[[3,9],[6,9]],[[2,2],[4,2]]],float)
        point,samples=a.pooled_bootstrap(values,self.config,'slopes')
        np.testing.assert_array_equal(samples[:,1,0],2*samples[:,0,0])
        np.testing.assert_array_equal(samples[:,1,1],samples[:,0,1])
        self.assertEqual(a.rate(point)[0],500.)
        self.assertNotEqual(a.rate(point)[0],np.mean(a.rate(values)[:,0]))
        _, repeat=a.pooled_bootstrap(values,self.config,'slopes')
        np.testing.assert_array_equal(samples,repeat)

    def test_nonpositive_difference_draws_not_discarded(self):
        y=np.array([[1,2,4],[1,0,4],[1,-1,4],[np.nan,2,4]])
        slopes=a.log_slope(y,[1,2,4]);self.assertAlmostEqual(slopes[0],1.)
        result=a.uncertainty(slopes)
        self.assertEqual(result['undefined_bootstrap_draws'],3)
        self.assertIsNone(result['interval95'])
        self.assertFalse(result['conditional_on_defined_draws'])

    def test_declared_linear_bracket_exact_and_fallback_rules(self):
        p=[.1,.2,.3]
        values=np.array([[1,5,9],[3,3,5],[2,2.5,2.9],[4,2,5],[np.nan,4,5]])
        chosen,valid,exact,bracket,finite=a.select(p,values)
        np.testing.assert_allclose(chosen,[.15,.1,.15,.2+1/30,.15])
        np.testing.assert_array_equal(valid,[True,True,False,True,False])
        np.testing.assert_array_equal(exact,[False,True,False,False,False])
        self.assertEqual((~bracket).sum(),3)
        self.assertFalse(finite[-1])

    def test_weekly_conservation_cross_phase_and_final_partial(self):
        def blocks(days):
            return [dict(day_start=i,day_end=i+1,days=1,acquired_colonizations=i%2,patient_days=30) for i in range(days)]
        warm,obs=blocks(90),blocks(56)
        def summary(items):return {k:sum(b[k] for b in items) for k in ('acquired_colonizations','patient_days')}
        row=dict(warmup_blocks=warm,warmup=summary(warm),arms=dict(baseline=dict(observation_blocks=obs,observation=summary(obs))))
        values=a.daily_values(row,self.config)
        self.assertEqual(values.shape,(21,2));self.assertEqual(values[-1,1],180)
        self.assertEqual(values[12,1],210)
        np.testing.assert_array_equal(values.sum(axis=0),[73,4380])
        obs[-1]['day_end']=57
        with self.assertRaises(ValueError):a.daily_values(row,self.config)


if __name__=='__main__':unittest.main()
