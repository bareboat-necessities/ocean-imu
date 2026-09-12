"""Magnetic reference write-generation regressions."""
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_reference_events as X
from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG


class Tests(unittest.TestCase):
    def test_startup_write_creates_generation_zero(self):
        model=MAG.Model((22,0,43),(1,1,1))
        out=X.first_write(X.WriteWitness('startup_provisional',model,'one-producer'))
        self.assertEqual(out.generation,0); self.assertEqual(out.root,'one-producer'); self.assertEqual(out.model,model)
        for reason in ('live_refinement','continuous_hard_iron'):
            with self.assertRaisesRegex(ValueError,'first active'):
                X.first_write(X.WriteWitness(reason,model,'one-producer'))

    def test_refinement_and_continuous_updates_increment_same_root(self):
        sigma=(1,1,1); a=X.first_write(X.WriteWitness('startup_provisional',MAG.Model((22,0,43),sigma),'root'))
        b=X.rewrite(a,X.WriteWitness('live_refinement',MAG.Model((21,0,44),sigma),'root'))
        c=X.rewrite(b,X.WriteWitness('continuous_hard_iron',MAG.Model((20,0,45),sigma),'root'))
        self.assertEqual((a.generation,b.generation,c.generation),(0,1,2))
        self.assertEqual(c.root,a.root)

    def test_reference_write_cannot_restart_root_or_rewrite_Rmag(self):
        a=X.first_write(X.WriteWitness('startup_provisional',MAG.Model((22,0,43),(1,2,3)),'root'))
        with self.assertRaisesRegex(ValueError,'root cannot restart'):
            X.rewrite(a,X.WriteWitness('live_refinement',MAG.Model((21,0,44),(1,2,3)),'other'))
        with self.assertRaisesRegex(ValueError,'cannot change constructor Rmag'):
            X.rewrite(a,X.WriteWitness('live_refinement',MAG.Model((21,0,44),(1,2,4)),'root'))
        with self.assertRaisesRegex(ValueError,'cannot restart'):
            X.rewrite(a,X.WriteWitness('startup_provisional',a.model,'root'))

    def test_no_write_is_literal_identity(self):
        a=X.first_write(X.WriteWitness('startup_provisional',MAG.Model((22,0,43),(1,1,1)),'root'))
        self.assertIs(X.no_write(a),a)

    def test_readiness_keeps_numeric_producers_open(self):
        r=X.readiness()
        self.assertTrue(r['all_mag_reference_writes_use_named_setMagWorldRef_funnel'])
        self.assertTrue(r['reference_write_preserves_constructor_Rmag_sigma'])
        self.assertFalse(r['MagAutoTuner_startup_numeric_output_attached'])
        self.assertFalse(r['continuous_hard_iron_numeric_output_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
