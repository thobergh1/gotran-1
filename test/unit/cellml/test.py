__author__ = "Johan Hake (hake.dev@gmail.com)"
__date__ = "2012-05-07 -- 2012-09-11"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU LGPL Version 3.0 or later"

import unittest
from gotran2.input.cellml import *

class CellMLTester(unittest.TestCase):

    def test_winslow_1999(self):
        ode = cellml2ode("winslow_rice_jafri_marban_ororke_1999.cellml")
        self.assertEqual(ode.num_states, 33)

    def test_terkildsen_2008(self):
        ode = cellml2ode("terkildsen_niederer_crampin_hunter_smith_2008.cellml",\
                         ["FVRT", "FVRT_Ca", "Ca_b"])
        self.assertEqual(ode.num_states, 22)

    def test_hinch_2004(self):
        ode = cellml2ode("Hinch_et_al_2004.cellml")
        self.assertEqual(ode.num_states, 6)

    def test_pandit_2001(self):
        ode = cellml2ode("Pandit_et_al_2001_endo.cellml")
        self.assertEqual(ode.num_states, 26)

    def test_pandit_niederer(self):
        ode = cellml2ode("Pandit_Hinch_Niederer.cellml")
        self.assertEqual(ode.num_states, 6)

    def test_niederer_2006(self):
        ode = cellml2ode("niederer_hunter_smith_2006.cellml", ["Ca_b", "Ca_i"])
        self.assertEqual(ode.num_states, 5)

    def test_iyer_2004(self):
        ode = cellml2ode("iyer_mazhari_winslow_2004.cellml", ["RT_over_F"])
        self.assertEqual(ode.num_states, 67)

    def test_shannon_2004(self):
        ode = cellml2ode("shannon_wang_puglisi_weber_bers_2004_b.cellml")
        self.assertEqual(ode.num_states, 45)

    def test_niederer_et_al(self):
        ode = cellml2ode("Niederer_et_al_2006.cellml", ["Ca_b", "Ca_i"])
        self.assertEqual(ode.num_states, 5)
        
if __name__ == "__main__":
    unittest.main()
