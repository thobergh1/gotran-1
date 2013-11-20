"""test for odeobjects module"""

import unittest

from modelparameters.logger import suppress_logging
from modelparameters.codegeneration import sympycode
from modelparameters.sympytools import symbols_from_expr

from gotran.common import GotranException

from gotran.model.odeobjects2 import *
from gotran.model.odeobjects2 import ODEObject, ODEValueObject
from gotran.model.utils import ode_primitives

from sympy import Symbol, Derivative


suppress_logging()

class TestODEObject(unittest.TestCase):

    def test_odeobjects(self):
        with self.assertRaises(TypeError) as cm:
            ODEObject(45)
        self.assertEqual(str(cm.exception), "expected a 'str' (got '45' which "\
                         "is 'int') as the first argument while instantiating "\
                         "'ODEObject'")

        with self.assertRaises(GotranException) as cm:
            ODEObject("_jada")
        self.assertEqual(str(cm.exception), "No ODEObject names can start "\
                         "with an underscore: '_jada'")
        
        obj0 = ODEObject("jada bada")
        self.assertEqual(str(obj0), "jada_bada")
        
        obj1 = ODEObject("jada_bada")

        self.assertEqual(obj0, obj1)
        
        obj0.rename("bada jada")
        self.assertEqual(str(obj0), "bada_jada")

        self.assertNotEqual(obj0, obj1)

    def test_odevalueobjects(self):
        with self.assertRaises(TypeError) as cm:
            ODEValueObject("jada", "bada")
        
        with self.assertRaises(GotranException) as cm:
            ODEValueObject("_jada", 45)
        self.assertEqual(str(cm.exception), "No ODEObject names can start "\
                         "with an underscore: '_jada'")
        
        obj = ODEValueObject("bada", 45)
        
        self.assertFalse(obj.is_field)
        self.assertEqual(Symbol(obj.name), obj.sym)
        self.assertEqual(45, obj.value)

    def test_state(self):
        t = Time("t")
        
        with self.assertRaises(TypeError) as cm:
            State("jada", "bada", t)
        
        with self.assertRaises(GotranException) as cm:
            State("_jada", 45, t)
        self.assertEqual(str(cm.exception), "No ODEObject names can start "\
                         "with an underscore: '_jada'")
        
        s = State("s", 45., t)
        self.assertFalse(s.is_field)
        s.toggle_field()
        self.assertTrue(s.is_field)
        
        a = State("a", 56., t)
        b = State("b", [40.], t)

        s_s = s.sym
        a_s = a.sym
        b_s = b.sym
        t_s = t.sym

        # Create expression from just states symbols
        self.assertEqual(ode_primitives(s_s**2*a_s + t_s*b_s*a_s, t_s), \
                         set([s_s, a_s, b_s, t_s]))

        # Create composite symbol
        sa_s = Symbol("sa")(s_s, a_s)
        
        self.assertEqual(symbols_from_expr(sa_s*a_s + t_s*b_s*a_s), \
                         set([sa_s, a_s, b_s, t_s]))

        # Test derivative
        self.assertEqual(s.derivative.sym, Symbol("s")(t_s).diff(t_s))
        #print sympycode(s.derivative.sym)
        #print sympycode(sa_s*a_s + t_s*b_s*a_s)
        #
        #print sympycode(Derivative(sa_s, s_s))
        #print sympycode((sa_s*a_s + t_s*b_s*a_s).diff(s_s))
        #print sympycode(Derivative(sa_s*a_s + t_s*b_s*a_s, s_s).doit())
        
    def test_param(self):
        
        with self.assertRaises(TypeError) as cm:
            Parameter("jada", "bada")
        
        with self.assertRaises(GotranException) as cm:
            Parameter("_jada", 45)
        self.assertEqual(str(cm.exception), "No ODEObject names can start "\
                         "with an underscore: '_jada'")
        
        s = Parameter("s", 45.)
        self.assertFalse(s.is_field)
        s.toggle_field()
        self.assertTrue(s.is_field)
        
        from gotran.model.utils import ode_primitives

        t = Time("t")
        a = State("a", 56., t)
        b = State("b", [40.], t)

        s_s = s.sym
        a_s = a.sym
        b_s = b.sym
        t_s = t.sym
        
        # Create expression from just states symbols
        self.assertEqual(ode_primitives(s_s**2*a_s + t_s*b_s*a_s, t_s), \
                         set([s_s, a_s, b_s, t_s]))

        # Create composite symbol
        sa_s = Symbol("sa")(s_s, a_s)
        
        self.assertEqual(symbols_from_expr(sa_s*a_s + t_s*b_s*a_s), \
                         set([sa_s, a_s, b_s, t_s]))


if __name__ == "__main__":
    unittest.main()