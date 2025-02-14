import pytest
import numpy as np
from polysolve.polysolve import quadratic

def test_quadratic():
    """tests the roots are zero"""
    params = [3., 0., -1.]
    roots = quadratic(*params)
    assert all(np.isclose(np.polyval(params, root), 0.) for root in roots)

#def test_quadratic_fails():
#    """Check bad quadratic raises error."""
#    with pytest.raises(RuntimeWarning, 
#                       match="invalid value encountered in scalar divide"):
#        # There are infinite roots on this flat line.
#        quadratic(1., 2., 3.)

