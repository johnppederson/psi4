from addons import *

@ctest_labeler("nbody;gradient")
def test_nbody_vmfc_gradient():
    ctest_runner(__file__)

