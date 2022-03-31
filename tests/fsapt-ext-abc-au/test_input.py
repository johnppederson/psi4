from pathlib import Path
import psi4
from addons import *

@ctest_labeler("quick;sapt;cart;extern")
def test_fsapt_ext_abc_au():
    fsaptpy_installed = (Path(psi4.executable) / ".." / ".." / "share" / "psi4" / "fsapt" / "fsapt.py").resolve()

    ctest_runner(__file__, [
        fsaptpy_installed,
    ])

