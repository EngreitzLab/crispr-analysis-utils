"""The compiled Rust core ships and loads.

If this fails, the wheel was built without the extension or the editable
install did not run maturin.
"""

import crispr_analysis_utils
from crispr_analysis_utils import _core


def test_core_is_a_compiled_extension():
    assert _core.__file__.endswith((".so", ".pyd", ".dylib"))


def test_core_version_matches_the_package():
    assert _core.version() == crispr_analysis_utils.__version__
    assert _core.__version__ == crispr_analysis_utils.__version__
