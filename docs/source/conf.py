from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

project = "QDD"
copyright = "QDD contributors"
author = "QDD contributors"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "breathe",
]

myst_enable_extensions = ["colon_fence"]
html_theme = "sphinx_rtd_theme"

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"

breathe_projects = {
    "qdd": str(ROOT / "docs" / "build" / "doxygen" / "xml"),
}
breathe_default_project = "qdd"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
if not (ROOT / "docs" / "build" / "doxygen" / "xml" / "index.xml").exists():
    exclude_patterns.append("api/cpp.md")
    suppress_warnings = ["toc.not_readable"]
