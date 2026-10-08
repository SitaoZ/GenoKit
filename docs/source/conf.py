"""Sphinx configuration: build the user guide without importing GenoKit."""
from pathlib import Path
import runpy
project = 'GenoKit'
author = 'Sitao Zhu'
copyright = '2026, Sitao Zhu'
version_file = Path(__file__).resolve().parents[2] / 'GenoKit/version.py'
release = runpy.run_path(str(version_file))['__version__'] if version_file.is_file() else '1.0.0'
version = release
language = 'en'
extensions = ['sphinx.ext.mathjax']
templates_path = ['_templates']
# Retain pre-existing apidoc stubs on disk but do not publish stale import-based API pages.
exclude_patterns = ['GenoKit/**', '_tools/**', '_downloads/**', '**/.DS_Store']
html_theme = 'alabaster'
html_static_path = ['_static']
html_title = f'GenoKit {release} documentation'
html_theme_options = {'description': 'Annotation, sequence extraction, design and visualization'}
