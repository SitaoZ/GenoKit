Build and maintain the documentation
====================================

Local HTML build
----------------

From the repository root, use Python 3.11 or newer for this Sphinx requirement:

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   bash sphinx.sh

Open ``docs/build/html/index.html``. The script works from any current directory,
uses the repository-relative source path and treats Sphinx warnings as build
failures. Set ``PYTHON`` to choose an interpreter; an optional first argument
selects an alternate output directory.

.. code-block:: bash

   PYTHON=/path/to/python bash sphinx.sh /tmp/genokit-docs-html

Do not repeatedly run ``sphinx-quickstart`` on an existing documentation tree.
That command initializes a project; it is not the normal rebuild step.
The old script built several output locations; the current entry point builds
one HTML tree by default.

Regenerate the command reference
--------------------------------

.. code-block:: bash

   python docs/source/_tools/update_cli_reference.py
   bash sphinx.sh

The generator reads parser definitions from ``GenoKit/genokit.py``, stops before
argument dispatch and writes ``command_reference.rst``. It requires only the
Python standard library and does not run extraction. The generated reference
omits parser epilog examples that can lag behind implementation. The user guide
still needs review when behavior, coordinates or algorithms change.

Sphinx itself does not import GenoKit's analysis modules. Pre-existing
``source/GenoKit`` apidoc stubs remain on disk but are excluded from this user
guide; they are not a verified public Python API reference. This prevents
import side effects and stale autodoc pages from breaking a basic HTML build.

Read the Docs configuration
---------------------------

The following is an example root-level ``.readthedocs.yaml`` for this static
user guide. Adapt it to an existing project configuration rather than creating
a competing file. Publishing still requires committing/pushing the documentation
to the repository connected to Read the Docs.

.. code-block:: yaml

   version: 2
   build:
     os: ubuntu-24.04
     tools:
       python: "3.11"
   sphinx:
     configuration: docs/source/conf.py
     fail_on_warning: true
   python:
     install:
       - requirements: docs/requirements.txt

Only documentation dependencies need installation because the HTML build uses
committed RST. The CLI reference should be regenerated locally when commands
change and committed with the other documentation. Check hosted build logs if
the published pages do not match the selected repository branch/version.

See the `Read the Docs configuration reference
<https://docs.readthedocs.com/platform/stable/config-file/v2.html>`_ for the
hosting settings. Local build success does not mean a site has been published.
