============
Installation
============

`effmass` is a Python 3 package and requires key packages from the `SciPy ecosystem <https://www.scipy.org/about.html>`_: SciPy, NumPy and Matplotlib. If you have not installed these packages before, it may be best to install them using your preferred package manager (eg: Homebrew). Note that together they will use >100MB of disk space. `effmass` can then be installed using the Python package manager `pip`:

.. code-block:: sh

   pip install effmass

To start the command line interface simply type

.. code-block:: sh

    effmass

Alternative Installation Methods
================================

If you use conda/anaconda, the safest thing to do is to create a new environment and then install `effmass` and all of its dependencies (which will use >300MB of disk space):

.. code-block:: sh

   conda create -n effmass python
   conda activate effmass
   pip install effmass

If you wish, you can install the very latest version of `effmass` from GitHub with the commands below. **Note**: The latest GitHub version may include more features and data format support that the latest release, but it is not a stable release, so may have more issues than usual. If you are unsure, use one of the above install methods instead.

.. code-block:: sh

   git clone https://github.com/lucydot/effmass.git
   # or git clone git@github.com:lucydot/effmass.git
   cd effmass
   pip install .
