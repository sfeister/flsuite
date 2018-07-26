# flsuite (Scott's FLASH Suite of Python Tools)

The **flsuite** is a set of Python tools related to FLASH simulations. They may be useful for certain stages of creating and analyzing the simulations.

The primary use case has been personal, but I'm sharing these tools in case others find them useful.

## Setup

### Dependencies
This module requires **Python 3.6+**. Installation requires **git**.

**OS X users:** Prior to installing dependencies, ensure an adequate Python installation by following [this guide](https://matplotlib.org/faq/installing_faq.html#osx-notes). The Python that ships with OS X may not work well with some required dependencies.

* [`numpy`](http://www.numpy.org/)
* [`scipy`](https://www.scipy.org/)
* [`matplotlib`](https://matplotlib.org/)
* [`h5py`](https://www.h5py.org/)
* [`yt`](https://yt-project.org/)

The dependencies may be installed according to the directions on 
their webpages, or with any Python
package manager that supports them. For example, one could use `pip` to install
them as
 ```bash
pip install numpy scipy matplotlib h5py yt
```
One could also use [Anaconda Python](https://anaconda.org/anaconda/python) to
install them as
```bash
conda install numpy scipy matplotlib h5py
conda install -c conda-forge yt
```

### Installation
After installing the required packages, we may install **flsuite**.

One way to install **flsuite** is via
```bash
pip install git+https://github.com/phyzicist/flsuite.git
```

Another way to install **flsuite** is via
```bash
git clone https://github.com/phyzicist/flsuite.git
cd flsuite/
python setup.py install
```

## Updating/Uninstalling
To update **flsuite** at a later date
```bash
pip install --upgrade git+https://github.com/phyzicist/flsuite.git
```

To uninstall **flsuite**
```shell
pip uninstall flsuite
```