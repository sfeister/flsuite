# flsuite (Scott's FLASH Suite)

**flsuite** is a set of Python tools (and bash scripts) related to FLASH simulations. They may be useful for certain stages of creating and analyzing the simulations.

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

**NOTE**: If you are on a cluster where you do not have write permissions to the python installation directory, you may need to add "--user" to your pip and setup calls here and below. E.g.
```bash
pip install --user numpy scipy matplotlib h5py yt
```

As an alternate to pip, one could also use [Anaconda Python](https://anaconda.org/anaconda/python) to
install the dependencies
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

To update **flsuite** at a later date
```bash
pip install --upgrade git+https://github.com/phyzicist/flsuite.git
```

An alternative way to install **flsuite** is via
```bash
git clone https://github.com/phyzicist/flsuite.git
cd flsuite/
python setup.py install
```

If you installed with the "--user" flag and "~/.local/bin" is not on your system path, add it to your path (e.g. by appending the line "PATH=$PATH:$HOME/.local/bin" to "~/.bashrc"). Documentation can be read by calling "ffpngs -h" from a bash terminal.

## Command-line tools
A useful bash script called "ffpngs" will be added to /bin upon installation (or ~/.local/bin, if --user flag was used). It is a wrapper for ffmpeg that turns a directory of PNG files into a cross-platform .mp4 movie file. This is convenient when trying to share between Powerpoint, MacOS, Linux, etc.
Requires bash shell (i.e. works natively with Linux or MacOS) and ffmpeg.

## flsuite Python package
Useful subpackages include:
* flsuite.parLaser (includes a class "parLasers" to programmatically define multiple lasers for a flash.par file)
* flsuite.parIO (includes general methods for reading and writing from flash.par files)
* flsuite.flyt.flyt (includes the function "get_simdata", a tool for extracting uniform-res arrays from FLASH sims using Yt)
* flsuite.sftools (includes the function "getH5Outs" which gets a sorted list of FLASH-generated HDF5 plotfiles or checkpoint files in a directory, and the function "subdir" which is a painless way to make a subdirectory while checking if it already exists)

Limited documentation for these modules appears in the source code.
 
## Uninstalling

To uninstall **flsuite**
```shell
pip uninstall flsuite
```