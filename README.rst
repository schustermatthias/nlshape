This is the code for obtaining the results of the paper *Shape optimization for interface identification in nonlocal models* 
by M. Schuster, C. Vollmann and V. Schulz that can be found on https://arxiv.org/abs/1909.08884.

Build and Install on Ubuntu
===========================
In order to clone the project do
::
  git clone https://github.com/schustermatthias/nlshape.git path/to/local_folder

| Since this code contains a customized version of **nlfem** the following **basic requirements** of nlfem are needed
| ``gcc, g++, python3-dev, python3-venv, libgmp-dev, libcgal-dev, metis, libmetis-dev, libarmadillo-dev``.
On Ubuntu this can be done via
::
  sudo apt-get install git gcc g++ libarmadillo-dev liblapack-dev libmetis-dev
  sudo apt-get install python3-venv python3-dev libgmp-dev libcgal-dev

| See https://gitlab.uni-trier.de/pde-opt/nonlocal-models/nlfem for more information.
| Moreover to run nlshape **legacy FEniCS(version 2019.1.0)** is required. In order to use FEniCS in a virtual environment, it may has to be 
installed globally and then inherited as a global site package. A virtual environment can be built and activated via
::
  mkdir venv
  python3 -m venv venv/
  source venv/bin/activate

Additionally the packages from the file **requirements.txt** are neccessary and can be installed by
::
  (venv) python3 -m pip install -r requirements.txt

The creation of the virtual environment and the installation of packages from requirements.txt can probably also be done via your IDE.
Finally, nlfem can be installed by
::
  (venv) python3 setup.py build --force install
  
Running the Examples from the Paper
===================================
The project contains four configuration files ``ex1_sym.py, ex2_sym.py, ex1_nonsym.py`` and ``ex2_nonsym.py`` that were used to obtain the results in Section 5 of the paper.
In order to test, e.g., ex1_sym.py we just need to change Line 3 of nlshape/main.py to 
::
  from ex1_sym.py import configuration
and run main.py.
  
Raw Data
========
The raw data which was the result of the tests in Section 5 and used to create the pictures in this section of the paper can be found in **nlshape/results/results_paper/**.

License
=======
nlshape is published under MIT license. Copyright (c) 2022 Matthias Schuster, Christian Vollmann

| Parts of the project are taken from **nlfem** and have been customized.
| nlfem is published under MIT license. Copyright (c) 2021 Manuel Klar, Christian Vollmann
  
