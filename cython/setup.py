from Cython.Distutils import build_ext
from distutils.core import setup
from distutils.extension import Extension
import os

# Environment variables
home = os.getenv("HOME")
os.environ['CC'] = 'gcc'
os.environ['CXX'] = 'g++'

# Project Name
name = "assemble"

ext_modules = [
    Extension(
        name=name,
        sources=["cython/assemble.pyx"],
        extra_link_args=['-fopenmp', '-llapack', '-lblas'],
        language="c++",
        include_dirs=["include", "cython"],
        library_dirs=[home+"/lib"],
        libraries=["Cassemble"]
    )
]

# include_dirs provides the path to the header (Cassemble.h) and source (Cassemble.cpp) files
# this is necessary if we want to "inline" our code into the Cython Code

# The above case only needs the libCassemble.a or libCassemble.so and the corresponding header. See also
# https://cython.readthedocs.io/en/latest/src/tutorial/clibraries.html

setup(
    name=name, ext_modules=ext_modules, cmdclass={"build_ext": build_ext}
)

