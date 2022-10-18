from setuptools import setup
from setuptools.extension import Extension
import numpy

try:
    from Cython.Build import cythonize
    cythonize("cython/nlfem.pyx", include_path=["include"])
except ModuleNotFoundError:
    print("\nWARNING: Cython was not found. Install Cython if you want that changes in cython/nlfem.pyx have any effect!\n")

# Project Name
name = "nlfem"
is_sequential = False

extra_link_args = ['-larmadillo', '-lgmp', '-lmpfr', '-lmetis', '-fopenmp']
extra_compile_args = ['-O3', '-DARMA_NO_DEBUG', '-fopenmp']

ext_modules = [
    Extension(
        name=name,
        sources=["cython/nlfem.cpp",
                 "src/Cassemble.cpp",
                 "./src/MeshTypes.cpp",
                 "./src/mathhelpers.cpp",
                 "./src/model.cpp",
                 "./src/integration.cpp"],
        extra_link_args=extra_link_args,
        extra_compile_args=extra_compile_args,
        language="c++",
        include_dirs=["include", "src", numpy.get_include()]
    )
]

setup(
    name=name,
    ext_modules=ext_modules,
    version="0.0.1",
    author="Manuel Klar, Christian Vollmann",
    author_email="klar@uni-trier.de, vollmann@uni-trier.de",
    description="This library provides a parallel assembly routine for a specific class of integral operators.",
    python_requires='>=3.6',
    include_package_data=True,
    install_requires=[['numpy'], ['scipy'], ['Cython']]
)
