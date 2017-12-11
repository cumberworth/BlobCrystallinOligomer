from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize

extensions = [
    Extension('*', ['*.pyx'],
            include_dirs=['../include'],
            extra_objects=['../lib/blobCrystallinOligomer.a'],
            language='c++',
            extra_compile_args=['-O3'])
]

setup(
    name = "blobalpha",
    ext_modules = cythonize(extensions, build_dir='build')
)
