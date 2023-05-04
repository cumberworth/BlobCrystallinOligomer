from distutils.core import setup
from distutils.extension import Extension

import numpy

from Cython.Build import cythonize

extensions = [
#    Extension('*', ['crystallinpy/*.pyx'],
    Extension('*', ['*.pyx'],
            include_dirs=['../include', numpy.get_include()],
            extra_objects=['../lib/blobCrystallinOligomer.a'],
            language='c++',
            extra_compile_args=['-O3'])
]

setup(
    name = "blobalpha",
#    ext_modules = cythonize(extensions, build_dir='build', gdb_debug=True)
    ext_modules = cythonize(extensions)
)
