from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
    name = "cpp_test_sum",
    ext_modules = cythonize(
        Extension("cpp_test_sum",
                  sources=["cpp_test_sum.pyx", "cpp_test_sum_.cpp"],
                  extra_compile_args=["-O3"],
                  language="c++",
                 )
    ),
    cmdclass = {'build_ext': build_ext},
)