from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
    name = "tradingrrl",
    ext_modules = cythonize(
        Extension("tradingrrl",
                  sources=["tradingrrl.pyx", "tradingrrl_.cpp"],
                  extra_compile_args=["-O3"],
                  language="c++",
                 )
    ),
    cmdclass = {'build_ext': build_ext},
)