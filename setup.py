from distutils.core import Extension, setup

swig_opts = ["-Wall", "-c++", "-outputtuple"]
libraries = ["bct", "gsl", "gslcblas"]

bct_gsl = Extension(
    "_bct_gsl",
    sources = ["bct_gsl.i"],
    swig_opts = swig_opts,
    libraries = libraries
)

bct_py = Extension(
    "_bct_py",
    sources = ["bct_py.i"],
    swig_opts = swig_opts,
    libraries = libraries
)

setup(
    name = "bct-cpp",
    version = "1.0",
    maintainer = "Steven Williams",
    maintainer_email = "stevencwilliams@gmail.com",
    url = "http://code.google.com/p/bct-cpp/",
    py_modules = ["bct_gsl", "bct_py"],
    ext_modules = [bct_gsl, bct_py]
)
