from distutils.core import setup, Extension
import numpy as np

module1 = Extension('tbcal', sources = ['cal.c', 'pyarray.c'])
module2 = Extension('tbfit', sources = ['fit.c', 'matrix.c', 'pyarray.c'])
# module3 = Extension('myanalyze', sources = ['analyze.c'])

setup(name = "TB_Module",
      version = "3.0",
      description = "intern_work",
      author = "taebum",
      author_email = "k.taebum@snu.ac.kr",
      include_dirs = [np.get_include()],
      ext_modules = [module1, module2]
      )

