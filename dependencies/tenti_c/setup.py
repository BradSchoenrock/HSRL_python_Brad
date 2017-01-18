from setuptools import setup, find_packages, Extension



classifiers = ""
version = '0.1.6'

module1 = Extension('_brill_c',
                   sources = ['brill_s6_py.i', 'brill_s6.c'],
                   include_package_data={'':['brill_s6_py.h']})

setup( name="tenti_s6",
       version=version, 
       py_modules = ['brill_c', 'tenti_s6'],
       include_package_data=True,
      zip_safe=False,
	  use_2to3 = True,
      description="Tenti S6 Python-swigged C code",
        ext_modules=[ module1],)
                                             
