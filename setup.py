import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, 'README.txt')).read()
CHANGES = open(os.path.join(here, 'CHANGES.txt')).read()

requires = [
    'scipy',
    'numpy',
    'matplotlib',
    'netCDF4',
    'ply',
    'bottleneck>=0.6.0',
    'dplkit>=0.3.0',
    'tenti_s6',
    ]

setup(name='hsrl',
      version='0.1',
      description='Data processing libraries for Lidar instruments',
      long_description=README + '\n\n' +  CHANGES,
      classifiers=[
        "Programming Language :: Python",
        ],
      author='E.Eloranta',
      author_email='ed.eloranta@ssec.wisc.edu',
      url='',
      keywords='lidar',
      packages=find_packages(),
      package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['config/*.json', 'config/*.cdl','config/*.txt'],
        'config': ['*.json', '*.cdl','*.txt'],
        'hsrl_config': ['*.json', '*.cdl', 'calvals_*.txt'],#this will go away
      },
      include_package_data=True,
      zip_safe=False,
      install_requires=requires,
      tests_require=requires,
      use_2to3 = True,
#      entry_points = """\
#      [paste.app_factory]
#      main = picnic:main
#      """,
      dependency_links = ['http://larch.ssec.wisc.edu/cgi-bin/repos.cgi']
      )

