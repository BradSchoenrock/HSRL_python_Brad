
This folder contains "easy_install"-able distributions of code maintained externally. In the future, these may be posted to a public location for easy retrieval.

dplkit-0.3.2.tar.gz   - DPL abstract classes for standardized data stream practices 
(under heavy development) - SSEC - ray.garcia@ssec.wisc.edu 

Most current version in
https://git.ucar.edu/eol_hsrl_dpl_tools.git/DplKit/python

tenti_s6-0.1.6.tar.gz - Tenti S6 spectral model. Uses standard C libraries, and should compile perfectly normal using pip or easy_install

The code for this is in dependencies/tenti_c (setup.py is used to construct the tarball), and is a derivative code base from an F2C-translated C source with all libf2c dependencies removed in favor of native C functions, libm, and complex values. The original fortran code is from the ahsrl_matlab codebase, and the main body fortran source is kept in tenti_c for proper documentation, and to facilitate any potential future evolution of the code.

virtual_radiosonde_source-0.1.tar.gz - Virtual Radiosonde Source, install using pip or easy_install

Virtual Radiosonde From NWP is designed to be an extension of current
gfs-based data providers. VRS processes the data to create a "virtual
radiosonde" with the properties of a radiosonde for use in place of
conventional radiosonde sources. GFS data is assumed to be provided by 
dreadnaught and if VRS is deployed on dreadnaught file downloads are 
not necessary. However, if VRS is deployed to an external machine, the 
needed grib files will be downloaded to a local cache for use offline.

Requires grib_api-devel openjpeg-devel libpng-devel RPMs
