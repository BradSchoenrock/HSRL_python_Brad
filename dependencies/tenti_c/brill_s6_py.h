#ifdef SWIG
%module brill_c
#endif
void brill_s6(double wavelen, double temp, double pres, double *freqs, unsigned long  int numspec, double * spect);
