%module brill_c

%include "carrays.i"
%array_class(double, doubleArray);

%{
#include "brill_s6_py.h"
%}

%include "brill_s6_py.h"
