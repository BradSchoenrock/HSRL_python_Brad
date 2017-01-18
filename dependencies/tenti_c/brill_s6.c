/* brill_s6.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include <stdio.h>

//#include "f2c.h"
#include "brill_s6_py.h"
#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>

/* Table of constant values */
#define integer long int
#define doublecomplex double complex
#define real double
#define ftnlen size_t
#define dabs fabs

static integer c__6 = 6;
static integer c__0 = 0;
static doublecomplex c_b22 = {.5,0.};
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__12 = 12;
static integer c__8 = 8;
static double c_b77 = 2.;
static real c_b78 = 1.f;
static real c_b85 = .5f;
static real c_b90 = -1.f;
static real c_b91 = 2.f;

#define cleararr(arr) for( i__1=0; i__1<sizeof(arr)/sizeof(arr[0]);i__1++) arr[i__1]=0

void w_(doublecomplex * ret_val, doublecomplex *z__);
void z_divxxxx(doublecomplex *c, doublecomplex *a, doublecomplex *b){
    *c=*a / *b;
}
void z_exp(doublecomplex *r, doublecomplex *z)
{
        double expx, zi = cimag(*z);

        expx = exp(creal(*z));
        *r = expx * cos(zi) + I*(expx * sin(zi));
}
void s_copy(char *d, const char *s, size_t dsize, size_t count){
    strncpy(d,s,dsize);
}
void z_div(doublecomplex *c, doublecomplex *a, doublecomplex *b)
{
        double ratio, den;
        double abr, abi;

        if( (abr = creal(*b)) < 0.)
                abr = - abr;
        if( (abi = cimag(*b)) < 0.)
                abi = - abi;
        if( abr <= abi )
                {
                if(abi == 0) {
                        if (cimag(*a) != 0 || creal(*a) != 0)
                                abi = 1.;
                        *c = (abi / abr) + I*(abi/abr);
                        return;
                        }
                ratio = creal(*b) / cimag(*b) ;
                den = cimag(*b) * (1 + ratio*ratio);
                *c = (creal(*a)*ratio + cimag(*a)) / den + \
                     I*(cimag(*a)*ratio - creal(*a)) / den;
                }

        else
                {
                ratio = cimag(*b) / creal(*b) ;
                den = creal(*b) * (1 + ratio*ratio);
                *c = (creal(*a) + cimag(*a)*ratio) / den + \
                     I*(cimag(*a) - creal(*a)*ratio) / den;
                }
}

/* Subroutine */ void brill_s6(double wavelen, double temp, 
	double pres, double *freqs, unsigned long int numspec, double *
	spect)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    ///* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    //double sqrt(double);

    /* Local variables */
    static integer i__;
    static real s, x, y, v0;
    static double kb;
    static real pi, cv;
    static double xp;
    static real axb, atm, qmag;
    static double spdc, mass;
    extern /* Subroutine */ int tent_(real *, real *, real *, real *, real *, 
	    real *);
    static char fname[32];
    static real etalam;
    static double snhalf;
    static real etazet;

/*      SUBROUTINE brill_s6(WAVENUM,TEMP,PRES,CENTSHIFT,DSIGMA, */
/*     & NUMSPEC,SPECT,SIZESP) */
/*     SUBROUTINE GETRB(GRID,NDWIDTH,RBVALUE) */
/*     THIS SUBROUTINE DETERMINES THE VALUE OF THE RB L@HAPE FUNC'NON WHICH */
/*     IS CALCULATED BY TENTI'S SUBROUTINE.  VALUES OF THE RB PROFELE ARE */
/*     TABULATED FOR FREQUENCIES FROM ZERO TO TOPLIM, IN INCREMENTS OF GRID */

/*     WAVELEN   = wavelength in meters ewe 10/28/03 */
/*     TEMP      = temperature in deg-Kelvin */
/*     PRES      = pressure in mb */
/*     FREQS     = frequencies in Hz */
/*       WIDTHRB = EFFECTIVE WIDTH OF THE RB PROFII-E (GHZ) */
/*       Y       = Y PARAMETER USED IN CALCULATING RB PROFILE */
/*       QMAG    = MAGNITUDE OF SCATTERING WAVEVECTOR OR K */
/*       V0 	= VELOCITY PARAMETER */

/*     I have attempted to express all quantities in MKS units */
/*     orignal code uses a mixture of cgs and MKS (ewe 10/29/03) */
/*      REAL*8 KB,DSIGMA,WAVENUM,SPDC,XP */
/*     new ewe 10/29/03 */
    /* fprintf(stderr, "Greetings from brill_s6 !!: wavelen = %g\n", wavelen); */
    /* fortran arrays  start with 1  */
    --spect;
    --freqs;

    /* Function Body */
    pi = M_PI;
    spdc = 2.998e8f;
    etalam = .198f;
    etazet = 1.407f;
    cv = 1.f;
    s_copy(fname, "S6_model_spect.out", (ftnlen)32, (ftnlen)18);
/*     Boltzman constant J/K */
/*     KB = 1.38062E-20    changed to MKS ewe 10/29/03 */
    kb = 1.38044e-23f;
/*      WAVE=1./WAVENUM     removed ewe 10/29/03 */
/*     molecule mass for dry air */
/*      MASS = 28.966/6.02297E23  changed from grams to kg ewe 10/29/03 */
    mass = 4.8092552345437547e-26f;

/* SNHALF = SIN(theta/2) Where theta is the scattering angle (180) */
    snhalf = 1.f;
/*                    ***** */
    axb = (temp + 110.4f) / (temp * temp);
/*      ATM=PRES/1013250  now uses press in mb ewe 10/29/03 */
    atm = pres / 1013.25f;
/*     Y is supposed to be the ratio of mean-free-path to wavelength */
    y = axb * .2308f * (atm / snhalf) * wavelen * 1e9f;
/*     notice that if wavelength is given in m this is just: */
/*     QMAG=4*pi/wavelength  with units of 1/m     ewe 10/29/03 */
/*     QMAG=1.2566E10*SNHALF/(WAVE/1.0E-9) new version follows ewe 10/29.03 */
    qmag = pi * 4 * snhalf / wavelen;
/*     V0 remains dimensionless after conversions units to MKS  ewe 10/29/03 */
    v0 = sqrt(kb * temp / mass);
    xp = pi * 2 * spdc / (sqrt(2.f) * qmag * v0);
/*        WRITE(*,*) */
/* 	WRITE(*,1111) Y,QMAG,V0,XP,NUMSPEC */
/* 1111	FORMAT('Y, QMAG, V0, XP, NUMSPEC',4E14.3,I6) */
/* 	WRITE(*,1111) FREQS(0),FREQS(1) */
/* 1111	FORMAT('freq0 freq1',4E14.3,I6) */

/*   Open file to write data */

/*      OPEN(UNIT=40,FILE=FNAME,STATUS='REPLACE') */
/* 109  FORMAT('% STEP (pm):  ',F10.3) */
/*      WRITE(40,109) DSIGMA*(WAVE**2)*1e12 */
/* 110  FORMAT('% TEMP (K):   ',F10.3) */
/*      WRITE(40,110) TEMP */
/* 111  FORMAT('% PRES (ATM): ',F10.3) */
/*      WRITE(40,111) ATM */
/* 112  FORMAT('% Y:          ',F10.3) */
/*      WRITE(40,112) Y */
/* 113  FORMAT('% K or QMAG   ',E14.3) */
/*      WRITE(40,113) QMAG */
/* 114  FORMAT('% V0:         ',F10.3) */
/*      WRITE(40,114) V0 */
/* 116  FORMAT('% SIZESP:    ',I5) */
/*      WRITE(40,116) NUMSPEC */
/* 117  FORMAT('% DX:       ',F10.3) */
/*      WRITE(40,117) XP*DSIGMA */

/* 221  FORMAT('%    X       SPECTRUM(I)     I') */
/* 222  FORMAT(F7.3,2X,F12.3,2X,I5) */
/*      WRITE(40,221) */

    s = 0.f;
    i__1 = numspec;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*         X=XP*(-CENTSHIFT+DSIGMA*(I-(NUMSPEC/2.0+1))) new version follows */
/*        ewe 10?29/03 */
	x = xp * freqs[i__] / spdc;
	r__1 = dabs(x);
	tent_(&r__1, &cv, &etalam, &etazet, &y, &s);
/*   FOLLOWING LINE ADDED BY J. FORKEY TO NORMALIZE S SO THAT INTEGRATING */
/*   ACROSS LINESHAPE IN GHZ WILL YIELD I */
	spect[i__] = pi * 2.f * 3e8f / sqrt(2.f) / (qmag * v0) * s;
/*         WRITE(40,222) X,S,I */
/* L300: */
    }
} /* brill_s6__ */

/* Computer Model Of Filtered Rayleigh Scattering Signals */
/* The following FORTRAN program, called FRSMODEL, evaluates */
/* equations 2.17, 2.18, */
/* 3.3, and 3.4, given the appropriate input information. */
/* The subroutine TENT, which is */
/* used to calculate the Rayleigh - Brillouin lineshape, */
/* was provided by Professor G. Tenti */
/* of the University of Waterloo.  Calculations in other subroutines */
/*  are based on expressions */
/* in chapters 2 and 3, and in appendices A and B of this dissertation. */
/* Subroutine */ int tent_(real *x, real *cv, real *etalam, real *etazet, 
	real *y, real *s)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Builtin functions */
    //double sqrt(double);

    /* Local variables */
    static integer i__;
    static real g2;
    static integer i6;
    static real s6, g10, g11, pi;
    static integer iy, nx, ny;
    static real g011, s6e, den, gp10, gp11;
    extern /* Subroutine */ int imn_(doublecomplex *, integer *, real *, real 
	    *);
    static real gbar;
    static doublecomplex intg[10];
    extern /* Subroutine */ int s6get_(real *, real *, real *, real *, real *,
	     real *, real *, doublecomplex *, integer *, real *, real *);
    static real dtbar;

    ny = 1;
    nx = 1;
    i6 = 1;

    pi = M_PI;

    cleararr(intg);
/* Computing 2nd power */
    r__1 = *cv / (*cv + 1.5f);
    g11 = .83333333333333337f - r__1 * r__1 * .55555555555555558f * *etazet;
/* Computing 2nd power */
    r__1 = *cv / (*cv + 1.5f);
    gp11 = sqrt(5.f / (*cv * 2.f)) * *etazet * (r__1 * r__1) / 3.f;
    gp11 = -gp11;
/* Computing 2nd power */
    r__1 = *cv / (*cv + 1.5f);
    g10 = 1.5f - *etazet * .66666666666666663f * (r__1 * r__1);
/* Computing 3rd power */
    r__1 = *cv;
/* Computing 2nd power */
    r__2 = *cv + 1.5f;
    gp10 = *etazet * sqrt(r__1 * (r__1 * r__1) * 2.f / 3.f) / (r__2 * r__2);
    gp10 = -gp10;
/* Computing 2nd power */
    r__1 = *cv + 1.5f;
    g2 = 1.5f - *etazet * *cv / (r__1 * r__1);
/* Computing 2nd power */
    r__1 = *cv + 1.5f;
/* Computing 2nd power */
    r__2 = *cv / (*cv + 1.5f) * *etazet;
    g011 = r__1 * r__1 * .4f + *cv * (*cv / 3.f + 1.f) * *etazet + r__2 * 
	    r__2 / (*etalam * 6.f);
/* Computing 2nd power */
    r__1 = *cv / (*cv + 1.5f);
    den = .26666666666666666f - *etalam + *etazet * .22222222222222221f * (
	    r__1 * r__1);
    g011 /= den;
/* Computing 2nd power */
    r__1 = *cv + 1.5f;
    g011 = g011 * 2.f * *cv * *etalam / (r__1 * r__1 * 3.f);
    g011 = 1.5f - g011;

    iy = 1;
    dtbar = 1.f / (*etalam * (*cv + 2.5f) * 2.f * *y);
    gbar = (1.f / *etazet + 1.3333333333333333f) / (*y * 2.f) + dtbar / (*cv 
	    + 1.5f);
/* 	Y2=1.5*Y */
/* 	Y3=H*Y */
    s6 = 0.f;
    imn_(intg, &c__6, x, y);
    s6get_(y, &g11, &gp11, &g10, &gp10, &g2, &g011, intg, &c__0, &s6, &s6e);
    *s = s6;
    return 0;
} /* tent_ */

/*       Subroutines for Tenti Program */
/*       from G. Tenti, University of Waterloo, Sept. 7, 1989 */
/*       from D. Seasholtz, NASA Lewis, March, 1994 */
/*       calls from J. Forkey's programs */
/*       CALL TENT(x,cv,etalam,etazet,s) */
/*       cv=internal specific/Boltzmann's constant */
/*       =0 for monoatomic gases */
/*       =0 for diatomic gases where rotational degrees of */
/* 	of freedom have fully kicked in (JL) */
/* 	etalam	reciprocal of Eucker ratio = viscosity*kb/thenn cond */
/* 	etazet	ratio of shear viscosity to bulk viscosity */
/* 	(bulk viscosity = 0 for monoatomic gases, so */
/* 	let etazet = 1000 for computational purposes) */
/* 	(for diatomic gas, the bulk viscosity depends on */
/* 	the Cv, the rotational relaxation time, and is */
/* 	linear in pressure.  For N2 at I atm, published */
/* 	values are between 1.367 and 1.407,( JL)) */
/*       .. revised for MS Fortran by J. Lock */
/*       .. revision date 12/9/90 */
/*       revised for MS Fortran by Jlock */
/*       SUB1.FOR */

/* Subroutine */ int imn_(doublecomplex *intg, integer *nmax, real *x, real *
	y)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    //double sqrt(double);
    //void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static double a;
    static integer k;
    //extern /* Double Complex */ VOID w_(doublecomplex *, doublecomplex *);
    static doublecomplex z__;
    static real ak;
    static double xx, yy;
    static doublecomplex sq2;
    static integer icg;
    static doublecomplex uni;
    static double rsq2;
    static doublecomplex acpl;
    static integer kmax;
    static doublecomplex zlam;

    /* Parameter adjustments */
    --intg;

    /* Function Body */
    uni = 0 + 1 * I;
    rsq2 = sqrt(2.);
    z__1 = rsq2 + 0.*I;
    sq2 = z__1;
    kmax = *nmax + 1;
    xx = (double) (*x);
    yy = (double) (*y);
    z__1 = xx + I * yy;
    z__ = z__1;
    z__2 = uni*sq2;//creal(uni) * creal(sq2) - cimag(uni) * cimag(sq2) + I*(creal(uni) * cimag(sq2) + cimag(uni) * creal(sq2));
    z__1 = z__2*z__;//creal(z__2) * creal(z__) - cimag(z__2) * cimag(z__) + I*(creal(z__2) * cimag(z__) + cimag(z__2) * creal(z__));
    zlam = z__1;
    k = 1;
    w_(&z__3, &z__);
    z__2 = uni*z__3;//creal(uni) * creal(z__3) - cimag(uni) * cimag(z__3) + I *(creal(uni) * cimag(z__3) + cimag(uni) * creal(z__3));
    z_div(&z__1, &z__2, &sq2);
    intg[1] = z__1;
    ak = 1.f;
    a = 1.;
    icg = 2;
    goto L30;
L5:
    switch (icg) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    ak *= k - 2;
    a = ak;
    icg = 2;
    goto L30;
L20:
    a = 0.f;
    icg = 1;
L30:
    z__1 = a + I*0;
    acpl = z__1;
    i__1 = k + 1;
    i__2 = k;
    z__3 = zlam*intg[i__2];// zlam.r * intg[i__2].r - zlam.i * intg[i__2].i, z__3.i = zlam.r * 
	    //intg[i__2].i + zlam.i * intg[i__2].r;
    z__2 = acpl - z__3;//, z__2.i = acpl.i - z__3.i;
    z__1 = z__2*uni;//z__2.r * uni.r - z__2.i * uni.i, z__1.i = z__2.r * uni.i + 
	   // z__2.i * uni.r;
    intg[i__1]=z__1;//.r = z__1.r, intg[i__1].i = z__1.i;
    ++k;
    if (k <= kmax) {
	goto L5;
    }
    return 0;
} /* imn_ */

/*       SUB2.FOR */

/* Double Complex */ void w_(doublecomplex * ret_val, doublecomplex *z__)
{
    /* Format strings */
    static char fmt_1001[] = "(\002 X AND Y ARE NOT IN THE FIRST QUADRANT"
	    "\002)";
    static char fmt_1003[] = "(\002CONTINUED FRACTION--N,XY\002,3x,i4,2x,e16"
	    ".5)";
    static char fmt_1004[] = "(\002POWER SERIES--N,XY\002,3x,i4,2x,e16.5)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    double d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    //double r_imag(complex *), sqrt(double);
    //integer s_wsfe(cilist *), e_wsfe(void);
    ///* Subroutine */ int s_stop(char *, ftnlen);
    //void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_exp(
	//    doublecomplex *, doublecomplex *);
    //integer do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublecomplex a, b;
    static integer n;
    static real x, y;
    static doublecomplex a0, a1, b0, b1, w1, w2;
    static double aa, pi;
    static complex cz;
    static real xy;
    static doublecomplex an1, bn1;
    static complex cz1;
    static real xy1;
    static integer icg;
    static doublecomplex sum, zsq, term, term1, sries, sqrpi;

    /* Fortran I/O blocks */
    //static cilist io___56 = { 0, 6, 0, fmt_1001, 0 };
    //static cilist io___77 = { 0, 6, 0, fmt_1003, 0 };
    //static cilist io___79 = { 0, 6, 0, fmt_1004, 0 };


    cz=*z__;
    x = creal(cz);
    y = cimag(cz);//r_imag(&cz);
    pi = M_PI;
    pi = sqrt(pi);
    z__1 = 0. + I*pi;
    sqrpi = z__1;
    assert (x >= 0.f && y >= 0.f) ;
L1:
    if (x <= 5.f && y <= 1.5f) {
	goto L50;
    }
    if (y >= 1.5f) {
	goto L30;
    }

    n = 1;
    z__1 = (*z__)*(*z__);// z__->r * z__->r - z__->i * z__->i, z__1.i = z__->r * z__->i + 
//	    z__->i * z__->r;
    zsq = z__1;//.r, zsq.i = z__1.i;
    z_div(&z__1, &c_b22, &zsq);
    term = z__1;//.r, term.i = z__1.i;
    sries = 1. + I*0.0;//, sries.i = 0.;
L5:
    ++n;
    z__1 = term + sries;//, z__1.i = term.i + sries.i;
    sries = z__1;//.r, sries.i = z__1.i;
    z__3 = term / 2.f;//, z__3.i = term.i / 2.f;
    z_div(&z__2, &z__3, &zsq);
    i__1 = (n << 1) - 1;
    d__1 = (double) i__1;
    z__1 = d__1 * z__2;//.r, z__1.i = d__1 * z__2.i;
    term1 = z__1;//.r, term1.i = z__1.i;
    cz = term;//.r, cz.i = term.i;
    cz1 = term1;//.r, cz1.i = term1.i;
    //xy = (r__1 = creal(cz), dabs(r__1)) + (r__2 = cimag(cz), dabs(r__2));
    //xy1 = (r__1 = creal(cz1), dabs(r__1)) + (r__2 = cimag(cz1), dabs(r__2));
    xy = dabs(creal(cz)) + dabs(cimag(cz));
    xy1 = dabs(creal(cz1)) + dabs(cimag(cz1));
    if (xy < 1e-15f || xy1 > xy) {
	goto L6;
    }
    term = term1;//.r, term.i = term1.i;
    goto L5;
L6:
    z__2 = -sries;//.r, z__2.i = -sries.i;
    z_div(&z__1, &z__2, z__);
    w1 = z__1;//.r, w1.i = z__1.i;
/* 	WRITE(*,1002) N,XY,XY1 */
/* 1002	FORMAT('#ASYMPTOTIC SERIES--N.XY,XY1',3X,I4,2X,2E16.5) */
     *ret_val = w1;//.r,  ret_val->i = w1.i;
    if (y != 0.f) {
	return ;
    }
/*       W=SQRPI*CDEXP(-ZXQ)+W1 'changed 12/9/90 */
    z__4 = -zsq;//.r, z__4.i = -zsq.i;
    z_exp(&z__3, &z__4);
    z__2 = sqrpi*z__3;//.r * z__3.r - sqrpi.i * z__3.i, z__2.i = sqrpi.r * z__3.i + 
	    //sqrpi.i * z__3.r;
    z__1 = z__2 + w1;//, z__1.i = z__2.i + w1.i;
     *ret_val = z__1;//.r,  ret_val->i = z__1.i;
    return ;

L30:
    n = 0;
    a0 = 1. + 0*I;//, a0.i = 0.;
    a1 = 0. + 0*I;//, a1.i = 0.;
    b0 = a1;//.r, b0.i = a1.i;
    b1 = a0;//.r, b1.i = a0.i;
    an1 = *z__;//->r, an1.i = z__->i;
    z__3 = -(*z__);//->r, z__3.i = -z__->i;
    z__2 = z__3*(*z__);//.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    //z__3.i * z__->r;
    z__1 = creal(z__2) + .5 + I*(cimag(z__2) + 0.);
    bn1 = z__1;//.r, bn1.i = z__1.i;
    icg = 1;
L35:
    z__2 = bn1*a1;//bn1.r * a1.r - bn1.i * a1.i, z__2.i = bn1.r * a1.i + bn1.i * 
	   // a1.r;
    z__3 = an1*a0;//an1.r * a0.r - an1.i * a0.i, z__3.i = an1.r * a0.i + an1.i * 
	  //  a0.r;
    z__1 = z__2 + z__3;//.r, z__1.i = z__2.i + z__3.i;
    a = z__1;//.r, a.i = z__1.i;
    z__2 = bn1*b1;//bn1.r * b1.r - bn1.i * b1.i, z__2.i = bn1.r * b1.i + bn1.i * 
	    //b1.r;
    z__3 = an1*b0;//.r * b0.r - an1.i * b0.i, z__3.i = an1.r * b0.i + an1.i * 
	 //   b0.r;
    z__1 = z__2 + z__3;//.r, z__1.i = z__2.i + z__3.i;
    b = z__1;//.r, b.i = z__1.i;
    switch (icg) {
	case 1:  goto L31;
	case 2:  goto L34;
    }
L31:
    an1 = 0.+ I*0.;//, an1.i = 0.;
    z_div(&z__1, &a, &b);
    w1 = z__1;//.r, w1.i = z__1.i;
    icg = 2;
    goto L33;
L34:
    z_div(&z__1, &a, &b);
    w2 = z__1;//.r, w2.i = z__1.i;
     *ret_val = w2;//.r,  ret_val->i = w2.i;
    z__1 = w2 - w1;//.r, z__1.i = w2.i - w1.i;
    cz = z__1;//.r, cz.i = z__1.i;
    //xy = (r__1 = creal(cz), dabs(r__1)) + (r__2 = cimag(cz), dabs(r__2));
    xy = dabs(creal(cz)) + dabs(cimag(cz));
    if (xy <= 1e-15f) {
	return ;
    }
    w1 = w2;//.r, w1.i = w2.i;
L33:
    a0 = a1;//.r, a0.i = a1.i;
    b0 = b1;//.r, b0.i = b1.i;
    a1 = a;//.r, a1.i = a.i;
    b1 = b;//.r, b1.i = b.i;
    z__1 = bn1 + 2.;//, z__1.i = bn1.i + 0.;
    bn1 = z__1;//.r, bn1.i = z__1.i;
    aa = n * -2.f - .5f;
    z__2 = aa;//, z__2.i = 0.;
    z__1 = an1 + z__2;//.r, z__1.i = an1.i + z__2.i;
    an1 = z__1;//.r, an1.i = z__1.i;
    ++n;
    if (n < 80) {
	goto L35;
    }
#if 0
    s_wsfe(&io___77);
    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&xy, (ftnlen)sizeof(real));
    e_wsfe();
#endif
    return ;

L50:
    n = 1;
    z__2 = (*z__)*(*z__);//->r * z__->r - z__->i * z__->i, z__2.i = z__->r * z__->i + 
	    //z__->i * z__->r;
    z__1 = z__2 * -2.0;//.r * -2. - z__2.i * 0., z__1.i = z__2.r * 0. + z__2.i * -2.;
    zsq = z__1;//.r, zsq.i = z__1.i;
    sum = 1.;//, sum.i = 0.;
    z__1 = zsq / 3.;//, z__1.i = zsq.i / 3.;
    term = z__1;//.r, term.i = z__1.i;
L55:
    z__1 = sum + term;//.r, z__1.i = sum.i + term.i;
    sum = z__1;//.r, sum.i = z__1.i;
    ++n;
    z__2 = term * zsq;//.r - term.i * zsq.i, z__2.i = term.r * zsq.i + 
	  //  term.i * zsq.r;
    i__1 = (n << 1) + 1;
    d__1 = (double) i__1;
    z__1 = z__2 / d__1; // + I* (cimag(z__2) / d__1);
    term = z__1;//.r, term.i = z__1.i;
    z_div(&z__1, &term, &sum);
    w1 = z__1;//.r, w1.i = z__1.i;
    cz = w1;//.r, cz.i = w1.i;
    //xy = (r__1 = creal(cz), dabs(r__1)) + (r__2 = cimag(cz), dabs(r__2));
    xy = dabs(creal(cz)) + dabs(cimag(cz));
    if (xy < 1e-15f) {
	goto L60;
    }
    if (n < 3000) {
	goto L55;
    }
    #if 0
    s_wsfe(&io___79);
    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&xy, (ftnlen)sizeof(real));
    e_wsfe();
    #endif
L60:
    z__2 = sum * -2.f;//, z__2.i = sum.i * -2.f;
    z__1 = z__2 * (*z__);//.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * z__->i + 
	  //  z__2.i * z__->r;
    sum = z__1;//.r, sum.i = z__1.i;
    z__4 = zsq / 2.f;//, z__4.i = zsq.i / 2.f;
    z_exp(&z__3, &z__4);
    z__2 = sqrpi * z__3;//.r - sqrpi.i * z__3.i, z__2.i = sqrpi.r * z__3.i + 
	    //sqrpi.i * z__3.r;
    z__1 = z__2 + sum;//.r, z__1.i = z__2.i + sum.i;
     *ret_val = z__1;//.r,  ret_val->i = z__1.i;
    return ;
} /* w_ */

/*       SUB3.FOR */

/* Subroutine */ int s6get_(real *y1, real *g11, real *gp11, real *g10, real *
	gp10, real *g2, real *g011, doublecomplex *intg, integer *ie, real *s,
	 real *se)
{
    /* Format strings */
    static char fmt_1005[] = "(\002 MATRIX BB6 IS SINGULAR\002)";
    static char fmt_1006[] = "(\002MATRIX BBE6 IS SINGULAR\002)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    //integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static doublecomplex b[36]	/* was [6][6] */, c__[6];
    static integer i__, j, k;
    static double bb[144]	/* was [12][12] */, cc[12];
    static doublecomplex be[16]	/* was [4][4] */, ce[4];
    static real pi;
    static integer ks;
    static double bbe[64]	/* was [8][8] */, cce[8];
    static doublecomplex gam[25]	/* was [5][5] */, one;
    extern /* Subroutine */ int simq_(double *, double *, integer *, 
	    integer *, integer *);
    static doublecomplex unit, zero;
    extern /* Subroutine */ int gamma_(doublecomplex *, doublecomplex *, 
	    integer *);

    /* Fortran I/O blocks */
    //static cilist io___97 = { 0, 1, 0, fmt_1005, 0 };
    //static cilist io___98 = { 0, 6, 0, fmt_1006, 0 };


    /* Parameter adjustments */
    --intg;

    /* Function Body */
    ks = 0;
    unit= 1*I;//r = 0., unit.i = 1.;
    one = 1.;//, one.i = 0.;
    zero = 0.;//, zero.i = 0.;
    pi = M_PI;
    cleararr(gam);
    cleararr(bbe);
    cleararr(cce);
    cleararr(bb);
    cleararr(cc);
    gamma_(gam, &intg[1], &c__4);
    for (i__ = 1; i__ <= 4; ++i__) {
	if (i__ <= 3) {
	    k = i__;
	}
	if (i__ == 4) {
	    k = 5;
	}
	i__1 = i__ - 1;
	i__2 = k - 1;
	c__[i__1] = gam[i__2];//.r, c__[i__1].i = gam[i__2].i;
	i__1 = i__ - 1;
	r__1 = -(*y1);
	i__2 = k - 1;
	z__1 = r__1 * gam[i__2];//.r, z__1.i = r__1 * gam[i__2].i;
	b[i__1] = z__1;//.r, b[i__1].i = z__1.i;
	i__1 = i__ + 5;
	r__1 = -(*y1);
	i__2 = k + 4;
	z__1 = r__1 * gam[i__2];//.r, z__1.i = r__1 * gam[i__2].i;
	b[i__1] = z__1;//.r, b[i__1].i = z__1.i;
	i__1 = i__ + 11;
	r__1 = -(*y1) * (*g11 - .5f);
	i__2 = k + 9;
	z__1 = r__1 * gam[i__2];//.r, z__1.i = r__1 * gam[i__2].i;
	b[i__1] = z__1;//.r, b[i__1].i = z__1.i;
	i__1 = i__ + 17;
	r__1 = -(*y1) * (*g10 - .5f);
	i__2 = k + 19;
	z__1 = r__1 * gam[i__2];//.r, z__1.i = r__1 * gam[i__2].i;
	b[i__1] = z__1;//.r, b[i__1].i = z__1.i;
	i__1 = i__ + 23;
	r__1 = -(*y1) * *gp10;
	i__2 = k + 19;
	z__1 = r__1 * gam[i__2];//.r, z__1.i = r__1 * gam[i__2].i;
	b[i__1] = z__1;//.r, b[i__1].i = z__1.i;
	i__1 = i__ + 29;
	r__1 = -(*y1) * *gp11;
	i__2 = k + 9;
	z__1 = r__1 * gam[i__2];//.r, z__1.i = r__1 * gam[i__2].i;
	b[i__1] = z__1;//.r, b[i__1].i = z__1.i;
	if (*ie == 0) {
	    goto L10;
	}
	i__1 = i__ - 1;
	i__2 = i__ - 1;
	ce[i__1] = c__[i__2];//.r, ce[i__1].i = c__[i__2].i;
	i__1 = i__ - 1;
	i__2 = i__ - 1;
	be[i__1] = b[i__2];//.r, be[i__1].i = b[i__2].i;
	i__1 = i__ + 3;
	i__2 = i__ + 5;
	be[i__1] = b[i__2];//.r, be[i__1].i = b[i__2].i;
	i__1 = i__ + 7;
	r__1 = -(*y1);
	i__2 = k + 9;
	z__2 = r__1 * gam[i__2];//.r, z__2.i = r__1 * gam[i__2].i;
	z__1 = z__2 / 3.f;//, z__1.i = z__2.i / 3.f;
	be[i__1] = z__1;//.r, be[i__1].i = z__1.i;
	i__1 = i__ + 11;
	r__1 = -(*y1);
	i__2 = k + 19;
	z__1 = r__1 * gam[i__2];//.r, z__1.i = r__1 * gam[i__2].i;
	be[i__1] = z__1;//.r, be[i__1].i = z__1.i;
L10:
	;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    if (i__ != j) {
		goto L15;
	    }
	    i__1 = i__ + j * 6 - 7;
	    i__2 = i__ + j * 6 - 7;
	    z__1 = b[i__2] - one;//.r, z__1.i = b[i__2].i - one.i;
	    b[i__1] = z__1;//.r, b[i__1].i = z__1.i;
	    if (*ie == 1) {
		i__1 = i__ + (j << 2) - 5;
		i__2 = i__ + (j << 2) - 5;
		z__1 = be[i__2] - one;//.r, z__1.i = be[i__2].i - one.i;
		be[i__1] = z__1;//.r, be[i__1].i = z__1.i;
	    }
L15:
	    ;
	}
/* L14: */
    }
    c__[4] = zero;//.r, c__[4].i = zero.i;
    c__[5] = zero;//.r, c__[5].i = zero.i;
    b[4] = zero;//.r, b[4].i = zero.i;
    b[10] = zero;//.r, b[10].i = zero.i;
    r__1 = -(*y1) * *gp11;
    z__1 = r__1 * gam[5];//.r, z__1.i = r__1 * gam[5].i;
    b[16] = z__1;//, b[16].i = z__1.i;
    r__1 = -(*y1) * *gp10;
    z__1 = r__1 * gam[0];//.r, z__1.i = r__1 * gam[0].i;
    b[22] = z__1;//.r, b[22].i = z__1.i;
    r__1 = -(*y1) * (*g2 - .5f);
    z__2 = r__1 * gam[0];//.r, z__2.i = r__1 * gam[0].i;
    z__1 = z__2 - one;//, z__1.i = z__2.i - one.i;
    b[28] = z__1;//.r, b[28].i = z__1.i;
    r__1 = -(*y1) * (*g011 - .5f);
    z__1 = r__1 * gam[5];//.r, z__1.i = r__1 * gam[5].i;
    b[34] = z__1;//, b[34].i = z__1.i;
    b[5] = zero;//, b[5].i = zero.i;
    b[11] = zero;//.r, b[11].i = zero.i;
    r__1 = -(*y1) * *gp11;
    z__1 = r__1 * gam[6];//.r, z__1.i = r__1 * gam[6].i;
    b[17] = z__1;//, b[17].i = z__1.i;
    r__1 = -(*y1) * *gp10;
    z__1 = r__1 * gam[1];//, z__1.i = r__1 * gam[1].i;
    b[23] = z__1;//, b[23].i = z__1.i;
    r__1 = -(*y1) * (*g2 - .5f);
    z__1 = r__1 * gam[1];//, z__1.i = r__1 * gam[1].i;
    b[29] = z__1;//.r, b[29].i = z__1.i;
    r__1 = -(*y1) * (*g011 - .5f);
    z__2 = r__1 * gam[6];//r, z__2.i = r__1 * gam[6].i;
    z__1 = z__2 - one;//.r, z__1.i = z__2.i - one.i;
    b[35] = z__1;//.r, b[35].i = z__1.i;
    for (i__ = 1; i__ <= 6; ++i__) {
	i__1 = i__ - 1;
	i__2 = i__ - 1;
	cc[i__1] = creal(c__[i__2]);
	i__1 = i__ + 5;
	z__2 = -unit;//.r, z__2.i = -unit.i;
	i__2 = i__ - 1;
	z__1 = z__2 * c__[i__2];//.r - z__2.i * c__[i__2].i, z__1.i = z__2.r 
	//	* c__[i__2].i + z__2.i * c__[i__2].r;
	cc[i__1] = creal(z__1);//.r;
	if (i__ > 4 || *ie == 0) {
	    goto L19;
	}
	cce[i__ - 1] = cc[i__ - 1];
/*       CCE(TT+4)=-UNIT*CE(I) 'changed 12/9/90 */
	i__1 = i__ + 3;
	z__2 = -unit;//.r, z__2.i = -unit.i;
	i__2 = i__ - 1;
	z__1 = z__2 * ce[i__2];//.r - z__2.i * ce[i__2].i, z__1.i = z__2.r * 
	//	ce[i__2].i + z__2.i * ce[i__2].r;
	cce[i__1] = creal(z__1);//.r;
L19:
	for (j = 1; j <= 6; ++j) {
	    i__1 = i__ + j * 12 - 13;
	    i__2 = i__ + j * 6 - 7;
	    bb[i__1] = creal(b[i__2]);//.r;
	    bb[i__ + 6 + (j + 6) * 12 - 13] = bb[i__ + j * 12 - 13];
	    i__1 = i__ + (j + 6) * 12 - 13;
	    i__2 = i__ + j * 6 - 7;
	    z__1 = unit * b[i__2];//.r - unit.i * b[i__2].i, z__1.i = unit.r 
		    //* b[i__2].i + unit.i * b[i__2].r;
	    bb[i__1] = creal(z__1);
	    bb[i__ + 6 + j * 12 - 13] = -bb[i__ + (j + 6) * 12 - 13];
	    if (i__ > 4 || j > 4 || *ie == 0) {
		goto L18;
	    }
	    i__1 = i__ + (j << 3) - 9;
	    i__2 = i__ + (j << 2) - 5;
	    bbe[i__1] = creal(be[i__2]);
	    bbe[i__ + 4 + (j + 4 << 3) - 9] = bbe[i__ + (j << 3) - 9];
	    i__1 = i__ + (j + 4 << 3) - 9;
	    i__2 = i__ + (j << 2) - 5;
	    z__1 = unit * be[i__2];//.r - unit.i * be[i__2].i, z__1.i = 
		    //unit.r * be[i__2].i + unit.i * be[i__2].r;
	    bbe[i__1] = creal(z__1);
	    bbe[i__ + 4 + (j << 3) - 9] = -bbe[i__ + (j + 4 << 3) - 9];
L18:
	    ;
	}
/* L17: */
    }
    simq_(bb, cc, &c__12, &ks, &c__12);
    if (ks == 1) {
        printf(fmt_1005);
	//s_wsfe(&io___97);
	//e_wsfe();
    }
    *s = cc[0] / pi;
    if (*ie == 1) {
	goto L32;
    }
    *se = 0.f;
    return 0;
L32:
    simq_(bbe, cce, &c__8, &ks, &c__8);
    if (ks == 1) {
        printf(fmt_1006);
	//s_wsfe(&io___98);
	//e_wsfe();
    }
    *se = cce[0] / pi;
    return 0;
} /* s6get_ */

/*       SUB4.FOR */

/* Subroutine */ int gamma_(doublecomplex *gam, doublecomplex *intg, integer *
	index)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    real r__1, r__2, r__3, r__4;
    double d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    //double pow_dd(double *, double *), sqrt(double), pow_ri(real *
	//    , integer *);

    /* Local variables */
    static integer j, k, l, jm, km, ir, jt, kt, lt, ime, jme;
    static integer lme[5] = {0,1,1,2,0};
    static double crl;
    static integer jmt, kmt, irt;
    static doublecomplex gama;
    static double bigg;
    static integer irme[5] = {0,0,1,0,1}; 
    extern /* Double Complex */ void intmn_(doublecomplex *, integer *, 
	    integer *, doublecomplex *);
    extern double barnes_(real *, integer *);

/* 	INTEGER IRME(5)/0,0,1,0,1/ */
/* 	INTEGER LME(5)/0,1,1,2,0/ */
    /* Parameter adjustments */
    --intg;
    gam -= 6;

    /* Function Body */
    for (ime = 1; ime <= 5; ++ime) {
	if (ime == *index) {
	    goto L100;
	}
	ir = irme[ime - 1];
	l = lme[ime - 1];
	for (jme = 1; jme <= 5; ++jme) {
	    if (jme == *index) {
		goto L200;
	    }
	    irt = irme[jme - 1];
	    lt = lme[jme - 1];
/* L10: */
	    d__1 = (double) ((l + lt) * 1.5f + 1.f);
	    crl = 1. / (double) pow(c_b77, d__1);
	    i__1 = l << 1;
/* Computing 2nd power */
	    d__1 = barnes_(&c_b78, &l);
	    i__2 = lt << 1;
/* Computing 2nd power */
	    d__2 = barnes_(&c_b78, &lt);
	    crl = crl * (barnes_(&c_b78, &i__1) / (d__1 * d__1)) * (barnes_(&
		    c_b78, &i__2) / (d__2 * d__2)) / sqrt(barnes_(&c_b78, &ir)
		     * barnes_(&c_b78, &irt));
/* L15: */
	    i__1 = l + ir + 1;
	    i__2 = lt + irt + 1;
	    crl *= sqrt(((l << 1) + 1) * ((lt << 1) + 1) / barnes_(&c_b85, &
		    i__1) / barnes_(&c_b85, &i__2));
	    km = l / 2 + 1;
/* L20: */
	    kmt = lt / 2 + 1;
	    jm = ir + 1;
	    jmt = irt + 1;
	    gama = 0.f + I* 0.f;
/* L30: */
	    i__1 = jm;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = jmt;
		for (jt = 1; jt <= i__2; ++jt) {
		    i__3 = km;
		    for (k = 1; k <= i__3; ++k) {
			i__4 = kmt;
			for (kt = 1; kt <= i__4; ++kt) {
			    i__5 = j + jt + k + kt;
			    i__6 = (k + kt << 1) - 4 + jm - j + jmt - jt;
			    bigg = (double) pow(c_b90, i__5) / (
				    double) pow(c_b91, i__6);
/* L40: */
			    r__1 = jm - j + 1.f;
			    i__5 = j - 1;
			    i__6 = j - 1;
			    r__2 = jmt - jt + 1.f;
			    i__7 = jt - 1;
			    i__8 = jt - 1;
			    bigg = bigg * barnes_(&r__1, &i__5) / barnes_(&
				    c_b78, &i__6) * barnes_(&r__2, &i__7) / 
				    barnes_(&c_b78, &i__8);
			    r__1 = irt + lt + 2.5f - jt;
			    i__5 = jt - 1;
			    r__2 = lt + 1.5f - kt;
			    i__6 = kt - 1;
			    r__3 = ir + l + 2.5f - j;
			    i__7 = j - 1;
			    r__4 = l + 1.5f - k;
			    i__8 = k - 1;
			    bigg = bigg * barnes_(&r__1, &i__5) / barnes_(&
				    r__2, &i__6) * barnes_(&r__3, &i__7) / 
				    barnes_(&r__4, &i__8);
			    r__1 = l - (k << 1) + 3.f;
			    i__5 = k - 1 << 1;
			    i__6 = k - 1;
			    r__2 = lt - (kt << 1) + 3.f;
			    i__7 = kt - 1 << 1;
			    i__8 = kt - 1;
			    bigg = bigg * barnes_(&r__1, &i__5) / barnes_(&
				    c_b78, &i__6) * barnes_(&r__2, &i__7) / 
				    barnes_(&c_b78, &i__8);
/* L42: */
			    i__5 = ir + irt - j - jt + k + kt;
			    i__6 = l + lt - (k + kt << 1) + 4;
			    intmn_(&z__3, &i__5, &i__6, &intg[1]);
			    z__2 = bigg * z__3;//.r, z__2.i = bigg * z__3.i;
			    z__1 = z__2 + gama;//.r, z__1.i = z__2.i + 
				   // gama.i;
			    gama = z__1;//.r, gama.i = z__1.i;
/* L70: */
			}
		    }
		}
	    }
/* L90: */
	    z__2 = crl * gama;//.r, z__2.i = crl * gama.i;
	    d__1 = sqrt(2.);
	    z__1 = d__1 * z__2;//.r, z__1.i = d__1 * z__2.i;
	    gama = z__1;//, gama.i = z__1.i;
	    i__4 = ime + jme * 5;
	    gam[i__4] = gama;//.r, gam[i__4].i = gama.i;
L200:
	    ;
	}
L100:
	;
    }
    return 0;
} /* gamma_ */

/* SUB5.FOR */

double barnes_(real *b, integer *m)
{
    /* System generated locals */
    integer i__1;
    double ret_val;

    /* Local variables */
    static integer n;

    ret_val = 1.;
    if (*m == 0) {
	return ret_val;
    }
    i__1 = *m;
    for (n = 1; n <= i__1; ++n) {
	ret_val *= *b + (double) (n - 1);
/* L10: */
    }
    return ret_val;
} /* barnes_ */

/* SUB6.FOR */

/* Double Complex */ void intmn_(doublecomplex * ret_val, integer *m, integer 
	*n, doublecomplex *intg)
{
    /* System generated locals */
    integer i__1, i__2;
    double d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    //double pow_di(double *, integer *);

    /* Local variables */
    static integer ip, iq;
    static double subs;
    static integer iqmx;
    extern double barnes_(real *, integer *);

    /* Parameter adjustments */
    --intg;

    /* Function Body */
/* L1: */
    ip = 0;
     *ret_val = 0.;//,  ret_val->i = 0.;
L5:
    iqmx = *m - ip;
/* L7: */
    subs = 0.;
/* L10: */
    iq = 0;
L15:
    i__1 = *m - ip - iq;
    i__2 = *m - ip - iq;
    subs = barnes_(&c_b85, &iq) / barnes_(&c_b78, &iq) * barnes_(&c_b85, &
	    i__1) / barnes_(&c_b78, &i__2) + subs;
    ++iq;
    if (iq - 1 < iqmx) {
	goto L15;
    }
    i__1 = *n + (ip << 1) + 1;
    z__4 = subs * intg[i__1];//.r, z__4.i = subs * intg[i__1].i;
    i__2 = *m - ip;
    d__1 = pow(c_b77, i__2);
    z__3 = d__1 * z__4;//.r, z__3.i = d__1 * z__4.i;
    d__2 = barnes_(&c_b78, &ip);
    z__2 = z__3 / d__2;//, z__2.i = z__3.i / d__2;
    z__1 = z__2 +  *ret_val;//->r, z__1.i = z__2.i +  ret_val->i;
     *ret_val = z__1;//.r,  ret_val->i = z__1.i;
    ++ip;
    if (ip - 1 < *m) {
	goto L5;
    }
    d__1 = barnes_(&c_b78, m);
    z__1 = d__1 *  (*ret_val);//->r, z__1.i = d__1 *  ret_val->i;
     *ret_val = z__1;//.r,  ret_val->i = z__1.i;
/* L30: */
    return ;
} /* intmn_ */

/* SUB7.FOR */

/* Subroutine */ int simq_(double *a, double *b, integer *n, integer *
	ks, integer *m)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    double d__1;

    /* Local variables */
    static integer i__, j, k, ib, ic, it, ix, jx, jy, ny;
    static double tol, biga, save;
    static integer imax;

/*       REAL*8 A(M,M),B(l),TOL,BIGA,SAVE */
    /* Parameter adjustments */
    --b;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    tol = 0.f;
    *ks = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jy = j + 1;
	biga = 0.;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {

	    if (fabs(biga) - (d__1 = a[i__ + j * a_dim1], fabs(d__1)) >= 0.) {
		goto L30;
	    } else {
		goto L20;
	    }
L20:
	    biga = a[i__ + j * a_dim1];
	    imax = i__;
L30:
	    ;
	}

	if (fabs(biga) - tol <= 0.) {
	    goto L35;
	} else {
	    goto L40;
	}
L35:
	*ks = 1;
	return 0;

L40:
	i__2 = *n;
	for (k = j; k <= i__2; ++k) {
	    save = a[j + k * a_dim1];
	    a[j + k * a_dim1] = a[imax + k * a_dim1];
	    a[imax + k * a_dim1] = save;

/* L50: */
	    a[j + k * a_dim1] /= biga;
	}
	save = b[imax];
	b[imax] = b[j];
	b[j] = save / biga;

	if (j - *n != 0) {
	    goto L55;
	} else {
	    goto L70;
	}
L55:
	i__2 = *n;
	for (ix = jy; ix <= i__2; ++ix) {
	    it = j - ix;
	    i__3 = *n;
	    for (jx = jy; jx <= i__3; ++jx) {
/* L60: */
		a[ix + jx * a_dim1] -= a[ix + (ix + it) * a_dim1] * a[ix + it 
			+ jx * a_dim1];
	    }
/* L65: */
	    b[ix] -= b[j] * a[ix + (ix + it) * a_dim1];
	}
    }

L70:
    ny = *n - 1;
    it = *n * *n;
    i__2 = ny;
    for (j = 1; j <= i__2; ++j) {
	ib = *n - j;
	ic = *n;
	i__1 = j;
	for (k = 1; k <= i__1; ++k) {
	    b[ib] -= a[ib + ic * a_dim1] * b[ic];
/* L80: */
	    --ic;
	}
    }
    return 0;
} /* simq_ */

