# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <iostream>

using namespace std;

# include "elliptic_integral.hpp"

int main ( );
void rc_test ( );
void rc_test2 ( );
void rd_test ( );
void rf_test ( );
void rj_test ( );
void elliptic_ea_test ( );
void elliptic_ek_test ( );
void elliptic_em_test ( );
void elliptic_fa_test ( );
void elliptic_fk_test ( );
void elliptic_fm_test ( );
void elliptic_pia_test ( );
void elliptic_pik_test ( );
void elliptic_pim_test ( );
void elliptic_inc_ea_test ( );
void elliptic_inc_ek_test ( );
void elliptic_inc_em_test ( );
void elliptic_inc_fa_test ( );
void elliptic_inc_fk_test ( );
void elliptic_inc_fm_test ( );
void elliptic_inc_pia_test ( );
void elliptic_inc_pik_test ( );
void elliptic_inc_pim_test ( );
void jacobi_cn_test ( );
void jacobi_dn_test ( );
void jacobi_sn_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ELLIPTIC_INTEGRAL_TEST.
//
//  Discussion:
//
//    ELLIPTIC_INTEGRAL_TEST tests ELLIPTIC_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ELLIPTIC_INTEGRAL_TEST\n";
  cout << "  C++ version\n";
  cout << "  ELLIPTIC_INTEGRAL evaluates elliptic integral functions\n";
  cout << "  using Carlson's elliptic functions.\n";
//
//  Carloson integrals.
//
  rc_test ( );
  rc_test2 ( );
  rd_test ( );
  rf_test ( );
  rj_test ( );
//
//  Complete elliptic integrals, first kind.
//
  elliptic_fa_test ( );
  elliptic_fk_test ( );
  elliptic_fm_test ( );
//
//  Complete elliptic integrals, second kind.
//
  elliptic_ea_test ( );
  elliptic_ek_test ( );
  elliptic_em_test ( );
//
//  Complete elliptic integrals, third kind.
//
  elliptic_pia_test ( );
  elliptic_pik_test ( );
  elliptic_pim_test ( );
//
//  Incomplete elliptic integrals, first kind.
//
  elliptic_inc_fa_test ( );
  elliptic_inc_fk_test ( );
  elliptic_inc_fm_test ( );
//
//  Incomplete elliptic integrals, second kind.
//
  elliptic_inc_ea_test ( );
  elliptic_inc_ek_test ( );
  elliptic_inc_em_test ( );
//
//  Incomplete elliptic integrals, third kind.
//
  elliptic_inc_pia_test ( );
  elliptic_inc_pik_test ( );
  elliptic_inc_pim_test ( );
//
//  Jacobi elliptic functions.
//
  jacobi_cn_test ( );
  jacobi_dn_test ( );
  jacobi_sn_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ELLIPTIC_INTEGRAL_TEST\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void rc_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RC_TEST tests RC.
//
//  Discussion:
//
//    This driver tests the function for the
//    integral RC(X,Y).  The first six sets of values of X and Y are
//    extreme points of the region of valid arguments defined by the
//    machine-dependent constants LOLIM and UPLIM.  The values of LOLIM,
//    UPLIM, X, Y, and ERRTOL (see comments in void) may be used on
//    most machines but provide a severe test of robustness only on the
//    ibm 360/370 series.  The seventh set tests the failure exit.  The
//    next three sets are check values: RC(0,0.25) = RC(0.0625,0.125) = PI
//    and RC(2.25,2) = LN(2).  The remaining sets show the dependence on X
//    when Y = 1.  Fixing Y entails no loss here because RC is homogeneous.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2018
//
//  Author:
//
//    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//    This C++ version by John Burkardt.
//
{
  double eliptc;
  double errtol;
  int i;
  int ierr;
  static double x[] = {
   1.51E-78, 
   3.01E-78, 
   0.00E+00, 
   0.99E+75, 
   0.00E+00, 
   0.99E+75, 
   0.00E+00, 
   0.00E+00, 
   6.25E-02, 
   2.25E+00, 
   0.01E+00, 
   0.02E+00, 
   0.05E+00, 
   0.10E+00, 
   0.20E+00, 
   0.40E+00, 
   0.60E+00, 
   0.80E+00, 
   1.00E+00, 
   1.20E+00, 
   1.50E+00, 
   2.00E+00, 
   3.00E+00, 
   4.00E+00, 
   5.00E+00, 
   1.00E+01, 
   2.00E+01, 
   5.00E+01, 
   1.00E+02, 
   1.00E+03, 
   1.00E+04, 
   1.00E+05, 
   1.00E+06, 
   1.00E+07, 
   1.00E+08, 
   1.00E+09, 
   1.00E+10, 
   1.00E+12, 
   1.00E+15, 
   1.00E+20, 
   1.00E+30, 
   1.00E+40, 
   1.00E+50 };

  static double y[] = { 
   1.51E-78, 
   0.55E-78, 
   3.01E-78, 
   0.55E-78, 
   0.99E+75, 
   0.99E+75, 
   2.00E-78, 
   2.50E-01, 
   1.25E-01, 
   2.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00 };

  cout << "\n";
  cout << "RC_TEST\n";
  cout << "  RC evaluates the elementary integral RC(X,Y)\n";

  cout << "\n";
  cout << "              X                          Y                         RC(X,Y)\n";
  cout << "\n";

  errtol = 1.0E-3;

  for ( i = 0; i < 43; i++ )
  {
    eliptc = rc ( x[i], y[i], errtol, ierr );
    cout << "  " << setw(27) << setprecision(16) << x[i]
         << "  " << setw(27) << setprecision(16) << y[i];
    if ( ierr == 0 )
    {
      cout << setw(27) << setprecision(16) << eliptc << "\n";
    }
    else
    {
      cout << "  ***Error***\n";
    }
  }

  return;
}
//****************************************************************************80

void rc_test2 ( )

//****************************************************************************80
//
//  Purpose:
//
//    RC_TEST2 checks RC by examining special values.
//
//  Discussion:
//
//    This driver compares values of (LOG X)/(X-1) and ARCTAN(X)
//    calculated on one hand from the void RC and on the other
//    from library LOG and ARCTAN routines.  to avoid over/underflows
//    for extreme values of X, we write (LOG X)/(X-1) = RC(Y,X/Y)/SQRT(Y),
//    where Y = (1+X)/2, and ARCTAN(X) = SQRT(X)*RC(Z,Z+X), where Z = 1/X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2018
//
//  Author:
//
//    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//    This C++ version by John Burkardt.
//
{
  double errtol;
  int i;
  double ibmarc;
  double ibmlog;
  int ierr;
  int ipower;
  int j;
  int m;
  double myarc;
  double mylog;
  double v;
  double w;
  double x;
  static double x_vec[] = {
   1.0E-75, 
   1.0E-15, 
   1.0E-03, 
   1.0E-01, 
   2.0E-01, 
   5.0E-01, 
   1.0E+00, 
   2.0E+00, 
   5.0E+00, 
   1.0E+01, 
   1.0E+03, 
   1.0E+15, 
   1.0E+75 };
  double y;
  double z;

  cout << "\n";
  cout << "RC_TEST2\n";
  cout << "  Compare LOG(X)/(X-1) and ARCTAN(X) with\n";
  cout << "  values based on RC.\n";

  cout << "\n";
  cout << "     X                From LOG                   From RC\n";
  cout << "\n";

  errtol = 1.0E-3;

  for ( j = 1; j <= 10; j++ )
  {
    x = 0.2 * ( double ) ( j );
    y = ( 1.0 + x ) / 2.0;
    v = x / y;
    mylog = rc ( y, v, errtol, ierr ) / sqrt ( y );
    cout << setw(9) << setprecision(1) << x << "     ";
    if ( j == 5 )
    {
      cout << "**** ZERO DIVIDE *****";
    }
    else
    {
      ibmlog = log ( x ) / ( x - 1.0 );
      cout << setw(27) << setprecision(16) << ibmlog;
    }
    cout << setw(27) << setprecision(16) << mylog << "\n";
  }

  cout << "\n";
  cout << "  Extreme values of X\n";
  cout << "\n";
  cout << "     X                From LOG                   From RC\n";
  cout << "\n";

  for ( i = 0; i < 16; i++ )
  {
    ipower = - 75 + 10 * i;
    x = pow ( 10.0, ipower );
    y = ( 1.0 + x ) / 2.0;
    v = x / y;
    mylog = rc ( y, v, errtol, ierr ) / sqrt ( y );
    ibmlog = log ( x ) / ( x - 1.0 );
    cout << setw(8) << setprecision(1) << x
         << setw(27) << setprecision(16) << ibmlog
         << setw(27) << setprecision(16) << mylog << "\n";
  }

  cout << "\n";
  cout << "     X              From ARCTAN                 From RC\n";
  cout << "\n";

  for ( m = 0; m < 13; m++ )
  {
    x = x_vec[m];
    z = 1.0 / x;
    w = z + x;
    myarc = sqrt ( x ) * rc ( z, w, errtol, ierr );
    ibmarc = atan ( x );
    cout << setw(8) << setprecision(1) << x
         << setw(27) << setprecision(16) << ibmarc
         << setw(27) << setprecision(16) << myarc << "\n";
  }

  return;
}
//****************************************************************************80

void rd_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RD_TEST tests RD.
//
//  Discussion:
//
//    This driver tests the function for the
//    integral RD(X,Y,Z), which is symmetric in X and Y.  The first
//    twelve sets of values of X, Y, Z are extreme points of the region of
//    valid arguments defined by the machine-dependent constants LOLIM
//    and UPLIM.  The values of LOLIM, UPLIM, X, Y, Z, and ERRTOL (see
//    comments in void) may be used on most machines but provide a
//    severe test of robustness only on the ibm 360/370 series.  The
//    thirteenth set tests the failure exit.  The fourteenth set is a
//    check value: RD(0,2,1) = 3B = 3(PI)/4A, where A and B are the
//    lemniscate constants.  The remaining sets show the dependence
//    on Z when Y = 1 (no loss of generality because of homogeneity)
//    and X = 0.5 (midway between the complete case X = 0 and the
//    degenerate case X = Y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2018
//
//  Author:
//
//    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//    This C++ version by John Burkardt.
//
{
  double eliptc;
  double errtol;
  int i;
  int ierr;
  static double x[] = {
   0.00E+00, 
   0.55E-78, 
   0.00E+00, 
   0.55E-78, 
   0.00E+00, 
   0.55E-78, 
   0.00E+00, 
   0.55E-78, 
   3.01E-51, 
   3.01E-51, 
   0.99E+48, 
   0.99E+48, 
   0.00E+00, 
   0.00E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00 };

  static double y[] = { 
   6.01E-51, 
   6.01E-51, 
   6.01E-51, 
   6.01E-51, 
   0.99E+48, 
   0.99E+48, 
   0.99E+48, 
   0.99E+48, 
   3.01E-51, 
   3.01E-51, 
   0.99E+48, 
   0.99E+48, 
   3.01E-51, 
   2.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00 };

  static double z[] = {
   6.01E-51, 
   6.01E-51, 
   0.99E+48, 
   0.99E+48, 
   6.01E-51, 
   6.01E-51, 
   0.99E+48, 
   0.99E+48, 
   6.01E-51, 
   0.99E+48, 
   6.01E-51, 
   0.99E+48, 
   1.00E+00, 
   1.00E+00, 
   1.00E-10, 
   1.00E-05, 
   1.00E-02, 
   1.00E-01, 
   2.00E-01, 
   5.00E-01, 
   1.00E+00, 
   2.00E+00, 
   5.00E+00, 
   1.00E+01, 
   1.00E+02, 
   1.00E+05, 
   1.00E+10 };
          
  cout << "\n";
  cout << "RD_TEST\n";
  cout << "  RD evaluates the Carlson elliptic integral\n";
  cout << "  of the second kind, RD(X,Y,Z)\n";
  cout << "\n";
  cout << "               X                          Y";
  cout << "                          Z                         RD(X,Y,Z)\n";
  cout << "\n";

  errtol = 1.0E-03;

  for ( i = 0; i < 27; i++ )
  {
    eliptc = rd ( x[i], y[i], z[i], errtol, ierr );
    cout << setw(27) << setprecision(16) << x[i]
         << setw(27) << setprecision(16) << y[i]
         << setw(27) << setprecision(16) << z[i];
    if (ierr == 0 )
    {
      cout << setw(27) << setprecision(16) << eliptc << "\n";
    }
    else
    {
      cout << "  ***Error***\n";
    }
  }

  return;
}
//****************************************************************************80

void rf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RF_TEST tests RF.
//
//  Discussion:
//
//    This driver tests the function for the
//    integral RF(X,Y,Z), which is symmetric in X, Y, Z.  The first nine
//    sets of values of X, Y, Z are extreme points of the region of valid
//    arguments defined by the machine-dependent constants LOLIM and
//    UPLIM.  The values of LOLIM, UPLIM, X, Y, Z, and ERRTOL (see
//    comments in void) may be used on most machines but provide a
//    severe test of robustness only on the ibm 360/370 series.  The
//    tenth set tests the failure exit.  The eleventh set is a check
//    value: RF(0,1,2) = A, where A is the first lemniscate constant.
//    The remaining sets show the dependence on Z when Y = 1 (no loss of
//    generality because of homogeneity) and X = 0.5 (midway between the
//    complete case X = 0 and the degenerate case X = Y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2018
//
//  Author:
//
//    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//    This C++ version by John Burkardt.
//
{
  double eliptc;
  double errtol;
  int i;
  int ierr;
  static double x[] = {
   1.51E-78, 
   1.51E-78, 
   0.00E+00, 
   0.00E+00, 
   0.00E+00, 
   0.99E+75, 
   0.55E-78, 
   0.55E-78, 
   0.55E-78, 
   0.00E+00, 
   0.00E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00 };

  static double y[] = {
   1.51E-78, 
   1.51E-78, 
   3.01E-78, 
   3.01E-78, 
   0.99E+75, 
   0.99E+75, 
   3.01E-78, 
   3.01E-78, 
   0.99E+75, 
   2.00E-78, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00 };

  static double z[] = {
   1.51E-78, 
   0.99E+75, 
   3.01E-78, 
   0.99E+75, 
   0.99E+75, 
   0.99E+75, 
   3.01E-78, 
   0.99E+75, 
   0.99E+75, 
   1.00E+00, 
   2.00E+00, 
   1.00E+00, 
   1.10E+00, 
   1.20E+00, 
   1.30E+00, 
   1.40E+00, 
   1.50E+00, 
   1.60E+00, 
   1.70E+00, 
   1.80E+00, 
   1.90E+00, 
   2.00E+00, 
   2.20E+00, 
   2.40E+00, 
   2.60E+00, 
   2.80E+00, 
   3.00E+00, 
   3.50E+00, 
   4.00E+00, 
   4.50E+00, 
   5.00E+00, 
   6.00E+00, 
   7.00E+00, 
   8.00E+00, 
   9.00E+00, 
   1.00E+01, 
   2.00E+01, 
   3.00E+01, 
   4.00E+01, 
   5.00E+01, 
   1.00E+02, 
   2.00E+02, 
   5.00E+02, 
   1.00E+03, 
   1.00E+04, 
   1.00E+05, 
   1.00E+06, 
   1.00E+08, 
   1.00E+10, 
   1.00E+12, 
   1.00E+15, 
   1.00E+20, 
   1.00E+30, 
   1.00E+40, 
   1.00E+50 };

  cout << "\n";
  cout << "RF_TEST\n";
  cout << "  RF evaluates the Carlson elliptic integral\n";
  cout << "  of the first kind, RF(X,Y,Z)\n";
  cout << "\n";
  cout << "               X                          Y";
  cout << "                          Z                         RF(X,Y,Z)\n";
  cout << " \n";

  errtol = 1.0E-3;

  for ( i = 0; i < 55; i++ )
  {
    eliptc = rf ( x[i], y[i], z[i], errtol, ierr );
    cout << setw(27) << setprecision(16) << x[i]
         << setw(27) << setprecision(16) << y[i]
         << setw(27) << setprecision(16) << z[i];
    if (ierr == 0 )
    {
      cout << setw(27) << setprecision(16) << eliptc << "\n";
    }
    else
    {
      cout << "  ***Error***\n";
    }
  }

  return;
}
//****************************************************************************80

void rj_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RJ_TEST tests RJ.
//
//  Discussion:
//
//    This driver tests the function for the
//    integral Rj(X,Y,Z,P), which is symmetric in X, Y, Z.  The first
//    twenty sets of values of X, Y, Z, P are extreme points of the region
//    of valid arguments defined by the machine-dependent constants
//    LOLIM and UPLIM.  The values of LOLIM, UPLIM, X, Y, Z, P, and
//    ERRTOL (see comments in void) may be used on most machines
//    but provide a severe test of robustness only on the ibm 360/370
//    series.  The twenty-first set tests the failure exit.  The twenty-
//    second set is a check value:
//      RJ(2,3,4,5) = 0.1429757966715675383323308.
//    The remaining sets show the dependence on Z and P
//    when Y = 1 (no loss of generality because of homogeneity) and
//    X = 0.5 (midway between the complete case x = 0 and the degenerate
//    case X = Y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2018
//
//  Author:
//
//    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//    This C++ version by John Burkardt.
//
{
  double eliptc;
  double errtol;
  int i;
  int ierr;

  static double p[] = { 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   1.00E+00, 
   5.00E+00, 
   0.25E+00, 
   0.75E+00, 
   1.00E+00, 
   2.00E+00, 
   0.25E+00, 
   0.75E+00, 
   1.50E+00, 
   4.00E+00, 
   0.25E+00, 
   0.75E+00, 
   3.00E+00, 
   1.00E+01, 
   0.25E+00, 
   0.75E+00, 
   5.00E+00, 
   2.00E+01, 
   0.25E+00, 
   0.75E+00, 
   5.00E+01, 
   2.00E+02 };

  static double x[]  = {
   1.01E-26, 
   1.01E-26, 
   0.00E+00, 
   0.00E+00, 
   0.00E+00, 
   2.99E+24, 
   0.55E-78, 
   0.55E-78, 
   0.55E-78, 
   2.01E-26, 
   1.01E-26, 
   1.01E-26, 
   0.00E+00, 
   0.00E+00, 
   0.00E+00, 
   2.99E+24, 
   0.55E-78, 
   0.55E-78, 
   0.55E-78, 
   2.01E-26, 
   0.00E+00, 
   2.00E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00, 
   0.50E+00 };

  static double y[] = {
   1.01E-26, 
   1.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.99E+24, 
   2.99E+24, 
   2.01E-26, 
   2.01E-26, 
   2.99E+24, 
   2.01E-26, 
   1.01E-26, 
   1.01E-26, 
   2.01E-26, 
   2.01E-26, 
   2.99E+24, 
   2.99E+24, 
   2.01E-26, 
   2.01E-26, 
   2.99E+24, 
   2.01E-26, 
   1.90E-26, 
   3.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00 };

  static double z[] = {
   1.01E-26, 
   2.99E+24, 
   2.01E-26, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.01E-26, 
   2.99E+24, 
   2.99E+24, 
   2.01E-26, 
   1.01E-26, 
   2.99E+24, 
   2.01E-26, 
   2.99E+24, 
   2.99E+24, 
   2.99E+24, 
   2.01E-26, 
   2.99E+24, 
   2.99E+24, 
   2.01E-26, 
   1.90E-26, 
   4.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   1.00E+00, 
   2.00E+00, 
   2.00E+00, 
   2.00E+00, 
   2.00E+00, 
   5.00E+00, 
   5.00E+00, 
   5.00E+00, 
   5.00E+00, 
   1.00E+01, 
   1.00E+01, 
   1.00E+01, 
   1.00E+01, 
   1.00E+02, 
   1.00E+02, 
   1.00E+02, 
   1.00E+02 };

  cout << "\n";
  cout << "RJ_TEST\n";
  cout << "  RJ evaluates the Carlson elliptic integral\n";
  cout << "  of the third kind, RJ(X,Y,Z,P)\n";
  cout << "\n";
  cout << "               X                          Y";
  cout << "                          Z                         P";
  cout << "                         RJ(X,Y,Z,P)\n";
  cout << "\n";

  errtol = 1.0E-3;

  for ( i = 0; i < 42; i++ )
  {
    eliptc = rj ( x[i], y[i], z[i], p[i], errtol, ierr );
    cout << setw(27) << setprecision(16) << x[i]
         << setw(27) << setprecision(16) << y[i]
         << setw(27) << setprecision(16) << z[i]
         << setw(27) << setprecision(16) << p[i];
    if (ierr == 0 )
    {
      cout << setw(27) << setprecision(16) << eliptc << "\n";
    }
    else
    {
      cout << "  ***Error***\n";
    }
  }

  return;
}
//****************************************************************************80

void elliptic_ea_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_EA_TEST tests ELLIPTIC_EA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;

  cout << "\n";
  cout << "ELLIPTIC_EA_TEST:\n";
  cout << "  ELLIPTIC_EA returns values of\n";
  cout << "  the complete elliptic integral of the\n";
  cout << "  second kind, with parameter angle A.\n";
  cout << "\n";
  cout << "      A       E(A)          E(A)\n";
  cout << "          Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_ea_values ( n_data, a, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_ea ( a );

    cout << "  " << setw(14) << setprecision(6) << a;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void elliptic_ek_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_EK_TEST tests ELLIPTIC_EK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double k;
  int n_data;

  cout << "\n";
  cout << "ELLIPTIC_EK_TEST:\n";
  cout << "  ELLIPTIC_EK returns values of\n";
  cout << "  the complete elliptic integral of the\n";
  cout << "  second kind, with parameter K.\n";
  cout << "\n";
  cout << "      K       E(K)          E(K)\n";
  cout << "          Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_ek_values ( n_data, k, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_ek ( k );

    cout << "  " << setw(14) << setprecision(6) << k;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_em_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_EM_TEST tests ELLIPTIC_EM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double m;
  int n_data;

  cout << "\n";
  cout << "ELLIPTIC_EM_TEST:\n";
  cout << "  ELLIPTIC_EM returns values of\n";
  cout << "  the complete elliptic integral of the\n";
  cout << "  second kind, with parameter M.\n";
  cout << "\n";
  cout << "      M       E(M)          E(M)\n";
  cout << "          Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_em_values ( n_data, m, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_em ( m );

    cout << "  " << setw(14) << setprecision(6) << m;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_fa_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_FA_TEST tests ELLIPTIC_FA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;

  cout << "\n";
  cout << "ELLIPTIC_FA_TEST:\n";
  cout << "  ELLIPTIC_FA returns values of\n";
  cout << "  the complete elliptic integral of the first\n";
  cout << "  kind, with parameter angle A.\n";
  cout << "\n";
  cout << "      A       F(A)          F(A)\n";
  cout << "          Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_fa_values ( n_data, a, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_fa ( a );

    cout << "  " << setw(14) << setprecision(6) << a;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_fk_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_FK_TEST tests ELLIPTIC_FK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double k;
  int n_data;

  cout << "\n";
  cout << "ELLIPTIC_FK_TEST:\n";
  cout << "  ELLIPTIC_FK returns values of\n";
  cout << "  the complete elliptic integral of the first\n";
  cout << "  kind, with parameter K.\n";
  cout << "\n";
  cout << "      K       F(K)          F(K)\n";
  cout << "          Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_fk_values ( n_data, k, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_fk ( k );

    cout << "  " << setw(14) << setprecision(6) << k;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_fm_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_FM_TEST tests ELLIPTIC_FM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double m;
  int n_data;

  cout << "\n";
  cout << "ELLIPTIC_FM_TEST:\n";
  cout << "  ELLIPTIC_FM returns values of\n";
  cout << "  the complete elliptic integral of the first\n";
  cout << "  kind, with parameter M.\n";
  cout << "\n";
  cout << "      M       F(M)          F(M)\n";
  cout << "          Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_fm_values ( n_data, m, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_fm ( m );

    cout << "  " << setw(14) << setprecision(6) << m;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_pia_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_PIA_TEST tests ELLIPTIC_PIA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double n;
  int n_data;
  double pia1;
  double pia2;

  cout << "\n";
  cout << "ELLIPTIC_PIA_TEST:\n";
  cout << "  ELLIPTIC_PIA returns values of\n";
  cout << "  the complete elliptic integral of the\n";
  cout << "  third kind, with parameter angle A.\n";
  cout << "\n";
  cout << "      N     A   Pi(N,A)      Pi(N,A)\n";
  cout << "                Tabulated    Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_pia_values ( n_data, n, a, pia1 );

    if ( n_data == 0 )
    {
      break;
    }

    pia2 = elliptic_pia ( n, a );

    cout << "  " << setw(14) << setprecision(6) << n;
    cout << "  " << setw(14) << setprecision(6) << a;
    cout << "  " << setw(24) << setprecision(16) << pia1;
    cout << "  " << setw(24) << setprecision(16) << pia2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_pik_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_PIK_TEST tests ELLIPTIC_PIK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double k;
  double n;
  int n_data;
  double pik1;
  double pik2;

  cout << "\n";
  cout << "ELLIPTIC_PIK_TEST:\n";
  cout << "  ELLIPTIC_PIK returns values of\n";
  cout << "  the complete elliptic integral of the\n";
  cout << "  third kind, with parameter K.\n";
  cout << "\n";
  cout << "      N     K    Pi(N,K)           Pi(N,K)\n";
  cout << "                 Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_pik_values ( n_data, n, k, pik1 );

    if ( n_data == 0 )
    {
      break;
    }

    pik2 = elliptic_pik ( n, k );

    cout << "  " << setw(14) << setprecision(6) << n;
    cout << "  " << setw(14) << setprecision(6) << k;
    cout << "  " << setw(24) << setprecision(16) << pik1;
    cout << "  " << setw(24) << setprecision(16) << pik2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_pim_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_PIM_TEST tests ELLIPTIC_PIM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double m;
  double n;
  int n_data;
  double pim1;
  double pim2;

  cout << "\n";
  cout << "ELLIPTIC_PIM_TEST:\n";
  cout << "  ELLIPTIC_PIM returns values of\n";
  cout << "  the complete elliptic integral of the\n";
  cout << "  third kind, with parameter modulus M.\n";
  cout << "\n";
  cout << "      N     M    Pi(N,M)           Pi(N,M)\n";
  cout << "                 Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_pim_values ( n_data, n, m, pim1 );

    if ( n_data == 0 )
    {
      break;
    }

    pim2 = elliptic_pim ( n, m );

    cout << "  " << setw(14) << setprecision(6) << n;
    cout << "  " << setw(14) << setprecision(6) << m;
    cout << "  " << setw(24) << setprecision(16) << pim1;
    cout << "  " << setw(24) << setprecision(16) << pim2 << "\n";
  }

  return;
}
//****************************************************************************80

void elliptic_inc_ea_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_EA_TEST tests ELLIPTIC_INC_EA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;
  double phi;

  cout << "\n";
  cout << "ELLIPTIC_INC_EA_TEST:\n";
  cout << "  ELLIPTIC_INC_EA returns values of\n";
  cout << "  the incomplete elliptic integral of the\n";
  cout << "  second kind, with parameters PHI, A.\n";
  cout << "\n";
  cout << "      PHI             A       E(PHI,A)          E(PHI,A)\n";
  cout << "                              Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_ea_values ( n_data, phi, a, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_inc_ea ( phi, a );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << a;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void elliptic_inc_ek_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_EK_TEST tests ELLIPTIC_INC_EK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double k;
  int n_data;
  double phi;

  cout << "\n";
  cout << "ELLIPTIC_INC_EK_TEST:\n";
  cout << "  ELLIPTIC_INC_EK returns values of\n";
  cout << "  the incomplete elliptic integral of the\n";
  cout << "  second kind, with parameters PHI, K.\n";
  cout << "\n";
  cout << "      PHI             K       E(PHI,K)          E(PHI,K)\n";
  cout << "                              Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_ek_values ( n_data, phi, k, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_inc_ek ( phi, k );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << k;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_inc_em_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_EM_TEST tests ELLIPTIC_INC_EM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double m;
  int n_data;
  double phi;

  cout << "\n";
  cout << "ELLIPTIC_INC_EM_TEST:\n";
  cout << "  ELLIPTIC_INC_EM returns values of\n";
  cout << "  the incomplete elliptic integral of the\n";
  cout << "  second kind, with parameters PHI, M.\n";
  cout << "\n";
  cout << "      PHI             M       E(PHI,M)          E(PHI,M)\n";
  cout << "                              Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_em_values ( n_data, phi, m, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_inc_em ( phi, m );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << m;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_inc_fa_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_FA_TEST tests ELLIPTIC_INC_FA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;
  double phi;

  cout << "\n";
  cout << "ELLIPTIC_INC_FA_TEST:\n";
  cout << "  ELLIPTIC_INC_FA returns values of\n";
  cout << "  the incomplete elliptic integral of the first\n";
  cout << "  kind, with parameters PHI, A.\n";
  cout << "\n";
  cout << "      PHI             A       F(PHI,A)          F(PHI,A)\n";
  cout << "                              Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_fa_values ( n_data, phi, a, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_inc_fa ( phi, a );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << a;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_inc_fk_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_FK_TEST tests ELLIPTIC_INC_FK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double k;
  int n_data;
  double phi;

  cout << "\n";
  cout << "ELLIPTIC_INC_FK_TEST:\n";
  cout << "  ELLIPTIC_INC_FK returns values of\n";
  cout << "  the incomplete elliptic integral of the first\n";
  cout << "  kind, with parameters PHI, K.\n";
  cout << "\n";
  cout << "      PHI             K       F(PHI,K)          F(PHI,K)\n";
  cout << "                              Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_fk_values ( n_data, phi, k, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_inc_fk ( phi, k );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << k;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_inc_fm_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_FM_TEST tests ELLIPTIC_INC_FM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double m;
  int n_data;
  double phi;

  cout << "\n";
  cout << "ELLIPTIC_INC_FM_TEST:\n";
  cout << "  ELLIPTIC_INC_FM returns values of\n";
  cout << "  the incomplete elliptic integral of the first\n";
  cout << "  kind, with parameters PHI, M.\n";
  cout << "\n";
  cout << "      PHI             M       F(PHI,M)          F(PHI,M)\n";
  cout << "                              Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_fm_values ( n_data, phi, m, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = elliptic_inc_fm ( phi, m );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << m;
    cout << "  " << setw(24) << setprecision(16) << fx;
    cout << "  " << setw(24) << setprecision(16) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_inc_pia_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_PIA_TEST tests ELLIPTIC_INC_PIA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double n;
  int n_data;
  double phi;
  double pia1;
  double pia2;

  cout << "\n";
  cout << "ELLIPTIC_INC_PIA_TEST:\n";
  cout << "  ELLIPTIC_INC_PIA returns values of\n";
  cout << "  the incomplete elliptic integral of the\n";
  cout << "  third kind, with parameters PHI, N, A.\n";
  cout << "\n";
  cout << "      PHI             N     A   Pi(PHI,N,A)  Pi(PHI,N,A)\n";
  cout << "                                Tabulated    Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_pia_values ( n_data, phi, n, a, pia1 );

    if ( n_data == 0 )
    {
      break;
    }

    pia2 = elliptic_inc_pia ( phi, n, a );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << n;
    cout << "  " << setw(14) << setprecision(6) << a;
    cout << "  " << setw(24) << setprecision(16) << pia1;
    cout << "  " << setw(24) << setprecision(16) << pia2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_inc_pik_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_PIK_TEST tests ELLIPTIC_INC_PIK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double k;
  double n;
  int n_data;
  double phi;
  double pik1;
  double pik2;

  cout << "\n";
  cout << "ELLIPTIC_INC_PIK_TEST:\n";
  cout << "  ELLIPTIC_INC_PIK returns values of\n";
  cout << "  the incomplete elliptic integral of the\n";
  cout << "  third kind, with parameters PHI, N, K.\n";
  cout << "\n";
  cout << "      PHI             N     K    Pi(PHI,N,K)       Pi(PHI,N,K)\n";
  cout << "                                 Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_pik_values ( n_data, phi, n, k, pik1 );

    if ( n_data == 0 )
    {
      break;
    }

    pik2 = elliptic_inc_pik ( phi, n, k );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << n;
    cout << "  " << setw(14) << setprecision(6) << k;
    cout << "  " << setw(24) << setprecision(16) << pik1;
    cout << "  " << setw(24) << setprecision(16) << pik2 << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_inc_pim_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_INC_PIM_TEST tests ELLIPTIC_INC_PIM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double m;
  double n;
  int n_data;
  double phi;
  double pim1;
  double pim2;

  cout << "\n";
  cout << "ELLIPTIC_INC_PIM_TEST:\n";
  cout << "  ELLIPTIC_INC_PIM returns values of\n";
  cout << "  the incomplete elliptic integral of the\n";
  cout << "  third kind, with parameters PHI, N, M.\n";
  cout << "\n";
  cout << "      PHI             N     M    Pi(PHI,N,M)       Pi(PHI,N,M)\n";
  cout << "                                 Tabulated         Calculated\n";
  cout << "\n";

  n_data = 0;
 
  while ( true )
  {
    elliptic_inc_pim_values ( n_data, phi, n, m, pim1 );

    if ( n_data == 0 )
    {
      break;
    }

    pim2 = elliptic_inc_pim ( phi, n, m );

    cout << "  " << setw(14) << setprecision(6) << phi;
    cout << "  " << setw(14) << setprecision(6) << n;
    cout << "  " << setw(14) << setprecision(6) << m;
    cout << "  " << setw(24) << setprecision(16) << pim1;
    cout << "  " << setw(24) << setprecision(16) << pim2 << "\n";
  }

  return;
}
//****************************************************************************80

void jacobi_cn_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_CN_TEST tests JACOBI_CN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double cn1;
  double cn2;
  double m;
  int n_data;
  double u;

  cout << "\n";
  cout << "JACOBI_CN_TEST:\n";
  cout << "  JACOBI_CN evaluates the Jacobi elliptic function CN.\n";
  cout << "\n";
  cout << "    U       M       Exact CN                CN(U,M)\n";
  cout << "\n";

  n_data = 0;

  while ( true )
  {
    jacobi_cn_values ( n_data, u, m, cn1 );

    if ( n_data == 0 )
    {
      break;
    }

    cn2 = jacobi_cn ( u, m );

    cout << setw(8) << setprecision(4) << u
         << setw(8) << setprecision(4) << m
         << setw(24) << setprecision(16) << cn1
         << setw(24) << setprecision(16) << cn2 << "\n";
  }

  return;
}
//****************************************************************************80

void jacobi_dn_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_DN_TEST tests JACOBI_DN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double dn1;
  double dn2;
  double m;
  int n_data;
  double u;

  cout << "\n";
  cout << "JACOBI_DN_TEST:\n";
  cout << "  JACOBI_DN evaluates the Jacobi elliptic function DN.\n";
  cout << "\n";
  cout << "    U       M       Exact DN                DN(U,M)\n";
  cout << "\n";

  n_data = 0;

  while ( true )
  {
    jacobi_dn_values ( n_data, u, m, dn1 );

    if ( n_data == 0 )
    {
      break;
    }

    dn2 = jacobi_dn ( u, m );

    cout << setw(8) << setprecision(4) << u
         << setw(8) << setprecision(4) << m
         << setw(24) << setprecision(16) << dn1
         << setw(24) << setprecision(16) << dn2 << "\n";
  }

  return;
}
//****************************************************************************80

void jacobi_sn_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SN_TEST tests JACOBI_SN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2018
//
//  Author:
//
//    John Burkardt
//
{
  double m;
  int n_data;
  double sn1;
  double sn2;
  double u;

  cout << "\n";
  cout << "JACOBI_SN_TEST:\n";
  cout << "  JACOBI_SN evaluates the Jacobi elliptic function SN.\n";
  cout << "\n";
  cout << "    U       M       Exact SN                SN(U,M)\n";
  cout << "\n";

  n_data = 0;

  while ( true )
  {
    jacobi_sn_values ( n_data, u, m, sn1 );

    if ( n_data == 0 )
    {
      break;
    }

    sn2 = jacobi_sn ( u, m );

    cout << setw(8) << setprecision(4) << u
         << setw(8) << setprecision(4) << m
         << setw(24) << setprecision(16) << sn1
         << setw(24) << setprecision(16) << sn2 << "\n";
  }

  return;
}
