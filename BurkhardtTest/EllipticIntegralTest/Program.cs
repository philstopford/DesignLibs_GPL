using System;
using Burkardt.Elliptic;

namespace EllipticIntegralTest;

using Integral = Integral;

internal class Program
{
    private static void Main(string[] args)
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
        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INTEGRAL_TEST");
        Console.WriteLine("  ELLIPTIC_INTEGRAL evaluates elliptic integral functions");
        Console.WriteLine("  using Carlson's elliptic functions.");
        //
        //  Carloson integrals.
        //
        rc_test();
        rc_test2();
        rd_test();
        rf_test();
        rj_test();
        //
        //  Complete elliptic integrals, first kind.
        //
        elliptic_fa_test();
        elliptic_fk_test();
        elliptic_fm_test();
        //
        //  Complete elliptic integrals, second kind.
        //
        elliptic_ea_test();
        elliptic_ek_test();
        elliptic_em_test();
        //
        //  Complete elliptic integrals, third kind.
        //
        elliptic_pia_test();
        elliptic_pik_test();
        elliptic_pim_test();
        //
        //  Incomplete elliptic integrals, first kind.
        //
        elliptic_inc_fa_test();
        elliptic_inc_fk_test();
        elliptic_inc_fm_test();
        //
        //  Incomplete elliptic integrals, second kind.
        //
        elliptic_inc_ea_test();
        elliptic_inc_ek_test();
        elliptic_inc_em_test();
        //
        //  Incomplete elliptic integrals, third kind.
        //
        elliptic_inc_pia_test();
        elliptic_inc_pik_test();
        elliptic_inc_pim_test();
        //
        //  Jacobi elliptic functions.
        //
        jacobi_cn_test();
        jacobi_dn_test();
        jacobi_sn_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INTEGRAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void rc_test()

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
        //    UPLIM, X, Y, and ERRTOL (see comments in static void) may be used on
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
        //    C++ version by John Burkardt.
        //
    {
        double eliptc;
        double errtol;
        int i;
        int ierr = 0;
        double[] x =
            {
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
                1.00E+50
            }
            ;

        double[] y =
            {
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
                1.00E+00
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("RC_TEST");
        Console.WriteLine("  RC evaluates the elementary integral RC(X,Y)");

        Console.WriteLine("");
        Console.WriteLine("              X                          Y                         RC(X,Y)");
        Console.WriteLine("");

        errtol = 1.0E-3;

        for (i = 0; i < 43; i++)
        {
            eliptc = Integral.rc(x[i], y[i], errtol, ref ierr);
            string cout = "  " + x[i].ToString("0.################").PadLeft(27)
                               + "  " + y[i].ToString("0.################").PadLeft(27);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + eliptc.ToString("0.################").PadLeft(27) + "");
                    break;
                default:
                    Console.WriteLine(cout + "  ***Error***");
                    break;
            }
        }
    }

    private static void rc_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RC_TEST2 checks RC by examining special values.
        //
        //  Discussion:
        //
        //    This driver compares values of (LOG X)/(X-1) and ARCTAN(X)
        //    calculated on one hand from the static void RC and on the other
        //    from library LOG and ARCTAN routines.  to astatic void over/underflows
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
        //    C++ version by John Burkardt.
        //
    {
        double errtol;
        int i;
        double ibmarc;
        double ibmlog;
        int ierr = 0;
        int ipower;
        int j;
        int m;
        double myarc;
        double mylog;
        double v;
        double w;
        double x;
        double[] x_vec =
            {
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
                1.0E+75
            }
            ;
        double y;
        double z;

        Console.WriteLine("");
        Console.WriteLine("RC_TEST2");
        Console.WriteLine("  Compare LOG(X)/(X-1) and ARCTAN(X) with");
        Console.WriteLine("  values based on RC.");

        Console.WriteLine("");
        Console.WriteLine("     X                From LOG                   From RC");
        Console.WriteLine("");

        errtol = 1.0E-3;

        for (j = 1; j <= 10; j++)
        {
            x = 0.2 * j;
            y = (1.0 + x) / 2.0;
            v = x / y;
            mylog = Integral.rc(y, v, errtol, ref ierr) / Math.Sqrt(y);
            string cout = x.ToString("0.#").PadLeft(9) + "     ";
            switch (j)
            {
                case 5:
                    cout += "**** ZERO DIVIDE *****";
                    break;
                default:
                    ibmlog = Math.Log(x) / (x - 1.0);
                    cout += ibmlog.ToString("0.################").PadLeft(27);
                    break;
            }

            Console.WriteLine(cout + mylog.ToString("0.################").PadLeft(27) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Extreme values of X");
        Console.WriteLine("");
        Console.WriteLine("     X                From LOG                   From RC");
        Console.WriteLine("");

        for (i = 0; i < 16; i++)
        {
            ipower = -75 + 10 * i;
            x = Math.Pow(10.0, ipower);
            y = (1.0 + x) / 2.0;
            v = x / y;
            mylog = Integral.rc(y, v, errtol, ref ierr) / Math.Sqrt(y);
            ibmlog = Math.Log(x) / (x - 1.0);
            Console.WriteLine(x.ToString("0.#").PadLeft(8)
                              + ibmlog.ToString("0.################").PadLeft(27)
                              + mylog.ToString("0.################").PadLeft(27) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("     X              From ARCTAN                 From RC");
        Console.WriteLine("");

        for (m = 0; m < 13; m++)
        {
            x = x_vec[m];
            z = 1.0 / x;
            w = z + x;
            myarc = Math.Sqrt(x) * Integral.rc(z, w, errtol, ref ierr);
            ibmarc = Math.Atan(x);
            Console.WriteLine(x.ToString("0.#").PadLeft(8)
                              + ibmarc.ToString("0.################").PadLeft(27)
                              + myarc.ToString("0.################").PadLeft(27) + "");
        }
    }

    private static void rd_test()

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
        //    comments in static void) may be used on most machines but provide a
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
        //    C++ version by John Burkardt.
        //
    {
        double eliptc;
        double errtol;
        int i;
        int ierr = 0;
        double[] x =
            {
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
                0.50E+00
            }
            ;

        double[] y =
            {
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
                1.00E+00
            }
            ;

        double[] z =
            {
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
                1.00E+10
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("RD_TEST");
        Console.WriteLine("  RD evaluates the Carlson elliptic integral");
        Console.WriteLine("  of the second kind, RD(X,Y,Z)");
        Console.WriteLine("");
        Console.WriteLine("               X                          Y" + 
                          "                          Z                         RD(X,Y,Z)");
        Console.WriteLine("");

        errtol = 1.0E-03;

        for (i = 0; i < 27; i++)
        {
            eliptc = Integral.rd(x[i], y[i], z[i], errtol, ref ierr);
            string cout = x[i].ToString("0.################").PadLeft(27)
                          + y[i].ToString("0.################").PadLeft(27)
                          + z[i].ToString("0.################").PadLeft(27);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + eliptc.ToString("0.################").PadLeft(27) + "");
                    break;
                default:
                    Console.WriteLine(cout + "  ***Error***");
                    break;
            }
        }
    }

    private static void rf_test()

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
        //    comments in static void) may be used on most machines but provide a
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
        //    C++ version by John Burkardt.
        //
    {
        double eliptc;
        double errtol;
        int i;
        int ierr = 0;
        double[] x =
            {
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
                0.50E+00
            }
            ;

        double[] y =
            {
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
                1.00E+00
            }
            ;

        double[] z =
            {
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
                1.00E+50
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("RF_TEST");
        Console.WriteLine("  RF evaluates the Carlson elliptic integral");
        Console.WriteLine("  of the first kind, RF(X,Y,Z)");
        Console.WriteLine("");
        Console.WriteLine("               X                          Y" +
                          "                          Z                         RF(X,Y,Z)");
        Console.WriteLine(" ");

        errtol = 1.0E-3;

        for (i = 0; i < 55; i++)
        {
            eliptc = Integral.rf(x[i], y[i], z[i], errtol, ref ierr);
            string cout = x[i].ToString("0.################").PadLeft(27)
                          + y[i].ToString("0.################").PadLeft(27)
                          + z[i].ToString("0.################").PadLeft(27);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + eliptc.ToString("0.################").PadLeft(27) + "");
                    break;
                default:
                    Console.WriteLine("  ***Error***");
                    break;
            }
        }
    }

    private static void rj_test()

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
        //    ERRTOL (see comments in static void) may be used on most machines
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
        //    C++ version by John Burkardt.
        //
    {
        double eliptc;
        double errtol;
        int i;
        int ierr = 0;

        double[] p =
            {
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
                2.00E+02
            }
            ;

        double[] x =
            {
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
                0.50E+00
            }
            ;

        double[] y =
            {
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
                1.00E+00
            }
            ;

        double[] z =
            {
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
                1.00E+02
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("RJ_TEST");
        Console.WriteLine("  RJ evaluates the Carlson elliptic integral");
        Console.WriteLine("  of the third kind, RJ(X,Y,Z,P)");
        Console.WriteLine("");
        Console.WriteLine("               X                          Y" +
                          "                          Z                         P" + 
                          "                         RJ(X,Y,Z,P)");
        Console.WriteLine("");

        errtol = 1.0E-3;

        for (i = 0; i < 42; i++)
        {
            eliptc = Integral.rj(x[i], y[i], z[i], p[i], errtol, ref ierr);
            string cout = x[i].ToString("0.################").PadLeft(27)
                          + y[i].ToString("0.################").PadLeft(27)
                          + z[i].ToString("0.################").PadLeft(27)
                          + p[i].ToString("0.################").PadLeft(27);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + eliptc.ToString("0.################").PadLeft(27) + "");
                    break;
                default:
                    Console.WriteLine(cout + "  ***Error***");
                    break;
            }
        }
    }

    private static void elliptic_ea_test()

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
        double a = 0;
        double fx = 0;
        double fx2;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_EA_TEST:");
        Console.WriteLine("  ELLIPTIC_EA returns values of");
        Console.WriteLine("  the complete elliptic integral of the");
        Console.WriteLine("  second kind, with parameter angle A.");
        Console.WriteLine("");
        Console.WriteLine("      A       E(A)          E(A)");
        Console.WriteLine("          Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            EA.values(ref n_data, ref a, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = EA.evaluate(a);

            Console.WriteLine("  " + a.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_ek_test()

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
        double fx = 0;
        double fx2;
        double k = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_EK_TEST:");
        Console.WriteLine("  ELLIPTIC_EK returns values of");
        Console.WriteLine("  the complete elliptic integral of the");
        Console.WriteLine("  second kind, with parameter K.");
        Console.WriteLine("");
        Console.WriteLine("      K       E(K)          E(K)");
        Console.WriteLine("          Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            EK.values(ref n_data, ref k, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = EK.evaluate(k);

            Console.WriteLine("  " + k.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_em_test()

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
        double fx = 0;
        double fx2;
        double m = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_EM_TEST:");
        Console.WriteLine("  ELLIPTIC_EM returns values of");
        Console.WriteLine("  the complete elliptic integral of the");
        Console.WriteLine("  second kind, with parameter M.");
        Console.WriteLine("");
        Console.WriteLine("      M       E(M)          E(M)");
        Console.WriteLine("          Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            EM.values(ref n_data, ref m, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = EM.evaluate(m);

            Console.WriteLine("  " + m.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_fa_test()

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
        double a = 0;
        double fx = 0;
        double fx2;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_FA_TEST:");
        Console.WriteLine("  ELLIPTIC_FA returns values of");
        Console.WriteLine("  the complete elliptic integral of the first");
        Console.WriteLine("  kind, with parameter angle A.");
        Console.WriteLine("");
        Console.WriteLine("      A       F(A)          F(A)");
        Console.WriteLine("          Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            FA.values(ref n_data, ref a, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FA.evaluate(a);

            Console.WriteLine("  " + a.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_fk_test()

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
        double fx = 0;
        double fx2;
        double k = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_FK_TEST:");
        Console.WriteLine("  ELLIPTIC_FK returns values of");
        Console.WriteLine("  the complete elliptic integral of the first");
        Console.WriteLine("  kind, with parameter K.");
        Console.WriteLine("");
        Console.WriteLine("      K       F(K)          F(K)");
        Console.WriteLine("          Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            FK.values(ref n_data, ref k, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FK.evaluate(k);

            Console.WriteLine("  " + k.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_fm_test()

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
        double fx = 0;
        double fx2;
        double m = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_FM_TEST:");
        Console.WriteLine("  ELLIPTIC_FM returns values of");
        Console.WriteLine("  the complete elliptic integral of the first");
        Console.WriteLine("  kind, with parameter M.");
        Console.WriteLine("");
        Console.WriteLine("      M       F(M)          F(M)");
        Console.WriteLine("          Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            FM.values(ref n_data, ref m, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FM.evaluate(m);

            Console.WriteLine("  " + m.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }
    //****************************************************************************80

    private static void elliptic_pia_test()

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
        double a = 0;
        double n = 0;
        int n_data;
        double pia1 = 0;
        double pia2;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_PIA_TEST:");
        Console.WriteLine("  ELLIPTIC_PIA returns values of");
        Console.WriteLine("  the complete elliptic integral of the");
        Console.WriteLine("  third kind, with parameter angle A.");
        Console.WriteLine("");
        Console.WriteLine("      N     A   Pi(N,A)      Pi(N,A)");
        Console.WriteLine("                Tabulated    Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            PIA.values(ref n_data, ref n, ref a, ref pia1);

            if (n_data == 0)
            {
                break;
            }

            pia2 = PIA.evaluate(n, a);

            Console.WriteLine("  " + n.ToString("0.######").PadLeft(14) + 
                              "  " + a.ToString("0.######").PadLeft(14) + 
                              "  " + pia1.ToString("0.################").PadLeft(24) + 
                              "  " + pia2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_pik_test()

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
        double k = 0;
        double n = 0;
        int n_data;
        double pik1 = 0;
        double pik2 = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_PIK_TEST:");
        Console.WriteLine("  ELLIPTIC_PIK returns values of");
        Console.WriteLine("  the complete elliptic integral of the");
        Console.WriteLine("  third kind, with parameter K.");
        Console.WriteLine("");
        Console.WriteLine("      N     K    Pi(N,K)           Pi(N,K)");
        Console.WriteLine("                 Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            PIK.values(ref n_data, ref n, ref k, ref pik1);

            if (n_data == 0)
            {
                break;
            }

            pik2 = PIK.evaluate(n, k);

            Console.WriteLine("  " + n.ToString("0.######").PadLeft(14) + 
                              "  " + k.ToString("0.######").PadLeft(14) + 
                              "  " + pik1.ToString("0.################").PadLeft(24) + 
                              "  " + pik2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_pim_test()

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
        double m = 0;
        double n = 0;
        int n_data;
        double pim1 = 0;
        double pim2 = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_PIM_TEST:");
        Console.WriteLine("  ELLIPTIC_PIM returns values of");
        Console.WriteLine("  the complete elliptic integral of the");
        Console.WriteLine("  third kind, with parameter modulus M.");
        Console.WriteLine("");
        Console.WriteLine("      N     M    Pi(N,M)           Pi(N,M)");
        Console.WriteLine("                 Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            PIM.values(ref n_data, ref n, ref m, ref pim1);

            if (n_data == 0)
            {
                break;
            }

            pim2 = PIM.evaluate(n, m);

            Console.WriteLine("  " + n.ToString("0.######").PadLeft(14) + 
                              "  " + m.ToString("0.######").PadLeft(14) + 
                              "  " + pim1.ToString("0.################").PadLeft(24) + 
                              "  " + pim2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_ea_test()

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
        double a = 0;
        double fx = 0;
        double fx2 = 0;
        int n_data;
        double phi = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_EA_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_EA returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the");
        Console.WriteLine("  second kind, with parameters PHI, A.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             A       E(PHI,A)          E(PHI,A)");
        Console.WriteLine("                              Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            EA_inc.values(ref n_data, ref phi, ref a, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = EA_inc.evaluate(phi, a);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + a.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_ek_test()

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
        double fx = 0;
        double fx2 = 0;
        double k = 0;
        int n_data;
        double phi = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_EK_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_EK returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the");
        Console.WriteLine("  second kind, with parameters PHI, K.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             K       E(PHI,K)          E(PHI,K)");
        Console.WriteLine("                              Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            EK_inc.values(ref n_data, ref phi, ref k, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = EK_inc.evaluate(phi, k);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + k.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_em_test()

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
        double fx = 0;
        double fx2 = 0;
        double m = 0;
        int n_data;
        double phi = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_EM_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_EM returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the");
        Console.WriteLine("  second kind, with parameters PHI, M.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             M       E(PHI,M)          E(PHI,M)");
        Console.WriteLine("                              Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            EM_inc.values(ref n_data, ref phi, ref m, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = EM_inc.evaluate(phi, m);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + m.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_fa_test()

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
        double a = 0;
        double fx = 0;
        double fx2 = 0;
        int n_data;
        double phi = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_FA_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_FA returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the first");
        Console.WriteLine("  kind, with parameters PHI, A.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             A       F(PHI,A)          F(PHI,A)");
        Console.WriteLine("                              Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            FA_inc.values(ref n_data, ref phi, ref a, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FA_inc.evaluate(phi, a);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + a.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_fk_test()

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
        double fx = 0;
        double fx2 = 0;
        double k = 0;
        int n_data;
        double phi = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_FK_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_FK returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the first");
        Console.WriteLine("  kind, with parameters PHI, K.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             K       F(PHI,K)          F(PHI,K)");
        Console.WriteLine("                              Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            FK_inc.values(ref n_data, ref phi, ref k, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FK_inc.evaluate(phi, k);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + k.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_fm_test()

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
        double fx = 0;
        double fx2 = 0;
        double m = 0;
        int n_data;
        double phi = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_FM_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_FM returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the first");
        Console.WriteLine("  kind, with parameters PHI, M.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             M       F(PHI,M)          F(PHI,M)");
        Console.WriteLine("                              Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            FM_inc.values(ref n_data, ref phi, ref m, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FM_inc.evaluate(phi, m);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + m.ToString("0.######").PadLeft(14) + 
                              "  " + fx.ToString("0.################").PadLeft(24) + 
                              "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_pia_test()

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
        double a = 0;
        double n = 0;
        int n_data;
        double phi = 0;
        double pia1 = 0;
        double pia2 = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_PIA_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_PIA returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the");
        Console.WriteLine("  third kind, with parameters PHI, N, A.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             N     A   Pi(PHI,N,A)  Pi(PHI,N,A)");
        Console.WriteLine("                                Tabulated    Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            PIA_inc.values(ref n_data, ref phi, ref n, ref a, ref pia1);

            if (n_data == 0)
            {
                break;
            }

            pia2 = PIA_inc.evaluate(phi, n, a);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + n.ToString("0.######").PadLeft(14) + 
                              "  " + a.ToString("0.######").PadLeft(14) + 
                              "  " + pia1.ToString("0.################").PadLeft(24) + 
                              "  " + pia2.ToString("0.################").PadLeft(24) + "");

        }
    }

    private static void elliptic_inc_pik_test()

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
        double k = 0;
        double n = 0;
        int n_data;
        double phi = 0;
        double pik1 = 0;
        double pik2 = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_PIK_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_PIK returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the");
        Console.WriteLine("  third kind, with parameters PHI, N, K.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             N     K    Pi(PHI,N,K)       Pi(PHI,N,K)");
        Console.WriteLine("                                 Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            PIK_inc.values(ref n_data, ref phi, ref n, ref k, ref pik1);

            if (n_data == 0)
            {
                break;
            }

            pik2 = PIK_inc.evaluate(phi, n, k);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + n.ToString("0.######").PadLeft(14) + 
                              "  " + k.ToString("0.######").PadLeft(14) + 
                              "  " + pik1.ToString("0.################").PadLeft(24) + 
                              "  " + pik2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void elliptic_inc_pim_test()

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
        double m = 0;
        double n = 0;
        int n_data;
        double phi = 0;
        double pim1 = 0;
        double pim2 = 0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPTIC_INC_PIM_TEST:");
        Console.WriteLine("  ELLIPTIC_INC_PIM returns values of");
        Console.WriteLine("  the incomplete elliptic integral of the");
        Console.WriteLine("  third kind, with parameters PHI, N, M.");
        Console.WriteLine("");
        Console.WriteLine("      PHI             N     M    Pi(PHI,N,M)       Pi(PHI,N,M)");
        Console.WriteLine("                                 Tabulated         Calculated");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            PIM_inc.values(ref n_data, ref phi, ref n, ref m, ref pim1);

            if (n_data == 0)
            {
                break;
            }

            pim2 = PIM_inc.evaluate(phi, n, m);

            Console.WriteLine("  " + phi.ToString("0.######").PadLeft(14) + 
                              "  " + n.ToString("0.######").PadLeft(14) + 
                              "  " + m.ToString("0.######").PadLeft(14) + 
                              "  " + pim1.ToString("0.################").PadLeft(24) + 
                              "  " + pim2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void jacobi_cn_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    jacobi_cn_test tests jacobi_cn().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 November 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double cn1 = 0;
        double cn2 = 0;
        double k = 0;
        double m = 0;
        int n_data;
        double u = 0;

        Console.WriteLine("");
        Console.WriteLine("jacobi_cn_test:");
        Console.WriteLine("  jacobi_cn() evaluates the Jacobi elliptic function CN.");
        Console.WriteLine("");
        Console.WriteLine("    U       M       Exact CN                CN(U,M)");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            Jacobi_cn.values(ref n_data, ref u, ref a, ref k, ref m, ref cn1);

            if (n_data == 0)
            {
                break;
            }

            cn2 = Jacobi_cn.evaluate(u, m);

            Console.WriteLine("  " + u.ToString("0.####").PadLeft(8) + 
                              "  " + m.ToString("0.####").PadLeft(8) + 
                              "  " + cn1.ToString("0.################").PadLeft(24) + 
                              "  " + cn2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void jacobi_dn_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    jacobi_dn_test tests jacobi_dn().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 November 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double dn1 = 0;
        double dn2 = 0;
        double k = 0;
        double m = 0;
        int n_data;
        double u = 0;

        Console.WriteLine("");
        Console.WriteLine("jacobi_dn_test:");
        Console.WriteLine("  jacobi_dn() evaluates the Jacobi elliptic function DN.");
        Console.WriteLine("");
        Console.WriteLine("    U       M       Exact DN                DN(U,M)");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            Jacobi_dn.values(ref n_data, ref u, ref a, ref k, ref m, ref dn1);

            if (n_data == 0)
            {
                break;
            }

            dn2 = Jacobi_dn.evaluate(u, m);

            Console.WriteLine("  " + u.ToString("0.####").PadLeft(8) + 
                              "  " + m.ToString("0.####").PadLeft(8) + 
                              "  " + dn1.ToString("0.################").PadLeft(24) + 
                              "  " + dn2.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void jacobi_sn_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    jacobi_sn_test tests jacobi_sn().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 November 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double k = 0;
        double m = 0;
        int n_data;
        double sn1 = 0;
        double sn2 = 0;
        double u = 0;

        Console.WriteLine("");
        Console.WriteLine("jacobi_sn_test:");
        Console.WriteLine("  jacobi_sn() evaluates the Jacobi elliptic function SN.");
        Console.WriteLine("");
        Console.WriteLine("    U       M       Exact SN                SN(U,M)");
        Console.WriteLine("");

        n_data = 0;

        while (true)
        {
            Jacobi_sn.values(ref n_data, ref u, ref a, ref k, ref m, ref sn1);

            if (n_data == 0)
            {
                break;
            }

            sn2 = Jacobi_sn.evaluate(u, m);

            Console.WriteLine("  " + u.ToString("0.####").PadLeft(8) + 
                              "  " + m.ToString("0.####").PadLeft(8) + 
                              "  " + sn1.ToString("0.################").PadLeft(24) + 
                              "  " + sn2.ToString("0.################").PadLeft(24) + "");
        }
    }
}