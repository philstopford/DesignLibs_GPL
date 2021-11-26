using System;
using System.Globalization;
using Burkardt.Elliptic;

namespace TOMS577Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS577_TEST.
        //
        //  Discussion:
        //
        //    TOMS577_TEST tests TOMS577.
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
        Console.WriteLine("");
        Console.WriteLine("TOMS577_TEST");
        Console.WriteLine("  TOMS577 evaluates Carlson's elliptic functions.");

        rc_test();
        rc_test2();
        rd_test();
        rf_test();
        rj_test();

        Console.WriteLine("");
        Console.WriteLine("TOMS577_TEST");
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
        //    C++ version by John Burkardt.
        //
    {
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
        };

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
        };

        Console.WriteLine("");
        Console.WriteLine("RC_TEST");
        Console.WriteLine("  RC evaluates the elementary integral RC(X,Y)");

        Console.WriteLine("");
        Console.WriteLine("              X                          Y                         RC(X,Y)");
        Console.WriteLine("");

        const double errtol = 1.0E-3;

        for (i = 0; i < 43; i++)
        {
            double eliptc = Integral.rc(x[i], y[i], errtol, ref ierr);
            string cout = "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                               + "  " + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(27);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + eliptc.ToString(CultureInfo.InvariantCulture).PadLeft(27) + "");
                    break;
                default:
                    Console.WriteLine("  ***Error***");
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
        //    C++ version by John Burkardt.
        //
    {
        int i;
        double ibmlog;
        int ierr = 0;
        int j;
        int m;
        double mylog;
        double v;
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
        };
        double y;

        Console.WriteLine("");
        Console.WriteLine("RC_TEST2");
        Console.WriteLine("  Compare LOG(X)/(X-1) and ARCTAN(X) with");
        Console.WriteLine("  values based on RC.");

        Console.WriteLine("");
        Console.WriteLine("     X                From LOG                   From RC");
        Console.WriteLine("");

        const double errtol = 1.0E-3;

        for (j = 1; j <= 10; j++)
        {
            string cout = "";
            x = 0.2 * j;
            y = (1.0 + x) / 2.0;
            v = x / y;
            mylog = Integral.rc(y, v, errtol, ref ierr) / Math.Sqrt(y);
            cout = x.ToString("0.#").PadLeft(9) + "     ";
            switch (j)
            {
                case 5:
                    cout = "**** ZERO DIVIDE *****";
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
            int ipower = -75 + 10 * i;
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
            double z = 1.0 / x;
            double w = z + x;
            double myarc = Math.Sqrt(x) * Integral.rc(z, w, errtol, ref ierr);
            double ibmarc = Math.Atan(x);
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
        //    C++ version by John Burkardt.
        //
    {
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
        };

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
        };

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
        };

        Console.WriteLine("");
        Console.WriteLine("RD_TEST");
        Console.WriteLine("  RD evaluates the Carlson elliptic integral");
        Console.WriteLine("  of the second kind, RD(X,Y,Z)");
        Console.WriteLine("");
        Console.WriteLine("               X                          Y" +
                          "                          Z                         RD(X,Y,Z)");
        Console.WriteLine("");

        const double errtol = 1.0E-03;

        for (i = 0; i < 27; i++)
        {
            double eliptc = Integral.rd(x[i], y[i], z[i], errtol, ref ierr);
            string cout = x[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                          + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                          + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(27);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + eliptc.ToString(CultureInfo.InvariantCulture).PadLeft(27) + "");
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
        //    C++ version by John Burkardt.
        //
    {
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
        };

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
        };

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
        };

        Console.WriteLine("");
        Console.WriteLine("RF_TEST");
        Console.WriteLine("  RF evaluates the Carlson elliptic integral");
        Console.WriteLine("  of the first kind, RF(X,Y,Z)");
        Console.WriteLine("");
        Console.WriteLine("               X                          Y" +
                          "                          Z                         RF(X,Y,Z)");
        Console.WriteLine(" ");

        const double errtol = 1.0E-3;

        for (i = 0; i < 55; i++)
        {
            string cout = x[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                          + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                          + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(27);
            double eliptc = Integral.rf(x[i], y[i], z[i], errtol, ref ierr);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + " " + eliptc.ToString(CultureInfo.InvariantCulture).PadLeft(27) + "");
                    break;
                default:
                    Console.WriteLine(cout + "  ***Error***");
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
        //    C++ version by John Burkardt.
        //
    {
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
        };

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
        };

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
        };

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
        };

        Console.WriteLine("");
        Console.WriteLine("RJ_TEST");
        Console.WriteLine("  RJ evaluates the Carlson elliptic integral");
        Console.WriteLine("  of the third kind, RJ(X,Y,Z,P)");
        Console.WriteLine("");
        Console.WriteLine("               X                          Y" +
                          "                          Z                         P" +
                          "                         RJ(X,Y,Z,P)");
        Console.WriteLine("");

        const double errtol = 1.0E-3;

        for (i = 0; i < 42; i++)
        {
            double eliptc = Integral.rj(x[i], y[i], z[i], p[i], errtol, ref ierr);
            string cout = x[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                          + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                          + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(27)
                          + p[i].ToString(CultureInfo.InvariantCulture).PadLeft(27);
            switch (ierr)
            {
                case 0:
                    Console.WriteLine(cout + eliptc.ToString(CultureInfo.InvariantCulture).PadLeft(27) + "");
                    break;
                default:
                    Console.WriteLine(cout + "  ***Error***");
                    break;
            }
        }

    }
}