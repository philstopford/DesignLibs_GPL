﻿using System;
using System.Globalization;
using Burkardt;
using Burkardt.SolveNS;

namespace BisectionRCTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BISECTION_RC_TEST.
        //
        //  Discussion:
        //
        //    BISECTION_RC_TEST tests the BISECTION_RC library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BISECTION_RC_TEST:");
        Console.WriteLine("  Test the BISECTION_RC library.");

        test01();
        test02();
        test03();
        test04();
        test05();

        Console.WriteLine("");
        Console.WriteLine("BISECTION_RC_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests BISECTION_RC, evaluating the function in a separate routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double x;
        BisectionRC.BisectionData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Demonstrate BISECTION_RC on a simple example.");
        Console.WriteLine("  The function is evaluated in a separate routine.");

        const double fx_tol = 1.0E-08;
        const double dx_tol = 1.0E-06;
        int it = 0;
        const int it_max = 30;

        double a = 0.0;
        double b = 1.0;
        double fx = 0.0;
        int job = 0;

        Console.WriteLine("");
        Console.WriteLine("     I      X               FX              DX");
        Console.WriteLine("");

        for (;;)
        {
            x = BisectionRC.bisection_rc(ref data, ref a, ref b, fx, ref job);

            if (job < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error return.");
                break;
            }

            it += 1;

            fx = f01(x);

            double dx = it switch
            {
                <= 2 => Math.Abs(b - a),
                _ => 0.5 * Math.Abs(b - a)
            };

            Console.WriteLine("  " + it.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + dx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (Math.Abs(fx) <= fx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Function is small.");
                break;
            }

            if (dx <= dx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Interval is tiny.");
                break;
            }

            if (it_max > it)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("  Reached iteration limit.");
            break;

        }

        Console.WriteLine("");
        Console.WriteLine("  A = " + a.ToString(CultureInfo.InvariantCulture).PadLeft(14) + " F(A) = " + f01(a) + "");
        Console.WriteLine("  X = " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + " F(X) = " + f01(x) + "");
        Console.WriteLine("  B = " + b.ToString(CultureInfo.InvariantCulture).PadLeft(14) + " F(B) = " + f01(b) + "");

    }

    private static double f01(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F01 evaluates the function f(x) = Math.Cos ( x ) - x which is zero around 0.74
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F01, the function value.
        //
    {
        double value = Math.Cos(x) - x;

        return value;
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests BISECTION_RC, evaluating the function within the routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double x;
        BisectionRC.BisectionData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Demonstrate BISECTION_RC on a simple example.");
        Console.WriteLine("  The function is evaluated within this routine.");

        const double fx_tol = 1.0E-09;
        const double dx_tol = 1.0E-09;
        int it = 0;
        const int it_max = 30;

        double a = 0.0;
        double b = 1.0;
        double fx = 0.0;
        int job = 0;

        Console.WriteLine("");
        Console.WriteLine("     I      X               FX              DX");
        Console.WriteLine("");

        for (;;)
        {
            x = BisectionRC.bisection_rc(ref data, ref a, ref b, fx, ref job);

            if (job < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error return.");
                break;
            }

            it += 1;

            fx = Math.Cos(100.0 * x) - 4.0 * Helpers.Erf(30.0 * x - 10.0);

            double dx = it switch
            {
                <= 2 => Math.Abs(b - a),
                _ => 0.5 * Math.Abs(b - a)
            };

            Console.WriteLine("  " + it.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + dx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (Math.Abs(fx) <= fx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Function is small.");
                break;
            }

            if (dx <= dx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Interval is tiny.");
                break;
            }

            if (it_max > it)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("  Reached iteration limit.");
            break;

        }

        Console.WriteLine("");
        fx = Math.Cos(100.0 * a) - 4.0 * Helpers.Erf(30.0 * a - 10.0);
        Console.WriteLine("  A = " + a.ToString(CultureInfo.InvariantCulture).PadLeft(14) + ", F(A) = " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        fx = Math.Cos(100.0 * x) - 4.0 * Helpers.Erf(30.0 * x - 10.0);
        Console.WriteLine("  X = " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + ", F(X) = " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        fx = Math.Cos(100.0 * b) - 4.0 * Helpers.Erf(30.0 * b - 10.0);
        Console.WriteLine("  B = " + b.ToString(CultureInfo.InvariantCulture).PadLeft(14) + ", F(B) = " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests BISECTION_RC, to invert the cardiod CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double alpha = 0.0;
        const double beta = 0.25;
        double cdf;
        double x;
        BisectionRC.BisectionData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Demonstrate BISECTION_RC on a probability example.");
        Console.WriteLine("");
        Console.WriteLine("  The cardioid probability density function has a");
        Console.WriteLine("  cumulative density function of the form:");
        Console.WriteLine("    CDF(X) = ( pi + x - alpha + 2 beta * Math.Sin ( x - alpha ) ) / ( 2 * pi )");
        Console.WriteLine("  where alpha and beta are parameters, and x is a value");
        Console.WriteLine("  in the range -pi <= x <= +pi.");
        Console.WriteLine("");
        Console.WriteLine("  CDF(X) is the probability that a random sample will have");
        Console.WriteLine("  a value less than or equal to X.");
        Console.WriteLine("");
        Console.WriteLine("  As X moves from -pi to +pi,");
        Console.WriteLine("  the CDF rises from 0 (no probability)");
        Console.WriteLine("  to 1 (certain probability).");
        Console.WriteLine("");
        Console.WriteLine("  Assuming that:");
        Console.WriteLine("  * ALPHA = " + alpha + "");
        Console.WriteLine("  * BETA =  " + beta + "");
        Console.WriteLine("  determine the value X where the Cardioid CDF is exactly 0.75.");

        const double fx_tol = 1.0E-05;
        const double dx_tol = 1.0E-08;
        int it = 0;
        const int it_max = 30;

        int job = 0;
        double a = -Math.PI;
        double b = +Math.PI;

        double fx = 0.0;

        Console.WriteLine("");
        Console.WriteLine("     I      X               FX              DX");
        Console.WriteLine("");

        for (;;)
        {
            x = BisectionRC.bisection_rc(ref data, ref a, ref b, fx, ref job);

            if (job < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error return.");
                break;
            }

            it += 1;

            cdf = (Math.PI + x - alpha + 2.0 * beta * Math.Sin(x - alpha)) / (2.0 * Math.PI);
            fx = cdf - 0.75;

            double dx = it switch
            {
                <= 2 => Math.Abs(b - a),
                _ => 0.5 * Math.Abs(b - a)
            };

            Console.WriteLine("  " + it.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + dx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (Math.Abs(fx) <= fx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Function is small.");
                break;
            }

            if (dx <= dx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Interval is tiny.");
                break;
            }

            if (it_max > it)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("  Reached iteration limit.");
            break;

        }

        Console.WriteLine("");
        cdf = (Math.PI + a - alpha + 2.0 * beta * Math.Sin(a - alpha)) / (2.0 * Math.PI);
        fx = cdf - 0.75;
        Console.WriteLine("  A = " + a.ToString(CultureInfo.InvariantCulture).PadLeft(14) + ", F(A) = " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        cdf = (Math.PI + x - alpha + 2.0 * beta * Math.Sin(x - alpha)) / (2.0 * Math.PI);
        fx = cdf - 0.75;
        Console.WriteLine("  X = " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + ", F(X) = " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        cdf = (Math.PI + b - alpha + 2.0 * beta * Math.Sin(b - alpha)) / (2.0 * Math.PI);
        fx = cdf - 0.75;
        Console.WriteLine("  B = " + b.ToString(CultureInfo.InvariantCulture).PadLeft(14) + ", F(B) = " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        Console.WriteLine("");
        Console.WriteLine("  Look at the actual cardioid CDF value now:");
        Console.WriteLine("");
        cdf = (Math.PI + x - alpha + 2.0 * beta * Math.Sin(x - alpha)) / (2.0 * Math.PI);
        Console.WriteLine("  Cardioid(" + x + ") = " + cdf + "");

    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests BISECTION_RC for the pipe freezing problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Cleve Moler,
        //    Numerical Computing with MATLAB,
        //    SIAM, 2004,
        //    ISBN13: 978-0-898716-60-3,
        //    LC: QA297.M625,
        //    ebook: http://www.mathworks.com/moler/chapters.html
        //
    {
        double x;
        BisectionRC.BisectionData data = new();

        Console.WriteLine("");
        Console.WriteLine("BISECTION_RC_TEST04");
        Console.WriteLine("  The freezing pipe problem.");
        Console.WriteLine("");
        Console.WriteLine("  At the beginning of a cold spell, the soil is at a uniform");
        Console.WriteLine("  temperature of Ti.  The cold spell applies a uniform air");
        Console.WriteLine("  temperature of Tc, which begins to cool the soil.");
        Console.WriteLine("  As a function of depth x and time t, the soil temperature");
        Console.WriteLine("  will now cool down as:");
        Console.WriteLine("    ( T(x,t) - Tc ) / ( Ti - Tc ) = Helpers.Erf ( 0.5 * x / Math.Sqrt ( alpha * t ) ).");
        Console.WriteLine("  where:");
        Console.WriteLine("    Ti =  20 degrees centigrade,");
        Console.WriteLine("    Tc = -15 degrees centigrade,");
        Console.WriteLine("    alpha = 0.000000138 meter^2 / second, thermal conductivity;");
        Console.WriteLine("    and erf() is the error function.");
        Console.WriteLine("  Water freezes at 0 degrees centigrade.");
        Console.WriteLine("");
        Console.WriteLine("  What depth x in meters must a water pipe be buried so that it will");
        Console.WriteLine("  not freeze even if this cold snap lasts for 60 days?");
        //
        //  Problem parameters.
        //
        const double ti = 20.0;
        const double tc = -15.0;
        const double t = 60.0 * 24.0 * 60.0 * 60.0;
        const double alpha = 0.000000138;
        //
        //  Iteration parameters.
        //
        const double fx_tol = 1.0E-09;
        const double dx_tol = 1.0E-09;
        int it = 0;
        const int it_max = 30;
        int job = 0;
        double fx = 0.0;
        //
        //  Initial guess for interval.
        //
        double a = 0.0;
        double b = 1000.0;

        Console.WriteLine("");
        Console.WriteLine("     I      X               FX              DX");
        Console.WriteLine("");

        for (;;)
        {
            x = BisectionRC.bisection_rc(ref data, ref a, ref b, fx, ref job);

            if (job < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error return.");
                break;
            }

            it += 1;

            fx = tc + (ti - tc) * Helpers.Erf(0.5 * x / Math.Sqrt(alpha * t));

            double dx = it switch
            {
                <= 2 => Math.Abs(b - a),
                _ => 0.5 * Math.Abs(b - a)
            };

            Console.WriteLine("  " + it.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + dx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (Math.Abs(fx) <= fx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Function is small.");
                break;
            }

            if (dx <= dx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Interval is tiny.");
                break;
            }

            if (it_max > it)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("  Reached iteration limit.");
            break;

        }

        Console.WriteLine("");
        fx = tc + (ti - tc) * Helpers.Erf(0.5 * a / Math.Sqrt(alpha * t));
        Console.WriteLine("  A = " + a + ", F(A) = " + fx + "");
        fx = tc + (ti - tc) * Helpers.Erf(0.5 * x / Math.Sqrt(alpha * t));
        Console.WriteLine("  X = " + x + ", F(X) = " + fx + "");
        fx = tc + (ti - tc) * Helpers.Erf(0.5 * b / Math.Sqrt(alpha * t));
        Console.WriteLine("  B = " + b + ", F(B) = " + fx + "");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests BISECTION_RC for Kepler's problem.
        //
        //  Discussion:
        //
        //    Kepler's equation has the form:
        //
        //      X = M + E * Math.Sin ( X )
        //
        //    X represents the eccentric anomaly of a planet, the angle between the
        //    perihelion (the point on the orbit nearest to the sun) through the sun 
        //    to the center of the ellipse, and the line from the center of the ellipse
        //    to the planet.
        //
        //    There are two parameters, E and M:
        //
        //    * E is the eccentricity of the orbit, which should be between 0 and 1.0;
        //
        //    * M is the angle from the perihelion made by a fictitious planet traveling
        //      on a circular orbit centered at the sun, and traveling at a constant
        //      angular velocity equal to the average angular velocity of the true
        //      planet.  M is usually between 0 and 180 degrees, but can have any value.
        //
        //    For convenience, X and M are measured in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Cleve Moler,
        //    Numerical Computing with MATLAB,
        //    SIAM, 2004,
        //    ISBN13: 978-0-898716-60-3,
        //    LC: QA297.M625,
        //    ebook: http://www.mathworks.com/moler/chapters.html
        //
    {
        double xr;
        BisectionRC.BisectionData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  The Kepler equation.");
        Console.WriteLine("");
        Console.WriteLine("  Kepler's equation has the form");
        Console.WriteLine("");
        Console.WriteLine("    X = M + E * Math.Sin ( X )");
        Console.WriteLine("");
        Console.WriteLine("  X represents the eccentric anomaly of a planet, the angle between the");
        Console.WriteLine("  perihelion (the point on the orbit nearest to the sun) through the sun");
        Console.WriteLine("  to the center of the ellipse, and the line from the center of the ellipse");
        Console.WriteLine("  to the planet.");
        Console.WriteLine("");
        Console.WriteLine("  There are two parameters, E and M:");
        Console.WriteLine("");
        Console.WriteLine("  * E is the eccentricity of the orbit, which should be between 0 and 1.0;");
        Console.WriteLine("");
        Console.WriteLine("  * M is the angle from the perihelion made by a fictitious planet traveling");
        Console.WriteLine("    on a circular orbit centered at the sun, and traveling at a constant");
        Console.WriteLine("    angular velocity equal to the average angular velocity of the true");
        Console.WriteLine("    planet.  M is usually between 0 and 180 degrees, but can have any value.");
        Console.WriteLine("");
        Console.WriteLine("  For convenience, X and M are measured in degrees.");
        //
        //  Problem parameters.
        //
        const double md = 24.851090;
        const double mr = md * Math.PI / 180.0;
        const double e = 0.1;

        Console.WriteLine("");
        Console.WriteLine("  Given eccentricity E = " + e + "");
        Console.WriteLine("  Given angle M = " + md + " (degrees)");
        Console.WriteLine("                = " + mr + " (radians)");
        Console.WriteLine("");
        Console.WriteLine("  Given E and M, find corresponding X.");
        //
        //  Iteration parameters.
        //
        const double fx_tol = 1.0E-09;
        const double dx_tol = 1.0E-09;
        int it = 0;
        const int it_max = 30;
        int job = 0;
        double fx = 0.0;
        //
        //  Initial guess for interval.
        //
        double ad = 0.0;
        double bd = 180.0;

        double ar = ad * Math.PI / 180.0;
        double br = bd * Math.PI / 180.0;

        Console.WriteLine("");
        Console.WriteLine("     I      X               FX              DX");
        Console.WriteLine("");

        for (;;)
        {
            xr = BisectionRC.bisection_rc(ref data, ref ar, ref br, fx, ref job);

            if (job < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error return.");
                break;
            }

            it += 1;

            fx = xr - mr - e * Math.Sin(xr);

            double dx = it switch
            {
                <= 2 => Math.Abs(br - ar),
                _ => 0.5 * Math.Abs(br - ar)
            };

            Console.WriteLine("  " + it.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + xr.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + dx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (Math.Abs(fx) <= fx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Function is small.");
                break;
            }

            if (dx <= dx_tol)
            {
                Console.WriteLine("");
                Console.WriteLine("  Interval is tiny.");
                break;
            }

            if (it_max > it)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("  Reached iteration limit.");
            break;

        }

        Console.WriteLine("");
        Console.WriteLine("  In Radians:");
        Console.WriteLine("");
        fx = ar - mr - e * Math.Sin(ar);
        Console.WriteLine("  A = " + ar + ", F(A) = " + fx + "");
        fx = xr - mr - e * Math.Sin(xr);
        Console.WriteLine("  X = " + xr + ", F(X) = " + fx + "");
        fx = br - mr - e * Math.Sin(br);
        Console.WriteLine("  B = " + br + ", F(B) = " + fx + "");

        ad = ar * 180.0 / Math.PI;
        double xd = xr * 180.0 / Math.PI;
        bd = br * 180.0 / Math.PI;

        Console.WriteLine("");
        Console.WriteLine("  In Degrees:");
        Console.WriteLine("");
        fx = (ad - md) * Math.PI / 180.0 - e * Math.Sin(ad * Math.PI / 180.0);
        Console.WriteLine("  A = " + ad + ", F(A) = " + fx + "");
        fx = (xd - md) * Math.PI / 180.0 - e * Math.Sin(xd * Math.PI / 180.0);
        Console.WriteLine("  X = " + xd + ", F(X) = " + fx + "");
        fx = (bd - md) * Math.PI / 180.0 - e * Math.Sin(bd * Math.PI / 180.0);
        Console.WriteLine("  B = " + bd + ", F(B) = " + fx + "");
    }
}