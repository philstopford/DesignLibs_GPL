using System;
using Burkardt.Bisection;

namespace Burkardt.BisectionRCTest
{
    class Program
    {
        static void Main(string[] args)
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
            //  LicenMath.Sing:
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

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests BISECTION_RC, evaluating the function in a separate routine.
            //
            //  LicenMath.Sing:
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
            double a;
            double b;
            double dx;
            double dx_tol;
            double fx;
            double fx_tol;
            int it;
            int it_max;
            int job;
            double x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Demonstrate BISECTION_RC on a simple example.");
            Console.WriteLine("  The function is evaluated in a separate routine.");

            fx_tol = 1.0E-08;
            dx_tol = 1.0E-06;
            it = 0;
            it_max = 30;

            a = 0.0;
            b = 1.0;
            fx = 0.0;
            job = 0;

            Console.WriteLine("");
            Console.WriteLine("     I      X               FX              DX");
            Console.WriteLine("");

            for (;;)
            {
                RC_data data = new RC_data();
                x = RC.bisection_rc(ref data, ref a, ref b, fx, ref job);

                if (job < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Error return.");
                    break;
                }

                it = it + 1;

                fx = f01(x);

                if (it <= 2)
                {
                    dx = Math.Abs(b - a);
                }
                else
                {
                    dx = 0.5 * Math.Abs(b - a);
                }

                Console.WriteLine("  " + it.ToString().PadLeft(4)
                    + "  " + x.ToString().PadLeft(14)
                    + "  " + fx.ToString().PadLeft(14)
                    + "  " + dx.ToString().PadLeft(14) + "");

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

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Reached iteration limit.");
                    break;
                }

            }

            Console.WriteLine("");
            Console.WriteLine("  A = " + a.ToString().PadLeft(14) + " F(A) = " + f01(a) + "");
            Console.WriteLine("  X = " + x.ToString().PadLeft(14) + " F(X) = " + f01(x) + "");
            Console.WriteLine("  B = " + b.ToString().PadLeft(14) + " F(B) = " + f01(b) + "");

            return;
        }

        static double f01(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F01 evaluates the function f(x) = Math.Cos ( x ) - x which is zero around 0.74
            //
            //  LicenMath.Sing:
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

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests BISECTION_RC, evaluating the function within the routine.
            //
            //  LicenMath.Sing:
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
            double a;
            double b;
            double dx;
            double dx_tol;
            double fx;
            double fx_tol;
            int it;
            int it_max;
            int job;
            double x;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Demonstrate BISECTION_RC on a simple example.");
            Console.WriteLine("  The function is evaluated within this routine.");

            fx_tol = 1.0E-09;
            dx_tol = 1.0E-09;
            it = 0;
            it_max = 30;

            a = 0.0;
            b = 1.0;
            fx = 0.0;
            job = 0;

            Console.WriteLine("");
            Console.WriteLine("     I      X               FX              DX");
            Console.WriteLine("");

            for (;;)
            {
                RC_data data = new RC_data();
                x = RC.bisection_rc(ref data, ref a, ref b, fx, ref job);
                
                if (job < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Error return.");
                    break;
                }

                it = it + 1;

                fx = Math.Cos(100.0 * x) - 4.0 * Helpers.Erf(30.0 * x - 10.0);

                if (it <= 2)
                {
                    dx = Math.Abs(b - a);
                }
                else
                {
                    dx = 0.5 * Math.Abs(b - a);
                }

                Console.WriteLine("  " + it.ToString().PadLeft(4)
                    + "  " + x.ToString().PadLeft(14)
                    + "  " + fx.ToString().PadLeft(14)
                    + "  " + dx.ToString().PadLeft(14) + "");

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

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Reached iteration limit.");
                    break;
                }

            }

            Console.WriteLine("");
            fx = Math.Cos(100.0 * a) - 4.0 * Helpers.Erf(30.0 * a - 10.0);
            Console.WriteLine("  A = " + a.ToString().PadLeft(14) + ", F(A) = " + fx.ToString().PadLeft(14) + "");
            fx = Math.Cos(100.0 * x) - 4.0 * Helpers.Erf(30.0 * x - 10.0);
            Console.WriteLine("  X = " + x.ToString().PadLeft(14) + ", F(X) = " + fx.ToString().PadLeft(14) + "");
            fx = Math.Cos(100.0 * b) - 4.0 * Helpers.Erf(30.0 * b - 10.0);
            Console.WriteLine("  B = " + b.ToString().PadLeft(14) + ", F(B) = " + fx.ToString().PadLeft(14) + "");
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests BISECTION_RC, to invert the cardiod CDF.
            //
            //  LicenMath.Sing:
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
            double a;
            double alpha = 0.0;
            double b;
            double beta = 0.25;
            double cdf;
            double dx;
            double dx_tol;
            double fx;
            double fx_tol;
            int it;
            int it_max;
            int job;
            double r8_pi = 3.141592653589793;
            double x;

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

            fx_tol = 1.0E-05;
            dx_tol = 1.0E-08;
            it = 0;
            it_max = 30;

            job = 0;
            a = -r8_pi;
            b = +r8_pi;

            fx = 0.0;

            Console.WriteLine("");
            Console.WriteLine("     I      X               FX              DX");
            Console.WriteLine("");

            for (;;)
            {
                RC_data data = new RC_data();
                x = RC.bisection_rc(ref data, ref a, ref b, fx, ref job);
                if (job < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Error return.");
                    break;
                }

                it = it + 1;

                cdf = (r8_pi + x - alpha + 2.0 * beta * Math.Sin(x - alpha)) / (2.0 * r8_pi);
                fx = cdf - 0.75;

                if (it <= 2)
                {
                    dx = Math.Abs(b - a);
                }
                else
                {
                    dx = 0.5 * Math.Abs(b - a);
                }

                Console.WriteLine("  " + it.ToString().PadLeft(4)
                    + "  " + x.ToString().PadLeft(14)
                    + "  " + fx.ToString().PadLeft(14)
                    + "  " + dx.ToString().PadLeft(14) + "");

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

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Reached iteration limit.");
                    break;
                }

            }

            Console.WriteLine("");
            cdf = (r8_pi + a - alpha + 2.0 * beta * Math.Sin(a - alpha)) / (2.0 * r8_pi);
            fx = cdf - 0.75;
            Console.WriteLine("  A = " + a.ToString().PadLeft(14) + ", F(A) = " + fx.ToString().PadLeft(14) + "");
            cdf = (r8_pi + x - alpha + 2.0 * beta * Math.Sin(x - alpha)) / (2.0 * r8_pi);
            fx = cdf - 0.75;
            Console.WriteLine("  X = " + x.ToString().PadLeft(14) + ", F(X) = " + fx.ToString().PadLeft(14) + "");
            cdf = (r8_pi + b - alpha + 2.0 * beta * Math.Sin(b - alpha)) / (2.0 * r8_pi);
            fx = cdf - 0.75;
            Console.WriteLine("  B = " + b.ToString().PadLeft(14) + ", F(B) = " + fx.ToString().PadLeft(14) + "");

            Console.WriteLine("");
            Console.WriteLine("  Look at the actual cardioid CDF value now:");
            Console.WriteLine("");
            cdf = (r8_pi + x - alpha + 2.0 * beta * Math.Sin(x - alpha)) / (2.0 * r8_pi);
            Console.WriteLine("  Cardioid(" + x + ") = " + cdf + "");

            return;
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests BISECTION_RC for the pipe freezing problem.
            //
            //  LicenMath.Sing:
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
            double a;
            double alpha;
            double b;
            double dx;
            double dx_tol;
            double fx;
            double fx_tol;
            int it;
            int it_max;
            int job;
            double t;
            double tc;
            double ti;
            double x;

            Console.WriteLine("");
            Console.WriteLine("BISECTION_RC_TEST04");
            Console.WriteLine("  The freezing pipe problem.");
            Console.WriteLine("");
            Console.WriteLine("  At the beginning of a cold spell, the soil is at a uniform");
            Console.WriteLine("  temperature of Ti.  The cold spell applies a uniform air");
            Console.WriteLine("  temperature of Tc, which begins to cool the soil.");
            Console.WriteLine("  As a function of depth x and time t, the soil temperature");
            Console.WriteLine("  will now cool down as:");
            Console.WriteLine("    ( T(x,t) - Tc ) / ( Ti - Tc ) = erf ( 0.5 * x / Math.Sqrt ( alpha * t ) ).");
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
            ti = 20.0;
            tc = -15.0;
            t = 60.0 * 24.0 * 60.0 * 60.0;
            alpha = 0.000000138;
            //
            //  Iteration parameters.
            //
            fx_tol = 1.0E-09;
            dx_tol = 1.0E-09;
            it = 0;
            it_max = 30;
            job = 0;
            fx = 0.0;
            //
            //  Initial guess for interval.
            //
            a = 0.0;
            b = 1000.0;

            Console.WriteLine("");
            Console.WriteLine("     I      X               FX              DX");
            Console.WriteLine("");

            for (;;)
            {
                RC_data data = new RC_data();
                x = RC.bisection_rc(ref data, ref a, ref b, fx, ref job);
                
                if (job < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Error return.");
                    break;
                }

                it = it + 1;

                fx = tc + (ti - tc) * Helpers.Erf(0.5 * x / Math.Sqrt(alpha * t));

                if (it <= 2)
                {
                    dx = Math.Abs(b - a);
                }
                else
                {
                    dx = 0.5 * Math.Abs(b - a);
                }

                Console.WriteLine("  " + it.ToString().PadLeft(4)
                    + "  " + x.ToString().PadLeft(14)
                    + "  " + fx.ToString().PadLeft(14)
                    + "  " + dx.ToString().PadLeft(14) + "");

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

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Reached iteration limit.");
                    break;
                }

            }

            Console.WriteLine("");
            fx = tc + (ti - tc) * Helpers.Erf(0.5 * a / Math.Sqrt(alpha * t));
            Console.WriteLine("  A = " + a + ", F(A) = " + fx + "");
            fx = tc + (ti - tc) * Helpers.Erf(0.5 * x / Math.Sqrt(alpha * t));
            Console.WriteLine("  X = " + x + ", F(X) = " + fx + "");
            fx = tc + (ti - tc) * Helpers.Erf(0.5 * b / Math.Sqrt(alpha * t));
            Console.WriteLine("  B = " + b + ", F(B) = " + fx + "");

            return;
        }

        static void test05()

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
            //  LicenMath.Sing:
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
            double ad;
            double ar;
            double bd;
            double br;
            double dx;
            double dx_tol;
            double e;
            double fx;
            double fx_tol;
            int it;
            int it_max;
            int job;
            double md;
            double mr;
            const double r8_pi = 3.141592653589793;
            double xd;
            double xr;

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
            md = 24.851090;
            mr = md * r8_pi / 180.0;
            e = 0.1;

            Console.WriteLine("");
            Console.WriteLine("  Given eccentricity E = " + e + "");
            Console.WriteLine("  Given angle M = " + md + " (degrees)");
            Console.WriteLine("                = " + mr + " (radians)");
            Console.WriteLine("");
            Console.WriteLine("  Given E and M, find corresponding X.");
            //
            //  Iteration parameters.
            //
            fx_tol = 1.0E-09;
            dx_tol = 1.0E-09;
            it = 0;
            ;
            it_max = 30;
            job = 0;
            fx = 0.0;
            //
            //  Initial guess for interval.
            //
            ad = 0.0;
            bd = 180.0;

            ar = ad * r8_pi / 180.0;
            br = bd * r8_pi / 180.0;

            Console.WriteLine("");
            Console.WriteLine("     I      X               FX              DX");
            Console.WriteLine("");

            for (;;)
            {
                RC_data data = new RC_data();
                xr = RC.bisection_rc(ref data, ref ar, ref br, fx, ref job);

                if (job < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Error return.");
                    break;
                }

                it = it + 1;

                fx = xr - mr - e * Math.Sin(xr);

                if (it <= 2)
                {
                    dx = Math.Abs(br - ar);
                }
                else
                {
                    dx = 0.5 * Math.Abs(br - ar);
                }

                Console.WriteLine("  " + it.ToString().PadLeft(4)
                    + "  " + xr.ToString().PadLeft(14)
                    + "  " + fx.ToString().PadLeft(14)
                    + "  " + dx.ToString().PadLeft(14) + "");

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

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Reached iteration limit.");
                    break;
                }

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

            ad = ar * 180.0 / r8_pi;
            xd = xr * 180.0 / r8_pi;
            bd = br * 180.0 / r8_pi;

            Console.WriteLine("");
            Console.WriteLine("  In Degrees:");
            Console.WriteLine("");
            fx = (ad - md) * r8_pi / 180.0 - e * Math.Sin(ad * r8_pi / 180.0);
            Console.WriteLine("  A = " + ad + ", F(A) = " + fx + "");
            fx = (xd - md) * r8_pi / 180.0 - e * Math.Sin(xd * r8_pi / 180.0);
            Console.WriteLine("  X = " + xd + ", F(X) = " + fx + "");
            fx = (bd - md) * r8_pi / 180.0 - e * Math.Sin(bd * r8_pi / 180.0);
            Console.WriteLine("  B = " + bd + ", F(B) = " + fx + "");
        }
    }
}