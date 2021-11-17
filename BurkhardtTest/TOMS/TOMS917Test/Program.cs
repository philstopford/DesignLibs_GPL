using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using Burkardt;
using Burkardt.Types;

namespace TOMS917Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS917_TEST.
        //
        //  Discussion:
        //
        //    WALSH_TEST tests the WALSH library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double pi = Math.PI;
        Complex z;

        Console.WriteLine("");
        Console.WriteLine("TOMS917_TEST");
        Console.WriteLine("  Test the TOMS917 library.");

        a = 0.0;
        b = 0.0;
        z = new Complex(a, b);
        driver(z);

        a = 1.0;
        b = 0.0;
        z = new Complex(a, b);
        driver(z);

        a = 1.0 + Math.Exp(1.0);
        b = 0.0;
        z = new Complex(a, b);
        driver(z);

        a = -1.0;
        b = pi;
        z = new Complex(a, b);
        driver(z);

        a = -1.0;
        b = -pi;
        z = new Complex(a, b);
        driver(z);

        a = -2.0 + Math.Log(2.0);
        b = pi;
        z = new Complex(a, b);
        driver(z);

        a = -2.0 + Math.Log(2.0);
        b = -pi;
        z = new Complex(a, b);
        driver(z);

        a = 0.0;
        b = 1.0 + pi / 2.0;
        z = new Complex(a, b);
        driver(z);

        a = 0.0;
        b = pi;
        z = new Complex(a, b);
        driver(z);

        a = 1.0;
        b = 1.0;
        z = new Complex(a, b);
        driver(z);
        //
        //  Test the function near the region boundaries.
        //
        test_boundary();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("TOMS917_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void driver(Complex z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DRIVER calls the simple and extended Wright Omega evaluators.
        //
        //  Modified:
        //
        //    14 May 2016
        //
        //  Author:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey
        //
        //  Reference:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey,
        //    Algorithm 917: Complex Double-Precision Evaluation of the Wright Omega 
        //    Function,
        //    ACM Transactions on Mathematical Software,
        //    Volume 38, Number 3, Article 20, April 2012, 17 pages.
        //
        //  Parameters:
        //
        //    Input, Complex Z, the argument of the Wright Omega function.
        //
    {
        Complex condest = new();
        Complex e = new();
        Complex r = new();
        Complex r_ult = new();
        Complex w = new();

        Console.WriteLine("");
        Console.WriteLine("DRIVER:");
        Console.WriteLine("  Demonstrate simple and extended Wright Omega evaluators.");
        //
        //  Simple evaluator.
        //
        w = WrightOmega.wrightomega(z);

        Console.WriteLine("");
        Console.WriteLine("  Calling:");
        Console.WriteLine("    w = wrightomega(z);");
        Console.WriteLine("  returns:");
        Console.WriteLine("    w = omega(" + z.Real
                                           + ", " + z.Imaginary
                                           + ") =  ( " + w.Real
                                           + ", " + w.Imaginary + ")");
        //
        //  Extended evaluator.
        //
        WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref condest);

        Console.WriteLine("");
        Console.WriteLine("  Calling:");
        Console.WriteLine("    wrightomega_ext ( z, w, e, r, condest );");
        Console.WriteLine("  returns:");
        Console.WriteLine("    w = omega(" + z.Real
                                           + ", " + z.Imaginary
                                           + ") =  ( " + w.Real
                                           + ", " + w.Imaginary + ")");
        Console.WriteLine("  e = last update step = ( " + e.Real
                                                        + ", " + e.Imaginary + ")");
        Console.WriteLine("  r = penultimate residual = ( " + r.Real
                                                            + ", " + r.Imaginary + ")");
        Console.WriteLine("  condest = condition number estimate = ( " + condest.Real
                                                                       + ", " + condest.Imaginary + ")");
        //
        //  Calculate and print ultimate residual.
        //
        r_ult = (2.0 * w * w - 8.0 * w - 1.0)
            / Complex.Pow(1.0 + w, 6.0) * r * r * r * r;
        Console.WriteLine("");
        Console.WriteLine("  ultimate residual = ( " + r_ult.Real
                                                     + ", " + r_ult.Imaginary + ")");

    }

    private static void test_boundary()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST_BOUNDARY tests wrightomega() near boundaries.
        //
        //  Discussion:
        //
        //    This is a test driver to evaluate the Wright Omega function along the
        //    boundaries of the different regions of the approximations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2016
        //
        //  Author:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey
        //
        //  Reference:
        //
        //    Piers Lawrence, Robert Corless, David Jeffrey,
        //    Algorithm 917: Complex Double-Precision Evaluation of the Wright Omega 
        //    Function,
        //    ACM Transactions on Mathematical Software,
        //    Volume 38, Number 3, Article 20, April 2012, 17 pages.
        //
    {
        double a;
        double b;
        Complex cond = new();
        Complex e = new();
        double exp_num = 160.0;
        string filename = "results.txt";
        List<string> fp = new();
        int i;
        int n = 100;
        double pi = Math.PI;
        Complex r = new();
        double td;
        Complex w = new();
        double[] x = new double[2];
        double[] y = new double[2];
        Complex z = new();

        Console.WriteLine("");
        Console.WriteLine("TEST_BOUNDARY:");
        Console.WriteLine("  Test wrightomega_ext() near approximation region boundaries.");
        Console.WriteLine("  Store results in a file for comparison with benchmark data.");

        //
        //  We want trailing zeros, to make comparison with benchmark easier.
        //
        //
        //  Region 1;
        //  x=(-2.0,1.0] ,y=2*pi
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        y[0] = NextAfterNS.NextAfter.nextafter(2.0 * pi, -1.0);
        td = (x[1] - x[0]) / n;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 2;
        //  x=1.0 ,y=(1.0,2*pi)
        //
        x[0] = 1.0;
        y[0] = NextAfterNS.NextAfter.nextafter(1.0, 2.0);
        y[1] = NextAfterNS.NextAfter.nextafter(2.0 * pi, 1.0);
        td = -(y[1] - y[0]) / n;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[1] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 3;
        //  x=(-2.0,1.0] ,y=1.0
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        td = -(x[1] - x[0]) / n;
        y[0] = NextAfterNS.NextAfter.nextafter(1.0, 2.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[1] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 4;
        //  x=-2.0 ,y=(1.0,2*pi)
        //
        y[0] = NextAfterNS.NextAfter.nextafter(1.0, 2.0);
        y[1] = NextAfterNS.NextAfter.nextafter(2.0 * pi, -1.0);
        td = (y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 5;
        //  x=(-2.0,1.0] ,y=-2*pi
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        y[0] = NextAfterNS.NextAfter.nextafter(-2.0 * pi, 1.0);
        td = (x[1] - x[0]) / n;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 6;
        //  x=1.0, y=(-2*pi,-1.0)
        //
        y[0] = NextAfterNS.NextAfter.nextafter(-2.0 * pi, 1.0);
        y[1] = NextAfterNS.NextAfter.nextafter(-1.0, -2.0);
        td = (y[1] - y[0]) / n;
        x[0] = 1.0;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 7;
        //  x=(-2.0,1.0] ,y=-1.0
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        td = (x[0] - x[1]) / n;
        y[0] = NextAfterNS.NextAfter.nextafter(-1.0, -2.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[1] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 8;
        //  x=-2.0 ,y=(-2*pi,-1.0)
        //
        y[0] = NextAfterNS.NextAfter.nextafter(-2.0 * pi, 1.0);
        y[1] = NextAfterNS.NextAfter.nextafter(-1.0, -2.0);
        td = (y[0] - y[1]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[1] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 9;
        //  x=-2.0 y=[-1.0,1.0]
        //
        y[0] = -1.0;
        y[1] = 1.0;
        td = (y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 10;
        //  x=(-2.0,1.0] y=1.0
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        td = (x[1] - x[0]) / n;
        y[0] = 1.0;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 11
        //  x=1.0 y=[1.0,pi]
        //
        y[0] = 1.0;
        y[1] = pi;
        td = (y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(1.0, 2.0);
        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 12
        //  (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
        //  (on inside)
        //
        td = pi / n;
        x[0] = pi / 2.0;

        for (i = 0; i < n; i++)
        {
            a = NextAfterNS.NextAfter.nextafter(pi, -1.0) * Math.Cos(x[0] - td * i)
                + NextAfterNS.NextAfter.nextafter(1.0, -1.0);
            b = NextAfterNS.NextAfter.nextafter(pi, -1.0) * Math.Sin(x[0] - td * i);
            z = new Complex(a, b);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 13
        //  x=1.0 y=[-pi,-1.0]
        //
        y[0] = -pi;
        y[1] = -1.0;
        td = (y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(1.0, 2.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 14
        //  x=(-2.0,1.0] y=-1.0
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        td = (x[1] - x[0]) / n;
        y[0] = -1.0;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 15
        //  x=(-inf,-2) y=pi^+
        //
        for (i = 0; i < n; i++)
        {
            x[0] = NextAfterNS.NextAfter.nextafter(-1.0 - Math.Exp((n - 1 - i) / exp_num),
                typeMethods.r8_huge());
            y[0] = NextAfterNS.NextAfter.nextafter(pi - 0.75 * (x[0] + 1.0), typeMethods.r8_huge());
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 16
        //
        y[0] = 0.75 + pi;
        y[1] = 2.0 * pi;
        td = (y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, -3.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 17
        //  x=(-2.0,1.0] ,y=2*pi
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        y[0] = 2.0 * pi;
        td = (x[1] - x[0]) / n;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        /*
        Region 18
        x=1.0 ,y=(pi,2*pi)
        */
        y[0] = NextAfterNS.NextAfter.nextafter(pi, 6.0);
        y[1] = NextAfterNS.NextAfter.nextafter(2.0 * pi, 1.0);
        td = -(y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(1.0, 2.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[1] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 19
        //  (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
        //  (on outside)
        //
        td = pi / (n - 1);
        y[0] = pi / 2.0;

        for (i = 0; i < n; i++)
        {
            y[1] = pi * Math.Sin(y[0] - td * i);
            x[0] = Math.Sqrt(pi * pi - y[1] * y[1]) + 1.0;
            z = y[1] switch
            {
                < 0 => new Complex(NextAfterNS.NextAfter.nextafter(x[0], typeMethods.r8_huge()),
                    NextAfterNS.NextAfter.nextafter(y[1], -typeMethods.r8_huge())),
                _ => new Complex(NextAfterNS.NextAfter.nextafter(x[0], typeMethods.r8_huge()),
                    NextAfterNS.NextAfter.nextafter(y[1], typeMethods.r8_huge()))
            };

            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 20;
        //  x=1.0 ,y=(-2*pi,-pi)
        //
        y[0] = NextAfterNS.NextAfter.nextafter(-2.0 * pi, 1.0);
        y[1] = NextAfterNS.NextAfter.nextafter(-pi, -6.0);
        td = -(y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(1.0, 2.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[1] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 21;
        //  x=(-2.0,1.0] ,y=-2*pi
        //
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, 1.0);
        x[1] = 1.0;
        y[0] = NextAfterNS.NextAfter.nextafter(-2.0 * pi, -7.0);
        td = -(x[1] - x[0]) / n;

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[1] + td * i, y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 22
        //
        y[0] = -0.75 - pi;
        y[1] = -2.0 * pi;
        td = -(y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, -3.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[1] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 23
        //  x=(-inf,-2) y=pi^+
        //
        for (i = 0; i < n; i++)
        {
            x[0] = NextAfterNS.NextAfter.nextafter(-1.0 - Math.Exp(i / exp_num), typeMethods.r8_huge());
            y[0] = NextAfterNS.NextAfter.nextafter(-pi + 0.75 * (x[0] + 1.0), -typeMethods.r8_huge());
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 24
        //  x=(-inf,-2) y=pi^+
        //
        for (i = 0; i < n; i++)
        {
            x[0] = NextAfterNS.NextAfter.nextafter(-1.0 - Math.Exp((n - 1 - i) / exp_num),
                -typeMethods.r8_huge());
            y[0] = NextAfterNS.NextAfter.nextafter(-pi + 0.75 * (x[0] + 1.0), typeMethods.r8_huge());
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 25
        //
        y[0] = -pi;
        y[1] = -0.75 - pi;
        td = -(y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, -3.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[1] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 26
        //  x=(-inf,-2) y=pi^+
        //
        y[0] = NextAfterNS.NextAfter.nextafter(-pi, -7.0);

        for (i = 0; i < n; i++)
        {
            x[0] = -1.0 - Math.Exp(i / exp_num);
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 27
        //  x=(-inf,-2) y=pi^+
        //
        y[0] = NextAfterNS.NextAfter.nextafter(-pi, 1.0);

        for (i = 0; i < n; i++)
        {
            x[0] = -1.0 - Math.Exp((n - 1 - i) / exp_num);
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 28
        //
        y[0] = NextAfterNS.NextAfter.nextafter(-pi, 1.0);
        y[1] = pi;
        td = (y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, -3.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 29
        //  x=(-inf,-2) y=pi^+
        //
        y[0] = NextAfterNS.NextAfter.nextafter(pi, 1.0);

        for (i = 0; i < n; i++)
        {
            x[0] = -1.0 - Math.Exp(i / exp_num);
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 30
        //  x=(-inf,-2) y=pi^+
        //
        y[0] = NextAfterNS.NextAfter.nextafter(pi, 7.0);

        for (i = 0; i < n; i++)
        {
            x[0] = -1.0 - Math.Exp((n - 1 - i) / exp_num);
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 31
        //
        y[0] = NextAfterNS.NextAfter.nextafter(pi, 7.0);
        y[1] = 0.75 + pi;
        td = (y[1] - y[0]) / n;
        x[0] = NextAfterNS.NextAfter.nextafter(-2.0, -3.0);

        for (i = 0; i < n; i++)
        {
            z = new Complex(x[0], y[0] + td * i);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        //
        //  Region 32
        //  x=(-inf,-2) y=pi^+
        //
        for (i = 0; i < n; i++)
        {
            x[0] = -1.0 - Math.Exp((n - 1 - i) / exp_num);
            y[0] = NextAfterNS.NextAfter.nextafter(pi - 0.75 * (x[0] + 1.0), 0.1);
            z = new Complex(x[0], y[0]);
            WrightOmega.wrightomega_ext(z, ref w, ref e, ref r, ref cond);
            fp.Add(z.Real + " "
                          + z.Imaginary + " "
                          + w.Real + " "
                          + w.Imaginary + "");
        }

        File.WriteAllLines(filename, fp);

        Console.WriteLine("");
        Console.WriteLine("TEST_BOUNDARY:");
        Console.WriteLine("  Results saved in file '" + filename + "'");

    }
}