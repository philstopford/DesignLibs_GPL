using System;
using Burkardt.Transform;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SineTransformTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SINE_TRANSFORM_TEST.
        //
        //  Discussion:
        //
        //    SINE_TRANSFORM_TEST tests SINE_TRANSFORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SINE_TRANSFORM_TEST");
        Console.WriteLine("  Test the SINE_TRANSFORM library.");

        sine_transform_test01();
        sine_transform_test02();
        sine_transform_test03();
        sine_transform_test04();

        Console.WriteLine("");
        Console.WriteLine("SINE_TRANSFORM_TEST");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void sine_transform_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINE_TRANSFORM_TEST01 demonstrates that the transform is its own inverse.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 10;
        int seed;
        double[] r;
        double[] s;
        double[] t;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("SINE_TRANSFORM_TEST01:");
        Console.WriteLine("  SINE_TRANSFORM_DATA does a sine transform of data");
        Console.WriteLine("  defined by a vector.");
        Console.WriteLine("");
        Console.WriteLine("  Demonstrate that the transform is its own inverse.");
        Console.WriteLine("  Let R be a random N vector.");
        Console.WriteLine("  Let S be the transform of D.");
        Console.WriteLine("  Let T be the transform of E.");
        Console.WriteLine("  Then R and T will be equal.");

        r = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        s = Sine.sine_transform_data(n, r);
        t = Sine.sine_transform_data(n, s);

        Console.WriteLine("");
        Console.WriteLine("     I      R(I)        S(I)        T(I)");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + r[i].ToString().PadLeft(10)
                                   + "  " + s[i].ToString().PadLeft(10)
                                   + "  " + t[i].ToString().PadLeft(10) + "");
        }
    }

    private static void sine_transform_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINE_TRANSFORM_TEST02 uses the functional form of the routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] f1;
        double[] f2;
        double fa;
        double fb;
        int i;
        int n = 9;
        double[] s;
        double[] x;

        a = 1.0;
        b = 3.0;
        //
        //  Evenly spaced points between A and B, but omitting
        //  A and B themselves.
        //
        x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i) * a
                    + (i + 1) * b)
                   / (n + 1);
        }

        Console.WriteLine("");
        Console.WriteLine("SINE_TRANSFORM_TEST02:");
        Console.WriteLine("  SINE_TRANSFORM_FUNCTION does a sine transform of data");
        Console.WriteLine("  defined by a function F(X) evaluated at equally spaced");
        Console.WriteLine("  points in an interval [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Demonstrate that the transform is its own inverse.");
        Console.WriteLine("  Let X(0:N+1) be N+2 equally spaced points in [A,B].");
        Console.WriteLine("  Let S be the transform of F(X(1:N)).");
        Console.WriteLine("  Let F1 be the linear interpolant of (A,F(A)), (B,F(B)).");
        Console.WriteLine("  Let F2 be the transform of S.");
        Console.WriteLine("  Then F(X(1:N)) = F1(X(1:N)) + F2(1:N).");

        s = Sine.sine_transform_function(n, a, b, poly5);

        fa = poly5(a);
        fb = poly5(b);

        f1 = new double[n];

        for (i = 0; i < n; i++)
        {
            f1[i] = ((b - x[i]) * fa
                     + (x[i] - a) * fb)
                    / (b - a);
        }

        f2 = Sine.sine_transform_data(n, s);

        Console.WriteLine("");
        Console.WriteLine("     I      X(I)      F(X(I))       S           F1          F2          F1+F2");
        Console.WriteLine("");
        Console.WriteLine("  " + 0.ToString().PadLeft(4)
                               + "  " + a.ToString().PadLeft(10)
                               + "  " + poly5(a).ToString().PadLeft(10)
                               + "  " + 0.0.ToString().PadLeft(10)
                               + "  " + fa.ToString().PadLeft(10)
                               + "  " + 0.0.ToString().PadLeft(10)
                               + "  " + fa.ToString().PadLeft(10) + "");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(10)
                                   + "  " + poly5(x[i]).ToString().PadLeft(10)
                                   + "  " + s[i].ToString().PadLeft(10)
                                   + "  " + f1[i].ToString().PadLeft(10)
                                   + "  " + f2[i].ToString().PadLeft(10)
                                   + "  " + f1[i] + f2[i].ToString().PadLeft(10) + "");
        }

        Console.WriteLine("  " + (n + 1).ToString().PadLeft(10)
                               + "  " + b.ToString().PadLeft(10)
                               + "  " + poly5(b).ToString().PadLeft(10)
                               + "  " + 0.0.ToString().PadLeft(10)
                               + "  " + fb.ToString().PadLeft(10)
                               + "  " + 0.0.ToString().PadLeft(10)
                               + "  " + fb.ToString().PadLeft(10) + "");
    }

    private static void sine_transform_test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINE_TRANSFORM_TEST03 evaluates the sine transform interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] f2;
        double fa;
        double fb;
        int i;
        int n = 9;
        int n2 = 1 + 2 * (n + 1);
        double[] s;
        double[] x;
        double[] x2;

        Console.WriteLine("");
        Console.WriteLine("SINE_TRANSFORM_TEST03:");
        Console.WriteLine("  SINE_TRANSFORM_FUNCTION does a sine transform of data");
        Console.WriteLine("  defined by a function F(X) evaluated at N equally spaced");
        Console.WriteLine("  points in an interval [A,B].");
        Console.WriteLine("  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.");
        Console.WriteLine("");
        Console.WriteLine("  The interpolant will be 0 at the 0th and (N+1)-th points.");
        Console.WriteLine("  It equals the function at points 1 through N.");
        Console.WriteLine("  In between, it can approximate smooth functions,");
        Console.WriteLine("  and the approximation improves with N.");
        //
        //  N determines the number of data points, indexed by 1 to N.  
        //  However, we essentially have N+2 data points, indexed 0 to N+1,
        //  with the data value being 0 at the first and last auxilliary points.
        //
        a = 1.0;
        b = 4.0;
        //
        //  Evenly spaced points between A and B, but omitting
        //  A and B themselves.
        //
        x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i) * a
                    + (i + 1) * b)
                   / (n + 1);
        }

        //
        //  Determine the interpolant coefficients.
        //
        s = Sine.sine_transform_function(n, a, b, poly5);

        Console.WriteLine("");
        Console.WriteLine("     I      X(I)      F(X(I))        S(I)");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(10)
                                   + "  " + poly5(x[i]).ToString().PadLeft(10)
                                   + "  " + s[i].ToString().PadLeft(10) + "");
        }

        //
        //  Evaluate the interpolant.
        //
        fa = poly5(a);
        fb = poly5(b);
        //
        //  Evenly spaced points between A and B, including A and B,
        //  and twice the density of the previous set of points.
        //
        x2 = new double[n2];

        for (i = 0; i < n2; i++)
        {
            x2[i] = ((n2 - i - 1) * a
                     + i * b)
                    / (n2 - 1);
        }

        f2 = Sine.sine_transform_interpolant(n, a, b, fa, fb, s, n2, x2);

        Console.WriteLine("");
        Console.WriteLine("     I      X            F(X)        FHAT(X)");
        Console.WriteLine("");

        for (i = 0; i < n2; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x2[i].ToString().PadLeft(10)
                                   + "  " + poly5(x2[i]).ToString().PadLeft(10)
                                   + "  " + f2[i].ToString().PadLeft(10) + "");
        }
    }

    private static void sine_transform_test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINE_TRANSFORM_TEST04 evaluates the sine transform interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] f2;
        double fa;
        double fb;
        int i;
        int n = 15;
        int n2 = 1 + 5 * (n + 1);
        double[] s;
        double[] x;
        double[] x2;

        Console.WriteLine("");
        Console.WriteLine("SINE_TRANSFORM_TEST04:");
        Console.WriteLine("  SINE_TRANSFORM_FUNCTION does a sine transform of data");
        Console.WriteLine("  defined by a function F(X) evaluated at N equally spaced");
        Console.WriteLine("  points in an interval [A,B].");
        Console.WriteLine("  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.");
        Console.WriteLine("");
        Console.WriteLine("  The interpolant will be 0 at the 0th and (N+1)-th points.");
        Console.WriteLine("  It equals the function at points 1 through N.");
        Console.WriteLine("  In between, it can approximate smooth functions,");
        Console.WriteLine("  and the approximation improves with N.");
        //
        //  N determines the number of data points, indexed by 1 to N.  
        //  However, we essentially have N+2 data points, indexed 0 to N+1,
        //  with the data value being 0 at the first and last auxilliary points.
        //
        a = 0.0;
        b = 7.0;
        //
        //  Evenly spaced points between A and B, but omitting
        //  A and B themselves.
        //
        x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i) * a
                    + (i + 1) * b)
                   / (n + 1);
        }

        //
        //  Determine the interpolant coefficients.
        //
        s = Sine.sine_transform_function(n, a, b, cosine_sum);
        //
        //  Evaluate the interpolant.
        //
        fa = cosine_sum(a);
        fb = cosine_sum(b);
        //
        //  Evenly spaced points between A and B, including A and B,
        //  and twice the density of the previous set of points.
        //
        x2 = new double[n2];

        for (i = 0; i < n2; i++)
        {
            x2[i] = ((n2 - i - 1) * a
                     + i * b)
                    / (n2 - 1);
        }

        f2 = Sine.sine_transform_interpolant(n, a, b, fa, fb, s, n2, x2);

        Console.WriteLine("");
        Console.WriteLine("  Expect exact agreement every 5th sample.");
        Console.WriteLine("");

        for (i = 0; i < n2; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x2[i].ToString().PadLeft(10)
                                   + "  " + cosine_sum(x2[i]).ToString().PadLeft(10)
                                   + "  " + f2[i].ToString().PadLeft(10) + "");
        }
    }

    private static double cosine_sum(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_SUM evaluates a function which is a cosine sum.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double COSINE_SUM, the value.
        //
    {
        double value = 0;

        value = Math.Cos(x)
                + 5.0 * Math.Cos(1.6 * x)
                - 2.0 * Math.Cos(2.0 * x)
                + 5.0 * Math.Cos(4.5 * x)
                + 7.0 * Math.Cos(9.0 * x);

        return value;
    }

    private static double poly5(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY5 evaluates a particular fifth-degree polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double POLY5, the value of the polynomial at X.
        //
    {
        double value = 0;

        value = (x - 0.1) *
                (x - 0.2) *
                (x - 0.4) *
                (x - 2.1) *
                (x - 3.0);

        return value;
    }
}