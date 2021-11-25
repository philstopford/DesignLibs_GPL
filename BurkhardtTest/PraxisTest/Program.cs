using System;
using System.Globalization;
using Burkardt.Praxis;
using Burkardt.Types;
using Burkardt.Uniform;

namespace PraxisTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PRAXIS_TEST.
        //
        //  Discussion:
        //
        //    PRAXIS_TEST tests the PRAXIS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("PRAXIS_TEST");
        Console.WriteLine("  Test the PRAXIS library.");
        //
        //  Minimization tests.
        //
        beale_test();
        box_test();
        chebyquad_test();
        cube_test();
        helix_test();
        hilbert_test();
        powell3d_test();
        rosenbrock_test();
        singular_test();
        tridiagonal_test();
        watson_test();
        wood_test();
        //
        //  Utility tests.
        //
        minfit_test();
        svsort_test();

        Console.WriteLine("");
        Console.WriteLine("PRAXIS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void beale_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BEALE_TEST calls PRAXIS for the Beale function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 2;
        double[] x = {0.1, 0.1};

        Console.WriteLine("");
        Console.WriteLine("BEALE_TEST");
        Console.WriteLine("  The Beale function.");

        const double t0 = 0.00001;
        const double h0 = 0.25;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + beale_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, beale_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + beale_f(x, n) + "");
    }

    private static double beale_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BEALE_F evaluates the Beale function.
        //
        //  Discussion:
        //
        //    The function is the sum of the squares of three functions.
        //
        //    This function has a valley approaching the line X(2) = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    E Beale,
        //    On an Iterative Method for Finding a Local Minimum of a Function
        //    of More than One Variable,
        //    Technical Report 25, Statistical Techniques Research Group,
        //    Princeton University, 1958.
        //
        //    Richard Brent,
        //    Algorithms for Finding Zeros and Extrema of Functions Without
        //    Calculating Derivatives,
        //    Stanford University Technical Report STAN-CS-71-198.
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double BEALE_F, the function value.
        //
    {
        const double c1 = 1.5;
        const double c2 = 2.25;
        const double c3 = 2.625;

        double fx1 = c1 - x[0] * (1.0 - x[1]);
        double fx2 = c2 - x[0] * (1.0 - Math.Pow(x[1], 2));
        double fx3 = c3 - x[0] * (1.0 - Math.Pow(x[1], 3));

        double value = fx1 * fx1 + fx2 * fx2 + fx3 * fx3;

        return value;
    }

    private static void box_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_TEST calls PRAXIS for the Box function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 3;
        double[] x = {0.0, 10.0, 20.0};

        Console.WriteLine("");
        Console.WriteLine("BOX_TEST");
        Console.WriteLine("  The Box function.");

        const double t0 = 0.00001;
        const double h0 = 20.0;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + box_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, box_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + box_f(x, n) + "");
    }

    private static double box_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_F evaluates the Box function.
        //
        //  Discussion:
        //
        //    The function is formed by the sum of squares of 10 separate terms.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double BOX_F, the function value.
        //
    {
        int i;

        double value = 0.0;

        for (i = 1; i <= 10; i++)
        {
            double c = -(double) i / 10.0;

            double fx = Math.Exp(c * x[0]) - Math.Exp(c * x[1]) - x[2] * (Math.Exp(c)
                                                                          - Math.Exp(10.0 * c));

            value += fx * fx;
        }

        return value;
    }

    private static void chebyquad_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYQUAD_TEST calls PRAXIS for the Chebyquad function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 8;
        double[] x = new double[8];

        Console.WriteLine("");
        Console.WriteLine("CHEBYQUAD_TEST");
        Console.WriteLine("  The Chebyquad function.");

        const double t0 = 0.00001;
        const double h0 = 0.1;
        const int prin = 0;

        for (i = 0; i < n; i++)
        {
            x[i] = (i + 1) / (double) (n + 1);
        }

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + chebyquad_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, chebyquad_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + chebyquad_f(x, n) + "");
    }

    private static double chebyquad_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYQUAD_F evaluates the Chebyquad function.
        //
        //  Discussion:
        //
        //    The function is formed by the sum of squares of N separate terms.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double CHEBYQUAD_F, the function value.
        //
    {
        int i;
        int j;

        double[] fvec = new double[n];

        for (i = 0; i < n; i++)
        {
            fvec[i] = 0.0;
        }

        for (j = 0; j < n; j++)
        {
            double t1 = 1.0;
            double t2 = 2.0 * x[j] - 1.0;
            double t = 2.0 * t2;

            for (i = 0; i < n; i++)
            {
                fvec[i] += t2;
                double th = t * t2 - t1;
                t1 = t2;
                t2 = th;
            }

        }

        for (i = 0; i < n; i++)
        {
            fvec[i] /= n;
            switch (i % 2)
            {
                case 1:
                    fvec[i] += 1.0 / (i * (i + 2));
                    break;
            }
        }

        //
        //  Compute F.
        //
        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            value += fvec[i] * fvec[i];
        }

        return value;
    }

    private static void cube_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_TEST calls PRAXIS for the Cube function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 2;
        double[] x = {-1.2, -1.0};

        Console.WriteLine("");
        Console.WriteLine("CUBE_TEST");
        Console.WriteLine("  The Cube function.");

        const double t0 = 0.00001;
        const double h0 = 1.0;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + cube_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, cube_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + cube_f(x, n) + "");
    }

    private static double cube_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_F evaluates the Cube function.
        //
        //  Discussion:
        //
        //    The function is the sum of the squares of two functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double CUBE_F, the function value.
        //
    {
        double fx1 = 10.0 * (x[1] - x[0] * x[0] * x[0]);
        double fx2 = 1.0 - x[0];

        double value = fx1 * fx1 + fx2 * fx2;

        return value;
    }

    private static void helix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HELIX_TEST calls PRAXIS for the Fletcher-Powell Helix function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 3;
        double[] x = {-1.0, 0.0, 0.0};

        Console.WriteLine("");
        Console.WriteLine("HELIX_TEST");
        Console.WriteLine("  The Fletcher-Powell Helix function.");

        const double t0 = 0.00001;
        const double h0 = 1.0;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + helix_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, helix_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + helix_f(x, n) + "");
    }

    private static double helix_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HELIX_F evaluates the Helix function.
        //
        //  Discussion:
        //
        //    The function is the sum of the squares of three functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double HELIX_F, the function value.
        //
    {
        double theta = 0;

        double r = Math.Sqrt(x[0] * x[0] + x[1] * x[1]);

        theta = x[0] switch
        {
            >= 0.0 => 0.5 * Math.Atan2(x[1], x[0]) / Math.PI,
            < 0.0 => 0.5 * (Math.Atan2(x[1], x[0]) + Math.PI) / Math.PI,
            _ => theta
        };

        double fx1 = 10.0 * (x[2] - 10.0 * theta);
        double fx2 = 10.0 * (r - 1.0);
        double fx3 = x[2];

        double value = fx1 * fx1 + fx2 * fx2 + fx3 * fx3;

        return value;
    }

    private static void hilbert_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HILBERT_TEST calls PRAXIS for the Hilbert function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 10;
        double[] x = new double[10];

        Console.WriteLine("");
        Console.WriteLine("HILBERT_TEST");
        Console.WriteLine("  The Hilbert function.");

        const double t0 = 0.00001;
        const double h0 = 10.0;
        const int prin = 0;

        for (i = 0; i < n; i++)
        {
            x[i] = 1.0;
        }

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + hilbert_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, hilbert_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + hilbert_f(x, n) + "");
    }

    private static double hilbert_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HILBERT_F evaluates the Hilbert function.
        //
        //  Discussion:
        //
        //    The function is a positive definite quadratic double of the form
        //
        //      f(x) = x" A x
        //
        //    where A is the Hilbert matrix, A(I,J) = 1/(I+J-1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double HILBERT_F, the function value.
        //
    {
        int i;

        double value = 0.0;

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                value += x[i] * x[j] / (i + j + 1);
            }
        }

        return value;
    }

    private static void minfit_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MINFIT_TEST tests MINFIT, which is a sort of SVD computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        const int n = 5;

        Console.WriteLine("");
        Console.WriteLine("MINFIT_TEST");
        Console.WriteLine("  MINFIT computes part of the SVD of a matrix A.");
        Console.WriteLine("    SVD: A = U * D * V'");
        Console.WriteLine("  MINFIT is given A, and returns the diagonal D");
        Console.WriteLine("  and the orthogonal matrix V.");

        double[] a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (j == i - 1 || j == i + 1)
                {
                    a[i + j * n] = -1.0;
                }
                else if (i == j)
                {
                    a[i + j * n] = 2.0;
                }
                else
                {
                    a[i + j * n] = 0.0;
                }
            }
        }

        typeMethods.r8mat_print(n, n, a, "  The matrix A:");

        double tol = Math.Sqrt(typeMethods.r8_epsilon());

        double[] d = new double[n];

        MINFIT.minfit(n, tol, ref a, ref d);

        typeMethods.r8mat_print(n, n, a, "  The vector V:");

        typeMethods.r8vec_print(n, d, "  The singular values D:");
        //
        //  Because A is positive definite symmetric, the "missing" matrix V = U.
        //
        Console.WriteLine("");
        Console.WriteLine("  Because A is positive definite symmetric,");
        Console.WriteLine("  we can reconstruct it as A = V * D * V'");

        double[] a2 = new double[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                a2[i + j * n] = 0.0;
                int k;
                for (k = 0; k < n; k++)
                {
                    a2[i + j * n] += a[i + k * n] * d[k] * a[j + k * n];
                }
            }
        }

        typeMethods.r8mat_print(n, n, a2, "  The product A2 = V * D * V'");
    }

    private static void powell3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWELL3D_TEST calls PRAXIS for the Powell 3D function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 3;
        double[] x = {0.0, 1.0, 2.0};

        Console.WriteLine("");
        Console.WriteLine("POWELL3D_TEST");
        Console.WriteLine("  The Powell 3D function.");

        const double t0 = 0.00001;
        const double h0 = 1.0;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + powell3d_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, powell3d_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + powell3d_f(x, n) + "");
    }

    private static double powell3d_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWELL3D_F evaluates the Powell 3D function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    M J D Powell,
        //    An Efficient Method for Finding the Minimum of a double of
        //    Several Variables Without Calculating Derivatives,
        //    Computer Journal, 
        //    Volume 7, Number 2, pages 155-162, 1964.
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double POWELL3D_F, the function value.
        //
    {
        double value = 0;

        value = 3.0 - 1.0 / (1.0 + Math.Pow(x[0] - x[1], 2))
                    - Math.Sin(0.5 * Math.PI * x[1] * x[2])
                    - Math.Exp(-Math.Pow((x[0] - 2.0 * x[1] + x[2]) / x[1], 2));

        return value;
    }

    private static void rosenbrock_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROSENBROCK_TEST calls PRAXIS for the Rosenbrock function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 2;
        double[] x = {-1.2, 1.0};

        Console.WriteLine("");
        Console.WriteLine("ROSENBROCK_TEST");
        Console.WriteLine("  The Rosenbrock function.");

        const double t0 = 0.00001;
        const double h0 = 1.0;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + rosenbrock_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, rosenbrock_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + rosenbrock_f(x, n) + "");
    }

    private static double rosenbrock_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROSENBROCK_F evaluates the Rosenbrock function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double ROSENBROCK_F, the function value.
        //
    {
        int j;
        double value = 0;

        value = 0.0;

        for (j = 0; j < n; j++)
        {
            value += (j % 2) switch
            {
                0 => Math.Pow(1.0 - x[j], 2),
                _ => 100.0 * Math.Pow(x[j] - x[j - 1] * x[j - 1], 2)
            };
        }

        return value;
    }

    private static void singular_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINGULAR_TEST calls PRAXIS for the Powell Singular function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 4;
        double[] x = {3.0, -1.0, 0.0, 1.0};

        Console.WriteLine("");
        Console.WriteLine("SINGULAR_TEST");
        Console.WriteLine("  The Powell Singular function.");

        const double t0 = 0.00001;
        const double h0 = 1.0;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + singular_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, singular_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + singular_f(x, n) + "");
    }

    private static double singular_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINGULAR_F evaluates the Powell Singular function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double SINGULAR_F, the function value.
        //
    {
        int j;

        double value = 0.0;

        for (j = 1; j <= n; j += 4)
        {
            double xjp1 = j + 1 <= n ? x[j] : 0.0;

            double xjp2 = j + 2 <= n ? x[j + 1] : 0.0;

            double xjp3 = j + 3 <= n ? x[j + 2] : 0.0;

            double f1 = x[j - 1] + 10.0 * xjp1;

            double f2;
            if (j + 1 <= n)
            {
                f2 = xjp2 - xjp3;
            }
            else
            {
                f2 = 0.0;
            }

            double f3;
            if (j + 2 <= n)
            {
                f3 = xjp1 - 2.0 * xjp2;
            }
            else
            {
                f3 = 0.0;
            }

            double f4;
            if (j + 3 <= n)
            {
                f4 = x[j - 1] - xjp3;
            }
            else
            {
                f4 = 0.0;
            }

            value = value
                    + Math.Pow(f1, 2)
                    + 5.0 * Math.Pow(f2, 2)
                    + Math.Pow(f3, 4)
                    + 10.0 * Math.Pow(f4, 4);
        }

        return value;
    }

    private static void tridiagonal_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIDIAGONAL_TEST calls PRAXIS for the Tridiagonal function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 4;
        double[] x = new double[4];

        Console.WriteLine("");
        Console.WriteLine("TRIDIAGONAL_TEST");
        Console.WriteLine("  The Tridiagonal function.");

        const double t0 = 0.00001;
        const double h0 = 8.0;
        const int prin = 0;

        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + tridiagonal_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, tridiagonal_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + tridiagonal_f(x, n) + "");
    }

    private static double tridiagonal_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIDIAGONAL_F evaluates the tridiagonal function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double TRIDIAGONAL_F, the function value.
        //
    {
        int i;
        double value = 0;

        value = x[0] * x[0];
        for (i = 1; i < n; i++)
        {
            value += 2.0 * x[i] * x[i];
        }

        for (i = 0; i < n - 1; i++)
        {
            value -= 2.0 * x[i] * x[i + 1];
        }

        value -= 2.0 * x[0];

        return value;
    }

    private static void watson_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WATSON_TEST calls PRAXIS for the Watson function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 6;
        double[] x = new double[6];

        Console.WriteLine("");
        Console.WriteLine("WATSON_TEST");
        Console.WriteLine("  The Watson function.");

        const double t0 = 0.00001;
        const double h0 = 1.0;
        const int prin = 0;

        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + watson_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, watson_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + watson_f(x, n) + "");
    }

    private static double watson_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WATSON_F evaluates the Watson function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double WATSON_F, the function value.
        //
    {
        int i;

        double value = 0.0;

        for (i = 1; i <= 29; i++)
        {
            double s1 = 0.0;
            double d = 1.0;
            int j;
            for (j = 1; j < n; j++)
            {
                s1 += j * d * x[j];
                d = d * i / 29.0;
            }

            double s2 = 0.0;
            d = 1.0;
            for (j = 0; j < n; j++)
            {
                s2 += d * x[j];
                d = d * i / 29.0;
            }

            value += Math.Pow(s1 - s2 * s2 - 1.0, 2);
        }

        value = value + x[0] * x[0] + Math.Pow(x[1] - x[0] * x[0] - 1.0, 2);

        return value;
    }

    private static void wood_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WOOD_TEST calls PRAXIS for the Wood function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 4;
        double[] x = {-3.0, -1.0, -3.0, -1.0};

        Console.WriteLine("");
        Console.WriteLine("WOOD_TEST");
        Console.WriteLine("  The Wood function.");

        const double t0 = 0.00001;
        const double h0 = 10.0;
        const int prin = 0;

        typeMethods.r8vec_print(n, x, "  Initial point:");

        Console.WriteLine("  Function value = " + wood_f(x, n) + "");

        PRAXIS.praxis(t0, h0, n, prin, ref x, wood_f);

        typeMethods.r8vec_print(n, x, "  Computed minimizer:");

        Console.WriteLine("  Function value = " + wood_f(x, n) + "");
    }

    private static double wood_f(double[] x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WOOD_F evaluates the Wood function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[N], the evaluation point.
        //
        //    Input, int N, the number of variables.
        //
        //    Output, double WOOD_F, the function value.
        //
    {
        double f1 = x[1] - x[0] * x[0];
        double f2 = 1.0 - x[0];
        double f3 = x[3] - x[2] * x[2];
        double f4 = 1.0 - x[2];
        double f5 = x[1] + x[3] - 2.0;
        double f6 = x[1] - x[3];

        double value = 100.0 * f1 * f1
                       + f2 * f2
                       + 90.0 * f3 * f3
                       + f4 * f4
                       + 10.0 * f5 * f5
                       + 0.1 * f6 * f6;

        return value;
    }

    private static void svsort_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SVSORT_TEST tests SVSORT, which sorts singular value information.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] d = new double[5];
        int i;
        int j;
        const int n = 5;
        double[] v = new double[5 * 5];

        Console.WriteLine("");
        Console.WriteLine("SVSORT_TEST");
        Console.WriteLine("  SVSORT sorts a vector D, and the corresponding columns");
        Console.WriteLine("  of a matrix V.");

        int seed = 123456789;

        for (i = 0; i < n; i++)
        {
            d[i] = UniformRNG.r8_uniform_01(ref seed);
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                v[i + j * n] = 10 * (i + 1) + j + 1;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  First row = entries of D.");
        Console.WriteLine("  Corresponding columns of V below.");
        Console.WriteLine("");
        string cout = "";
        for (j = 0; j < n; j++)
        {
            cout += d[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            cout = "";

            for (j = 0; j < n; j++)
            {
                cout += v[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        SVSORT.svsort(n, ref d, ref v);

        Console.WriteLine("");
        Console.WriteLine("  After sorting D and rearranging V:");
        Console.WriteLine("");
        cout = "";
        for (j = 0; j < n; j++)
        {
            cout += d[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            cout = "";
            for (j = 0; j < n; j++)
            {
                cout += v[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

    }
}