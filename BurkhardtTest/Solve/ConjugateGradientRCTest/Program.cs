using System;
using Burkardt.MatrixNS;
using Burkardt.SolveNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ConjugateGradientRCTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CG_RC_TEST.
        //
        //  Discussion:
        //
        //    CG_RC_TEST tests the CG_RC library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CG_RC_TEST:");

        Console.WriteLine("  Test the CG_RC library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("CG_RC_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses CG_RC for the simple 1, -2, 1 matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double angle;
        double[] b;
        double bnrm2;
        double err;
        int i;
        int it;
        int it_max;
        int job;
        int n = 21;
        double[] p;
        double pi = 3.141592653589793;
        double[] q;
        double[] r;
        double rnrm2;
        double t;
        double tol;
        double[] x;
        double[] x_exact;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use CG_RC on the 1, -2, 1 matrix.");
        //
        //  In order to specify the right hand side, pick an exact solution,
        //  and multiply by the matrix.
        //
        x_exact = new double[n];
        for (i = 0; i < n; i++)
        {
            angle = 2.0 * pi * i / (n - 1);
            x_exact[i] = Math.Sin(angle);
        }

        b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = -2.0 * x_exact[i];
        }

        for (i = 0; i < n - 1; i++)
        {
            b[i] += x_exact[i + 1];
        }

        for (i = 1; i < n; i++)
        {
            b[i] += x_exact[i - 1];
        }

        //
        //  Here is the initial guess for the solution.
        //
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        //
        //  Parameters for the stopping test.
        //
        it = 0;
        it_max = 30;
        tol = 1.0E-05;
        bnrm2 = 0.0;
        for (i = 0; i < n; i++)
        {
            bnrm2 += b[i] * b[i];
        }

        bnrm2 = Math.Sqrt(bnrm2);
        //
        //  Set parameters for the CG_RC code.
        //
        r = new double[n];
        z = new double[n];
        p = new double[n];
        q = new double[n];

        job = 1;
        //
        //  Repeatedly call the CG_RC code, and on return, do what JOB tells you.
        //
        ConjugateGradientData data = new();
        for (;;)
        {
            job = ConjugateGradientRC.cg_rc(ref data, n, b, ref x, ref r, ref z, ref p, ref q, ref job);
            //
            //  Compute q = A * p.
            //
            if (job == 1)
            {
                for (i = 0; i < n; i++)
                {
                    q[i] = -2.0 * p[i];
                }

                for (i = 0; i < n - 1; i++)
                {
                    q[i] += p[i + 1];
                }

                for (i = 1; i < n; i++)
                {
                    q[i] += p[i - 1];
                }
            }
            //
            //  Solve M * z = r.
            //
            else if (job == 2)
            {
                for (i = 0; i < n; i++)
                {
                    z[i] = r[i] / -2.0;
                }
            }
            //
            //  Compute r = r - A * x.
            //
            else if (job == 3)
            {
                for (i = 0; i < n; i++)
                {
                    r[i] += 2.0 * x[i];
                }

                for (i = 0; i < n - 1; i++)
                {
                    r[i] -= x[i + 1];
                }

                for (i = 1; i < n; i++)
                {
                    r[i] -= x[i - 1];
                }
            }
            //
            //  Stopping test on R.
            //
            else if (job == 4)
            {
                rnrm2 = 0.0;
                for (i = 0; i < n; i++)
                {
                    rnrm2 += r[i] * r[i];
                }

                rnrm2 = Math.Sqrt(rnrm2);

                if (bnrm2 == 0.0)
                {
                    if (rnrm2 <= tol)
                    {
                        break;
                    }
                }
                else
                {
                    if (rnrm2 <= tol * bnrm2)
                    {
                        break;
                    }
                }

                it += 1;

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Iteration limit exceeded.");
                    Console.WriteLine("  Terminating early.");
                    break;
                }
            }

            job = 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations was " + it + "");
        Console.WriteLine("  Estimated error is " + rnrm2 + "");
        err = 0.0;
        for (i = 0; i < n; i++)
        {
            t = Math.Abs(x_exact[i] - x[i]);
            if (err < t)
            {
                err = t;
            }
        }

        Console.WriteLine("  Loo error is " + err + "");

        Console.WriteLine("");
        Console.WriteLine("     I      X(I)         X_EXACT(I)        B(I)");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(14)
                                   + "  " + x_exact[i].ToString().PadLeft(14)
                                   + "  " + b[i].ToString().PadLeft(14) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests CG_RC with the Wathen matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] ax;
        double[] b;
        double bnrm2;
        double err;
        int i;
        int it;
        int it_max;
        int job;
        int n;
        int nx;
        int ny;
        double[] p;
        double[] q;
        double[] r;
        double rnrm2;
        int seed;
        double t;
        double tol;
        double[] x;
        double[] x_exact;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Use CG_RC to solve a linear system");
        Console.WriteLine("  involving the Wathen matrix.");

        nx = 5;
        ny = 4;

        Console.WriteLine("");
        Console.WriteLine("  NX = " + nx + "");
        Console.WriteLine("  NY = " + ny + "");

        n = WathenMatrix.wathen_order(nx, ny);

        Console.WriteLine("  N  = " + n + "");

        a = WathenMatrix.wathen(nx, ny, n);

        seed = 123456789;
        x_exact = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        b = new double[n];
        typeMethods.r8mat_mv(n, n, a, x_exact, ref b);
        //
        //  Here is the initial guess for the solution.
        //
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        ax = new double[n];
        //
        //  Parameters for the stopping test.
        //
        it = 0;
        it_max = 30;
        tol = 1.0E-05;
        bnrm2 = 0.0;
        for (i = 0; i < n; i++)
        {
            bnrm2 += b[i] * b[i];
        }

        bnrm2 = Math.Sqrt(bnrm2);
        //
        //  Set parameters for the CG_RC code.
        //
        r = new double[n];
        z = new double[n];
        p = new double[n];
        q = new double[n];
        job = 1;
        //
        //  Repeatedly call the CG_RC code, and on return, do what JOB tells you.
        //
        ConjugateGradientData data = new();
        for (;;)
        {
            job = ConjugateGradientRC.cg_rc(ref data, n, b, ref x, ref r, ref z, ref p, ref q, ref job);
            //
            //  Compute q = A * p.
            //
            if (job == 1)
            {
                typeMethods.r8mat_mv(n, n, a, p, ref q);
            }
            //
            //  Solve M * z = r.
            //
            else if (job == 2)
            {
                for (i = 0; i < n; i++)
                {
                    z[i] = r[i] / a[i + i * n];
                }
            }
            //
            //  Compute r = r - A * x.
            //
            else if (job == 3)
            {
                typeMethods.r8mat_mv(n, n, a, x, ref ax);
                for (i = 0; i < n; i++)
                {
                    r[i] -= ax[i];
                }
            }
            //
            //  Stopping test.
            //
            else if (job == 4)
            {
                rnrm2 = 0.0;
                for (i = 0; i < n; i++)
                {
                    rnrm2 += r[i] * r[i];
                }

                rnrm2 = Math.Sqrt(rnrm2);

                if (bnrm2 == 0.0)
                {
                    if (rnrm2 <= tol)
                    {
                        break;
                    }
                }
                else
                {
                    if (rnrm2 <= tol * bnrm2)
                    {
                        break;
                    }
                }

                it += 1;

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Iteration limit exceeded.");
                    Console.WriteLine("  Terminating early.");
                    break;
                }
            }

            job = 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations was " + it + "");
        Console.WriteLine("  Estimated error is " + rnrm2 + "");
        err = 0.0;
        for (i = 0; i < n; i++)
        {
            t = Math.Abs(x_exact[i] - x[i]);
            if (err < t)
            {
                err = t;
            }
        }

        Console.WriteLine("  Loo error is " + err + "");

        Console.WriteLine("");
        Console.WriteLine("     I      X(I)         X_EXACT(I)        B(I)");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(14)
                                   + "  " + x_exact[i].ToString().PadLeft(14)
                                   + "  " + b[i].ToString().PadLeft(14) + "");
        }
    }
}