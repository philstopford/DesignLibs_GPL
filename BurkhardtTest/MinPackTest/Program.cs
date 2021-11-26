using System;
using System.Globalization;
using Burkardt.MinpackNS;
using Burkardt.SolveNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace MinPackTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MINPACK_TEST tests MINPACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("MINPACK_TEST");
        Console.WriteLine("  Test minpack().");

        chkder_test();
        hybrd1_test();
        qform_test();
        Console.WriteLine("");
        Console.WriteLine("MINPACK_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void chkder_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHKDER_TEST tests CHKDER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ido;

        const int m = 5;
        const int n = 5;
        double[] err = new double[m];
        double[] fjac = new double[n * n];
        double[] fvec = new double[m];
        double[] fvecp = new double[m];
        double[] x = new double[n];
        double[] xp = new double[n];

        Console.WriteLine("");
        Console.WriteLine("CHKDER_TEST");
        Console.WriteLine("  CHKDER compares a user supplied jacobian");
        Console.WriteLine("  and a finite difference approximation to it");
        Console.WriteLine("  and judges whether the jacobian is correct.");

        for (ido = 1; ido <= 2; ido++)
        {
            int seed = 123456789;
            switch (ido)
            {
                case 1:
                    Console.WriteLine("");
                    Console.WriteLine("  On the first test, use a correct jacobian.");
                    break;
                case 2:
                    Console.WriteLine("");
                    Console.WriteLine("  Repeat the test, but use a bad jacobian");
                    Console.WriteLine("  and see if the routine notices.");
                    break;
            }

            //
            //  Set the point at which the test is to be made:
            //
            int i;
            for (i = 0; i < n; i++)
            {
                x[i] = UniformRNG.r8_uniform_01(ref seed);
            }

            Console.WriteLine("");
            Console.WriteLine("  Evaluation point X:");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            int mode = 1;
            Minpack.chkder(m, n, x, fvec, fjac, n, xp, fvecp, mode, ref err);

            int iflag = 1;
            chkder_f(n, x, fvec, fjac, n, ref iflag);
            chkder_f(n, xp, fvecp, fjac, n, ref iflag);

            Console.WriteLine("");
            Console.WriteLine("  Sampled function values F(X) and F(XP)");
            Console.WriteLine("");
            for (i = 0; i < m; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + fvec[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + fvecp[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            iflag = 2;
            chkder_f(n, x, fvec, fjac, n, ref iflag);
            switch (ido)
            {
                //
                //  Here's where we put a mistake into the jacobian, on purpose.
                //
                case 2:
                    fjac[0 + 0 * n] = 1.01 * fjac[0 + 0 * n];
                    fjac[1 + 2 * n] = -fjac[1 + 2 * n];
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  Computed jacobian:");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < m; i++)
            {
                int j;
                for (j = 0; j < n; j++)
                {
                    cout += "  " + fjac[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(12);
                }

                Console.WriteLine(cout);
            }

            mode = 2;
            Minpack.chkder(m, n, x, fvec, fjac, n, xp, fvecp, mode, ref err);

            Console.WriteLine("");
            Console.WriteLine("  CHKDER gradient component error estimates:");
            Console.WriteLine("     > 0.5, the component is probably correct.");
            Console.WriteLine("     < 0.5, the component is probably incorrect.");
            Console.WriteLine("");
            for (i = 0; i < m; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + err[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void chkder_f(int n, double[] x, double[] fvec, double[] fjac, int ldfjac,
            ref int iflag )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHKDER_F is a function/jacobian routine for the CHKDER test.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of variables.
        //
        //    Input, double X[N], the variable values.
        //
        //    Output, double FVEC[N], the function values at X, if IFLAG = 1.
        //
        //    Output, double FJAC[LDFJAC*N], the N by N jacobian at X, if IFLAG = 2.
        //
        //    Input, int LDFJAC, the leading dimension of FJAC, which must
        //    be at least N.
        //
        //    Input, int &IFLAG:
        //    1, please compute F(I) (X).
        //    2, please compute FJAC(I,J) (X).
        //
    {
        int i;
        double x_prod;
        switch (iflag)
        {
            //
            //  If IFLAG is 1, we are supposed to evaluate F(X).
            //
            case 1:
            {
                double x_sum = 0.0;
                x_prod = 1.0;
                for (i = 0; i < n; i++)
                {
                    x_sum += x[i];
                    x_prod *= x[i];
                }

                for (i = 0; i < n - 1; i++)
                {
                    fvec[i] = x[i] - (n + 1) + x_sum;
                }

                fvec[n - 1] = x_prod - 1.0;
                break;
            }
            //
            //  If IFLAG is 2, we are supposed to evaluate FJAC(I,J) = d F(I)/d X(J)
            //
            case 2:
            {
                int j;
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < n - 1; i++)
                    {
                        fjac[i + j * ldfjac] = 1.0;
                    }
                }

                for (i = 0; i < n - 1; i++)
                {
                    fjac[i + i * ldfjac] = 2.0;
                }

                x_prod = 1.0;
                for (i = 0; i < n; i++)
                {
                    x_prod *= x[i];
                }

                for (j = 0; j < n; j++)
                {
                    fjac[n - 1 + j * ldfjac] = x_prod / x[j];
                }

                break;
            }
        }
    }

    private static void hybrd1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYBRD1_TEST tests HYBRD1.
        //
        //  Discussion:
        //
        //    This is an example of what your main program would look
        //    like if you wanted to use MINPACK to solve N nonlinear equations
        //    in N unknowns.  In this version, we avoid computing the jacobian
        //    matrix, and request that MINPACK approximate it for us.
        //
        //    The set of nonlinear equations is:
        //
        //      x1 * x1 - 10 * x1 + x2 * x2 + 8 = 0
        //      x1 * x2 * x2 + x1 - 10 * x2 + 8 = 0
        //
        //    with solution x1 = x2 = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 2;
        const double tol = 0.00001;

        const int lwa = n * (3 * n + 13) / 2;
        double[] fvec = new double[n];
        double[] wa = new double[lwa];
        double[] x = new double[n];

        Console.WriteLine("");
        Console.WriteLine("HYBRD1_TEST");
        Console.WriteLine("  HYBRD1 solves a nonlinear system of equations.");

        x[0] = 3.0;
        x[1] = 0.0;
        typeMethods.r8vec_print(n, x, "  Initial X");
        const int iflag = 1;
        fcnData res = hybrd1_f(n, x, fvec, iflag, 0, 0);
        fvec = res.fvec;

        typeMethods.r8vec_print(n, fvec, "  F(X)");

        int info = Minpack.hybrd1(hybrd1_f, n, x, fvec, tol, wa, lwa);

        Console.WriteLine("");
        Console.WriteLine("  Returned value of INFO = " + info + "");
        typeMethods.r8vec_print(n, x, "  X");
        typeMethods.r8vec_print(n, fvec, "  F(X)");
    }

    private static fcnData hybrd1_f(int n, double[] x, double[] fvec, int iflag, int arrayIndex, int array2Index )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYBRD1_F is a function routine for use with HYBRD1_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        fcnData ret = new()
        {
            iflag = iflag,
            fvec = fvec
        };
        fvec[0] = x[0] * x[0] - 10.0 * x[0] + x[1] * x[1] + 8.0;
        fvec[1] = x[0] * x[1] * x[1] + x[0] - 10.0 * x[1] + 8.0;

        return ret;
    }

    private static void qform_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QFORM_TEST tests QFORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        const int lda = 5;
        const int m = 5;
        const int n = 7;

        Console.WriteLine("");
        Console.WriteLine("QFORM_TEST:");
        Console.WriteLine("  QFORM constructs the Q factor explicitly");
        Console.WriteLine("  after the use of QRFAC.");
        //
        //  Set the matrix A.
        //
        double[] a = new double[m * n];

        int seed = 123456789;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        typeMethods.r8mat_print(m, n, a, "  Matrix A:");
        //
        //  Compute the QR factors.
        //
        const bool pivot = false;
        int[] ipivot = new int[n];
        int lipvt = n;
        double[] rdiag = new double[n];
        double[] acnorm = new double[n];

        QRSolve.qrfac(m, n, ref a, lda, pivot, ref ipivot, ref lipvt, ref rdiag, ref acnorm);
        //
        //  Extract the R factor.
        //
        double[] r = new double[m * n];
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                r[i + j * m] = 0.0;
            }
        }

        for (k = 0; k < Math.Min(m, n); k++)
        {
            r[k + k * m] = rdiag[k];
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < Math.Min(j, m); i++)
            {
                r[i + j * m] = a[i + j * m];
            }
        }

        typeMethods.r8mat_print(m, n, r, "  Matrix R:");
        //
        //  Call QRFORM to form the Q factor.
        //
        double[] q = new double[m * m];
        for (j = 0; j < m; j++)
        {
            for (i = 0; i < m; i++)
            {
                q[i + j * m] = 0.0;
            }
        }

        for (j = 0; j < Math.Min(m, n); j++)
        {
            for (i = 0; i < m; i++)
            {
                q[i + j * m] = a[i + j * m];
            }
        }

        QRSolve.qform(m, n, ref q, m);

        typeMethods.r8mat_print(m, m, q, "  Matrix Q:");
        //
        //  Compute Q*R.
        //
        double[] a2 = typeMethods.r8mat_mm_new(m, m, n, q, r);
        //
        //  Compare Q*R to A.
        //
        typeMethods.r8mat_print(m, n, a2, "  Matrix A2 = Q * R:");
    }
}