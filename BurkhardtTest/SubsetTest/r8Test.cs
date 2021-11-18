using System;
using Burkardt;
using Burkardt.RationalNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTestNS;

public static class r8Test
{
    public static void r8_agm_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_AGM_TEST tests R8_AGM;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int seed;
        double x;
        double y;
        double z;

        Console.WriteLine("");
        Console.WriteLine("R8_AGM_TEST");
        Console.WriteLine("  R8_AGM computes the arithmetic-geometric mean (AGM)");
        Console.WriteLine("  of two nonnegative real numbers.");

        Console.WriteLine("");
        Console.WriteLine("    X        Y    R8_AGM(X,Y)");
        Console.WriteLine("");

        seed = 123456789;

        for (i = 1; i <= 10; i++)
        {
            x = UniformRNG.i4_uniform_ab(1, 10, ref seed);

            y = UniformRNG.i4_uniform_ab(1, 10, ref seed);

            z = typeMethods.r8_agm(x, y);

            Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                       + y.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                       + z.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    public static void r8_choose_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHOOSE_TEST tests R8_CHOOSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double cnk;
        int k;
        int n;

        Console.WriteLine("");
        Console.WriteLine("R8_CHOOSE_TEST");
        Console.WriteLine("  R8_CHOOSE evaluates C(N,K).");
        Console.WriteLine("");
        Console.WriteLine("     N     K    CNK");
        Console.WriteLine("");

        for (n = 0; n <= 4; n++)
        {
            for (k = 0; k <= n; k++)
            {
                cnk = typeMethods.r8_choose(n, k);

                Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                          + k.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                          + cnk.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
            }
        }
    }

    public static void r8_epsilon_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_EPSILON_TEST tests R8_EPSILON
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double r;
        double s;

        Console.WriteLine("");
        Console.WriteLine("R8_EPSILON_TEST");
        Console.WriteLine("  R8_EPSILON produces the floating point machine precision.");
        Console.WriteLine("");

        r = typeMethods.r8_epsilon();
        Console.WriteLine("  R = R8_EPSILON() = " + r.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        s = 1.0 + r - 1.0;
        Console.WriteLine("  ( 1 + R ) - 1 = " + s.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        s = 1.0 + r / 2.0 - 1.0;
        Console.WriteLine("  ( 1 + (R/2) ) - 1 = " + s.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

    }

    public static void r8_fall_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FALL_TEST tests R8_FALL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double f1 = 0;
        double f2 = 0;
        int n = 0;
        int n_data = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_FALL_TEST");
        Console.WriteLine("  R8_FALL evaluates the falling factorial Fall(X,N).");
        Console.WriteLine("");
        Console.WriteLine("    X          N                Exact                  Computed");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            typeMethods.r8_fall_values(ref n_data, ref x, ref n, ref f1);

            if (n_data == 0)
            {
                break;
            }

            f2 = typeMethods.r8_fall(x, n);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + f1.ToString("0.################").PadLeft(24) + "  "
                              + f2.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void r8_rise_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RISE_TEST tests R8_RISE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double f1 = 0;
        double f2 = 0;
        int n = 0;
        int n_data = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_RISE_TEST");
        Console.WriteLine("  R8_RISE evaluates the rising factorial Fall(X,N).");
        Console.WriteLine("");
        Console.WriteLine("    X          N                Exact                  Computed");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            typeMethods.r8_rise_values(ref n_data, ref x, ref n, ref f1);

            if (n_data == 0)
            {
                break;
            }

            f2 = typeMethods.r8_rise(x, n);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + f1.ToString("0.################").PadLeft(24) + "  "
                              + f2.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void r8_to_cfrac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_CFRAC_TEST tests R8_TO_CFRAC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        int[] a = new int[N + 1];
        double error;
        int i;
        int[] p = new int[N + 2];
        int[] q = new int[N + 2];
        double r;
        double temp;

        Console.WriteLine("");
        Console.WriteLine("R8_TO_CFRAC_TEST");
        Console.WriteLine("  R8_TO_CFRAC converts a double precision number to");
        Console.WriteLine("  a sequence of continued fraction convergents.");

        r = 2.0 * Math.PI;

        Console.WriteLine("");
        Console.WriteLine("  Use the real number R = " + r + "");

        typeMethods.r8_to_cfrac(r, N, ref a, ref p, ref q);

        Console.WriteLine("");

        for (i = 0; i <= N; i++)
        {
            temp = p[i + 1] / (double)q[i + 1];

            error = r - temp;

            Console.WriteLine("  "
                              + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + p[i + 1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + q[i + 1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + temp.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                              + error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void r8_to_dec_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_DEC_TEST tests R8_TO_DEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a = 0;
        int b = 0;
        int dec_digit;
        int i;
        double r;
        double r2;
        double r8_hi;
        double r8_lo;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("R8_TO_DEC_TEST");
        Console.WriteLine("  R8_TO_DEC converts a real number to a decimal;");

        dec_digit = 5;

        Console.WriteLine("");
        Console.WriteLine("  The maximum number of digits allowed is " + dec_digit + "");

        r8_lo = -10.0;
        r8_hi = +10.0;
        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("     R   =>  A * 10^B  =>  R2");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            r = UniformRNG.r8_uniform_ab(r8_lo, r8_hi, ref seed);

            typeMethods.r8_to_dec(r, dec_digit, ref a, ref b);
            r2 = typeMethods.dec_to_r8(a, b);

            Console.WriteLine("  "
                              + r.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + r2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    public static void r8_to_rat_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TO_RAT_TEST tests R8_TO_RAT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a = 0;
        int b = 0;
        int i;
        int ndig = 4;
        double r;
        double r2;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("R8_TO_RAT_TEST");
        Console.WriteLine("  R8_TO_RAT converts a real number to a rational;");
        Console.WriteLine("");
        Console.WriteLine("  The maximum number of digits allowed is " + ndig + "");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("     R   =>  A / B  =>  R2");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            r = UniformRNG.r8_uniform_01(ref seed);
            r = 10.0 * (r - 0.25);

            typeMethods.r8_to_rat(r, ndig, ref a, ref b);
            r2 = Rational.rat_to_r8(a, b);

            Console.WriteLine("  "
                              + r.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + r2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void r8mat_det_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DET_TEST tests R8MAT_DET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N3 = 3;
        int N4 = 4;

        double[] a3 = new double[N3 * N3];
        double[] a4 = new double[N4 * N4];
        double det;
        int i;
        int j;
        int k;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_DET_TEST");
        Console.WriteLine("  R8MAT_DET: determinant of a real matrix.");
        Console.WriteLine("");

        k = 0;
        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                k += 1;
                a3[i + j * N3] = k;
            }
        }

        typeMethods.r8mat_print(N3, N3, a3, "  The 123/456/789 matrix:");

        det = typeMethods.r8mat_det(N3, a3);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the 123/456/789 matrix is " + det + "");

        for (i = 0; i < N4; i++)
        {
            for (j = 0; j < N4; j++)
            {
                a4[i + j * N4] = 1.0 / (2 + i + j);
            }
        }

        typeMethods.r8mat_print(N4, N4, a4, "  The Hilbert matrix:");

        det = typeMethods.r8mat_det(N4, a4);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the Hilbert matrix is " + det + "");

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                if (i == j)
                {
                    a3[i + j * N3] = 2.0;
                }
                else if (i == j + 1 || i == j - 1)
                {
                    a3[i + j * N3] = -1.0;
                }
                else
                {
                    a3[i + j * N3] = 0.0;
                }
            }
        }

        typeMethods.r8mat_print(N3, N3, a3, "  The -1,2,-1 matrix:");

        det = typeMethods.r8mat_det(N3, a3);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the -1,2,-1 matrix is " + det + "");
    }

    public static void r8mat_perm0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PERM0_TEST tests R8MAT_PERM0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 9;

        double[] a = new double[N * N];
        int i;
        int[] p = { 1, 2, 8, 5, 6, 7, 4, 3, 0 };
        int j;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_PERM0_TEST");
        Console.WriteLine("  R8MAT_PERM0 reorders a real matrix in place.");
        Console.WriteLine("  The rows and columns use the same permutation.");

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                a[i + j * N] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print(N, N, a, "  The original matrix");

        Permutation.perm0_print(N, p, "  The row and column permutation:");

        typeMethods.r8mat_perm0(N, ref a, p);

        typeMethods.r8mat_print(N, N, a, "  The permuted matrix");

    }

    public static void r8mat_2perm0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_2PERM0_TEST tests R8MAT_2PERM0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 9;
        int N = 7;

        double[] a = new double[M * N];
        int i;
        int j;
        int[] p = { 1, 2, 8, 5, 6, 7, 4, 3, 0 };
        int[] q = { 2, 3, 4, 5, 6, 0, 1 };

        Console.WriteLine("");
        Console.WriteLine("R8MAT_2PERM0_TEST");
        Console.WriteLine("  R8MAT_2PERM0 reorders a real matrix in place.");
        Console.WriteLine("  Rows and columns use different permutations.");

        for (i = 0; i < M; i++)
        {
            for (j = 0; j < N; j++)
            {
                a[i + j * M] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print(M, N, a, "  The original matrix");

        Permutation.perm0_print(M, p, "  The row permutation:");

        Permutation.perm0_print(N, q, "  The column permutation:");

        typeMethods.r8mat_2perm0(M, N, ref a, p, q);

        typeMethods.r8mat_print(M, N, a, "  The permuted matrix");
    }

    public static void r8mat_permanent_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PERMANENT_TEST tests R8MAT_PERMANENT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int i;
        int j;
        int n;
        double perm;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_PERMANENT_TEST");
        Console.WriteLine("  R8MAT_PERMANENT: the matrix permanent function.");
        Console.WriteLine("  We will analyze matrices with 0 diagonal and");
        Console.WriteLine("  1 on all offdiagonals.");
        Console.WriteLine("");
        Console.WriteLine("  Order	    Permanent.");
        Console.WriteLine("");
        Console.WriteLine("DEBUG");
        for (n = 2; n <= 12; n++)
        {
            a = new double[n * n];

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        a[i + j * n] = 0.0;
                    }
                    else
                    {
                        a[i + j * n] = 1.0;
                    }
                }
            }

            perm = typeMethods.r8mat_permanent(n, a);

            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + perm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            Console.WriteLine("DEBUG, N = " + n + "");
        }
    }

    public static void r8poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_TEST test R8POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        int i;
        double[] a = new double[N];
        int iopt = 0;
        int test = 0;
        double val = 0;
        double x0 = 0;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_TEST");
        Console.WriteLine("  R8POLY converts between power sum, factorial");
        Console.WriteLine("  and Taylor forms, and can evaluate a polynomial");
        Console.WriteLine("");

        for (test = 1; test <= 6; test++)
        {
            switch (test)
            {
                case 1:
                    iopt = -3;
                    break;
                case 2:
                    iopt = -2;
                    break;
                case 3:
                    iopt = -1;
                    x0 = 2.0;
                    break;
                case 4:
                    iopt = 0;
                    x0 = 2.0;
                    break;
                case 5:
                    iopt = 6;
                    x0 = 2.0;
                    break;
                case 6:
                    iopt = 6;
                    x0 = -2.0;
                    break;
            }

            for (i = 0; i < N - 1; i++)
            {
                a[i] = 0.0;
            }

            a[N - 1] = 1.0;

            switch (test)
            {
                case 1:
                    typeMethods.r8vec_print(N, a, "  All calls  have input A as follows:");
                    break;
            }

            typeMethods.r8poly(N, ref a, x0, iopt, ref val);

            Console.WriteLine("");
            Console.WriteLine("  Option IOPT = " + iopt + "");

            switch (iopt)
            {
                case >= -1:
                    Console.WriteLine("  X0 = " + x0 + "");
                    break;
            }

            switch (iopt)
            {
                case -3:
                case -2:
                case > 0:
                    typeMethods.r8vec_print(N, a, "  Output array:");
                    break;
                case -1:
                case 0:
                    Console.WriteLine("  Value = " + val + "");
                    break;
            }
        }
            
    }

    public static void r8poly_f2p_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_F2P_TEST tests R8POLY_F2P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        double[] a;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_F2P_TEST");
        Console.WriteLine("  R8POLY_F2P: factorial => power sum.");

        a = typeMethods.r8vec_indicator1_new(N);

        typeMethods.r8poly_print(N - 1, a, "  The power sum polynomial:");

        typeMethods.r8poly_p2f(N, ref a);

        typeMethods.r8vec_print(N, a, "  The factorial coefficients:");

        typeMethods.r8poly_f2p(N, ref a);

        typeMethods.r8poly_print(N - 1, a, "  The recovered power sum polynomial:");
    }

    public static void r8poly_fval_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_FVAL_TEST tests R8POLY_FVAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        double[] a;
        double val;
        double x;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_FVAL_TEST");
        Console.WriteLine("  R8POLY_FVAL evaluates a polynomial in factorial form.");

        a = typeMethods.r8vec_indicator1_new(N);

        typeMethods.r8vec_print(N, a, "  The factorial coefficients:");

        x = 2.0;

        val = typeMethods.r8poly_fval(N, a, x);

        Console.WriteLine("");
        Console.WriteLine("  RPOLY (" + x + ") = " + val + "");
        Console.WriteLine("  The correct value is 11.");
    }

    public static void r8poly_n2p_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_N2P_TEST tests R8POLY_N2P and R8POLY_P2N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        int i;
        double[] a;
        double[] a2 = new double[N];

        Console.WriteLine("");
        Console.WriteLine("R8POLY_N2P_TEST");
        Console.WriteLine("  R8POLY_N2P: Newton => power sum;");

        a = typeMethods.r8vec_indicator1_new(N);

        for (i = 0; i < N; i++)
        {
            a2[i] = 2.0 * a[i];
        }

        typeMethods.r8poly_print(N - 1, a, "  The power sum polynomial:");

        typeMethods.r8poly_p2n(N, ref a, a2);

        typeMethods.r8vec_print(N, a, "  Derived Newton form coefficients:");

        typeMethods.r8vec_print(N, a2, "  Newton form abscissas:");

        typeMethods.r8poly_n2p(N, ref a, ref a2);

        typeMethods.r8poly_print(N - 1, a, "  The recovered power sum polynomial:");

    }

    public static void r8poly_nval_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_NVAL_TEST tests R8POLY_NVAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        double[] a;
        double[] a2 = new double[N];
        int i;
        double val;
        double x;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_NVAL_TEST");
        Console.WriteLine("  R8POLY_NVAL evaluates a polynomial in Newton form.");

        a = typeMethods.r8vec_indicator1_new(N);

        for (i = 0; i < N; i++)
        {
            a2[i] = a[i] - 1.0;
        }

        typeMethods.r8vec_print(N, a, "  Newton polynomial coefficients:");

        typeMethods.r8vec_print(N, a2, "  Newton polynomial abscissas:");

        x = 2.0;

        val = typeMethods.r8poly_nval(N, a, a2, x);

        Console.WriteLine("");
        Console.WriteLine("  RPOLY (" + x + ") = " + val + "");
        Console.WriteLine("  The correct value is 11.");

    }

    public static void r8poly_nx_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_NX_TEST tests R8POLY_NX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int n;
        double x;
        double[] xarray;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_NXL_TEST");
        Console.WriteLine("  R8POLY_NX replaces one abscissa in a Newton polynomial.");

        n = 3;

        a = typeMethods.r8vec_indicator1_new(n);
        xarray = typeMethods.r8vec_indicator1_new(n);

        typeMethods.r8vec_print(n, a, "  Newton polynomial coefficients:");
        typeMethods.r8vec_print(n, xarray, "  Newton polynomial abscissas:");
        /*
        Shift the X array by inserting X=0.
        */
        x = 0.0;
        Console.WriteLine("");
        Console.WriteLine("  Replace one abscissa by X = " + x + "");

        typeMethods.r8poly_nx(n, ref a, ref xarray, x);
        /*
        Report the new polynomial form.
        */
        typeMethods.r8vec_print(n, a, "  Revised Newton polynomial coefficients:");
        typeMethods.r8vec_print(n, xarray, "  Revised Newton polynomial abscissas:");
    }

    public static void r8poly_p2f_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_P2F_TEST tests R8POLY_P2F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        double[] a;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_P2F_TEST");
        Console.WriteLine("  R8POLY_P2F: power sum => factorial;");

        a = typeMethods.r8vec_indicator1_new(N);

        typeMethods.r8poly_print(N - 1, a, "  The power sum polynomial:");

        typeMethods.r8poly_p2f(N, ref a);

        typeMethods.r8vec_print(N, a, "  The factorial coefficients:");

        typeMethods.r8poly_f2p(N, ref a);

        typeMethods.r8poly_print(N - 1, a, "  The recovered power sum polynomial:");

    }

    public static void r8poly_p2n_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_P2N_TEST tests R8POLY_P2N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        int i;
        double[] a;
        double[] a2 = new double[N];

        Console.WriteLine("");
        Console.WriteLine("R8POLY_P2N_TEST");
        Console.WriteLine("  R8POLY_P2N: Power sum => Newton.");

        a = typeMethods.r8vec_indicator1_new(N);

        for (i = 0; i < N; i++)
        {
            a2[i] = 2.0 * a[i];
        }

        typeMethods.r8poly_print(N - 1, a, "  The power sum polynomial:");

        typeMethods.r8poly_p2n(N, ref a, a2);

        typeMethods.r8vec_print(N, a, "  Derived Newton form coefficients:");

        typeMethods.r8vec_print(N, a2, "  Newton form abscissas:");

        typeMethods.r8poly_n2p(N, ref a, ref a2);

        typeMethods.r8poly_print(N - 1, a, "  The recovered power sum polynomial:");
    }

    public static void r8poly_p2t_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_P2T_TEST tests R8POLY_P2T.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        double[] a;
        double x = 2.0;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_P2T_TEST");
        ;
        Console.WriteLine("  R8POLY_P2T: Power sum => Taylor.");
        Console.WriteLine("  The Taylor form uses the base point X0 = " + x + "");

        a = typeMethods.r8vec_indicator1_new(N + 1);

        typeMethods.r8vec_print(N, a, "  Initial Taylor sum form:");

        typeMethods.r8poly_t2p(N, ref a, x);

        typeMethods.r8poly_print(N, a, "  Power sum form:");

        typeMethods.r8poly_p2t(N, ref a, x);

        typeMethods.r8vec_print(N, a, "  Recovered Taylor sum form:");

    }

    public static void r8poly_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_PRINT_TEST tests R8POLY_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a = { -2.0, 5.1, 2.2, 3.3, 1.4 };
        int n = 4;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_PRINT_TEST");
        Console.WriteLine("  R8POLY_PRINT prints an R8POLY.");

        typeMethods.r8poly_print(n, a, "  The polynomial:");
    }

    public static void r8poly_pval_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_PVAL_TEST tests R8POLY_PVAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        double[] a;
        double val;
        double x;

        a = typeMethods.r8vec_indicator1_new(N + 1);

        Console.WriteLine("");
        Console.WriteLine("R8POLY_PVAL_TEST");
        Console.WriteLine("  R8POLY_PVAL evaluates a polynomial");
        Console.WriteLine("  in power sum form.");

        typeMethods.r8poly_print(N, a, "  The polynomial to be evaluated:");

        x = 2.0;

        val = typeMethods.r8poly_pval(N, a, x);

        Console.WriteLine("  At X = " + x + "");
        Console.WriteLine("  Computed polynomial value is " + val + "");
        Console.WriteLine("  Correct value is 129.");
    }

    public static void r8poly_t2p_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R89POLY_T2P_TEST tests R8POLY_T2P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        double[] a;
        double x = 2.0;

        Console.WriteLine("");
        Console.WriteLine("R8POLY_T2P_TEST");
        Console.WriteLine("  R8POLY_T2P: Taylor => Power sum;");
        Console.WriteLine("  The Taylor form uses the base point X0 = " + x + "");

        a = typeMethods.r8vec_indicator1_new(N + 1);

        typeMethods.r8vec_print(N, a, "  Initial Taylor sum form:");

        typeMethods.r8poly_t2p(N, ref a, x);

        typeMethods.r8poly_print(N, a, "  Power sum form:");

        typeMethods.r8poly_p2t(N, ref a, x);

        typeMethods.r8vec_print(N, a, "  Recovered Taylor sum form:");

    }

    public static void r8vec_backtrack_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BACKTRACK_TEST tests R8VEC_BACKTRACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int found_num;
        int i;
        int indx;
        int k;
        int n = 8;
        int maxstack = 100;
        int[] ncan = new int[8];
        int nstack;
        double[] stacks = new double[100];
        double t;
        double total;
        double[] w = { 15.0, 22.0, 14.0, 26.0, 32.0, 9.0, 16.0, 8.0 };
        double[] x = new double[8];

        Console.WriteLine("");
        Console.WriteLine("R8VEC_BACKTRACK_TEST");
        Console.WriteLine("  I4VEC_BACKTRACK uses backtracking, seeking a vector X of");
        Console.WriteLine("  N values which satisfies some condition.");

        Console.WriteLine("");
        Console.WriteLine("  In this demonstration, we have 8 values W(I).");
        Console.WriteLine("  We seek all subsets that sum to 53.0.");
        Console.WriteLine("  X(I) is 0.0 or 1.0 if the entry is skipped or used.");
        Console.WriteLine("");

        t = 53.0;

        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        indx = 0;
        k = 0;
        nstack = 0;
        for (i = 0; i < n; i++)
        {
            ncan[i] = 0;
        }

        found_num = 0;

        for (;;)
        {
            typeMethods.r8vec_backtrack(n, maxstack, stacks, ref x, ref indx, ref k, ref nstack, ref ncan);

            if (indx == 1)
            {
                found_num += 1;
                string cout = "  " + found_num.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "   ";

                total = typeMethods.r8vec_dot_product(n, w, x);
                cout += "  " + total.ToString(CultureInfo.InvariantCulture).PadLeft(8) + ":  ";

                for (i = 0; i < n; i++)
                {
                    switch (x[i])
                    {
                        case 1.0:
                            cout += "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                            break;
                    }
                }

                Console.WriteLine(cout);
            }
            //
            //  Given that we've chose X(1:K-1), what are our choices for X(K)?
            //
            //     if T < TOTAL, 
            //       no choices
            //     if T = TOTAL, 
            //       X(K) = 0
            //     if T > TOTAL and K < N, 
            //       X(k) = 0
            //       If ( W(K)+TOTAL <= T ) X(K) = 1
            //     If T > TOTAL and K = N,
            //       If ( W(K) + TOTAL) = T ) X(K) = 1
            //
            else if (indx == 2)
            {
                total = typeMethods.r8vec_dot_product(k - 1, w, x);

                if (t < total)
                {
                    ncan[k - 1] = 0;
                }
                else if (t == total)
                {
                    ncan[k - 1] += 1;
                    stacks[nstack] = 0.0;
                    nstack += 1;
                }
                else if (total < t && k < n)
                {
                    ncan[k - 1] += 1;
                    stacks[nstack] = 0.0;
                    nstack += 1;

                    if (total + w[k - 1] <= t)
                    {
                        ncan[k - 1] += 1;
                        stacks[nstack] = 1.0;
                        nstack += 1;
                    }
                }
                else if (total < t && k == n)
                {
                    if (total + w[k - 1] == t)
                    {
                        ncan[k - 1] += 1;
                        stacks[nstack] = 1.0;
                        nstack += 1;
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  Done!");
                break;
            }
        }
    }

    public static void r8vec_frac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_FRAC_TEST tests R8VEC_FRAC;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        double[] a = new double[N];
        double ahi = 10.0;
        double alo = 0.0;
        int k;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_FRAC_TEST");
        Console.WriteLine("  R8VEC_FRAC: K-th smallest real vector entry;");

        UniformRNG.r8vec_uniform(N, alo, ahi, ref seed, ref a);

        typeMethods.r8vec_print(N, a, "  The real array to search: ");

        Console.WriteLine("");
        Console.WriteLine("Frac     Value");
        Console.WriteLine("");

        for (k = 1; k <= N; k++)
        {
            Console.WriteLine(k.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + typeMethods.r8vec_frac(N, ref a, k).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }
    }

    public static void r8vec_mirror_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MIRROR_NEXT_TEST tests R8VEC_MIRROR_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        double[] a = new double[N];
        bool done;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_MIRROR_NEXT_TEST");
        Console.WriteLine("  R8VEC_MIRROR_NEXT generates all sign variations");
        Console.WriteLine("  of a real vector.");

        a[0] = 1.0;
        a[1] = 2.0;
        a[2] = 3.0;

        for (;;)
        {
            typeMethods.r8vec_print(N, a, "  Next vector:");

            done = typeMethods.r8vec_mirror_next(N, ref a);

            if (done)
            {
                Console.WriteLine("");
                Console.WriteLine("  Done.");
                break;
            }
        }

        a[0] = 1.0;
        a[1] = 0.0;
        a[2] = 3.0;

        for (;;)
        {
            typeMethods.r8vec_print(N, a, "  Next vector:");

            done = typeMethods.r8vec_mirror_next(N, ref a);

            if (done)
            {
                Console.WriteLine("");
                Console.WriteLine("  Done.");
                break;
            }
        }
    }

}