﻿using System;
using Burkardt.Sampling;
using Burkardt.Sequence;
using Burkardt.Types;

namespace WishartTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for WISHART_TEST.
        //
        //  Discussion:
        //
        //    WISHART_TEST tests the WISHART library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST");
        Console.WriteLine("  Test the WISHART library.");

        wishart_unit_sample_test();
        bartlett_unit_sample_test();
        wishart_test03();
        wishart_test04();
        wishart_test05();
        wishart_test06();
        wishart_test07();
        wishart_test08();
        wishart_test09();

        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void wishart_unit_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_UNIT_SAMPLE_TEST demonstrates the unit Wishart sampling function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        int rot_num = 0;

        Console.WriteLine("");
        Console.WriteLine("WISHART_UNIT_SAMPLE_TEST:");
        Console.WriteLine("  WISHART_UNIT_SAMPLE samples unit Wishart matrices by:");
        Console.WriteLine("  W = Wishart.wishart_unit_sample ( n, df );");
        //
        //  Set the parameters and call.
        //
        int n = 5;
        int df = 8;
        double[] w = Wishart.wishart_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_unit_sample ( 5, 8 ):");
        //
        //  Calling again yields a new matrix.
        //
        w = Wishart.wishart_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_unit_sample ( 5, 8 ):");
        //
        //  Reduce DF
        //
        n = 5;
        df = 5;
        w = Wishart.wishart_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_unit_sample ( 5, 5 ):");
        //
        //  Try a smaller matrix.
        //
        n = 3;
        df = 5;
        w = Wishart.wishart_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_unit_sample ( 3, 5 ):");
        //
        //  What is the eigendecomposition of the matrix?
        //
        int it_max = 50;
        double[] v = new double[n * n];
        double[] lambda = new double[n];

        Jacobi.jacobi_eigenvalue(n, w, it_max, ref v, ref lambda, ref it_num, ref rot_num);
        typeMethods.r8mat_print(n, n, v, "  Eigenvectors of previous matrix:");
        typeMethods.r8vec_print(n, lambda, "  Eigenvalues of previous matrix:");
    }

    private static void bartlett_unit_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BARTLETT_UNIT_SAMPLE_TEST demonstrates the unit Bartlett sampling function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        int rot_num = 0;

        Console.WriteLine("");
        Console.WriteLine("BARTLETT_UNIT_SAMPLE_TEST:");
        Console.WriteLine("  BARTLETT_UNIT_SAMPLE samples unit Bartlett matrices by:");
        Console.WriteLine("  T = Bartlett.bartlett_unit_sample ( n, df );");
        //
        //   Set the parameters and call.
        //
        int n = 5;
        int df = 8;
        double[] t = Bartlett.bartlett_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_unit_sample ( 5, 8 ):");
        //
        //   Calling again yields a new matrix.
        //
        t = Bartlett.bartlett_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_unit_sample ( 5, 8 ):");
        //
        //   Reduce DF.
        //
        n = 5;
        df = 5;
        t = Bartlett.bartlett_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_unit_sample ( 5, 5 ):");
        //
        //   Try a smaller matrix.
        //
        n = 3;
        df = 5;
        t = Bartlett.bartlett_unit_sample(n, df);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_unit_sample ( 3, 5 ):");
        //
        //   What is the eigendecomposition of the matrix T' * T?
        //
        double[] w = typeMethods.r8mat_mtm_new(n, n, n, t, t);

        int it_max = 50;
        double[] v = new double[n * n];
        double[] lambda = new double[n];

        Jacobi.jacobi_eigenvalue(n, w, it_max, ref v, ref lambda, ref it_num, ref rot_num);
        typeMethods.r8mat_print(n, n, v, "  Eigenvectors of previous matrix:");
        typeMethods.r8vec_print(n, lambda, "  Eigenvalues of previous matrix:");
    }

    private static void wishart_test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_TEST03 compares the unit Wishart and Bartlett sample matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST03:");
        Console.WriteLine("  Verify that, if using the same set of random numbers,");
        Console.WriteLine("    W = T' * T,");
        Console.WriteLine("  where");
        Console.WriteLine("    W = Wishart.wishart_unit_sample ( n, df );");
        Console.WriteLine("    T = Bartlett.bartlett_unit_sample ( n, df );");
        //
        //   Set the parameters.
        //
        int n = 5;
        int df = 8;
        double[] w = Wishart.wishart_unit_sample(n, df);
        double[] t = Bartlett.bartlett_unit_sample(n, df);
        //
        //   Compute T' * T.
        //
        double[] tt = typeMethods.r8mat_mtm_new(n, n, n, t, t);
        //
        //   Compare T'T to W.
        //
        double diff = typeMethods.r8mat_norm_fro_affine(n, n, w, tt);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of error is " + diff + "");
    }

    private static void wishart_test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_TEST04 demonstrates the Wishart sampling function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        int rot_num = 0;
        //
        //  Note that R is an upper triangular matrix,
        //  whose entries here are listed in column major order.
        //
        double[] r =  {
                5.0, 0.0, 0.0,
                1.0, 4.0, 0.0,
                3.0, 2.0, 6.0
            }
            ;
        double[] sigma_diag =  {
                1.0, 2.0, 3.0, 4.0, 5.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST04:");
        Console.WriteLine("  We can compute sample Wishart matrices by:");
        Console.WriteLine("    W = Wishart.wishart_sample ( n, df, sigma );");
        //
        //   Set the parameters and call.
        //
        int n = 5;
        int df = 8;
        double[] sigma = typeMethods.r8mat_identity_new(n);
        double[] w = Wishart.wishart_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_sample ( 5, 8, Identity ):");
        //
        //   Calling again yields a new matrix.
        //
        w = Wishart.wishart_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_sample ( 5, 8, Identity ):");
        //
        //   Try a diagonal matrix.
        //
        sigma = typeMethods.r8mat_diagonal_new(n, sigma_diag);
        w = Wishart.wishart_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_sample ( 5, 8, diag(1,2,3,4,5) ):");
        //
        //   Try a smaller matrix.  Sigma must be positive definite symmetric.
        //
        n = 3;
        df = 3;
        sigma = typeMethods.r8mat_mtm_new(n, n, n, r, r);
        typeMethods.r8mat_print(n, n, sigma, "  Set covariance SIGMA:");
        w = Wishart.wishart_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, w, "  Wishart.wishart_sample ( 3, 3, sigma ):");
        //
        //   What is the eigendecomposition of this matrix?
        //
        int it_max = 50;
        double[] v = new double[n * n];
        double[] lambda = new double[n];

        Jacobi.jacobi_eigenvalue(n, w, it_max, ref v, ref lambda, ref it_num, ref rot_num);
        typeMethods.r8mat_print(n, n, v, "  Eigenvectors of previous matrix:");
        typeMethods.r8vec_print(n, lambda, "  Eigenvalues of previous matrix:");
    }

    private static void wishart_test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_TEST05 demonstrates the Bartlett sampling function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        //
        //  Note that R is an upper triangular matrix,
        //  whose entries here are listed in column major order.
        //
        double[] r =  {
                5.0, 0.0, 0.0,
                1.0, 4.0, 0.0,
                3.0, 2.0, 6.0
            }
            ;
        int rot_num = 0;
        double[] sigma_diag =  {
                1.0, 2.0, 3.0, 4.0, 5.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST05:");
        Console.WriteLine("  We can compute sample Bartlett matrices by:");
        Console.WriteLine("    T = Bartlett.bartlett_sample ( n, df, sigma );");
        //
        //   Set the parameters and call.
        //
        int n = 5;
        int df = 8;
        double[] sigma = typeMethods.r8mat_identity_new(n);
        double[] t = Bartlett.bartlett_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_sample ( 5, 8, Identity ):");
        //
        //   Calling again yields a new matrix.
        //
        t = Bartlett.bartlett_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_sample ( 5, 8, Identity ):");
        //
        //   Try a diagonal matrix.
        //
        sigma = typeMethods.r8mat_diagonal_new(n, sigma_diag);
        t = Bartlett.bartlett_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_sample ( 5, 8, diag(1,2,3,4,5) ):");
        //
        //   Try a smaller matrix.
        //
        n = 3;
        df = 3;
        sigma = typeMethods.r8mat_mtm_new(n, n, n, r, r);
        typeMethods.r8mat_print(n, n, sigma, "  Set covariance SIGMA:");
        t = Bartlett.bartlett_sample(n, df, sigma);
        typeMethods.r8mat_print(n, n, t, "  Bartlett.bartlett_sample ( 3, 3, sigma ):");
        //
        //   What is the eigendecomposition of T' * T?
        //
        double[] w = typeMethods.r8mat_mtm_new(n, n, n, t, t);
        int it_max = 50;
        double[] v = new double[n * n];
        double[] lambda = new double[n];

        Jacobi.jacobi_eigenvalue(n, w, it_max, ref v, ref lambda, ref it_num, ref rot_num);
        typeMethods.r8mat_print(n, n, v, "  Eigenvectors of previous matrix:");
        typeMethods.r8vec_print(n, lambda, "  Eigenvalues of previous matrix:");
    }

    private static void wishart_test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_TEST06 compares the Wishart and Bartlett sample matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        //
        //  Note that R is an upper triangular matrix,
        //  whose entries here are listed in column major order.
        //
        double[] r =  {
                5.0, 0.0, 0.0,
                1.0, 4.0, 0.0,
                3.0, 2.0, 6.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST06:");
        Console.WriteLine("  Verify that, if using the same set of random numbers,");
        Console.WriteLine("    W = T'' * T,");
        Console.WriteLine("  where");
        Console.WriteLine("    W = Wishart.wishart_sample ( n, df, sigma );");
        Console.WriteLine("    T = Bartlett.bartlett_sample ( n, df, sigma );");
        //
        //   Set the parameters.
        //
        int n = 3;
        int df = 5;
        double[] sigma = typeMethods.r8mat_mtm_new(n, n, n, r, r);
        typeMethods.r8mat_print(n, n, sigma, "  Covariance SIGMA:");
        //
        //   Initialize the random number package and compute W.
        //
        double[] w = Wishart.wishart_sample(n, df, sigma);
        //
        //   Initialize the random number package again, and compute T.
        //
        double[] t = Bartlett.bartlett_sample(n, df, sigma);
        //
        //   Compute T' * T.
        //
        double[] tt = typeMethods.r8mat_mtm_new(n, n, n, t, t);
        //
        //   Compare T'T to W.
        //
        double diff = typeMethods.r8mat_norm_fro_affine(n, n, w, tt);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of error is " + diff + "");
    }

    private static void wishart_test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_TEST07 demonstrates a property of the Wishart distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        //
        //  Note that R is an upper triangular matrix,
        //  whose entries here are listed in column major order.
        //
        double[] r =  {
                5.0, 0.0, 0.0,
                1.0, 4.0, 0.0,
                3.0, 2.0, 6.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST07:");
        Console.WriteLine("  For given values of N, DF, SIGMA, the random");
        Console.WriteLine("  matrices from the Wishart distribution:");
        Console.WriteLine("    W = Wishart.wishart_sample ( n, df, sigma );");
        Console.WriteLine("  should have mean DF * SIGMA.");
        //
        //   Set the parameters.
        //
        int n = 3;
        Console.WriteLine("  Fix N = " + n + "");
        int df = 5;
        Console.WriteLine("  Fix DF = " + df + "");
        double[] sigma = typeMethods.r8mat_mtm_new(n, n, n, r, r);
        typeMethods.r8mat_print(n, n, sigma, "  Fix covariance SIGMA:");
        //
        //   Sample many times and average.
        //
        int sample_num = 1000;
        double[] w_average = typeMethods.r8mat_zero_new(n, n);
        for (i = 1; i <= sample_num; i++)
        {
            double[] w = Wishart.wishart_sample(n, df, sigma);
            typeMethods.r8mat_add(n, n, w, ref w_average);
        }

        double divisor = sample_num;
        typeMethods.r8mat_divide(n, n, divisor, ref w_average);
        //
        //   Compare SIGMA and W_SAMPLE / DF.
        //
        divisor = df;
        typeMethods.r8mat_divide(n, n, divisor, ref w_average);

        typeMethods.r8mat_print(n, n, w_average, "  W_Average / DF: ");

        double diff = typeMethods.r8mat_norm_fro_affine(n, n, sigma, w_average);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of SIGMA-W_average/DF = " + diff + "");
    }

    private static void wishart_test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_TEST08 samples the unit Wishart and unit Wishart inverse matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST08:");
        Console.WriteLine("  Verify that, if using the same set of random numbers,");
        Console.WriteLine("    inverse(W) = M,");
        Console.WriteLine("  where");
        Console.WriteLine("    W = Wishart.wishart_unit_sample ( n, df );");
        Console.WriteLine("    M = Wishart.wishart_unit_sample_inverse ( n, df );");
        //
        //   Set the parameters.
        //
        int n = 5;
        int df = 8;
        double[] w = Wishart.wishart_unit_sample(n, df);
        double[] m = Wishart.wishart_unit_sample_inverse(n, df);
        //
        //   Compute W * M.
        //
        double[] wm = typeMethods.r8mat_mm_new(n, n, n, w, m);
        //
        //   Compare M * W to I.
        //
        double[] ident = typeMethods.r8mat_identity_new(n);
        double diff = typeMethods.r8mat_norm_fro_affine(n, n, wm, ident);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of error is " + diff + "");
    }

    private static void wishart_test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_TEST09 samples the Wishart and Wishart inverse matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        //
        //  Note that R is an upper triangular matrix,
        //  whose entries here are listed in column major order.
        //
        double[] r =  {
                3.0, 0.0, 0.0, 0.0, 0.0,
                1.0, 7.0, 0.0, 0.0, 0.0,
                1.0, 1.0, 5.0, 0.0, 0.0,
                1.0, 2.0, 1.0, 4.0, 0.0,
                1.0, 3.0, 3.0, 2.0, 6.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("WISHART_TEST09:");
        Console.WriteLine("  Verify that, if using the same set of random numbers,");
        Console.WriteLine("    inverse(W) = M,");
        Console.WriteLine("  where");
        Console.WriteLine("    W = Wishart.wishart_sample ( n, df, sigma );");
        Console.WriteLine("    M = Wishart.wishart_sample_inverse ( n, df, sigma );");
        //
        //   Set the parameters.
        //
        int n = 5;
        int df = 8;
        double[] sigma = typeMethods.r8mat_mtm_new(n, n, n, r, r);

        double[] w = Wishart.wishart_sample(n, df, sigma);
        double[] m = Wishart.wishart_sample_inverse(n, df, sigma);
        //
        //   Compute W * M.
        //
        double[] wm = typeMethods.r8mat_mm_new(n, n, n, w, m);
        //
        //   Compare M * W to I.
        //
        double[] ident = typeMethods.r8mat_identity_new(n);
        double diff = typeMethods.r8mat_norm_fro_affine(n, n, wm, ident);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of error is " + diff + "");
    }
}