﻿using System;
using Burkardt.MatrixNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace VandermondeTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        /*
        Purpose:
        
        MAIN is the main program for VANDERMONDE_TEST.
        
        Discussion:
        
        VANDERMONDE_TEST tests the VANDERMONDE library.
        
        Licensing:
        
        This code is distributed under the GNU LGPL license.
        
        Modified:
        
        01 May 2014
        
        Author:
        
        John Burkardt
        */
    {
        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_TEST");
        Console.WriteLine("  Test the VANDERMONDE library.");

        bivand1_test();
        bivand2_test();
        dvand_test();
        dvandprg_test();
        pvand_test();
        pvandprg_test();

        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void bivand1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIVAND1_TEST tests BIVAND1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 3;

        double[] alpha = {1.0, 2.0, 3.0};
        double[] beta = {10.0, 20.0, 30.0};

        Console.WriteLine("");
        Console.WriteLine("BIVAND1_TEST:");
        Console.WriteLine("  Compute a bidimensional Vandermonde matrix");
        Console.WriteLine("  associated with polynomials of");
        Console.WriteLine("  total degree less than N.");

        typeMethods.r8vec_print(N, alpha, "  Vandermonde vector ALPHA:");
        typeMethods.r8vec_print(N, beta, "  Vandermonde vector BETA:");

        double[] a = VandermondeMatrix.bivand1(N, alpha, beta);

        const int n2 = N * (N + 1) / 2;
        typeMethods.r8mat_print(n2, n2, a, "  Bidimensional Vandermonde matrix:");

    }

    private static void bivand2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIVAND2_TEST tests BIVAND2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 3;

        double[] alpha = {1.0, 2.0, 3.0};
        double[] beta = {10.0, 20.0, 30.0};

        Console.WriteLine("");
        Console.WriteLine("BIVAND2_TEST:");
        Console.WriteLine("  Compute a bidimensional Vandermonde matrix");
        Console.WriteLine("  associated with polynomials of maximum degree less than N.");

        typeMethods.r8vec_print(N, alpha, "  Vandermonde vector ALPHA:");
        typeMethods.r8vec_print(N, beta, "  Vandermonde vector BETA:");

        double[] a = VandermondeMatrix.bivand2(N, alpha, beta);

        int n2 = N * N;
        typeMethods.r8mat_print(n2, n2, a, "  Bidimensional Vandermonde matrix:");

    }

    private static void dvand_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVAND_TEST tests DVAND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        double[] alpha = null;
        double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
        int seed = 12345;
        int test;
        double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

        Console.WriteLine("");
        Console.WriteLine("DVAND_TEST:");
        Console.WriteLine("  Solve a Vandermonde linear system A'*x=b");

        for (test = 1; test <= 2; test++)
        {
            alpha = test switch
            {
                1 => typeMethods.r8vec_copy_new(N, alpha1),
                2 => UniformRNG.r8vec_uniform_01_new(N, ref seed),
                _ => alpha
            };

            typeMethods.r8vec_print(N, alpha, "  Vandermonde vector ALPHA:");

            double[] a = VandermondeMatrix.vand1(N, alpha);

            double[] x = typeMethods.r8vec_copy_new(N, x1);
            double[] b = typeMethods.r8mat_mtv_new(N, N, a, x);
            typeMethods.r8vec_print(N, b, "  Right hand side B:");

            x = VandermondeMatrix.dvand(N, alpha, b);
            typeMethods.r8vec_print(N, x, "  Solution X:");
        }
    }

    private static void dvandprg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVANDPRG_TEST tests DVANDPRG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        double[] alpha = null;
        double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
        int test;
        double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

        Console.WriteLine("");
        Console.WriteLine("DVANDPRG_TEST:");
        Console.WriteLine("  Solve a Vandermonde linear system A'*x=b");
        Console.WriteLine("  progressively.");
        Console.WriteLine("  First we use ALPHA = 0, 1, 2, 3, 4.");
        Console.WriteLine("  Then we choose ALPHA at random.");

        for (test = 1; test <= 2; test++)
        {
            switch (test)
            {
                case 1:
                    alpha = typeMethods.r8vec_copy_new(N, alpha1);
                    break;
                case 2:
                    int seed = 123456789;
                    alpha = UniformRNG.r8vec_uniform_01_new(N, ref seed);
                    break;
            }

            typeMethods.r8vec_print(N, alpha, "  Vandermonde vector ALPHA:");

            double[] a = VandermondeMatrix.vand1(N, alpha);

            double[] x = typeMethods.r8vec_copy_new(N, x1);
            double[] b = typeMethods.r8mat_mtv_new(N, N, a, x);
            typeMethods.r8vec_print(N, b, "  Right hand side B:");

            x = new double[N];
            double[] c = new double[N];
            double[] m = new double[N];

            int nsub;
            for (nsub = 1; nsub <= N; nsub++)
            {
                VandermondeMatrix.dvandprg(nsub, alpha, b, ref x, ref c, ref m);
                typeMethods.r8vec_print(nsub, x, "  Solution X:");
            }
        }
    }

    private static void pvand_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PVAND_TEST tests PVAND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        double[] alpha = null;
        double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
        int test;
        double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

        Console.WriteLine("");
        Console.WriteLine("PVAND_TEST:");
        Console.WriteLine("  Solve a Vandermonde linear system A*x=b");

        for (test = 1; test <= 2; test++)
        {
            switch (test)
            {
                case 1:
                    alpha = typeMethods.r8vec_copy_new(N, alpha1);
                    break;
                case 2:
                    int seed = 123456789;
                    alpha = UniformRNG.r8vec_uniform_01_new(N, ref seed);
                    break;
            }

            typeMethods.r8vec_print(N, alpha, "  Vandermonde vector ALPHA:");

            double[] a = VandermondeMatrix.vand1(N, alpha);

            double[] x = typeMethods.r8vec_copy_new(N, x1);
            double[] b = typeMethods.r8mat_mv_new(N, N, a, x);
            typeMethods.r8vec_print(N, b, "  Right hand side B:");

            x = VandermondeMatrix.pvand(N, alpha, b);
            typeMethods.r8vec_print(N, x, "  Solution X:");

        }
    }

    private static void pvandprg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PVANDPRG_TEST tests PVANDPRG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        double[] alpha = null;
        double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
        int test;
        double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

        Console.WriteLine("");
        Console.WriteLine("PVANDPRG_TEST:");
        Console.WriteLine("  Solve a Vandermonde linear system A*x=b");

        for (test = 1; test <= 2; test++)
        {
            switch (test)
            {
                case 1:
                    alpha = typeMethods.r8vec_copy_new(N, alpha1);
                    break;
                case 2:
                    int seed = 123456789;
                    alpha = UniformRNG.r8vec_uniform_01_new(N, ref seed);
                    break;
            }

            typeMethods.r8vec_print(N, alpha, "  Vandermonde vector ALPHA:");

            double[] a = VandermondeMatrix.vand1(N, alpha);

            double[] x = typeMethods.r8vec_copy_new(N, x1);
            double[] b = typeMethods.r8mat_mv_new(N, N, a, x);
            typeMethods.r8vec_print(N, b, "  Right hand side B:");

            x = new double[N];
            double[] c = new double[N];
            double[] m = new double[N];

            int nsub;
            for (nsub = 1; nsub <= N; nsub++)
            {
                VandermondeMatrix.pvandprg(nsub, alpha, b, ref x, ref c, ref m);
                typeMethods.r8vec_print(nsub, x, "  Solution X:");
            }
        }

    }
}