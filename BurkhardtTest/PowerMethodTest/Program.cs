﻿using System;
using System.Numerics;
using Burkardt.MatrixNS;
using Burkardt.PowerMethodNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace PowerMethodTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POWER_METHOD_TEST.
        //
        //  Discussion:
        //
        //    POWER_METHOD_TEST tests the POWER_METHOD library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("POWER_METHOD_TEST");
        Console.WriteLine("  Test the POWER_METHOD library.");

        test01();
        test02();
        test03();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("POWER_METHOD_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses POWER_METHOD on the Fibonacci2 matrix.
        //
        //  Discussion:
        //
        //    This matrix, despite having a single dominant eigenvalue, will generally
        //    converge only very slowly under the power method.  This has to do with
        //    the fact that the matrix has only 3 eigenvectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int it_num = 0;
        double lambda = 0;
        const int n = 50;
        double[] x = new double[n];

        double[] a = Fibonacci.fibonacci2(n);

        int seed = 123456789;
        UniformRNG.r8vec_uniform_01(n, ref seed, ref x);

        const int it_max = 300;
        const double tol = 0.000001;

        double phi = (1.0 + Math.Sqrt(5.0)) / 2.0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use POWER_METHOD on the Fibonacci2 matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N       = " + n + "");
        Console.WriteLine("  Maximum iterations   = " + it_max + "");
        Console.WriteLine("  Error tolerance      = " + tol + "");

        DateTime ctime1 = DateTime.Now;

        PowerMethod.power_method(n, a, x, it_max, tol, ref lambda, ref it_num);

        DateTime ctime2 = DateTime.Now;
        double ctime = (ctime2 - ctime1).Seconds;

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + it_num + "");
        Console.WriteLine("  CPU time             = " + ctime + "");
        Console.WriteLine("  Estimated eigenvalue = " + lambda.ToString("0.##############") + "");
        Console.WriteLine("  Correct value        = " + phi.ToString("0.##############") + "");
        Console.WriteLine("  Error                = " + Math.Abs(lambda - phi) + "");
        //
        //  X2 is the exact eigenvector.
        //
        double[] x2 = new double[n];

        x2[0] = 1.0;
        for (i = 1; i < n; i++)
        {
            x2[i] = phi * x2[i - 1];
        }

        double norm = typeMethods.r8vec_norm_l2(n, x2);
        for (i = 0; i < n; i++)
        {
            x2[i] /= norm;
        }

        //
        //  The sine of the angle between X and X2 is a measure of error.
        //
        double cos_x1x2 = typeMethods.r8vec_dot(n, x, x2);
        double sin_x1x2 = Math.Sqrt((1.0 - cos_x1x2) * (1.0 + cos_x1x2));

        Console.WriteLine("");
        Console.WriteLine("  Sine of angle between true and estimated vectors = " + sin_x1x2 + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses POWER_METHOD2 on the Fibonacci2 matrix.
        //
        //  Discussion:
        //
        //    This matrix, despite having a single dominant eigenvalue, will generally
        //    converge only very slowly under the power method.  This has to do with
        //    the fact that the matrix has only 3 eigenvectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        Complex lambda = new();
        const int n = 50;
        double[] x = new double[n];

        double[] a = Fibonacci.fibonacci2(n);
        Complex[] v = new Complex [n];

        int seed = 123456789;
        UniformRNG.r8vec_uniform_01(n, ref seed, ref x);

        const int it_max = 300;
        const double tol = 0.000001;

        double phi = (1.0 + Math.Sqrt(5.0)) / 2.0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Use POWER_METHOD2 on the Fibonacci2 matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N       = " + n + "");
        Console.WriteLine("  Maximum iterations   = " + it_max + "");
        Console.WriteLine("  Error tolerance      = " + tol + "");

        DateTime ctime1 = DateTime.Now;

        PowerMethod.power_method2(n, a, x, it_max, tol, ref lambda, v, ref it_num);

        DateTime ctime2 = DateTime.Now;
        double ctime = (ctime2 - ctime1).Seconds;

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + it_num + "");
        Console.WriteLine("  CPU time             = " + ctime + "");
        Console.WriteLine("  Estimated eigenvalue = "
                          + "  " + lambda.Real.ToString("0.##############")
                          + "  " + lambda.Imaginary.ToString("0.##############") + "");
        Console.WriteLine("  Correct value        = " + phi.ToString("0.##############") + "");
        Console.WriteLine("  Error                = " + Complex.Abs(lambda - phi) + "");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses POWER_METHOD2 on the TRIS matrix.
        //
        //  Discussion:
        //
        //    This matrix, despite having a single dominant eigenvalue, will generally
        //    converge only very slowly under the power method.  This has to do with
        //    the fact that the matrix has only 3 eigenvectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int it_num = 0;
        Complex lambda = new();
        const int n = 50;
        double[] x = new double[n];

        const double alpha = -1.0;
        const double beta = 10.0;
        const double gamma = 8.0;

        double[] a = Tridiagonal.tris(n, n, alpha, beta, gamma);

        Complex[] v = new Complex[n];

        int seed = 123456789;
        UniformRNG.r8vec_uniform_01(n, ref seed, ref x);

        const int it_max = 4000;
        const double tol = 0.000001;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Use POWER_METHOD2 on the TRIS (tridiagonal scalar) matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N         = " + n + "");
        Console.WriteLine("  Maximum iterations     = " + it_max + "");
        Console.WriteLine("  Error tolerance        = " + tol + "");

        DateTime ctime1 = DateTime.Now;

        PowerMethod.power_method2(n, a, x, it_max, tol, ref lambda, v, ref it_num);

        DateTime ctime2 = DateTime.Now;
        double ctime = (ctime2 - ctime1).Seconds;

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations   = " + it_num + "");
        Console.WriteLine("  CPU time               = " + ctime + "");
        Console.WriteLine("  Estimated eigenvalue   = "
                          + lambda.Real.ToString("0.##############")
                          + "  " + lambda.Imaginary.ToString("0.##############") + "");

        Complex[] lambda_vec = Tridiagonal.tris_eigenvalues(n, alpha, beta, gamma);

        Complex lambda_max = lambda_vec[0];
        for (i = 1; i < n; i++)
        {
            if (Complex.Abs(lambda_max) < Complex.Abs(lambda_vec[i]))
            {
                lambda_max = lambda_vec[i];
            }
        }

        Console.WriteLine("  Correct max eigenvalue = "
                          + lambda_max.Real.ToString("0.##############")
                          + "  " + lambda_max.Imaginary.ToString("0.##############") + "");

        Console.WriteLine("  Error                  = " + Complex.Abs(lambda - lambda_max) + "");
    }
}