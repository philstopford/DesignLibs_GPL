using System;
using Burkardt.Probability;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace QRSolveTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for QR_SOLVE_TEST.
        //
        //  Discussion:
        //
        //    QR_SOLVE_TEST tests the QR_SOLVE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("QR_SOLVE_TEST");
        Console.WriteLine("  Test the QR_SOLVE library.");
        Console.WriteLine("  QR_SOLVE needs the R8LIB library.");
        Console.WriteLine("  This test also needs the TEST_LS library.");

        normal_solve_test();
        qr_solve_test();
        svd_solve_test();
        dqrls_test();

        Console.WriteLine("");
        Console.WriteLine("QR_SOLVE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void normal_solve_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_SOLVE_TEST tests NORMAL_SOLVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b;
        double b_norm;
        int flag = 0;
        int i;
        int m;
        int n;
        int prob;
        int prob_num;
        double[] r1;
        double r1_norm;
        double[] r2;
        double r2_norm;
        double x_diff_norm;
        double[] x1;
        double x1_norm;
        double[] x2;
        double x2_norm;

        Console.WriteLine("");
        Console.WriteLine("NORMAL_SOLVE_TEST");
        Console.WriteLine("  NORMAL_SOLVE is a function with a simple interface which");
        Console.WriteLine("  solves a linear system A*x = b in the least squares sense.");
        Console.WriteLine("  Compare a tabulated solution X1 to the NORMAL_SOLVE result X2.");
        Console.WriteLine("");
        Console.WriteLine("  NORMAL_SOLVE cannot be applied when N < M,");
        Console.WriteLine("  or if the matrix does not have full column rank.");

        prob_num = ProbabilityFunctions.p00_prob_num();

        Console.WriteLine("");
        Console.WriteLine("  Number of problems = " + prob_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||");
        Console.WriteLine("");

        for (prob = 1; prob <= prob_num; prob++)
        {
            //
            //  Get problem size.
            //
            m = ProbabilityFunctions.p00_m(prob);
            n = ProbabilityFunctions.p00_n(prob);
            //
            //  Retrieve problem data.
            //
            a = ProbabilityFunctions.p00_a(prob, m, n);
            b = ProbabilityFunctions.p00_b(prob, m);
            x1 = ProbabilityFunctions.p00_x(prob, n);

            b_norm = typeMethods.r8vec_norm(m, b);
            x1_norm = typeMethods.r8vec_norm(n, x1);
            r1 = typeMethods.r8mat_mv_new(m, n, a, x1);
            for (i = 0; i < m; i++)
            {
                r1[i] -= b[i];
            }

            r1_norm = typeMethods.r8vec_norm(m, r1);
            //
            //  Use NORMAL_SOLVE on the problem.
            //
            x2 = NormalSolve.normal_solve(m, n, a, b, ref flag);

            if (flag != 0)
            {
                Console.WriteLine("  " + prob.ToString().PadLeft(5)
                                       + "  " + m.ToString().PadLeft(4)
                                       + "  " + n.ToString().PadLeft(4)
                                       + "  " + b_norm.ToString().PadLeft(12)
                                       + "  " + "------------"
                                       + "  " + x1_norm.ToString().PadLeft(12)
                                       + "  " + "------------"
                                       + "  " + r1_norm.ToString().PadLeft(12)
                                       + "  " + "------------" + "");
            }
            else
            {
                x2_norm = typeMethods.r8vec_norm(n, x2);
                r2 = typeMethods.r8mat_mv_new(m, n, a, x2);
                for (i = 0; i < m; i++)
                {
                    r2[i] -= b[i];
                }

                r2_norm = typeMethods.r8vec_norm(m, r2);
                //
                //  Compare tabulated and computed solutions.
                //
                x_diff_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
                //
                //  Report results for this problem.
                //
                Console.WriteLine("  " + prob.ToString().PadLeft(5)
                                       + "  " + m.ToString().PadLeft(4)
                                       + "  " + n.ToString().PadLeft(4)
                                       + "  " + b_norm.ToString().PadLeft(12)
                                       + "  " + x_diff_norm.ToString().PadLeft(12)
                                       + "  " + x1_norm.ToString().PadLeft(12)
                                       + "  " + x2_norm.ToString().PadLeft(12)
                                       + "  " + r1_norm.ToString().PadLeft(12)
                                       + "  " + r2_norm.ToString().PadLeft(12) + "");

            }

        }

    }

    private static void qr_solve_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QR_SOLVE_TEST tests QR_SOLVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b;
        double b_norm;
        int i;
        int m;
        int n;
        int prob;
        int prob_num;
        double[] r1;
        double r1_norm;
        double[] r2;
        double r2_norm;
        double x_diff_norm;
        double[] x1;
        double x1_norm;
        double[] x2;
        double x2_norm;

        Console.WriteLine("");
        Console.WriteLine("QR_SOLVE_TEST");
        Console.WriteLine("  QR_SOLVE is a function with a simple interface which");
        Console.WriteLine("  solves a linear system A*x = b in the least squares sense.");
        Console.WriteLine("  Compare a tabulated solution X1 to the QR_SOLVE result X2.");

        prob_num = ProbabilityFunctions.p00_prob_num();

        Console.WriteLine("");
        Console.WriteLine("  Number of problems = " + prob_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||");
        Console.WriteLine("");

        for (prob = 1; prob <= prob_num; prob++)
        {
            //
            //  Get problem size.
            //
            m = ProbabilityFunctions.p00_m(prob);
            n = ProbabilityFunctions.p00_n(prob);
            //
            //  Retrieve problem data.
            //
            a = ProbabilityFunctions.p00_a(prob, m, n);
            b = ProbabilityFunctions.p00_b(prob, m);
            x1 = ProbabilityFunctions.p00_x(prob, n);

            b_norm = typeMethods.r8vec_norm(m, b);
            x1_norm = typeMethods.r8vec_norm(n, x1);
            r1 = typeMethods.r8mat_mv_new(m, n, a, x1);
            for (i = 0; i < m; i++)
            {
                r1[i] -= b[i];
            }

            r1_norm = typeMethods.r8vec_norm(m, r1);
            //
            //  Use QR_SOLVE on the problem.
            //
            x2 = QRSolve.qr_solve(m, n, a, b);

            x2_norm = typeMethods.r8vec_norm(n, x2);
            r2 = typeMethods.r8mat_mv_new(m, n, a, x2);
            for (i = 0; i < m; i++)
            {
                r2[i] -= b[i];
            }

            r2_norm = typeMethods.r8vec_norm(m, r2);
            //
            //  Compare tabulated and computed solutions.
            //
            x_diff_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
            //
            //  Report results for this problem.
            //
            Console.WriteLine("  " + prob.ToString().PadLeft(5)
                                   + "  " + m.ToString().PadLeft(4)
                                   + "  " + n.ToString().PadLeft(4)
                                   + "  " + b_norm.ToString().PadLeft(12)
                                   + "  " + x_diff_norm.ToString().PadLeft(12)
                                   + "  " + x1_norm.ToString().PadLeft(12)
                                   + "  " + x2_norm.ToString().PadLeft(12)
                                   + "  " + r1_norm.ToString().PadLeft(12)
                                   + "  " + r2_norm.ToString().PadLeft(12) + "");
        }
    }

    private static void svd_solve_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SVD_SOLVE_TEST tests SVD_SOLVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b;
        double b_norm;
        int i;
        int m;
        int n;
        int prob;
        int prob_num;
        double[] r1;
        double r1_norm;
        double[] r2;
        double r2_norm;
        double x_diff_norm;
        double[] x1;
        double x1_norm;
        double[] x2;
        double x2_norm;

        Console.WriteLine("");
        Console.WriteLine("SVD_SOLVE_TEST");
        Console.WriteLine("  SVD_SOLVE is a function with a simple interface which");
        Console.WriteLine("  solves a linear system A*x = b in the least squares sense.");
        Console.WriteLine("  Compare a tabulated solution X1 to the SVD_SOLVE result X2.");

        prob_num = ProbabilityFunctions.p00_prob_num();

        Console.WriteLine("");
        Console.WriteLine("  Number of problems = " + prob_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||");
        Console.WriteLine("");

        for (prob = 1; prob <= prob_num; prob++)
        {
            //
            //  Get problem size.
            //
            m = ProbabilityFunctions.p00_m(prob);
            n = ProbabilityFunctions.p00_n(prob);
            //
            //  Retrieve problem data.
            //
            a = ProbabilityFunctions.p00_a(prob, m, n);
            b = ProbabilityFunctions.p00_b(prob, m);
            x1 = ProbabilityFunctions.p00_x(prob, n);

            b_norm = typeMethods.r8vec_norm(m, b);
            x1_norm = typeMethods.r8vec_norm(n, x1);
            r1 = typeMethods.r8mat_mv_new(m, n, a, x1);
            for (i = 0; i < m; i++)
            {
                r1[i] -= b[i];
            }

            r1_norm = typeMethods.r8vec_norm(m, r1);
            //
            //  Use SVD_SOLVE on the problem.
            //
            x2 = SVDSolve.svd_solve(m, n, a, b);

            x2_norm = typeMethods.r8vec_norm(n, x2);
            r2 = typeMethods.r8mat_mv_new(m, n, a, x2);
            for (i = 0; i < m; i++)
            {
                r2[i] -= b[i];
            }

            r2_norm = typeMethods.r8vec_norm(m, r2);
            //
            //  Compare tabulated and computed solutions.
            //
            x_diff_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
            //
            //  Report results for this problem.
            //
            Console.WriteLine("  " + prob.ToString().PadLeft(5)
                                   + "  " + m.ToString().PadLeft(4)
                                   + "  " + n.ToString().PadLeft(4)
                                   + "  " + b_norm.ToString().PadLeft(12)
                                   + "  " + x_diff_norm.ToString().PadLeft(12)
                                   + "  " + x1_norm.ToString().PadLeft(12)
                                   + "  " + x2_norm.ToString().PadLeft(12)
                                   + "  " + r1_norm.ToString().PadLeft(12)
                                   + "  " + r2_norm.ToString().PadLeft(12) + "");
        }
    }

    private static void dqrls_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DQRLS_TEST tests DQRLS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b =  {
                1.0, 2.3, 4.6, 3.1, 1.2
            }
            ;
        int i;
        int ind;
        int itask;
        int j;
        int[] jpvt;
        int kr = 0;
        int m = 5;
        int n = 3;
        double[] qraux;
        double tol;
        double[] x;

        a = new double[m * n];
        jpvt = new int[n];
        qraux = new double[n];
        x = new double[n];
        //
        //  Set up least-squares problem
        //  quadratic model, equally-spaced points
        //
        Console.WriteLine("");
        Console.WriteLine("DQRLS_TEST");
        Console.WriteLine("  DQRLS solves a linear system A*x = b in the least squares sense.");

        for (i = 0; i < m; i++)
        {
            a[i + 0 * m] = 1.0;
            for (j = 1; j < n; j++)
            {
                a[i + j * m] = a[i + (j - 1) * m] * (i + 1);
            }
        }

        tol = 1.0E-06;

        typeMethods.r8mat_print(m, n, a, "  Coefficient matrix A:");

        typeMethods.r8vec_print(m, b, "  Right hand side b:");
        //
        //  Solve least-squares problem
        //
        itask = 1;
        ind = QRSolve.dqrls(ref a, m, m, n, tol, ref kr, b, ref x, ref b, jpvt, qraux, itask);
        //
        //  Print results
        //
        Console.WriteLine("");
        Console.WriteLine("  Error code = " + ind + "");
        Console.WriteLine("  Estimated matrix rank = " + kr + "");

        typeMethods.r8vec_print(n, x, "  Least squares solution x:");

        typeMethods.r8vec_print(m, b, "  Residuals A*x-b");
    }
}