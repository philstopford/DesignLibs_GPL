﻿using System;
using System.Globalization;
using Burkardt.MatrixNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace WathenMatrixTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for wathen_test.
        //
        //  Discussion:
        //
        //    wathen_test tests wathen.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 February 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("wathen_test");
        Console.WriteLine("  Test wathen.");

        test01();
        test02();
        test05();
        test06();
        test07();
        test08();
        test10();
        test11();
        test115();
        wathen_xy_test();

        Console.WriteLine("");
        Console.WriteLine("wathen_test");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 assembles, factor and solve using WATHEN_GE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Assemble, factor and solve a Wathen system");
        Console.WriteLine("  defined by WATHEN_GE.");
        Console.WriteLine("");

        const int nx = 4;
        const int ny = 4;
        Console.WriteLine("  Elements in X direction NX = " + nx + "");
        Console.WriteLine("  Elements in Y direction NY = " + ny + "");
        Console.WriteLine("  Number of elements = " + nx * ny + "");
        //
        //  Compute the number of unknowns.
        //
        int n = WathenMatrix.wathen_order(nx, ny);
        Console.WriteLine("  Number of nodes N = " + n + "");
        //
        //  Set up a random solution X.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the matrix.
        //
        seed = 123456789;
        double[] a = WathenMatrix.wathen_ge(nx, ny, n, ref seed);
        //
        //  Compute the corresponding right hand side B.
        //
        double[] b = MatbyVector.mv_ge(n, n, a, x1);
        //
        //  Solve the linear system.
        //
        int[] ipvt = new int[n];
        Matrix.dgefa(ref a, n, n, ref ipvt);

        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = b[i];
        }

        int job = 0;
        Matrix.dgesl(a, n, n, ipvt, ref x2, job);
        //
        //  Compute the maximum solution error.
        //
        double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
        Console.WriteLine("  Maximum solution error is " + e + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 assembles, factors and solves using WATHEN_GB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int md = 0;
        int ml = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Assemble, factor and solve a Wathen system");
        Console.WriteLine("  using WATHEN_GB.");
        Console.WriteLine("");

        int nx = 4;
        int ny = 4;
        Console.WriteLine("  Elements in X direction NX = " + nx + "");
        Console.WriteLine("  Elements in Y direction NY = " + ny + "");
        Console.WriteLine("  Number of elements = " + nx * ny + "");
        //
        //  Compute the number of unknowns.
        //
        int n = WathenMatrix.wathen_order(nx, ny);
        Console.WriteLine("  Number of nodes N = " + n + "");
        //
        //  Compute the bandwidth.
        //
        WathenMatrix.wathen_bandwidth(nx, ny, ref ml, ref md, ref mu);
        Console.WriteLine("  Lower bandwidth ML = " + ml + "");
        Console.WriteLine("  Upper bandwidth MU = " + mu + "");
        //
        //  Set up a random solution X1.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the matrix.
        //
        seed = 123456789;
        double[] a = WathenMatrix.wathen_gb(nx, ny, n, ref seed);
        //
        //  Compute the corresponding right hand side B.
        //
        double[] b = MatbyVector.mv_gb(n, n, ml, mu, a, x1);
        //
        //  Solve the linear system.
        //
        int lda = 2 * ml + mu + 1;
        int[] ipvt = new int[n];
        Matrix.dgbfa(ref a, lda, n, ml, mu, ref ipvt);

        double[] x2 = new double[n];
        for (int i = 0; i < n; i++)
        {
            x2[i] = b[i];
        }

        int job = 0;
        Matrix.dgbsl(a, lda, n, ml, mu, ipvt, ref x2, job);
        //
        //  Compute the maximum solution error.
        //
        double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
        Console.WriteLine("  Maximum solution error is " + e + "");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 measures the storage needed for the Wathen system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int bd1 = 0;
        int bd2 = 0;
        int bl1 = 0;
        int bl2 = 0;
        int bu1 = 0;
        int bu2 = 0;
        int bw2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  For various problem sizes and storage schemes,");
        Console.WriteLine("  measure the storage used for the Wathen system.");
        Console.WriteLine("");
        Console.WriteLine("                                   Predicted  Observed");
        Console.WriteLine("                              GE        Band      Band      "
                          + "Band    Sparse");
        Console.WriteLine("    NX  Elements   Nodes   storage     width     width   "
                          + "  storage   storage");
        Console.WriteLine("");

        int nx = 1;
        int ny = 1;

        for (int test = 1; test <= 6; test++)
        {
            //
            //  Compute the number of unknowns.
            //
            int n = WathenMatrix.wathen_order(nx, ny);
            //
            //  Predict the bandwidth.
            //
            WathenMatrix.wathen_bandwidth(nx, ny, ref bl1, ref bd1, ref bu1);
            int bw1 = bl1 + bd1 + bu1;
            //
            //  Compute the matrix.
            //
            int seed = 123456789;
            double[] a = WathenMatrix.wathen_ge(nx, ny, n, ref seed);

            int storage_ge = n * n;

            Matrix.bandwidth(n, n, a, ref bw2, ref bl2, ref bd2, ref bu2);
            int storage_gb = (2 * bl2 + 1 + bu2) * n;

            int storage_sparse = Matrix.nonzeros(n, n, a);
            //
            //  Report.
            //
            Console.WriteLine(nx.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "      "
                                                       + (nx * ny).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                       + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                       + storage_ge.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                       + bw1.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                       + bw2.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                       + storage_gb.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                       + storage_sparse.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
            //
            //  Ready for next iteration.
            //
            nx *= 2;
            ny *= 2;
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 times WATHEN_GE assembly and solution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  For various problem sizes,");
        Console.WriteLine("  time the assembly and factorization of a Wathen system");
        Console.WriteLine("  using the WATHEN_GE function.");
        Console.WriteLine("");
        Console.WriteLine("    NX  Elements   Nodes   Storage    Assembly      Factor      Error");
        Console.WriteLine("");

        int nx = 1;
        int ny = 1;

        for (test = 1; test <= 6; test++)
        {
            //
            //  Compute the number of unknowns.
            //
            int n = WathenMatrix.wathen_order(nx, ny);
            int storage_ge = n * n;
            //
            //  Set up a random solution X1.
            //
            int seed = 123456789;
            double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            //
            //  Compute the matrix, and measure the storage required.
            //
            seed = 123456789;

            DateTime t0 = DateTime.Now;
            double[] a = WathenMatrix.wathen_ge(nx, ny, n, ref seed);
            TimeSpan t1 = DateTime.Now - t0;
            //
            //  Compute the corresponding right hand side B.
            //
            double[] b = MatbyVector.mv_ge(n, n, a, x1);
            //
            //  Solve the system.
            //
            int[] ipvt = new int[n];
            double[] x2 = new double[n];
            int i;
            for (i = 0; i < n; i++)
            {
                x2[i] = b[i];
            }

            int job = 0;

            t0 = DateTime.Now;
            Matrix.dgefa(ref a, n, n, ref ipvt);
            Matrix.dgesl(a, n, n, ipvt, ref x2, job);
            TimeSpan t2 = DateTime.Now - t0;
            //
            //  Compute the maximum solution error.
            //
            double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
            //
            //  Report.
            //
            Console.WriteLine(nx.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "      "
                                                       + (nx * ny).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                       + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                       + storage_ge.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                       + t1.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                       + t2.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                       + e.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            //
            //  Ready for next iteration.
            //
            nx *= 2;
            ny *= 2;
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 times WATHEN_GB assembly and solution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int md = 0;
        int ml = 0;
        int mu = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  For various problem sizes,");
        Console.WriteLine("  time the assembly and factorization of a Wathen system");
        Console.WriteLine("  using the WATHEN_GB function.");
        Console.WriteLine("");
        Console.WriteLine("    NX  Elements   Nodes   Storage    Assembly      Factor      Error");
        Console.WriteLine("");

        int nx = 1;
        int ny = 1;

        for (test = 1; test <= 6; test++)
        {
            //
            //  Compute the number of unknowns.
            //
            int n = WathenMatrix.wathen_order(nx, ny);
            //
            //  Compute the bandwidth.
            //
            WathenMatrix.wathen_bandwidth(nx, ny, ref ml, ref md, ref mu);
            int 
            
            
            
                storage_gb = (2 * ml + mu + 1) * n;
            //
            //  Set up a random solution X1.
            //
            int seed = 123456789;
            double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            //
            //  Compute the matrix.
            //
            seed = 123456789;
            DateTime t0 = DateTime.Now;
            double[] a = WathenMatrix.wathen_gb(nx, ny, n, ref seed);
            TimeSpan t1 = DateTime.Now - t0;
            //
            //  Compute the corresponding right hand side B.
            //
            double[] b = MatbyVector.mv_gb(n, n, ml, mu, a, x1);
            //
            //  Solve the system.
            //
            int lda = 2 * ml + mu + 1;
            int[] ipvt = new int[n];
            double[] x2 = new double[n];
            int i;
            for (i = 0; i < n; i++)
            {
                x2[i] = b[i];
            }

            int job = 0;

            t0 = DateTime.Now;
            Matrix.dgbfa(ref a, lda, n, ml, mu, ref ipvt);
            Matrix.dgbsl(a, lda, n, ml, mu, ipvt, ref x2, job);
            TimeSpan t2 = DateTime.Now - t0;
            //
            //  Compute the maximum solution error.
            //
            double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
            //
            //  Report.
            //
            Console.WriteLine(nx.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "      "
                                                       + (nx * ny).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                       + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                       + storage_gb.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                       + t1.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                       + t2.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                       + e.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            //
            //  Ready for next iteration.
            //
            nx *= 2;
            ny *= 2;
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 times WATHEN_GE/WATHEN_GB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int md = 0;
        int ml = 0;
        int mu = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  For various problem sizes,");
        Console.WriteLine("  time the assembly and factorization of a Wathen system");
        Console.WriteLine("  WATHEN_GE/WATHEN_GB");
        Console.WriteLine("");
        Console.WriteLine("                   NX  Elements   Nodes   Storage    "
                          + "  Assembly      Factor      Error");

        int nx = 1;
        int ny = 1;

        for (test = 1; test <= 6; test++)
        {
            //
            //  Compute the number of unknowns.
            //
            int n = WathenMatrix.wathen_order(nx, ny);
            int storage_ge = n * n;
            //
            //  Set up a random solution X1.
            //
            int seed = 123456789;
            double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            //
            //  Compute the matrix.
            //
            seed = 123456789;
            DateTime t0 = DateTime.Now;
            double[] a = WathenMatrix.wathen_ge(nx, ny, n, ref seed);
            TimeSpan t1 = DateTime.Now - t0;
            //
            //  Compute the corresponding right hand side B.
            //
            double[] b = MatbyVector.mv_ge(n, n, a, x1);
            //
            //  Solve the system.
            //
            int[] ipvt = new int[n];
            double[] x2 = new double[n];
            int i;
            for (i = 0; i < n; i++)
            {
                x2[i] = b[i];
            }

            int job = 0;

            t0 = DateTime.Now;
            Matrix.dgefa(ref a, n, n, ref ipvt);
            Matrix.dgesl(a, n, n, ipvt, ref x2, job);
            TimeSpan t2 = DateTime.Now - t0;
            //
            //  Compute the maximum solution error.
            //
            double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
            //
            //  Report.
            //
            Console.WriteLine("");
            Console.WriteLine("  WATHEN_GE      "
                              + nx.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "      "
                              + (nx * ny).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + storage_ge.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + t1.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + t2.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + e.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

            //
            //  Compute the bandwidth.
            //
            WathenMatrix.wathen_bandwidth(nx, ny, ref ml, ref md, ref mu);
            int 
                storage_gb = (2 * ml + mu + 1) * n;
            //
            //  Set up a random solution X1.
            //
            seed = 123456789;
            x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            //
            //  Compute the matrix.
            //
            seed = 123456789;
            t0 = DateTime.Now;
            a = WathenMatrix.wathen_gb(nx, ny, n, ref seed);
            t1 = DateTime.Now - t0;
            //
            //  Compute the corresponding right hand side B.
            //
            b = MatbyVector.mv_gb(n, n, ml, mu, a, x1);
            //
            //  Solve the system.
            //
            int lda = 2 * ml + mu + 1;
            ipvt = new int[n];
            x2 = new double[n];
            for (i = 0; i < n; i++)
            {
                x2[i] = b[i];
            }

            job = 0;

            t0 = DateTime.Now;
            Matrix.dgbfa(ref a, lda, n, ml, mu, ref ipvt);
            Matrix.dgbsl(a, lda, n, ml, mu, ipvt, ref x2, job);
            t2 = DateTime.Now - t0;
            //
            //  Compute the maximum solution error.
            //
            e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
            //
            //  Report.
            //
            Console.WriteLine("  WATHEN_GB      "
                              + nx.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "      "
                              + (nx * ny).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + storage_gb.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + t1.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + t2.TotalSeconds.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + e.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            //
            //  Ready for next iteration.
            //
            nx *= 2;
            ny *= 2;
        }
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 assembles, factor and solve using WATHEN_GE and CG_GE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  Assemble, factor and solve a Wathen system");
        Console.WriteLine("  defined by WATHEN_GE and CG_GE.");
        Console.WriteLine("");

        const int nx = 1;
        const int ny = 1;
        Console.WriteLine("  Elements in X direction NX = " + nx + "");
        Console.WriteLine("  Elements in Y direction NY = " + ny + "");
        Console.WriteLine("  Number of elements = " + nx * ny + "");
        //
        //  Compute the number of unknowns.
        //
        int n = WathenMatrix.wathen_order(nx, ny);
        Console.WriteLine("  Number of nodes N = " + n + "");
        //
        //  Set up a random solution X.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the matrix.
        //
        seed = 123456789;
        double[] a = WathenMatrix.wathen_ge(nx, ny, n, ref seed);
        //
        //  Compute the corresponding right hand side B.
        //

        double[] b = MatbyVector.mv_ge(n, n, a, x1);
        //
        //  Solve the linear system.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        ConjugateGradient.cg_ge(n, a, b, ref x2);
        //
        //  Compute the maximum solution error.
        //
        double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
        Console.WriteLine("  Maximum solution error is " + e + "");
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 assemble, factor and solve using WATHEN_ST + CG_ST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  Assemble, factor and solve a Wathen system");
        Console.WriteLine("  defined by WATHEN_ST and CG_ST.");
        Console.WriteLine("");

        const int nx = 1;
        const int ny = 1;
        Console.WriteLine("  Elements in X direction NX = " + nx + "");
        Console.WriteLine("  Elements in Y direction NY = " + ny + "");
        Console.WriteLine("  Number of elements = " + nx * ny + "");
        //
        //  Compute the number of unknowns.
        //
        int n = WathenMatrix.wathen_order(nx, ny);
        Console.WriteLine("  Number of nodes N = " + n + "");
        //
        //  Set up a random solution X1.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the matrix size.
        //
        int nz_num = WathenMatrix.wathen_st_size(nx, ny);
        Console.WriteLine("  Number of nonzeros NZ_NUM = " + nz_num + "");
        //
        //  Compute the matrix.
        //
        seed = 123456789;
        int[] row = new int[nz_num];
        int[] col = new int[nz_num];
        double[] a = WathenMatrix.wathen_st(nx, ny, nz_num, ref seed, ref row, ref col);
        //
        //  Compute the corresponding right hand side B.
        //
        double[] b = MatbyVector.mv_st(n, n, nz_num, row, col, a, x1);
        //
        //  Solve the linear system.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        ConjugateGradient.cg_st(n, nz_num, row, col, a, b, ref x2);
        //
        //  Compute the maximum solution error.
        //
        double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
        Console.WriteLine("  Maximum solution error is " + e + "");
    }

    private static void test115()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST115 assembles, factors and solves using WATHEN_GB and CG_GB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int md = 0;
        int ml = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST115");
        Console.WriteLine("  Assemble, factor and solve a Wathen system");
        Console.WriteLine("  using WATHEN_GB and CG_GB.");
        Console.WriteLine("");

        const int nx = 4;
        const int ny = 4;
        Console.WriteLine("  Elements in X direction NX = " + nx + "");
        Console.WriteLine("  Elements in Y direction NY = " + ny + "");
        Console.WriteLine("  Number of elements = " + nx * ny + "");
        //
        //  Compute the number of unknowns.
        //
        int n = WathenMatrix.wathen_order(nx, ny);
        Console.WriteLine("  Number of nodes N = " + n + "");
        //
        //  Compute the bandwidth.
        //
        WathenMatrix.wathen_bandwidth(nx, ny, ref ml, ref md, ref mu);
        Console.WriteLine("  Lower bandwidth ML = " + ml + "");
        Console.WriteLine("  Upper bandwidth MU = " + mu + "");
        //
        //  Set up a random solution X1.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the matrix.
        //
        seed = 123456789;
        double[] a = WathenMatrix.wathen_gb(nx, ny, n, ref seed);
        //
        //  Compute the corresponding right hand side B.
        //
        double[] b = MatbyVector.mv_gb(n, n, ml, mu, a, x1);
        //
        //  Solve the linear system.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        ConjugateGradient.cg_gb(n, ml, mu, a, b, ref x2);
        //
        //  Compute the maximum solution error.
        //
        double e = typeMethods.r8vec_diff_norm_li(n, x1, x2);
        Console.WriteLine("  Maximum solution error is " + e + "");
    }

    private static void wathen_xy_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    wathen_xy_test tests wathen_xy.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 February 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine("wathen_xy_test");
        Console.WriteLine("  wathen_xy returns the (X,Y) coordinates of nodes in the");
        Console.WriteLine("  Wathen finite element system.");

        const int nx = 3;
        const int ny = 3;
        int n = WathenMatrix.wathen_order(nx, ny);
        double[] xy = WathenMatrix.wathen_xy(nx, ny, n);

        Console.WriteLine("");
        Console.WriteLine("   k   i   j         x          y");
        Console.WriteLine("");

        int k = 0;
        for (j = 0; j <= 2 * ny; j++)
        {
            int i;
            switch (j % 2)
            {
                case 0:
                {
                    for (i = 0; i <= 2 * ny; i++)
                    {
                        Console.WriteLine("  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                               + "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                               + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                               + "  " + xy[k].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + "  " + xy[k + n].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                        k += 1;
                    }

                    break;
                }
                default:
                {
                    for (i = 0; i <= ny; i++)
                    {
                        Console.WriteLine("  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                               + "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                               + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                               + "  " + xy[k].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + "  " + xy[k + n].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                        k += 1;
                    }

                    break;
                }
            }
        }
    }
}