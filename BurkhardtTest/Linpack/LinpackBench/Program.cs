﻿using System;
using System.Globalization;
using Burkardt.Linpack;
using Burkardt.Types;

namespace LinpackBench;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LINPACK_BENCH.
        //
        //  Discussion:
        //
        //    LINPACK_BENCH drives the double precision LINPACK benchmark program.
        //
        //  Modified:
        //
        //    30 November 2006
        //
        //  Parameters:
        //
        //    N is the problem size.
        //
    {
        const int N = 1000;
        const int LDA = N + 1;

        const double cray = 0.056;
        int i;
        int j;
        double[] time = new double[6];

        Console.WriteLine("");
        Console.WriteLine("LINPACK_BENCH");
        Console.WriteLine("");
        Console.WriteLine("  The LINPACK benchmark.");
        Console.WriteLine("  Datatype: Double precision");
        Console.WriteLine("  Matrix order N               = " + N + "");
        Console.WriteLine("  Leading matrix dimension LDA = " + LDA + "");

        double ops = 2 * N * N * N / 3.0 + 2.0 * (N * N);
        //
        //  Allocate space for arrays.
        //
        double[] a = typeMethods.r8mat_gen(LDA, N);
        double[] b = new double[N];
        int[] ipvt = new int[N];
        double[] resid = new double[N];
        double[] rhs = new double[N];
        double[] x = new double[N];

        double a_max = 0.0;
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                a_max = Math.Max(a_max, a[i + j * LDA]);
            }
        }

        for (i = 0; i < N; i++)
        {
            x[i] = 1.0;
        }

        for (i = 0; i < N; i++)
        {
            b[i] = 0.0;
            for (j = 0; j < N; j++)
            {
                b[i] += a[i + j * LDA] * x[j];
            }
        }

        DateTime t1 = DateTime.Now;

        int info = DGEFA.dgefa(ref a, LDA, N, ref ipvt);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("LINPACK_BENCH - Fatal error!");
            Console.WriteLine("  The matrix A is apparently singular.");
            Console.WriteLine("  Abnormal end of execution.");
            return;
        }

        DateTime t2 = DateTime.Now;
        time[0] = (t2 - t1).TotalSeconds;

        t1 = DateTime.Now;

        int job = 0;
        DGESL.dgesl(a, LDA, N, ipvt, ref b, job);

        t2 = DateTime.Now;
        time[1] = (t2 - t1).TotalSeconds;

        double total = time[0] + time[1];

        //
        //  Compute a residual to verify results.
        //
        a = typeMethods.r8mat_gen(LDA, N);

        for (i = 0; i < N; i++)
        {
            x[i] = 1.0;
        }

        for (i = 0; i < N; i++)
        {
            rhs[i] = 0.0;
            for (j = 0; j < N; j++)
            {
                rhs[i] += a[i + j * LDA] * x[j];
            }
        }

        for (i = 0; i < N; i++)
        {
            resid[i] = -rhs[i];
            for (j = 0; j < N; j++)
            {
                resid[i] += a[i + j * LDA] * b[j];
            }
        }

        double resid_max = 0.0;
        for (i = 0; i < N; i++)
        {
            resid_max = Math.Max(resid_max, Math.Abs(resid[i]));
        }

        double b_max = 0.0;
        for (i = 0; i < N; i++)
        {
            b_max = Math.Max(b_max, Math.Abs(b[i]));
        }

        double eps = typeMethods.r8_epsilon();

        double residn = resid_max / (N * a_max * b_max * eps);

        time[2] = total;
        time[3] = total switch
        {
            > 0.0 => ops / (1.0E+06 * total),
            _ => -1.0
        };

        time[4] = 2.0 / time[3];
        time[5] = total / cray;

        Console.WriteLine("");
        Console.WriteLine("     Norm. Resid      Resid           MACHEP         X[1]          X[N]");
        Console.WriteLine("");
        Console.WriteLine(residn.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                        + resid_max.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                        + eps.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                        + b[0].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                        + b[N - 1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("");
        Console.WriteLine("      Factor     Solve      Total     MFLOPS       Unit      Cray-Ratio");
        Console.WriteLine("");
        Console.WriteLine(time[0].ToString(CultureInfo.InvariantCulture).PadLeft(9) + "  "
                                                        + time[1].ToString(CultureInfo.InvariantCulture).PadLeft(9) + "  "
                                                        + time[2].ToString(CultureInfo.InvariantCulture).PadLeft(9) + "  "
                                                        + time[3].ToString(CultureInfo.InvariantCulture).PadLeft(9) + "  "
                                                        + time[4].ToString(CultureInfo.InvariantCulture).PadLeft(9) + "  "
                                                        + time[5].ToString(CultureInfo.InvariantCulture).PadLeft(9) + "");

        Console.WriteLine("");
        Console.WriteLine("LINPACK_BENCH");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}