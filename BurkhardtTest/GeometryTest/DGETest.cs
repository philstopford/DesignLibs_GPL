using System;
using Burkardt.MatrixNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest;

public static class DGETest
{
    public static void test1787()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1787 tests DGE_FA and DGE_SL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        double[] a;
        double[] alu = new double[N * N];
        double[] b = new double[N];
        int i;
        int info;
        int j;
        int job;
        int[] pivot = new int[N];
        int seed;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST1787");
        Console.WriteLine("  DGE_FA factors a general linear system,");
        Console.WriteLine("  DGE_SL solves a factored system.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + N + "");
        //
        //  Set the matrix.
        //
        seed = 123456789;
        a = UniformRNG.r8mat_uniform_01_new(N, N, ref seed);
        typeMethods.r8mat_print(N, N, a, "  Matrix A:");
        //
        //  Set the desired solution.
        //
        for (i = 0; i < N; i++)
        {
            x[i] = i + 1;
        }

        typeMethods.r8vec_print(N, x, "  Desired solution vector:");
        //
        //  Compute the corresponding right hand side.
        //
        for (i = 0; i < N; i++)
        {
            b[i] = 0.0;
            for (j = 0; j < N; j++)
            {
                b[i] += a[i + j * N] * x[j];
            }
        }

        typeMethods.r8vec_print(N, b, "  Right hand side vector:");
        //
        //  Make a copy of the matrix.
        //
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                alu[i + j * N] = a[i + j * N];
            }
        }

        //
        //  Factor the matrix.
        //
        info = Matrix.dge_fa(N, ref alu, ref pivot);

        typeMethods.r8mat_print(N, N, alu, "  Factored matrix ALU:");

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  Fatal error!");
            Console.WriteLine("  DGE_FA declares the matrix is singular!");
            Console.WriteLine("  The value of INFO is " + info + "");
            return;
        }

        //
        //  Solve the linear system.
        //
        job = 0;
        Matrix.dge_sl(N, alu, pivot, ref b, job);

        typeMethods.r8vec_print(N, b, "  Solution: (Should be 1, 2, 3,...)");
        //
        //  Set another the desired solution.
        //
        for (i = 0; i < N; i++)
        {
            x[i] = 1.0;
        }

        //
        //  Compute the corresponding right hand side.
        //
        for (i = 0; i < N; i++)
        {
            b[i] = 0.0;
            for (j = 0; j < N; j++)
            {
                b[i] += a[i + j * N] * x[j];
            }
        }

        //
        //  Solve the system
        //
        job = 0;
        Matrix.dge_sl(N, alu, pivot, ref b, job);

        typeMethods.r8vec_print(N, b, "  Solution: (Should be 1, 1, 1,...)");
        //
        //  Set the desired solution to a problem involving the transposed matrix.
        //
        for (i = 0; i < N; i++)
        {
            x[i] = i + 1;
        }

        //
        //  Compute the corresponding right hand side using the transposed matrix.
        //
        for (i = 0; i < N; i++)
        {
            b[i] = 0.0;
            for (j = 0; j < N; j++)
            {
                b[i] += a[j + i * N] * x[j];
            }
        }

        //
        //  Solve the transposed system.
        //
        job = 1;
        Matrix.dge_sl(N, alu, pivot, ref b, job);

        typeMethods.r8vec_print(N, b,
            "  Solution of transposed system: (Should be 1, 2, 3,...)");

    }


}