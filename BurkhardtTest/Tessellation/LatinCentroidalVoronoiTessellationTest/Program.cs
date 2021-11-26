﻿using System;
using Burkardt.Tessellation;
using Burkardt.Types;

namespace LatinCentroidalVoronoiTessellationTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LCVT_TEST.
        //
        //  Discussion:
        //
        //    LCVT_TEST tests LCVT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("LCVT_TEST");
        Console.WriteLine("  Test the LCVT library.");

        for (i = -1; i <= 2; i++)
        {
            int sample_function_cvt = i;
            test01(sample_function_cvt);
        }

        test02();

        Console.WriteLine("");
        Console.WriteLine("LCVT_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int sample_function_cvt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CVT, R8MAT_LATINIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 2;
        const int N = 25;
        LCVData ldata = new(M);

        double[] generator = new double[M * N];
        int i;
        const int latin_steps = 3;
        const int sample_function_init = 0;
        const int sample_num_cvt = 100000;
        const int sample_num_steps = 50;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  R8MAT_LATINIZE makes it a Latin Hypersquare.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we vary the sampling used during the");
        Console.WriteLine("  CVT Latin iteration.");
        //
        //  GET_SEED can be used to produce a different seed on each run.
        //  But using a fixed seed is useful for debugging.
        //
        int seed = entropyRNG.RNG.nextint();

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M =        " + M + "");
        Console.WriteLine("  Number of generators =       " + N + "");
        Console.WriteLine("  Initial random number seed = " + seed + "");
        ;
        Console.WriteLine("");

        switch (sample_function_init)
        {
            case -1:
                Console.WriteLine("  Initialize using RANDOM_NUMBER (C++ STDLIB intrinsic).");
                break;
            case 0:
                Console.WriteLine("  Initialize using UNIFORM.");
                break;
            case 1:
                Console.WriteLine("  Initialize using HALTON.");
                break;
            case 2:
                Console.WriteLine("  Initialize using GRID.");
                break;
            case 3:
                Console.WriteLine("  USER will initialize data.");
                break;
        }

        switch (sample_function_cvt)
        {
            case -1:
                Console.WriteLine("  Sample using RANDOM_NUMBER (C++ STDLIB intrinsic).");
                break;
            case 0:
                Console.WriteLine("  Sample using UNIFORM.");
                break;
            case 1:
                Console.WriteLine("  Sample using HALTON.");
                break;
            case 2:
                Console.WriteLine("  Sample using GRID.");
                break;
        }

        Console.WriteLine("  Number of sample points = " + sample_num_cvt + "");
        Console.WriteLine("  Number of sample steps =  " + sample_num_steps + "");

        for (i = 1; i <= latin_steps; i++)
        {
            LatinCentroidalVoronoi.cvt(ref ldata, M, N, sample_function_init, sample_function_cvt,
                sample_num_cvt, sample_num_steps, ref seed, ref generator);

            typeMethods.r8mat_transpose_print(M, N, generator, "  After CVT steps:");

            typeMethods.r8mat_latinize(M, N, ref generator);

            typeMethods.r8mat_transpose_print(M, N, generator, "  After Latin step:");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests CVT. R8MAT_LATINIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 2;
        const int N = 25;
        LCVData ldata = new(M);

        double[] generator = new double[M * N];
        int i;
        const int latin_steps = 3;
        int rank;
        const int sample_function_cvt = 0;
        const int sample_function_init = 3;
        const int sample_num_cvt = 100000;
        const int sample_num_steps = 50;
        int[] tuple = new int[M];

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  R8MAT_LATINIZE makes it a Latin Hypersquare.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we initialize the generators to");
        Console.WriteLine("  grid points; this is an unstable CVT solution.");
        //
        //  GET_SEED can be used to produce a different seed on each run.
        //  But using a fixed seed is useful for debugging.
        //
        int seed = entropyRNG.RNG.nextint();

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M =        " + M + "");
        Console.WriteLine("  Number of generators =       " + N + "");
        Console.WriteLine("  Initial random number seed = " + seed + "");
        ;
        Console.WriteLine("");

        switch (sample_function_init)
        {
            case -1:
                Console.WriteLine("  Initialize using RANDOM_NUMBER (C++ STDLIB intrinsic).");
                break;
            case 0:
                Console.WriteLine("  Initialize using UNIFORM.");
                break;
            case 1:
                Console.WriteLine("  Initialize using HALTON.");
                break;
            case 2:
                Console.WriteLine("  Initialize using GRID.");
                break;
            case 3:
                Console.WriteLine("  USER will initialize data.");
                break;
        }

        switch (sample_function_cvt)
        {
            case -1:
                Console.WriteLine("  Sample using RANDOM_NUMBER (C++ STDLIB intrinsic).");
                break;
            case 0:
                Console.WriteLine("  Sample using UNIFORM.");
                break;
            case 1:
                Console.WriteLine("  Sample using HALTON.");
                break;
            case 2:
                Console.WriteLine("  Sample using GRID.");
                break;
        }

        Console.WriteLine("  Number of sample points = " + sample_num_cvt + "");
        Console.WriteLine("  Number of sample steps =  " + sample_num_steps + "");

        int ngrid = 5;

        for (rank = 0; rank <= N - 1; rank++)
        {
            BTuple.tuple_next_fast(ref ldata.data.tdata, ngrid, M, rank, ref tuple);
            for (i = 0; i < M; i++)
            {
                generator[i + rank * M] = (2 * tuple[i] - 1)
                                          / (double) (2 * ngrid);
            }
        }

        typeMethods.r8mat_transpose_print(M, N, generator, "  Initial generators (rows):");

        for (i = 1; i <= latin_steps; i++)
        {
            LatinCentroidalVoronoi.cvt(ref ldata, M, N, sample_function_init, sample_function_cvt,
                sample_num_cvt, sample_num_steps, ref seed, ref generator);

            typeMethods.r8mat_transpose_print(M, N, generator, "  After CVT steps:");

            typeMethods.r8mat_latinize(M, N, ref generator);

            typeMethods.r8mat_transpose_print(M, N, generator, "  After Latin step:");

        }
    }
}