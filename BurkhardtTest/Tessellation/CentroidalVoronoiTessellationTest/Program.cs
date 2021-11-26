using System;
using Burkardt.Tessellation;
using Burkardt.Types;

namespace CentroidalVoronoiTessellationTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CVT_TEST.
        //
        //  Discussion:
        //
        //    CVT_TEST tests the CVT library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CVT_TEST");
            
        Console.WriteLine("  Test the CVT library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test12();
        test13();
        test14();

        Console.WriteLine("");
        Console.WriteLine("CVT_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CVT with uniform initialization and uniform sampling.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 10000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 repeats test 1, but uses twice as many iterations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1, but with twice the number of iterations.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 80;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 10000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 repeats test 1 but uses 100 times as many sample points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1, but with 100 times the sample points.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 1000000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 repeats test 1 with uniform initialization and Halton sampling.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1, but with Halton sampling.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 1;
        const int sample_num = 10000;
        const string sample_string = "halton";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 repeats test 1 with uniform initialization and grid sampling.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1, but with grid sampling.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 2;
        const int sample_num = 10000;
        const string sample_string = "grid";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 repeats test 1 with uniform initialization and C++ RANDOM sampling.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1, but with C++ RANDOM sampling.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = -1;
        const int sample_num = 10000;
        const string sample_string = "random";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests CVT with a different seed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1 with a different seed.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 10000;
        const string sample_string = "uniform";
        int seed = 987654321;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 repeats test 1 with a different batch size.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1 with a different batch size.");

        const int batch = 5;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 10000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 repeats test 1, but with IT_FIXED = IT_MAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Repeat test 1, but with IT_FIXED = IT_MAX.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int sample = 0;
        const int sample_num = 10000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_max,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_max + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 generates 100 points in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 100;
        const int DIM_NUM = 3;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  Compute 100 points in 3D.");

        const int batch = 1000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 10000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print_some(DIM_NUM, N, r, 1, 1, DIM_NUM, 10,
            "  First 10 Generators (rows):");

    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests CVT.
        //
        //  Discussion:
        //
        //    In this test, we initialize the generators to grid points; this is 
        //    an unstable CVT solution.  The data would "prefer" to be in a
        //    different form.  However, even if we take 2000 steps of CVT iteration,
        //    the data is still only slowly progressing towards that other 
        //    configuration.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 16;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        int j;
        double[] r = new double[DIM_NUM * N];
        int[] tuple = new int[DIM_NUM];

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we initialize the generators to");
        Console.WriteLine("  grid points; this is an unstable CVT solution.");

        const int batch = 1000;
        const int init = 4;
        const string init_string = "user initialization";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 1000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;
        //
        //  Initialize the tuple generator.
        //
        int rank = -1;
        const int ngrid = 4;

        CVTHaltonData data = new(DIM_NUM);

        BTuple.tuple_next_fast(ref data.tupledata, ngrid, DIM_NUM, rank, ref tuple);
        //
        //  Pick points on a grid.
        //
        for (j = 0; j < N; j++)
        {
            rank = j;
            BTuple.tuple_next_fast(ref data.tupledata, ngrid, DIM_NUM, rank, ref tuple);
            int i;
            for (i = 0; i < DIM_NUM; i++)
            {
                r[i + j * DIM_NUM] = (2 * tuple[i] - 1)
                                     / (double) (2 * ngrid);
            }
        }

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Initial generators (rows):");

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");

    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests CVT with 'RANDOM' initialization.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 2;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  The \"random\" initialization option calls the");
        Console.WriteLine("  system random number generator.  There is some");
        Console.WriteLine("  question about whether this works correctly.");
        Console.WriteLine("");
        Console.WriteLine("  The test is as follows:");
        Console.WriteLine("");
        Console.WriteLine("  CVT call #1:");
        Console.WriteLine("");
        Console.WriteLine("    DIM_NUM      =      2");
        Console.WriteLine("    N         =     10");
        Console.WriteLine("    INIT      =     -1");
        Console.WriteLine("    IT_MAX    =      0");
        Console.WriteLine("    SEED      = 100000");
        Console.WriteLine("");
        Console.WriteLine("    Print output values of SEED and R #1.");
        Console.WriteLine("");
        Console.WriteLine("  CVT call #2: (jump SEED)");
        Console.WriteLine("");
        Console.WriteLine("    DIM_NUM      =      2");
        Console.WriteLine("    N         =     10");
        Console.WriteLine("    INIT      =     -1");
        Console.WriteLine("    IT_MAX    =      0");
        Console.WriteLine("    SEED      = 200000.");
        Console.WriteLine("");
        Console.WriteLine("    Print output values of SEED and R #2.");
        Console.WriteLine("");
        Console.WriteLine("  CVT call #3: (restore SEED)");
        Console.WriteLine("");
        Console.WriteLine("    DIM_NUM      =      2");
        Console.WriteLine("    N         =     10");
        Console.WriteLine("    INIT      =     -1");
        Console.WriteLine("    IT_MAX    =      0");
        Console.WriteLine("    SEED_INIT = 100000");
        Console.WriteLine("");
        Console.WriteLine("    Print output values of SEED and R #3.");
        Console.WriteLine("");
        Console.WriteLine("  We expect that:");
        Console.WriteLine("  * the values of R #1 and R #2 differ;");
        Console.WriteLine("  AND");
        Console.WriteLine("  * the values of R #1 and R #3 agree.");
        //
        //  Run #1.
        //
        int batch = 1000;
        int init = -1;
        string init_string = "random";
        int it_max = 0;
        int it_fixed = 1;
        int sample = 0;
        int sample_num = 10000;
        string sample_string = "uniform";
        int seed = 100000;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
        //
        //  Run #2.
        //
        batch = 1000;
        init = -1;
        init_string = "random";
        it_max = 0;
        it_fixed = 1;
        sample = 0;
        sample_num = 10000;
        sample_string = "uniform";
        seed = 200000;

        seed_init = seed;

        // data = new CVTHaltonData(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
        //
        //  Run #3.
        //
        batch = 1000;
        init = -1;
        init_string = "random";
        it_max = 0;
        it_fixed = 1;
        sample = 0;
        sample_num = 10000;
        sample_string = "uniform";
        seed = 100000;

        seed_init = seed;

        // data = new CVTHaltonData(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests CVT with the "user" routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 100;
        const int DIM_NUM = 2;

        double energy = 0;
        const string file_out_name = "cvt_circle.txt";
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[DIM_NUM * N];

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  CVT computes a Centroidal Voronoi Tessellation.");
        Console.WriteLine("  In this example, we call the \"USER\" routine,");
        Console.WriteLine("  which allows the user to define the geometry and");
        Console.WriteLine("  density implicitly, by returning sample points.");

        const int batch = 1000;
        const int init = 3;
        string init_string = "user";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 3;
        const int sample_num = 10000;
        const string sample_string = "user";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_write(file_out_name, DIM_NUM, N, r);
    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 generates a CVT in the interval [0,1] using 10 points..
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        const int DIM_NUM = 1;

        double energy = 0;
        double it_diff = 0;
        int it_num = 0;
        double[] r = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  Generate a CVT in the interval [0,1] using 10 points.");

        const int batch = 10000;
        const int init = 0;
        const string init_string = "uniform";
        const int it_max = 40;
        const int it_fixed = 1;
        const int sample = 0;
        const int sample_num = 10000;
        const string sample_string = "uniform";
        int seed = 123456789;

        int seed_init = seed;

        CVTHaltonData data = new(DIM_NUM);

        CentroidalVoronoi.cvt(ref data, DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
            ref seed, ref r, ref it_num, ref it_diff, ref energy);

        Console.WriteLine("");
        Console.WriteLine("  Dimension DIM_NUM =        " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =       " + N + "");
        Console.WriteLine("  Initial SEED =             " + seed_init + "");
        Console.WriteLine("  Current SEED =             " + seed + "");
        Console.WriteLine("  INIT =                    \"" + init_string + "\".");
        Console.WriteLine("  Max iterations IT_MAX =    " + it_max + "");
        Console.WriteLine("  IT_FIXED (fixed samples) = " + it_fixed + "");
        Console.WriteLine("  Iterations IT_NUM =        " + it_num + "");
        Console.WriteLine("  Difference IT_DIFF =       " + it_diff + "");
        Console.WriteLine("  CVT ENERGY =               " + energy + "");
        Console.WriteLine("  SAMPLE =                  \"" + sample_string + "\".");
        Console.WriteLine("  Samples SAMPLE_NUM    =    " + sample_num + "");
        Console.WriteLine("  Sampling BATCH size =      " + batch + "");
        Console.WriteLine("  EPSILON (unit roundoff) =  " + typeMethods.r8_epsilon() + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, r, "  Generators (rows):");
    }
}