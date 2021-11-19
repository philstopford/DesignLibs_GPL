using System;
using Burkardt.Grid;
using Burkardt.RandomNS;
using Burkardt.Sequence;
using Burkardt.Types;
using Burkardt.Uniform;
using Normal = Burkardt.RandomNS.Normal;
using Polygon = Burkardt.Uniform.Polygon;
using Triangle = Burkardt.Uniform.Triangle;
using Sphere = Burkardt.Uniform.Sphere;
using Tetrahedron = Burkardt.Uniform.Tetrahedron;

namespace RandomDataTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for RANDOM_DATA_TEST.
        //
        //  Discussion:
        //
        //    RANDOM_DATA_TEST tests the RANDOM_DATA library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("RANDOM_DATA_TEST");
        Console.WriteLine("  Test the RANDOM_DATA library.");

        test005();
        test01();
        r8_normal_01_test();
        r8_uniform_01_test();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test115();
        test12();
        test125();
        test13();
        test14();
        test15();
        test16();
        test17();
        test18();
        test19();

        test20();
        test205();
        uniform_in_triangle_map1_test();
        test22();
        test23();
        test235();
        test24();
        test245();
        test25();
        test26();
        test264();
        test265();
        test267();
        test27();

        Console.WriteLine("");
        Console.WriteLine("RANDOM_DATA_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test005()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST005 tests BAD_IN_SIMPLEX01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        int n;
        string output_filename = "dummy.txt";
        int seed;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST005:");
        Console.WriteLine("  BAD_IN_SIMPLEX01 is a \"bad\" sampling technique");
        Console.WriteLine("  for the unit simplex.");


        for (dim_num = 2; dim_num <= 3; dim_num++)
        {
            seed = 123456789;
            n = 10000;
            output_filename = dim_num switch
            {
                2 => "bad_in_triangle.txt",
                3 => "bad_in_tetrahedron.txt",
                _ => output_filename
            };

            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension DIM_NUM =  " + dim_num + "");
            Console.WriteLine("  Number of points N =         " + n + "");
            Console.WriteLine("  Initial random number SEED = " + seed + "");

            x = BRandom.bad_in_simplex01(dim_num, n, ref seed);

            typeMethods.r8mat_write(output_filename, dim_num, n, x);

            Console.WriteLine("");
            Console.WriteLine("  Data written to file \"" + output_filename + "\".");
        }
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests BROWNIAN
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 100;

        string output_filename = "brownian.txt";
        int seed = 123456789;
        double[] x;
        typeMethods.r8NormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  BROWNIAN generates Brownian motion points.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = BRandom.brownian(DIM_NUM, N, ref data, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        Scale.scale_to_block01(DIM_NUM, N, ref x);

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");

    }

    private static void r8_normal_01_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_NORMAL_01_TEST tests R8_NORMAL_01
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int seed = 123456789;
        int seed_in;
        double x;
        typeMethods.r8NormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("R8_NORMAL_01_TEST:");
        Console.WriteLine("  R8_NORMAL_01 generates a single normal ");
        Console.WriteLine("  pseudorandom value.");
        Console.WriteLine("");
        Console.WriteLine("     Seed          Seed       D_NORMAL_01");
        Console.WriteLine("    (Input)       (Output)");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed_in = seed;
            x = typeMethods.r8_normal_01(ref data, ref seed);

            Console.WriteLine("  "
                              + seed_in.ToString().PadLeft(12) + "  "
                              + seed.ToString().PadLeft(12) + "  "
                              + x.ToString().PadLeft(12) + "");
        }
    }

    private static void r8_uniform_01_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_UNIFORM_01_TEST tests R8_UNIFORM_01
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int seed = 123456789;
        int seed_in;
        double x;

        Console.WriteLine("");
        Console.WriteLine("R8_UNIFORM_01_TEST:");
        Console.WriteLine("  R8_UNIFORM_01 generates a single uniform ");
        Console.WriteLine("  pseudorandom value.");
        Console.WriteLine("");
        Console.WriteLine("     Seed          Seed       D_UNIFORM_01");
        Console.WriteLine("    (Input)       (Output)");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            seed_in = seed;
            x = UniformRNG.r8_uniform_01(ref seed);

            Console.WriteLine("  "
                              + seed_in.ToString().PadLeft(12) + "  "
                              + seed.ToString().PadLeft(12) + "  "
                              + x.ToString().PadLeft(12) + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests GRID_IN_CUBE01
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 85;

        int center = 1;
        string output_filename = "grid_in_cube01.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  GRID_IN_CUBE01 generates grid points in the unit hypercube.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  CENTER option =              " + center + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Grid.grid_in_cube01(DIM_NUM, N, center, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");

    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests HALTON_IN_CIRCLE01_ACCEPT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 400;

        string output_filename = "halton_in_circle01_accept.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  HALTON_IN_CIRCLE01_ACCEPT accepts ");
        Console.WriteLine("  Halton points in the unit circle.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Halton.halton_in_circle01_accept(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests HALTON_IN_CIRCLE01_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 400;

        string output_filename = "halton_in_circle01_map.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  HALTON_IN_CIRCLE01_MAP maps ");
        Console.WriteLine("  Halton points into the unit circle.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Halton.halton_in_circle01_map(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests HALTON_IN_CUBE01
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 100;

        string output_filename = "halton_in_cube01.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  HALTON_IN_CUBE01 generates Halton points");
        Console.WriteLine("  in the unit hypercube.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Halton.halton_in_cube01(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests HAMMERSLEY_IN_CUBE01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 100;

        string output_filename = "hammersley_in_cube01.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  HAMMERSLEY_IN_CUBE01 generates Hammersley points");
        Console.WriteLine("  in the unit hypercube.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Hammersley.hammersley_in_cube01(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests NORMAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "normal.txt";
        int i;
        int info;
        int j;
        double[] mu =  {
                6.0, 100.0
            }
            ;
        double[] r = new double[DIM_NUM * DIM_NUM];
        int seed = 123456789;
        double[] v =  {
                5.0, 2.0, 2.0, 1.0
            }
            ;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  NORMAL generates normal points");
        Console.WriteLine("    in M dimensions, using a nonzero mean, and with");
        Console.WriteLine("    user-specified variance-covariance matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        typeMethods.r8vec_print(DIM_NUM, mu, "  Mean vector MU:");

        typeMethods.r8mat_print(DIM_NUM, DIM_NUM, v, "  Variance-covariance matrix V:");

        for (i = 0; i < DIM_NUM; i++)
        {
            for (j = 0; j < DIM_NUM; j++)
            {
                r[i + j * DIM_NUM] = v[i + j * DIM_NUM];
            }
        }

        info = typeMethods.r8po_fa(ref r, DIM_NUM, DIM_NUM);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST04 - Fatal error!");
            Console.WriteLine("  Variance-covariance matrix factorization failed.");
            return;
        }

        typeMethods.r8mat_print(DIM_NUM, DIM_NUM, r, "  Cholesky factor R:");

        x = Normal.normal(DIM_NUM, N, r, mu, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        Scale.scale_to_block01(DIM_NUM, N, ref x);

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests NORMAL_CIRCULAR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 2000;

        string output_filename = "normal_circular.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  NORMAL_CIRCULAR generates points in 2D");
        Console.WriteLine("    distributed according to a circular normal.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Normal.normal_circular(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        Scale.scale_to_block01(DIM_NUM, N, ref x);

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests NORMAL_SIMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "normal_simple.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  NORMAL_SIMPLE generates normal points");
        Console.WriteLine("    in M dimensions, using a zero mean, and with");
        Console.WriteLine("    the identity as the variance-covariance matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Normal.normal_simple(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        Scale.scale_to_block01(DIM_NUM, N, ref x);

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test115()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST115 tests UNIFORM_IN_ANNULUS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 400;

        string output_filename = "uniform_in_annulus.txt";
        double[] pc =  {
                10.0, 5.0
            }
            ;
        double r1 = 1.0;
        double r2 = 3.0;
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST115");
        Console.WriteLine("  UNIFORM_IN_ANNULUS generates uniform");
        Console.WriteLine("  points in an annulus by mapping.");
        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Center PC(1:2) =             " + pc[0] + "  " + pc[1] + "");
        Console.WriteLine("  Inner radius is R1 =         " + r1 + "");
        Console.WriteLine("  Outer radius is R2 =         " + r2 + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Annulus.uniform_in_annulus(pc, r1, r2, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests UNIFORM_IN_ANNULUS_ACCEPT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 400;

        string output_filename = "uniform_in_annulus_accept.txt";
        double[] pc =  {
                10.0, 5.0
            }
            ;
        double r1 = 1.0;
        double r2 = 3.0;
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  UNIFORM_IN_ANNULUS_ACCEPT generates uniform");
        Console.WriteLine("  points in an annulus by acceptance/rejection.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Center PC(1:2) =             " + pc[0] + "  " + pc[1] + "");
        Console.WriteLine("  Inner radius is R1 =         " + r1 + "");
        Console.WriteLine("  Outer radius is R2 =         " + r2 + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Annulus.uniform_in_annulus_accept(pc, r1, r2, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test125()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST125 tests UNIFORM_IN_ANNULUS_SECTOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 400;

        string output_filename = "uniform_in_annulus_sector.txt";
        double[] pc =  {
                10.0, 5.0
            }
            ;
        double r1 = 1.0;
        double r2 = 3.0;
        int seed = 123456789;
        double theta1 = 0.0;
        double theta2 = 1.5707964;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST125");
        Console.WriteLine("  UNIFORM_IN_ANNULUS_SECTOR generates uniform ");
        Console.WriteLine("  points in an annular sector by mapping.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Center PC(1:2) =             " + pc[0] + "  " + pc[1] + "");
        Console.WriteLine("  Inner radius is R1 =         " + r1 + "");
        Console.WriteLine("  Outer radius is R2 =         " + r2 + "");
        Console.WriteLine("  THETA1 =                     " + theta1 + "");
        Console.WriteLine("  THETA2 =                     " + theta2 + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Annulus.uniform_in_annulus_sector(pc, r1, r2, theta1, theta2, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests UNIFORM_IN_CIRCLE01_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 400;

        string output_filename = "uniform_in_circle01_map.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  UNIFORM_IN_CIRCLE01_MAP maps uniform ");
        Console.WriteLine("  points into the unit circle.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Circle.uniform_in_circle01_map(N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests UNIFORM_IN_CUBE01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "uniform_in_cube01.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  UNIFORM_IN_CUBE01 generates uniform");
        Console.WriteLine("  points in the unit hypercube.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Cube.uniform_in_cube01(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests UNIFORM_IN_ELLIPSOID_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        double[] a =  {
                3.0, 1.0, 1.0, 2.0
            }
            ;
        int fail_num;
        string output_filename = "uniform_in_ellipsoid_map.txt";
        int i;
        int j;
        int k;
        double r = 1.0;
        double r2;
        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();
        int success_num;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  UNIFORM_IN_ELLIPSOID_MAP maps uniform");
        Console.WriteLine("  points into an ellipsoid.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Ellipsoid.uniform_in_ellipsoid_map(DIM_NUM, N, a, r, ref data, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
        //
        //  Test the data.
        //
        fail_num = 0;
        success_num = 0;

        for (j = 0; j < N; j++)
        {

            r2 = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                for (k = 0; k < DIM_NUM; k++)
                {
                    r2 += x[i + j * DIM_NUM] * a[i + k * DIM_NUM] * x[k + j * DIM_NUM];
                }
            }

            r2 = Math.Sqrt(r2);

            if (r < r2)
            {
                fail_num += 1;
            }
            else
            {
                success_num += 1;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  " + fail_num
                               + " points failed the ellipsoid test.");

        Console.WriteLine("  " + success_num
                               + " points satisfy the ellipsoid test.");
    }

    private static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests UNIFORM_IN_PARALLELOGRAM_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "uniform_in_parallelogram_map.txt";
        int seed = 123456789;
        double[] v1 =  {
                0.75E+00, 0.90E+00
            }
            ;
        double[] v2 =  {
                0.00E+00, 0.20E+00
            }
            ;
        double[] v3 =  {
                1.10E+00, 0.65E+00
            }
            ;
        double[] v4 = new double[DIM_NUM];
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  UNIFORM_IN_PARALLELOGRAM_MAP maps uniform");
        Console.WriteLine("  points into a parallelogram.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        v4[0] = v3[0] + v2[0] - v1[0];
        v4[1] = v3[1] + v2[1] - v1[1];

        Console.WriteLine("");
        Console.WriteLine("  V1 = " + v1[0] + ", " + v1[1] + "");
        Console.WriteLine("  V2 = " + v2[0] + ", " + v2[1] + "");
        Console.WriteLine("  V3 = " + v3[0] + ", " + v3[1] + "");
        Console.WriteLine("  V4 = " + v4[0] + ", " + v4[1] + "");

        x = Parallelogram.uniform_in_parallelogram_map(v1, v2, v3, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        Scale.scale_to_block01(DIM_NUM, N, ref x);

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test17()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17 tests UNIFORM_IN_POLYGON_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;
        int NV = 10;

        string output_filename = "uniform_in_polygon_map.txt";
        int seed = 123456789;
        double[] v =  {
                0.0E+00, 0.0E+00,
                0.5E+00, 0.3E+00,
                1.0E+00, 0.0E+00,
                0.7E+00, 0.4E+00,
                1.0E+00, 0.6E+00,
                0.6E+00, 0.6E+00,
                0.5E+00, 1.0E+00,
                0.4E+00, 0.6E+00,
                0.0E+00, 0.6E+00,
                0.3E+00, 0.4E+00
            }
            ;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST17");
        Console.WriteLine("  UNIFORM_IN_POLYGON_MAP maps uniform");
        Console.WriteLine("  points into a polygon.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        typeMethods.r8mat_print(DIM_NUM, NV, v, "  Polygonal vertices:");

        x = Polygon.uniform_in_polygon_map(NV, v, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write("polygon_vertices.txt", DIM_NUM, NV, v);

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test18()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests UNIFORM_IN_SECTOR_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 300;

        string output_filename = "uniform_in_sector_map.txt";
        double r1 = 1.0;
        double r2 = 2.0;
        int seed = 123456789;
        double t1 = 0.78;
        double t2 = 2.35;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST18");
        Console.WriteLine("  UNIFORM_IN_SECTOR_MAP maps uniform");
        Console.WriteLine("  points into a circular sector.");
        Console.WriteLine("");
        Console.WriteLine("  R1 = " + r1 + "");
        Console.WriteLine("  R2 = " + r2 + "");
        Console.WriteLine("  T1 = " + t1 + "");
        Console.WriteLine("  T2 = " + t2 + "");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Sector.uniform_in_sector_map(r1, r2, t1, t2, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test19()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST19 tests UNIFORM_IN_SIMPLEX01_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "uniform_in_simplex01_map.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST19");
        Console.WriteLine("  UNIFORM_IN_SIMPLEX01_MAP maps uniform");
        Console.WriteLine("  points into the unit simplex");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Simplex.uniform_in_simplex01_map(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test20()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20 tests UNIFORM_IN_SPHERE01_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "uniform_in_sphere01_map.txt";
        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST20");
        Console.WriteLine("  UNIFORM_IN_SPHERE01_MAP maps uniform");
        Console.WriteLine("  points into the unit sphere.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Sphere.uniform_in_sphere01_map(DIM_NUM, N, ref data, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test205()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST205 tests UNIFORM_IN_TETRAHEDRON.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int N = 1000;

        string output_filename = "uniform_in_tetrahedron.txt";
        int seed = 123456789;
        double[] v =  {
                1.0, 2.0, 3.0,
                4.0, 1.0, 2.0,
                2.0, 4.0, 4.0,
                3.0, 2.0, 5.0
            }
            ;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST205");
        Console.WriteLine("  UNIFORM_IN_TETRAHEDRON returns uniform");
        Console.WriteLine("  points from a tetrahedron.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        typeMethods.r8mat_print(3, 4, v, "  Tetrahedron vertices:");

        x = Tetrahedron.uniform_in_tetrahedron(v, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void uniform_in_triangle_map1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_TRIANGLE_MAP1_TEST tests UNIFORM_IN_TRIANGLE_MAP1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "uniform_in_triangle_map1.txt";
        int seed = 123456789;
        double[] v1 =  {
                0.75E+00, 0.90E+00
            }
            ;
        double[] v2 =  {
                0.00E+00, 0.20E+00
            }
            ;
        double[] v3 =  {
                0.95E+00, 0.65E+00
            }
            ;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("UNIFORM_IN_TRIANGLE_MAP1_TEST");
        Console.WriteLine("  UNIFORM_IN_TRIANGLE_MAP1 maps uniform");
        Console.WriteLine("  points into a triangle, by Turk 1 mapping.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  V1 = " + v1[0] + ", " + v1[1] + "");
        Console.WriteLine("  V2 = " + v2[0] + ", " + v2[1] + "");
        Console.WriteLine("  V3 = " + v3[0] + ", " + v3[1] + "");

        x = Triangle.uniform_in_triangle_map1(v1, v2, v3, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test22()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST22 tests UNIFORM_IN_TRIANGLE_MAP2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 1000;

        string output_filename = "uniform_in_triangle_map2.txt";
        int seed = 123456789;
        double[] v1 =  {
                0.75E+00, 0.90E+00
            }
            ;
        double[] v2 =  {
                0.00E+00, 0.20E+00
            }
            ;
        double[] v3 =  {
                0.95E+00, 0.65E+00
            }
            ;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST22");
        Console.WriteLine("  UNIFORM_IN_TRIANGLE_MAP2 maps uniform");
        Console.WriteLine("  points into a triangle, by Turk 2 mapping.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  V1 = " + v1[0] + ", " + v1[1] + "");
        Console.WriteLine("  V2 = " + v2[0] + ", " + v2[1] + "");
        Console.WriteLine("  V3 = " + v3[0] + ", " + v3[1] + "");

        x = Triangle.uniform_in_triangle_map2(v1, v2, v3, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test23()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST23 tests UNIFORM_IN_TRIANGLE01_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 2000;

        string output_filename = "uniform_in_triangle01_map.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST23");
        Console.WriteLine("  UNIFORM_IN_TRIANGLE01_MAP maps uniform");
        Console.WriteLine("  points into the unit triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Triangle.uniform_in_triangle01_map(N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test235()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST235 tests UNIFORM_ON_CUBE01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 April 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 2;
        int N = 200;

        string output_filename = "uniform_on_cube01.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST235");
        Console.WriteLine("  UNIFORM_ON_CUBE01 samples N uniform points on");
        Console.WriteLine("  the surface of the unit M-dimensional cube.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M =        " + M + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Cube.uniform_on_cube01(M, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, M, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test24()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST24 tests UNIFORM_ON_ELLIPSOID_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 200;

        double[] a =  {
                3.0, 1.0, 1.0, 2.0
            }
            ;
        string output_filename = "uniform_on_ellipsoid_map.txt";
        double r = 1.0;
        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST24");
        Console.WriteLine("  UNIFORM_ON_ELLIPSOID_MAP maps uniform");
        Console.WriteLine("  points onto an ellipsoid.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Ellipsoid.uniform_on_ellipsoid_map(DIM_NUM, N, a, r, ref data, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test245()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST245 tests UNIFORM_ON_HEMISPHERE01_PHONG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int N = 50;

        string output_filename = "uniform_on_hemisphere01_phong.txt";
        int m = 2;
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST245");
        Console.WriteLine("  UNIFORM_ON_HEMISPHERE01_PHONG maps uniform");
        Console.WriteLine("  points onto the unit hemisphere with Phong density.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Phong exponent M =           " + m + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Sphere.uniform_on_hemisphere01_phong(N, m, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test25()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST25 tests UNIFORM_ON_SIMPLEX01_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 50;

        string output_filename = "uniform_on_simplex01_map.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST25");
        Console.WriteLine("  UNIFORM_ON_SIMPLEX01_MAP maps uniform ");
        Console.WriteLine("  points onto the unit simplex.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Simplex.uniform_on_simplex01_map(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test26()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST26 tests UNIFORM_ON_SPHERE01_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 50;

        string output_filename = "uniform_on_sphere01_map.txt";
        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST26");
        Console.WriteLine("  UNIFORM_ON_SPHERE01_MAP maps uniform");
        Console.WriteLine("  points onto the unit sphere.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Sphere.uniform_on_sphere01_map(DIM_NUM, N, ref data, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test264()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST264 tests UNIFORM_ON_SPHERE01_PATCH_TP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 50;

        string output_filename = "uniform_on_sphere01_patch_tp.txt";
        double phi1;
        double phi2;
        int seed = 123456789;
        double theta1;
        double theta2;
        double[] tp;

        phi1 = 75.0 * (Math.PI / 180.0);
        phi2 = 90.0 * (Math.PI / 180.0);
        theta1 = 0.0 * (Math.PI / 360.0);
        theta2 = 30.0 * (Math.PI / 360.0);

        Console.WriteLine("");
        Console.WriteLine("TEST264");
        Console.WriteLine("  UNIFORM_ON_SPHERE01_PATCH_TP maps uniform");
        Console.WriteLine("  points onto a TP (THETA,PHI) patch of the unit sphere.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension =          " + 3 + "");
        Console.WriteLine("  Data dimension =             " + 2 + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Latitudinal angle PHI1 =     " + phi1 + "");
        Console.WriteLine("  Latitudinal angle PHI2 =     " + phi2 + "");
        Console.WriteLine("  Longitudinal angle THETA1 =  " + theta1 + "");
        Console.WriteLine("  Longitudinal angle THETA2 =  " + theta2 + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        tp = Sphere.uniform_on_sphere01_patch_tp(N, phi1, phi2, theta1, theta2, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, 2, N, tp);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test265()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST265 tests UNIFORM_ON_SPHERE01_PATCH_XYZ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int N = 50;

        string output_filename = "uniform_on_sphere01_patch_xyz.txt";
        double phi1;
        double phi2;
        int seed = 123456789;
        double theta1;
        double theta2;
        double[] x;

        phi1 = 75.0 * (Math.PI / 180.0);
        phi2 = 90.0 * (Math.PI / 180.0);
        theta1 = 0.0 * (Math.PI / 360.0);
        theta2 = 30.0 * (Math.PI / 360.0);

        Console.WriteLine("");
        Console.WriteLine("TEST265");
        Console.WriteLine("  UNIFORM_ON_SPHERE01_PATCH_XYZ maps uniform");
        Console.WriteLine("  points onto an XYZ patch of the unit sphere.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Latitudinal angle PHI1 =     " + phi1 + "");
        Console.WriteLine("  Latitudinal angle PHI2 =     " + phi2 + "");
        Console.WriteLine("  Longitudinal angle THETA1 =  " + theta1 + "");
        Console.WriteLine("  Longitudinal angle THETA2 =  " + theta2 + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Sphere.uniform_on_sphere01_patch_xyz(N, phi1, phi2, theta1, theta2, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);
            
        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test267()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST267 tests UNIFORM_ON_SPHERE01_TRIANGLE_XYZ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int N = 500;

        string output_filename = "uniform_on_sphere01_triangle_xyz.txt";
        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();
        double[] v1;
        double[] v2;
        double[] v3;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST267");
        Console.WriteLine("  UNIFORM_ON_SPHERE01_TRIANGLE_XYZ maps uniform");
        Console.WriteLine("  points onto a spherical triangle using XYZ coordinates.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + 3 + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        switch (true)
        {
            case true:
                v1 = Sphere.uniform_on_sphere01_map(3, 1, ref data, ref seed);
                v2 = Sphere.uniform_on_sphere01_map(3, 1, ref data, ref seed);
                v3 = Sphere.uniform_on_sphere01_map(3, 1, ref data, ref seed);
                break;
            default:
                v1 = typeMethods.r8vec_zero_new(3);
                v1[0] = 1.0;
                v2 = typeMethods.r8vec_zero_new(3);
                v2[1] = 1.0;
                v3 = typeMethods.r8vec_zero_new(3);
                v3[2] = 1.0;
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Vertices of spherical triangle:");
        Console.WriteLine("");
        Console.WriteLine("  V1: (" + v1[0] + ", " + v1[1] + ", " + v1[2] + ")");
        Console.WriteLine("  V2: (" + v2[0] + ", " + v2[1] + ", " + v2[2] + ")");
        Console.WriteLine("  V3: (" + v3[0] + ", " + v3[1] + ", " + v3[2] + ")");

        x = Sphere.uniform_on_sphere01_triangle_xyz(N, v1, v2, v3, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Final random number SEED =   " + seed + "");

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }

    private static void test27()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST27 tests UNIFORM_WALK
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 400;

        string output_filename = "uniform_walk.txt";
        int seed = 123456789;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST27:");
        Console.WriteLine("  UNIFORM_WALK generates uniform random walk points.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + DIM_NUM + "");
        Console.WriteLine("  Number of points N =         " + N + "");
        Console.WriteLine("  Initial random number SEED = " + seed + "");

        x = Walk.uniform_walk(DIM_NUM, N, ref seed);

        Console.WriteLine("  Final random number SEED =   " + seed + "");

        Scale.scale_to_block01(DIM_NUM, N, ref x);

        typeMethods.r8mat_write(output_filename, DIM_NUM, N, x);

        Console.WriteLine("");
        Console.WriteLine("  Data written to file \"" + output_filename + "\".");
    }
}