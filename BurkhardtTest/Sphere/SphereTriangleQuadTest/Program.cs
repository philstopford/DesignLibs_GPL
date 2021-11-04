using System;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace SphereTriangleQuadTest
{
    class Program
    {
        static int[] e_save = new int[3];

        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SPHERE_TRIANGLE_QUAD_TEST.
            //
            //  Discussion:
            //
            //    SPHERE_TRIANGLE_QUAD_TEST tests the SPHERE_TRIANGLE_QUAD library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_TRIANGLE_QUAD_TEST");
            Console.WriteLine("  Test the SPHERE_TRIANGLE_QUAD library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();

            Console.WriteLine("");
            Console.WriteLine("SPHERE_TRIANGLE_QUAD_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests SPHERE01_TRIANGLE_QUAD_01, 02, 03.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] e = new int[3];
            int i;
            double result_01;
            double result_02;
            double result_03;
            int seed;
            double[] v1;
            double[] v2;
            double[] v3;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Approximate the integral of a function on a random spherical triangle.");
            Console.WriteLine("");
            Console.WriteLine("  QUAD_01 uses centroids of spherical triangles.");
            Console.WriteLine("  QUAD_02 uses vertices of spherical triangles.");
            Console.WriteLine("  QUAD_03 uses midsides of spherical triangles.");
            //
            //  Choose three points at random to define a spherical triangle.
            //
            v1 = MonteCarlo.sphere01_sample(1, ref seed);
            v2 = MonteCarlo.sphere01_sample(1, ref seed);
            v3 = MonteCarlo.sphere01_sample(1, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle:");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(3, v1, "  V1:");
            typeMethods.r8vec_transpose_print(3, v2, "  V2:");
            typeMethods.r8vec_transpose_print(3, v3, "  V3:");

            Console.WriteLine("");
            Console.WriteLine("QUAD_01      QUAD_02      QUAD_03");

            for (i = 1; i <= 17; i++)
            {
                if (i == 1)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 2)
                {
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 3)
                {
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                }
                else if (i == 4)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                }
                else if (i == 5)
                {
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 6)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 7)
                {
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 8)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 9)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                }
                else if (i == 10)
                {
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 11)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                }
                else if (i == 12)
                {
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                }
                else if (i == 13)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                }
                else if (i == 14)
                {
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                }
                else if (i == 15)
                {
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                }
                else if (i == 16)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                }
                else if (i == 17)
                {
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                }

                polyterm_exponent("SET", e);

                polyterm_exponent("PRINT", e);

                result_01 = TriangleQuad.sphere01_triangle_quad_01(v1, v2, v3, polyterm_value_3d);

                result_02 = TriangleQuad.sphere01_triangle_quad_02(v1, v2, v3, polyterm_value_3d);

                result_03 = TriangleQuad.sphere01_triangle_quad_03(v1, v2, v3, polyterm_value_3d);

                Console.WriteLine("  " + result_01.ToString().PadLeft(14)
                                       + "  " + result_02.ToString().PadLeft(14)
                                       + "  " + result_03.ToString().PadLeft(14) + "");

            }

            return;
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests SPHERE01_TRIANGLE_QUAD_00.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] e = new int[3];
            int i;
            int n_mc1 = 1000;
            int n_mc2 = 10000;
            int n_mc3 = 100000;
            double result_01;
            double result_02;
            double result_03;
            int seed;
            double[] v1;
            double[] v2;
            double[] v3;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Approximate the integral of a function on a random spherical triangle.");
            Console.WriteLine("");
            Console.WriteLine("  SPHERE01_TRIANGLE_QUAD_00 uses the Monte Carlo method.");
            Console.WriteLine("  QUAD_MC1 uses a Monte Carlo method with " + n_mc1 + " points.");
            Console.WriteLine("  QUAD_MC2 uses a Monte Carlo method with " + n_mc2 + " points.");
            Console.WriteLine("  QUAD_MC3 uses a Monte Carlo method with " + n_mc3 + " points.");
            //
            //  Choose three points at random to define a spherical triangle.
            //
            v1 = MonteCarlo.sphere01_sample(1, ref seed);
            v2 = MonteCarlo.sphere01_sample(1, ref seed);
            v3 = MonteCarlo.sphere01_sample(1, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle:");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(3, v1, "  V1:");
            typeMethods.r8vec_transpose_print(3, v2, "  V2:");
            typeMethods.r8vec_transpose_print(3, v3, "  V3:");

            Console.WriteLine("");
            Console.WriteLine("QUAD_MC1     QUAD_MC2     QUAD_MC3");

            for (i = 1; i <= 17; i++)
            {
                if (i == 1)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 2)
                {
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 3)
                {
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                }
                else if (i == 4)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                }
                else if (i == 5)
                {
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 6)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 7)
                {
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 8)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 9)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                }
                else if (i == 10)
                {
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 11)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                }
                else if (i == 12)
                {
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                }
                else if (i == 13)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                }
                else if (i == 14)
                {
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                }
                else if (i == 15)
                {
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                }
                else if (i == 16)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                }
                else if (i == 17)
                {
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                }

                polyterm_exponent("SET", e);

                polyterm_exponent("PRINT", e);

                result_01 = TriangleQuad.sphere01_triangle_quad_00(n_mc1, v1, v2, v3,
                    polyterm_value_3d, ref seed);

                result_02 = TriangleQuad.sphere01_triangle_quad_00(n_mc2, v1, v2, v3,
                    polyterm_value_3d, ref seed);

                result_03 = TriangleQuad.sphere01_triangle_quad_00(n_mc3, v1, v2, v3,
                    polyterm_value_3d, ref seed);

                Console.WriteLine("  " + result_01.ToString().PadLeft(14)
                                       + "  " + result_02.ToString().PadLeft(14)
                                       + "  " + result_03.ToString().PadLeft(14) + "");

            }

        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests SPHERE01_TRIANGLE_QUAD_ICOS1C.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double best;
            int[] e = new int[3];
            double error;
            int factor;
            int factor_log;
            int i;
            int node_num = 0;
            double result;
            int seed;
            double[] v1;
            double[] v2;
            double[] v3;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  Approximate the integral of a function on a random spherical triangle.");
            Console.WriteLine("");
            Console.WriteLine("  SPHERE01_TRIANGLE_QUAD_ICOS1C approximates the");
            Console.WriteLine("  integral of a function over a spherical triangle on");
            Console.WriteLine("  the surface of the unit sphere using a centroid rule.");
            Console.WriteLine("");
            Console.WriteLine("  We do not have an exact result, so we compare each");
            Console.WriteLine("  estimate to the final one.");
            //
            //  Choose three points at random to define a spherical triangle.
            //
            v1 = MonteCarlo.sphere01_sample(1, ref seed);
            v2 = MonteCarlo.sphere01_sample(1, ref seed);
            v3 = MonteCarlo.sphere01_sample(1, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle:");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(3, v1, "  V1:");
            typeMethods.r8vec_transpose_print(3, v2, "  V2:");
            typeMethods.r8vec_transpose_print(3, v3, "  V3:");

            Console.WriteLine("");
            Console.WriteLine("FACTOR   N   RESULT");

            for (i = 1; i <= 17; i++)
            {
                if (i == 1)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 2)
                {
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 3)
                {
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                }
                else if (i == 4)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                }
                else if (i == 5)
                {
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 6)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 7)
                {
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 8)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 9)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                }
                else if (i == 10)
                {
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 11)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                }
                else if (i == 12)
                {
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                }
                else if (i == 13)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                }
                else if (i == 14)
                {
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                }
                else if (i == 15)
                {
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                }
                else if (i == 16)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                }
                else if (i == 17)
                {
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                }

                polyterm_exponent("SET", e);

                polyterm_exponent("PRINT", e);

                factor = (int) Math.Pow(2, 9);

                best = TriangleQuad.sphere01_triangle_quad_icos1c(v1, v2, v3, factor,
                    polyterm_value_3d, ref node_num);

                factor = 1;

                for (factor_log = 0; factor_log <= 9; factor_log++)
                {
                    result = TriangleQuad.sphere01_triangle_quad_icos1c(v1, v2, v3, factor,
                        polyterm_value_3d, ref node_num);

                    error = Math.Abs(result - best);

                    Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                           + "  " + node_num.ToString().PadLeft(8)
                                           + "  " + result.ToString().PadLeft(16)
                                           + "  " + error.ToString().PadLeft(10) + "");

                    factor = factor * 2;
                }
            }

            return;
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests SPHERE01_TRIANGLE_QUAD_ICOS1M.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double best;
            int[] e = new int[3];
            double error;
            int factor;
            int factor_log;
            int i;
            int node_num = 0;
            double result;
            int seed;
            double[] v1;
            double[] v2;
            double[] v3;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  Approximate the integral of a function on a random spherical triangle.");
            Console.WriteLine("");
            Console.WriteLine("  SPHERE01_TRIANGLE_QUAD_ICOS1M approximates the");
            Console.WriteLine("  integral of a function over a spherical triangle on");
            Console.WriteLine("  the surface of the unit sphere using a midside rule.");
            Console.WriteLine("");
            Console.WriteLine("  We do not have an exact result, so we compare each");
            Console.WriteLine("  estimate to the final one.");
            //
            //  Choose three points at random to define a spherical triangle.
            //
            v1 = MonteCarlo.sphere01_sample(1, ref seed);
            v2 = MonteCarlo.sphere01_sample(1, ref seed);
            v3 = MonteCarlo.sphere01_sample(1, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle:");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(3, v1, "  V1:");
            typeMethods.r8vec_transpose_print(3, v2, "  V2:");
            typeMethods.r8vec_transpose_print(3, v3, "  V3:");

            Console.WriteLine("");
            Console.WriteLine("FACTOR   N   RESULT");

            for (i = 1; i <= 17; i++)
            {
                if (i == 1)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 2)
                {
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 3)
                {
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                }
                else if (i == 4)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                }
                else if (i == 5)
                {
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 6)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 7)
                {
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 8)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 9)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                }
                else if (i == 10)
                {
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 11)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                }
                else if (i == 12)
                {
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                }
                else if (i == 13)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                }
                else if (i == 14)
                {
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                }
                else if (i == 15)
                {
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                }
                else if (i == 16)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                }
                else if (i == 17)
                {
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                }

                polyterm_exponent("SET", e);

                polyterm_exponent("PRINT", e);

                factor = (int) Math.Pow(2, 9);

                best = TriangleQuad.sphere01_triangle_quad_icos1m(v1, v2, v3, factor,
                    polyterm_value_3d, ref node_num);

                factor = 1;

                for (factor_log = 0; factor_log <= 9; factor_log++)
                {
                    result = TriangleQuad.sphere01_triangle_quad_icos1m(v1, v2, v3, factor,
                        polyterm_value_3d, ref node_num);

                    error = Math.Abs(result - best);

                    Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                           + "  " + node_num.ToString().PadLeft(8)
                                           + "  " + result.ToString().PadLeft(16)
                                           + "  " + error.ToString().PadLeft(10) + "");

                    factor = factor * 2;
                }
            }

        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests SPHERE01_TRIANGLE_QUAD_ICOS1V.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double best;
            int[] e = new int[3];
            double error;
            int factor;
            int factor_log;
            int i;
            int node_num = 0;
            double result;
            int seed;
            double[] v1;
            double[] v2;
            double[] v3;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Approximate the integral of a function on a random spherical triangle.");
            Console.WriteLine("");
            Console.WriteLine("  SPHERE01_TRIANGLE_QUAD_ICOS1V approximates the");
            Console.WriteLine("  integral of a function over a spherical triangle on");
            Console.WriteLine("  the surface of the unit sphere using a vertex rule.");
            Console.WriteLine("");
            Console.WriteLine("  We do not have an exact result, so we compare each");
            Console.WriteLine("  estimate to the final one.");
            //
            //  Choose three points at random to define a spherical triangle.
            //
            v1 = MonteCarlo.sphere01_sample(1, ref seed);
            v2 = MonteCarlo.sphere01_sample(1, ref seed);
            v3 = MonteCarlo.sphere01_sample(1, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle:");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(3, v1, "  V1:");
            typeMethods.r8vec_transpose_print(3, v2, "  V2:");
            typeMethods.r8vec_transpose_print(3, v3, "  V3:");

            Console.WriteLine("");
            Console.WriteLine("FACTOR   N   RESULT");

            for (i = 1; i <= 17; i++)
            {
                if (i == 1)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 2)
                {
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 3)
                {
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                }
                else if (i == 4)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                }
                else if (i == 5)
                {
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 6)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 7)
                {
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 8)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 9)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                }
                else if (i == 10)
                {
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 11)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                }
                else if (i == 12)
                {
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                }
                else if (i == 13)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                }
                else if (i == 14)
                {
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                }
                else if (i == 15)
                {
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                }
                else if (i == 16)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                }
                else if (i == 17)
                {
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                }

                polyterm_exponent("SET", e);

                polyterm_exponent("PRINT", e);

                factor = (int) Math.Pow(2, 9);

                best = TriangleQuad.sphere01_triangle_quad_icos1v(v1, v2, v3, factor,
                    polyterm_value_3d, ref node_num);

                factor = 1;

                for (factor_log = 0; factor_log <= 9; factor_log++)
                {
                    result = TriangleQuad.sphere01_triangle_quad_icos1v(v1, v2, v3, factor,
                        polyterm_value_3d, ref node_num);

                    error = Math.Abs(result - best);

                    Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                           + "  " + node_num.ToString().PadLeft(8)
                                           + "  " + result.ToString().PadLeft(16)
                                           + "  " + error.ToString().PadLeft(10) + "");

                    factor = factor * 2;
                }
            }

            return;
        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests SPHERE01_TRIANGLE_QUAD_ICOS2V.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double best;
            int[] e = new int[3];
            double error;
            int factor;
            int factor_log;
            int i;
            int node_num = 0;
            double result;
            int seed;
            double[] v1;
            double[] v2;
            double[] v3;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  Approximate the integral of a function on a random spherical triangle.");
            Console.WriteLine("");
            Console.WriteLine("  SPHERE01_TRIANGLE_QUAD_ICOS2V approximates the");
            Console.WriteLine("  integral of a function over a spherical triangle on");
            Console.WriteLine("  the surface of the unit sphere using a vertex rule.");
            Console.WriteLine("");
            Console.WriteLine("  We do not have an exact result, so we compare each");
            Console.WriteLine("  estimate to the final one.");
            //
            //  Choose three points at random to define a spherical triangle.
            //
            v1 = MonteCarlo.sphere01_sample(1, ref seed);
            v2 = MonteCarlo.sphere01_sample(1, ref seed);
            v3 = MonteCarlo.sphere01_sample(1, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle:");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(3, v1, "  V1:");
            typeMethods.r8vec_transpose_print(3, v2, "  V2:");
            typeMethods.r8vec_transpose_print(3, v3, "  V3:");

            Console.WriteLine("");
            Console.WriteLine("FACTOR   N   RESULT");

            for (i = 1; i <= 17; i++)
            {
                if (i == 1)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 2)
                {
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 3)
                {
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                }
                else if (i == 4)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                }
                else if (i == 5)
                {
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                }
                else if (i == 6)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 7)
                {
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                }
                else if (i == 8)
                {
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 9)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                }
                else if (i == 10)
                {
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                }
                else if (i == 11)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                }
                else if (i == 12)
                {
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                }
                else if (i == 13)
                {
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                }
                else if (i == 14)
                {
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                }
                else if (i == 15)
                {
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                }
                else if (i == 16)
                {
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                }
                else if (i == 17)
                {
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                }

                polyterm_exponent("SET", e);

                polyterm_exponent("PRINT", e);

                factor = (int) Math.Pow(2, 9);

                best = TriangleQuad.sphere01_triangle_quad_icos2v(v1, v2, v3, factor,
                    polyterm_value_3d, ref node_num);

                factor = 1;

                for (factor_log = 0; factor_log <= 9; factor_log++)
                {
                    result = TriangleQuad.sphere01_triangle_quad_icos2v(v1, v2, v3, factor,
                        polyterm_value_3d, ref node_num);

                    error = Math.Abs(result - best);

                    Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                           + "  " + node_num.ToString().PadLeft(8)
                                           + "  " + result.ToString().PadLeft(16)
                                           + "  " + error.ToString().PadLeft(10) + "");

                    factor = factor * 2;
                }
            }

            return;
        }

        static void polyterm_exponent(string action, int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string ACTION.
            //    'GET' asks the routine to return the current values in E.
            //    'SET' asks the routine to set the current values to E.
            //
            //    Input/output, int E[3], storage used to set or get values.
            //
        {
            int i;

            if (action[0] == 'G')
            {
                for (i = 0; i < 3; i++)
                {
                    e[i] = e_save[i];
                }
            }
            else if (action[0] == 'P')
            {
                Console.WriteLine("");

                if (e_save[0] == 0 && e_save[1] == 0 && e_save[2] == 0)
                {
                    Console.WriteLine("P(X,Y,Z) = 1");
                }
                else
                {
                    string cout = "P(X,Y,Z) = ";

                    if (e_save[0] == 0)
                    {
                    }
                    else if (e_save[0] == 1)
                    {
                        cout += " X";
                    }
                    else
                    {
                        cout += " X^" + e_save[0];
                    }

                    if (e_save[1] == 0)
                    {
                    }
                    else if (e_save[1] == 1)
                    {
                        cout += " Y";
                    }
                    else
                    {
                        cout += " Y^" + e_save[1];
                    }

                    if (e_save[2] == 0)
                    {
                    }
                    else if (e_save[2] == 1)
                    {
                        cout += " Z";
                    }
                    else
                    {
                        cout += " Z^" + e_save[2];
                    }

                    Console.WriteLine(cout);
                }
            }
            else if (action[0] == 'S')
            {
                for (i = 0; i < 3; i++)
                {
                    e_save[i] = e[i];
                }
            }

        }

        static double polyterm_value_3d(double[] x, int index)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
            //
            //  Discussion:
            //
            //    The polynomial term has the form:
            //
            //      F(X) = X(1)^E(1) * X(2)^E(2) * X(3)^E(3)
            //
            //    The exponents E(1:3) are set by calling POLYTERM_EXPONENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X[3], the points where the polynomial term 
            //    is to be evaluated.
            //
            //    Output, double POLYTERM_VALUE_3D, the value of the polynomial term.
            //
        {
            int[] e = new int[3];
            int i;
            double value;

            polyterm_exponent("GET", e);

            value = 1.0;
            for (i = 0; i < 3; i++)
            {
                if (e[i] != 0)
                {
                    value = value * Math.Pow(x[(i + index) % x.Length], e[i]);
                }
            }

            return value;
        }
    }
}