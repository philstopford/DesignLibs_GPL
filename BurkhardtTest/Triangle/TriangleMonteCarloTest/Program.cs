using System;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace TriangleMonteCarloTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    TRIANGLE_MONTE_CARLO_TEST tests the TRIANGLE_MONTE_CARLO library.
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
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the TRIANGLE_MONTE_CARLO library.");
        //
        //  Try each sampler on the unit triangle, integrating X^2, X*Y, Y^2.
        //
        test01();
        test02();
        test03();
        test04();
        //
        //  Try each sampler on a general triangle, integrating a selection of functions.
        //
        test05();
        test06();
        test07();
        test08();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses TRIANGLE_SAMPLE_01 with an increasing number of points.
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
        int f_num = 3;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_01");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_03");
        Console.WriteLine("  Integration region is the unit triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a \"bad\" sampler.");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("     P_NUM      X^2             X*Y             Y^2");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_01,
                MonteCarlo.triangle_integrand_03, ref seed);

            Console.WriteLine("  " + p_num.ToString().PadLeft(8)
                                   + "  " + result[0].ToString().PadLeft(14)
                                   + "  " + result[1].ToString().PadLeft(14)
                                   + "  " + result[2].ToString().PadLeft(14) + "");

            p_num = 2 * p_num;
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses TRIANGLE_SAMPLE_02 with an increasing number of points.
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
        int f_num = 3;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_02");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_03");
        Console.WriteLine("  Integration region is the unit triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a good sampler.");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("     P_NUM      X^2             X*Y             Y^2");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_02,
                MonteCarlo.triangle_integrand_03, ref seed);

            Console.WriteLine("  " + p_num.ToString().PadLeft(8)
                                   + "  " + result[0].ToString().PadLeft(14)
                                   + "  " + result[1].ToString().PadLeft(14)
                                   + "  " + result[2].ToString().PadLeft(14) + "");

            p_num = 2 * p_num;
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses TRIANGLE_SAMPLE_03 with an increasing number of points.
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
        int f_num = 3;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_03");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_03");
        Console.WriteLine("  Integration region is the unit triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a good sampler.");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("     P_NUM      X^2             X*Y             Y^2");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_03,
                MonteCarlo.triangle_integrand_03, ref seed);

            Console.WriteLine("  " + p_num.ToString().PadLeft(8)
                                   + "  " + result[0].ToString().PadLeft(14)
                                   + "  " + result[1].ToString().PadLeft(14)
                                   + "  " + result[2].ToString().PadLeft(14) + "");

            p_num = 2 * p_num;

        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 uses TRIANGLE_SAMPLE_04 with an increasing number of points.
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
        int f_num = 3;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_04");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_03");
        Console.WriteLine("  Integration region is the unit triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a good sampler.");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("     P_NUM      X^2             X*Y             Y^2");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_04,
                MonteCarlo.triangle_integrand_03, ref seed);

            Console.WriteLine("  " + p_num.ToString().PadLeft(8)
                                   + "  " + result[0].ToString().PadLeft(14)
                                   + "  " + result[1].ToString().PadLeft(14)
                                   + "  " + result[2].ToString().PadLeft(14) + "");

            p_num = 2 * p_num;
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 uses TRIANGLE_SAMPLE_01 on a general triangle.
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
        int f_num = 8;
        int i;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            4.0, 1.0,
            8.0, 3.0,
            0.0, 9.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_01");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_USER");
        Console.WriteLine("  Integration region is over a general triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a \"bad\" sampler.");

        seed = 123456789;

        typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("     P_NUM");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_01,
                triangle_integrand_user, ref seed);

            string tmp = "  " + p_num.ToString().PadLeft(8);
            for (i = 0; i < f_num; i++)
            {
                tmp += "  " + result[i].ToString().PadLeft(12);
            }

            Console.WriteLine(tmp);

            p_num = 2 * p_num;
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 uses TRIANGLE_SAMPLE_02 on a general triangle.
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
        int f_num = 8;
        int i;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            4.0, 1.0,
            8.0, 3.0,
            0.0, 9.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_02");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_USER");
        Console.WriteLine("  Integration region is over a general triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a \"good\" sampler.");

        seed = 123456789;

        typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("     P_NUM");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_02,
                triangle_integrand_user, ref seed);

            string tmp = "  " + p_num.ToString().PadLeft(8);
            for (i = 0; i < f_num; i++)
            {
                tmp += "  " + result[i].ToString().PadLeft(12);
            }

            Console.WriteLine(tmp);

            p_num = 2 * p_num;
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 uses TRIANGLE_SAMPLE_03 on a general triangle.
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
        int f_num = 8;
        int i;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            4.0, 1.0,
            8.0, 3.0,
            0.0, 9.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_03");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_USER");
        Console.WriteLine("  Integration region is over a general triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a \"good\" sampler.");

        seed = 123456789;

        typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("     P_NUM");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_03,
                triangle_integrand_user, ref seed);

            string tmp = "  " + p_num.ToString().PadLeft(8);
            for (i = 0; i < f_num; i++)
            {
                tmp += "  " + result[i].ToString().PadLeft(12);
            }

            Console.WriteLine(tmp);

            p_num = 2 * p_num;
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 uses TRIANGLE_SAMPLE_04 on a general triangle.
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
        int f_num = 8;
        int i;
        int p_num;
        double[] result;
        int seed;
        double[] t =
        {
            4.0, 1.0,
            8.0, 3.0,
            0.0, 9.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  Sample using TRIANGLE_UNIT_SAMPLE_04");
        Console.WriteLine("  Integrate TRIANGLE_UNIT_INTEGRAND_USER");
        Console.WriteLine("  Integration region is over a general triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Use an increasing number of points P_NUM.");
        Console.WriteLine("  Note that the sample routine is a \"good\" sampler.");

        seed = 123456789;

        typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("     P_NUM");
        Console.WriteLine("");

        p_num = 1;

        while (p_num <= 65536)
        {
            result = MonteCarlo.triangle_monte_carlo(t, p_num, f_num, MonteCarlo.triangle_unit_sample_04,
                triangle_integrand_user, ref seed);

            string tmp = "  " + p_num.ToString().PadLeft(8);
            for (i = 0; i < f_num; i++)
            {
                tmp += "  " + result[i].ToString().PadLeft(12);
            }

            Console.WriteLine(tmp);

            p_num = 2 * p_num;
        }
    }

    private static double[] triangle_integrand_user(int p_num, double[] p, int f_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTEGRAND_USER evaluates 8 integrand functions defined by the user.
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
        //  Parameters:
        //
        //    Input, int P_NUM, the number of points.
        //
        //    Input, double P(2,P_NUM), the evaluation points.
        //
        //    Input, int F_NUM, the number of integrands.
        //
        //    Output, double FP(F_NUM,P_NUM), the integrand values.
        //
    {
        double[] fp;
        int j;

        fp = new double[f_num * p_num];

        for (j = 0; j < p_num; j++)
        {
            fp[0 + j * f_num] = 1.0;
            fp[1 + j * f_num] = p[0 + j * 2];
            fp[2 + j * f_num] = p[1 + j * 2];
            fp[3 + j * f_num] = p[0 + j * 2] * p[0 + j * 2];
            fp[4 + j * f_num] = p[0 + j * 2] * p[1 + j * 2];
            fp[5 + j * f_num] = p[1 + j * 2] * p[1 + j * 2];
            fp[6 + j * f_num] = p[0 + j * 2] * p[0 + j * 2] * p[1 + j * 2];
            fp[7 + j * f_num] = p[0 + j * 2] * p[0 + j * 2] * p[1 + j * 2] * p[1 + j * 2];
        }

        return fp;
    }
}