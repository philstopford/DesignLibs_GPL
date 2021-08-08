using System;
using Burkardt.Types;

namespace TetrahedronMonteCarloTest
{
    using MonteCarlo = Burkardt.TetrahedronNS.MonteCarlo;
    using Integrand = Burkardt.TetrahedronNS.Integrand;

    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TETRAHEDRON_MONTE_CARLO_TEST.
            //
            //  Discussion:
            //
            //    TETRAHEDRON_MONTE_CARLO_TEST tests the TETRAHEDRON_MONTE_CARLO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_MONTE_CARLO_TEST");
            Console.WriteLine("  Test the TETRAHEDRON_MONTE_CARLO library.");
            //
            //  Try each sampler on the unit tetrahedron, integrating quadratics.
            //
            test01();
            test02();
            test03();
            test04();
            //
            //  Try each sampler on a general tetrahedron, integrating a selection of functions.
            //
            test05();
            test06();
            test07();
            test08();
            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_MONTE_CARLO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 uses TETRAHEDRON_SAMPLE_01 with an increasing number of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                0.0, 0.0, 0.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_01");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_03");
            Console.WriteLine("  Integration region is the unit tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a \"bad\" sampler.");

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("     P_NUM      X^2             X*Y             X*Z" + 
                              "             Y^2             Y*Z             Z^2");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_01, Integrand.tetrahedron_integrand_03, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;

            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 uses TETRAHEDRON_SAMPLE_02 with an increasing number of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                0.0, 0.0, 0.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_02");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_03");
            Console.WriteLine("  Integration region is the unit tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a good sampler.");

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("     P_NUM      X^2             X*Y             X*Z" + 
                              "             Y^2             Y*Z             Z^2");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_02, Integrand.tetrahedron_integrand_03, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;
            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 uses TETRAHEDRON_SAMPLE_03 with an increasing number of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                0.0, 0.0, 0.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_03");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_03");
            Console.WriteLine("  Integration region is the unit tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a good sampler.");

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("     P_NUM      X^2             X*Y             X*Z" + 
                              "             Y^2             Y*Z             Z^2");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_03, Integrand.tetrahedron_integrand_03, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;
            }
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 uses TETRAHEDRON_SAMPLE_04 with an increasing number of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                0.0, 0.0, 0.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_04");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_03");
            Console.WriteLine("  Integration region is the unit tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a good sampler.");

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("     P_NUM      X^2             X*Y             X*Z" + 
                              "             Y^2             Y*Z             Z^2");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_04, Integrand.tetrahedron_integrand_03, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;
            }
        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 uses TETRAHEDRON_SAMPLE_01 on a general tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 2.0, 3.0,
                4.0, 1.0, 2.0,
                2.0, 4.0, 4.0,
                3.0, 2.0, 5.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_01");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER");
            Console.WriteLine("  Integration region is over a general tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a bad sampler.");

            seed = 123456789;

            typeMethods.r8mat_transpose_print(3, 4, t, "  Tetrahedron vertices:");

            Console.WriteLine("");
            Console.WriteLine("     P_NUM");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_01, tetrahedron_integrand_user, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;

            }
        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 uses TETRAHEDRON_SAMPLE_02 on a general tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 2.0, 3.0,
                4.0, 1.0, 2.0,
                2.0, 4.0, 4.0,
                3.0, 2.0, 5.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_02");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER");
            Console.WriteLine("  Integration region is over a general tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a good sampler.");

            seed = 123456789;

            typeMethods.r8mat_transpose_print(3, 4, t, "  Tetrahedron vertices:");

            Console.WriteLine("");
            Console.WriteLine("     P_NUM");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_02, tetrahedron_integrand_user, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;

            }
        }

        static void test07()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST07 uses TETRAHEDRON_SAMPLE_03 on a general tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 2.0, 3.0,
                4.0, 1.0, 2.0,
                2.0, 4.0, 4.0,
                3.0, 2.0, 5.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_03");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER");
            Console.WriteLine("  Integration region is over a general tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a good sampler.");

            seed = 123456789;

            typeMethods.r8mat_transpose_print(3, 4, t, "  Tetrahedron vertices:");

            Console.WriteLine("");
            Console.WriteLine("     P_NUM");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_03, tetrahedron_integrand_user, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;

            }
        }

        static void test08()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST08 uses TETRAHEDRON_SAMPLE_04 on a general tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int f_num = 6;
            int i;
            int p_num;
            double[] result;
            int seed;
            double[] t =
            {
                1.0, 2.0, 3.0,
                4.0, 1.0, 2.0,
                2.0, 4.0, 4.0,
                3.0, 2.0, 5.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  Sample using TETRAHEDRON_UNIT_SAMPLE_04");
            Console.WriteLine("  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER");
            Console.WriteLine("  Integration region is over a general tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Use an increasing number of points P_NUM.");
            Console.WriteLine("  Note that the sample routine is a good sampler.");

            seed = 123456789;

            typeMethods.r8mat_transpose_print(3, 4, t, "  Tetrahedron vertices:");

            Console.WriteLine("");
            Console.WriteLine("     P_NUM");
            Console.WriteLine("");

            p_num = 1;

            while (p_num <= 65536)
            {
                MonteCarlo.TetrahedronSampleResult tmp = MonteCarlo.tetrahedron_monte_carlo(t, p_num, f_num,
                    MonteCarlo.tetrahedron_unit_sample_04, tetrahedron_integrand_user, ref seed);
                result = tmp.result;
                seed = tmp.seed;

                string cout = "  " + p_num.ToString().PadLeft(8);
                for (i = 0; i < f_num; i++)
                {
                    cout += "  " + result[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                p_num = 2 * p_num;

            }
        }

        static double[] tetrahedron_integrand_user(int p_num, double[] p, int f_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_INTEGRAND_USER evaluates 6 integrand functions defined by the user.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int P_NUM, the number of points.
            //
            //    Input, double P[3*P_NUM], the evaluation points.
            //
            //    Input, int F_NUM, the number of integrands.
            //
            //    Output, double TETRAHEDRON_INTEGRAND_USER[F_NUM*P_NUM], the integrand values.
            //
        {
            int j;
            double[] fp;

            fp = new double[f_num * p_num];

            for (j = 0; j < p_num; j++)
            {
                fp[0 + j * f_num] = 1.0;
                fp[1 + j * f_num] = p[0 + j * 3];
                fp[2 + j * f_num] = Math.Pow(p[1 + j * 3], 2);
                fp[3 + j * f_num] = Math.Pow(p[2 + j * 3], 3);
                fp[4 + j * f_num] = p[0 + j * 3] * p[1 + j * 3] * Math.Pow(p[2 + j * 3], 2);
                fp[5 + j * f_num] = Math.Pow(p[0 + j * 3], 2) * Math.Pow(p[1 + j * 3], 2) * Math.Pow(p[2 + j * 3], 2);
            }

            return fp;
        }
    }
}