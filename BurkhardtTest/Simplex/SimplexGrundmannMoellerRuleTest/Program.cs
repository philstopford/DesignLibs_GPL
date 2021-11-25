using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.SimplexNS;
using Burkardt.Types;

namespace SimplexGrundmannMoellerRuleTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SIMPLEX_GM_RULE_TEST.
        //
        //  Discussion:
        //
        //    SIMPLEX_GM_RULE_TEST tests the SIMPLEX_GM_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_GM_RULE_TEST");
        Console.WriteLine("  Test the SIMPLEX_GM_RULE library.");

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

        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_GM_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SIMPLEX_UNIT_TO_GENERAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int m = 2;
        int n = 10;
        int seed = 123456789;
        double[] t =
        {
            1.0, 1.0,
            3.0, 1.0,
            2.0, 5.0
        };
        double[] t_unit =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0
        };
        string cout = "";
        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  SIMPLEX_UNIT_TO_GENERAL");
        Console.WriteLine("  maps points in the unit simplex to a general simplex.");
        Console.WriteLine("");
        Console.WriteLine("  Here we consider a simplex in 2D, a triangle.");
        Console.WriteLine("");
        Console.WriteLine("  The vertices of the general triangle are:");
        Console.WriteLine("");
        for (j = 0; j < m + 1; j++)
        {
            cout = "";
            for (i = 0; i < m; i++)
            {
                cout += "  " + t[i + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("   (  XSI     ETA )   ( X       Y  )");
        Console.WriteLine("");

        double[] phy_unit = new double[m * (m + 1)];

        Simplex.simplex_unit_to_general(m, m + 1, t, t_unit, ref phy_unit);

        for (j = 0; j < m + 1; j++)
        {
            cout = "";
            for (i = 0; i < m; i++)
            {
                cout += "  " + t_unit[i + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            for (i = 0; i < m; i++)
            {
                cout += "  " + phy_unit[i + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            Console.WriteLine(cout);
        }

        double[] ref_ = Simplex.simplex_unit_sample(m, n, ref seed);
        double[] phy = new double[m * n];

        Simplex.simplex_unit_to_general(m, n, t, ref_, ref phy);

        for (j = 0; j < n; j++)
        {
            cout = "";
            for (i = 0; i < m; i++)
            {
                cout += "  " + ref_[i + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            for (i = 0; i < m; i++)
            {
                cout += "  " + phy[i + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests SIMPLEX_UNIT_TO_GENERAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim;
        int j;
        const int m = 3;
        const int n = 10;
        int seed = 123456789;
        double[] t =
        {
            1.0, 1.0, 1.0,
            3.0, 1.0, 1.0,
            1.0, 4.0, 1.0,
            1.0, 1.0, 5.0
        };
        double[] t_unit =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };
        const int vertex_num = 3 + 1;
        string cout = "";
        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  SIMPLEX_UNIT_TO_GENERAL");
        Console.WriteLine("  maps points in the unit simplex to a general simplex.");
        Console.WriteLine("");
        Console.WriteLine("  Here we consider a simplex in 3D, a tetrahedron.");
        Console.WriteLine("");
        Console.WriteLine("  The vertices of the general tetrahedron are:");
        Console.WriteLine("");
        for (j = 0; j < vertex_num; j++)
        {
            cout = "";
            for (dim = 0; dim < m; dim++)
            {
                cout += "  " + t[dim + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("   (  XSI     ETA     MU )    ( X       Y       Z )");
        Console.WriteLine("");

        double[] phy_unit = new double[m * (m + 1)];

        Simplex.simplex_unit_to_general(m, m + 1, t, t_unit, ref phy_unit);

        for (j = 0; j < m + 1; j++)
        {
            cout = "";
            for (dim = 0; dim < m; dim++)
            {
                cout += "  " + t_unit[dim + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            for (dim = 0; dim < m; dim++)
            {
                cout += "  " + phy_unit[dim + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            Console.WriteLine(cout);
        }

        double[] ref_ = Simplex.simplex_unit_sample(m, n, ref seed);
        double[] phy = new double[m * n];

        Simplex.simplex_unit_to_general(m, n, t, ref_, ref phy);

        for (j = 0; j < n; j++)
        {
            cout = "";
            for (dim = 0; dim < m; dim++)
            {
                cout += "  " + ref_[dim + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            for (dim = 0; dim < m; dim++)
            {
                cout += "  " + phy[dim + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(9);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests GM_RULE_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 4;

        int[] m_test = { 2, 3, 5, 10 };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  GM_RULE_SIZE returns N, the number of points");
        Console.WriteLine("  associated with a Grundmann-Moeller quadrature rule");
        Console.WriteLine("  for the unit simplex of dimension M");
        Console.WriteLine("  with rule index RULE");
        Console.WriteLine("  and degree of exactness DEGREE = 2*RULE+1.");

        Console.WriteLine("");
        Console.WriteLine("   M      RULE    DEGREE N");

        for (test = 0; test < TEST_NUM; test++)
        {
            int m = m_test[test];

            Console.WriteLine("");

            int rule;
            for (rule = 0; rule <= 5; rule++)
            {
                int n = GrundmannMoellerRule.gm_rule_size(rule, m);
                int degree = 2 * rule + 1;

                Console.WriteLine("  " + m.ToString().PadLeft(8)
                                       + "  " + rule.ToString().PadLeft(8)
                                       + "  " + degree.ToString().PadLeft(8)
                                       + "  " + n.ToString().PadLeft(8) + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests GM_UNIT_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int point;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  GM_UNIT_RULE_SET determines the weights and abscissas");
        Console.WriteLine("  of a Grundmann-Moeller quadrature rule for");
        Console.WriteLine("  the M dimensional unit simplex,");
        Console.WriteLine("  using a rule of index RULE,");
        Console.WriteLine("  which will have degree of exactness 2*RULE+1.");

        const int m = 3;
        const int rule = 2;

        Console.WriteLine("");
        Console.WriteLine("  Here we use M = " + m + "");
        Console.WriteLine("  RULE = " + rule + "");
        Console.WriteLine("  DEGREE = " + 2 * rule + 1 + "");

        int n = GrundmannMoellerRule.gm_rule_size(rule, m);

        double[] w = new double[n];
        double[] x = new double[m * n];

        GrundmannMoellerRule.gm_unit_rule_set(rule, m, n, ref w, ref x);

        Console.WriteLine("");
        Console.WriteLine("     POINT        W             X             Y             Z");
        Console.WriteLine("");

        for (point = 0; point < n; point++)
        {
            string cout = "  " + (point + 1).ToString().PadLeft(8)
                               + "  " + w[point].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            int dim;
            for (dim = 0; dim < m; dim++)
            {
                cout += "  " + x[dim + point * m].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests GM_UNIT_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 4;

        int[] m_test = { 2, 3, 5, 10 };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  GM_UNIT_RULE_SET determines the weights and abscissas");
        Console.WriteLine("  of a Grundmann-Moeller quadrature rule for");
        Console.WriteLine("  the M dimensional unit simplex,");
        Console.WriteLine("  using a rule of index RULE,");
        Console.WriteLine("  which will have degree of exactness 2*RULE+1.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we compute various rules, and simply");
        Console.WriteLine("  report the number of points, and the sum of weights.");

        Console.WriteLine("");
        Console.WriteLine("   M      RULE    N  WEIGHT SUM");

        for (test = 0; test < TEST_NUM; test++)
        {
            int m = m_test[test];

            Console.WriteLine("");

            int rule;
            for (rule = 0; rule <= 5; rule++)
            {
                int n = GrundmannMoellerRule.gm_rule_size(rule, m);

                double[] w = new double[n];
                double[] x = new double[m * n];

                GrundmannMoellerRule.gm_unit_rule_set(rule, m, n, ref w, ref x);

                double w_sum = 0.0;
                int point;
                for (point = 0; point < n; point++)
                {
                    w_sum += w[point];
                }

                Console.WriteLine("  " + m.ToString().PadLeft(8)
                                       + "  " + rule.ToString().PadLeft(8)
                                       + "  " + n.ToString().PadLeft(8)
                                       + "  " + w_sum.ToString("0.################").PadLeft(24) + "");
            }
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests GM_UNIT_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int point;
        List<string> w_unit = new();
        List<string> x_unit = new();

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  GM_UNIT_RULE_SET determines the weights and abscissas");
        Console.WriteLine("  of a Grundmann-Moeller quadrature rule for");
        Console.WriteLine("  the M dimensional unit simplex,");
        Console.WriteLine("  using a rule of index RULE,");
        Console.WriteLine("  which will have degree of exactness 2*RULE+1.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we write a rule to a file.");

        const int m = 3;
        const int rule = 2;

        Console.WriteLine("");
        Console.WriteLine("  Here we use M = " + m + "");
        Console.WriteLine("  RULE = " + rule + "");
        Console.WriteLine("  DEGREE = " + 2 * rule + 1 + "");

        int n = GrundmannMoellerRule.gm_rule_size(rule, m);

        double[] w = new double[n];
        double[] x = new double[m * n];

        GrundmannMoellerRule.gm_unit_rule_set(rule, m, n, ref w, ref x);

        string w_file = "gm" + rule + "_" + m + "d_w.txt";

        for (point = 0; point < n; point++)
        {
            w_unit.Add(w[point].ToString("0.################").PadLeft(20) + "");
        }

        File.WriteAllLines(w_file, w_unit);

        string x_file = "gm" + rule + "_" + m + "d_x.txt";

        for (point = 0; point < n; point++)
        {
            string tmp = "";
            int dim;
            for (dim = 0; dim < m; dim++)
            {
                tmp += x[dim + point * m].ToString("0.################").PadLeft(20);
            }

            x_unit.Add(tmp);
        }

        File.WriteAllLines(x_file, x_unit);

        Console.WriteLine("");
        Console.WriteLine("  Wrote rule " + rule
                                          + " to \"" + w_file
                                          + "\" and \"" + x_file + "");

    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests GM_UNIT_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int degree;
        const int degree_max = 4;
        int h = 0;
        const int m = 5;
        const int rule_max = 3;
        int t = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  GM_UNIT_RULE_SET determines the weights and abscissas");
        Console.WriteLine("  of a Grundmann-Moeller quadrature rule for");
        Console.WriteLine("  the M dimensional unit simplex,");
        Console.WriteLine("  using a rule of index RULE,");
        Console.WriteLine("  which will have degree of exactness 2*RULE+1.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, look at all the monomials up to");
        Console.WriteLine("  some maximum degree, choose a few low order rules");
        Console.WriteLine("  and determine the quadrature error for each.");
        Console.WriteLine("");
        Console.WriteLine("  Here we use M = " + m + "");

        Console.WriteLine("");
        Console.WriteLine("      Rule     Order     Quad_Error");
        Console.WriteLine("");

        int[] expon = new int[m];

        for (degree = 0; degree <= degree_max; degree++)
        {
            bool more = false;

            for (;;)
            {
                Comp.comp_next(degree, m, ref expon, ref more, ref h, ref t);

                Console.WriteLine("");
                Console.WriteLine("  F(X) = X1^" + expon[0]
                                                 + " * X2^" + expon[1]
                                                 + " * X3^" + expon[2]
                                                 + " * X4^" + expon[3]
                                                 + " * X5^" + expon[4] + "");
                Console.WriteLine("");

                int rule;
                for (rule = 0; rule <= rule_max; rule++)
                {
                    int n = GrundmannMoellerRule.gm_rule_size(rule, m);

                    double[] w = new double[n];
                    double[] x = new double[m * n];

                    GrundmannMoellerRule.gm_unit_rule_set(rule, m, n, ref w, ref x);

                    double quad_error = Simplex.simplex_unit_monomial_quadrature(m, expon,
                        n, x, w);

                    Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + quad_error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }

                if (!more)
                {
                    break;
                }
            }
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests GM_GENERAL_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        double[] t =
        {
            1.0, 0.0, 0.0,
            2.0, 0.0, 0.0,
            1.0, 2.0, 0.0,
            1.0, 0.0, 3.0
        };
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  GM_GENERAL_RULE_SET determines the weights and abscissas");
        Console.WriteLine("  of a Grundmann-Moeller quadrature rule for");
        Console.WriteLine("  the M dimensional general simplex,");
        Console.WriteLine("  using a rule of index RULE,");
        Console.WriteLine("  which will have degree of exactness 2*RULE+1.");

        const int m = 3;
        const int rule = 2;

        Console.WriteLine("");
        Console.WriteLine("  Here we use M = " + m + "");
        Console.WriteLine("  RULE = " + rule + "");
        Console.WriteLine("  DEGREE = " + 2 * rule + 1 + "");

        Console.WriteLine("");
        Console.WriteLine("  Simplex vertices:");
        Console.WriteLine("");
        for (j = 0; j < 4; j++)
        {
            cout = "";
            for (i = 0; i < 3; i++)
            {
                cout += "  " + t[i + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        int n = GrundmannMoellerRule.gm_rule_size(rule, m);

        double[] w = new double[n];
        double[] x = new double[m * n];

        GrundmannMoellerRule.gm_general_rule_set(rule, m, n, t, ref w, ref x);

        Console.WriteLine("");
        Console.WriteLine("     POINT        W             X             Y             Z");
        Console.WriteLine("");

        for (j = 0; j < n; j++)
        {
            cout = "  " + j.ToString().PadLeft(8)
                        + "  " + w[j].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            for (i = 0; i < m; i++)
            {
                cout += "  " + x[i + j * m].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    private static void test09()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    TEST09 tests GM_UNIT_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        int[] e_test =
        {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            2, 0, 0,
            1, 1, 0,
            1, 0, 1,
            0, 2, 0,
            0, 1, 1,
            0, 0, 2
        };
        const int m = 3;
        int rule;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  GM_UNIT_RULE_SET determines the weights and abscissas");
        Console.WriteLine("  of a Grundmann-Moeller quadrature rule for");
        Console.WriteLine("  the M dimensional unit simplex,");
        Console.WriteLine("  using a rule of index RULE,");
        Console.WriteLine("  which will have degree of exactness 2*RULE+1.");

        Console.WriteLine("");
        Console.WriteLine("  In this test, look at all the monomials up to");
        Console.WriteLine("  some maximum degree, choose a few low order rules");
        Console.WriteLine("  and determine the quadrature error for each.");

        double volume = Simplex.simplex_unit_volume(m);
        Console.WriteLine("");
        Console.WriteLine("  Simplex volume = " + volume + "");

        Console.WriteLine("");
        Console.WriteLine("  Here we use M = " + m + "");

        Console.WriteLine("");
        Console.WriteLine("         N        1               X               Y " +
                          "              Z               X^2              XY             XZ" +
                          "              Y^2             YZ               Z^2");
        Console.WriteLine("");

        for (rule = 0; rule <= 5; rule++)
        {
            int n = GrundmannMoellerRule.gm_rule_size(rule, m);

            double[] w = new double[n];
            double[] x = new double[m * n];

            GrundmannMoellerRule.gm_unit_rule_set(rule, m, n, ref w, ref x);

            string cout = "  " + n.ToString().PadLeft(8);

            int k;
            for (k = 0; k < 10; k++)
            {
                int i;
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + k * m];
                }

                double[] value = Monomial.monomial_value(m, n, e, x);

                double result = typeMethods.r8vec_dot_product(n, w, value);

                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

            }

            Console.WriteLine(cout);

        }
    }

    private static void test10()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    TEST10 tests GM_GENERAL_RULE_SET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        int[] e_test =
        {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            2, 0, 0,
            1, 1, 0,
            1, 0, 1,
            0, 2, 0,
            0, 1, 1,
            0, 0, 2
        };
        int i;
        int j;
        const int m = 3;
        int rule;
        double[] t =
        {
            1.0, 0.0, 0.0,
            2.0, 0.0, 0.0,
            1.0, 2.0, 0.0,
            1.0, 0.0, 3.0
        };
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  GM_GENERAL_RULE_SET determines the weights and abscissas");
        Console.WriteLine("  of a Grundmann-Moeller quadrature rule for");
        Console.WriteLine("  the M dimensional general simplex,");
        Console.WriteLine("  using a rule of index RULE,");
        Console.WriteLine("  which will have degree of exactness 2*RULE+1.");

        Console.WriteLine("");
        Console.WriteLine("  In this test, look at all the monomials up to");
        Console.WriteLine("  some maximum degree, choose a few low order rules");
        Console.WriteLine("  and determine the quadrature error for each.");

        Console.WriteLine("");
        Console.WriteLine("  Simplex vertices:");
        Console.WriteLine("");
        for (j = 0; j < 4; j++)
        {
            cout = "";
            for (i = 0; i < 3; i++)
            {
                cout += "  " + t[i + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        double volume = Simplex.simplex_general_volume(m, t);
        Console.WriteLine("");
        Console.WriteLine("  Simplex volume = " + volume + "");

        Console.WriteLine("");
        Console.WriteLine("  Here we use M = " + m + "");

        Console.WriteLine("");
        Console.WriteLine("         N        1               X               Y " +
                          "              Z               X^2              XY             XZ" +
                          "              Y^2             YZ               Z^2");
        Console.WriteLine("");

        for (rule = 0; rule <= 5; rule++)
        {
            int n = GrundmannMoellerRule.gm_rule_size(rule, m);

            double[] w = new double[n];
            double[] x = new double[m * n];

            GrundmannMoellerRule.gm_general_rule_set(rule, m, n, t, ref w, ref x);

            cout = "  " + n.ToString().PadLeft(8);

            int k;
            for (k = 0; k < 10; k++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + k * m];
                }

                double[] value = Monomial.monomial_value(m, n, e, x);

                double result = typeMethods.r8vec_dot_product(n, w, value);

                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

            }

            Console.WriteLine(cout);
        }
    }
}