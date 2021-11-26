using System;
using Burkardt.MonomialNS;
using Burkardt.Uniform;

namespace LegendreProductPolynomialTest;

public static class monoTest
{
    public static void mono_next_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_NEXT_GRLEX_TEST tests MONO_NEXT_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 4;
        int i;

        Console.WriteLine("");
        Console.WriteLine("MONO_NEXT_GRLEX_TEST");
        Console.WriteLine("  MONO_NEXT_GRLEX computes the next monomial");
        Console.WriteLine("  in M variables, in grlex order.");
        Console.WriteLine("");
        Console.WriteLine("  Let M =  " + m + "");

        const int a = 0;
        const int b = 3;
        int seed = 123456789;

        for (i = 1; i <= 10; i++)
        {
            int[] x = UniformRNG.i4vec_uniform_ab_new(m, a, b, ref seed);
            Console.WriteLine("");
            string cout = "  ";
            int k;
            for (k = 0; k < m; k++)
            {
                cout += x[k].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);

            int j;
            for (j = 1; j <= 5; j++)
            {
                Monomial.mono_next_grlex(m, ref x);
                cout = "  ";
                for (k = 0; k < m; k++)
                {
                    cout += x[k].ToString().PadLeft(2);
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static void mono_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_PRINT_TEST tests MONO_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] f1 = {5};
        int[] f2 = {-5};
        int[] f3 = {2, 1, 0, 3};
        int[] f4 = {17, -3, 199};

        Console.WriteLine("");
        Console.WriteLine("MONO_PRINT_TEST");
        Console.WriteLine("  MONO_PRINT can print out a monomial.");
        Console.WriteLine("");

        int m = 1;
        Monomial.mono_print(m, f1, "  Monomial [5]:");

        m = 1;
        Monomial.mono_print(m, f2, "  Monomial [5]:");

        m = 4;
        Monomial.mono_print(m, f3, "  Monomial [2,1,0,3]:");

        m = 3;
        Monomial.mono_print(m, f4, "  Monomial [17,-3,199]:");

    }

    public static void mono_rank_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_RANK_GRLEX_TEST tests MONO_RANK_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 3;
        int j;
        int test;
        const int test_num = 8;
        int[] x = new int[3];
        int[] x_test =
        {
            0, 0, 0,
            1, 0, 0,
            0, 0, 1,
            0, 2, 0,
            1, 0, 2,
            0, 3, 1,
            3, 2, 1,
            5, 2, 1
        };

        Console.WriteLine("");
        Console.WriteLine("MONO_RANK_GRLEX_TEST");
        Console.WriteLine("  MONO_RANK_GRLEX returns the rank of a monomial in the sequence");
        Console.WriteLine("  of all monomials in M dimensions, in grlex order.");

        Console.WriteLine("");
        Console.WriteLine("  Print a monomial sequence with ranks assigned.");

        const int n = 4;

        Console.WriteLine("");
        Console.WriteLine("  Let M = " + m + "");
        Console.WriteLine("      N = " + n + "");
        Console.WriteLine("");

        x[0] = 0;
        x[1] = 0;
        x[2] = 0;

        int i = 1;

        for (;;)
        {
            string cout = "  " + i.ToString().PadLeft(3) + "    ";
            for (j = 0; j < m; j++)
            {
                cout += x[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);

            if (x[0] == n)
            {
                break;
            }

            Monomial.mono_upto_next_grlex(m, n, ref x);
            i += 1;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now, given a monomial, retrieve its rank in the sequence:");
        Console.WriteLine("");

        for (test = 0; test < test_num; test++)
        {
            for (j = 0; j < m; j++)
            {
                x[j] = x_test[j + test * m];
            }

            int rank = Monomial.mono_rank_grlex(m, x);

            string cout = "  " + rank.ToString().PadLeft(3) + "    ";
            for (j = 0; j < m; j++)
            {
                cout += x[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

    }

    public static void mono_unrank_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_UNRANK_GRLEX_TEST tests MONO_UNRANK_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 3;
        int i;
        int j;
        int test;

        Console.WriteLine("");
        Console.WriteLine("MONO_UNRANK_GRLEX_TEST");
        Console.WriteLine("  MONO_UNRANK_GRLEX is given a rank, and returns the corresponding");
        Console.WriteLine("  monomial in the sequence of all monomials in M dimensions");
        Console.WriteLine("  in grlex order.");

        Console.WriteLine("");
        Console.WriteLine("  For reference, print a monomial sequence with ranks.");

        const int n = 4;
        int rank_max = Monomial.mono_upto_enum(m, n);

        Console.WriteLine("");
        Console.WriteLine("  Let M = " + m + "");
        Console.WriteLine("      N = " + n + "");
        Console.WriteLine("");

        int[] x = new int[n];
        for (i = 0; i < m; i++)
        {
            x[i] = 0;
        }

        i = 1;

        for (;;)
        {
            string cout = "  " + i.ToString().PadLeft(3) + "    ";
            for (j = 0; j < m; j++)
            {
                cout += x[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);

            if (x[0] == n)
            {
                break;
            }

            Monomial.mono_upto_next_grlex(m, n, ref x);
            i += 1;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now choose random ranks between 1 and " + rank_max + "");
        Console.WriteLine("");

        int seed = 123456789;
        const int test_num = 5;

        for (test = 1; test <= test_num; test++)
        {
            int rank = UniformRNG.i4_uniform_ab(1, rank_max, ref seed);
            x = Monomial.mono_unrank_grlex(m, rank);
            string cout = "  " + rank.ToString().PadLeft(3) + "    ";
            for (j = 0; j < m; j++)
            {
                cout += x[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

    }

    public static void mono_upto_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_UPTO_ENUM_TEST tests MONO_UPTO_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;

        Console.WriteLine("");
        Console.WriteLine("MONO_UPTO_ENUM_TEST");
        Console.WriteLine("  MONO_UPTO_ENUM can enumerate the number of monomials");
        Console.WriteLine("  in M variables, of total degree 0 up to N.");

        Console.WriteLine("");
        Console.WriteLine("    N:");
        string cout = "";
        for (n = 0; n <= 8; n++)
        {
            cout += "  " + n.ToString().PadLeft(4);
        }

        Console.WriteLine(cout);
        cout = "";
        Console.WriteLine("   m +------------------------------------------------------");
        for (m = 1; m <= 8; m++)
        {
            cout += "  " + m.ToString().PadLeft(2) + "  |";
            for (n = 0; n <= 8; n++)
            {
                int v = Monomial.mono_upto_enum(m, n);
                cout += " " + v.ToString().PadLeft(5);
            }

            Console.WriteLine(cout);
            cout = "";
        }
    }

    public static void mono_upto_next_grlex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_UPTO_NEXT_GRLEX_TEST tests MONO_UPTO_NEXT_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 3;
        int[] x = new int[3];

        Console.WriteLine("");
        Console.WriteLine("MONO_UPTO_NEXT_GRLEX_TEST");
        Console.WriteLine("  MONO_UPTO_NEXT_GRLEX can list the monomials");
        Console.WriteLine("  in M variables, of total degree up to N,");
        Console.WriteLine("  in grlex order, one at a time.");
        Console.WriteLine("");
        Console.WriteLine("  We start the process with (0,0,...,0,0).");
        Console.WriteLine("  The process ends with (N,0,...,0,0)");

        const int n = 4;

        Console.WriteLine("");
        Console.WriteLine("  Let M = " + m + "");
        Console.WriteLine("      N = " + n + "");
        Console.WriteLine("");

        x[0] = 0;
        x[1] = 0;
        x[2] = 0;

        int i = 1;

        for (;;)
        {
            string cout = "  " + i.ToString().PadLeft(2) + "    ";
            int j;
            for (j = 0; j < m; j++)
            {
                cout += x[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);

            if (x[0] == n)
            {
                break;
            }

            Monomial.mono_upto_next_grlex(m, n, ref x);
            i += 1;
        }

    }

    public static void mono_upto_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_UPTO_RANDOM_TEST tests MONO_UPTO_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 3;
        int rank = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("MONO_UPTO_RANDOM_TEST");
        Console.WriteLine("  MONO_UPTO_RANDOM selects at random a monomial");
        Console.WriteLine("  in M dimensions of total degree no greater than N.");

        const int n = 4;

        Console.WriteLine("");
        Console.WriteLine("  Let M = " + m + "");
        Console.WriteLine("      N = " + n + "");
        Console.WriteLine("");

        int seed = 123456789;
        const int test_num = 5;

        for (test = 1; test <= test_num; test++)
        {
            int[] x = Monomial.mono_upto_random(m, n, ref seed, ref rank);
            string cout = "  " + rank.ToString().PadLeft(3) + "    ";
            int j;
            for (j = 0; j < m; j++)
            {
                cout += x[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

    }

    public static void mono_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONO_VALUE_TEST tests MONO_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 3;
        const int nx = 2;
        int rank = 0;
        int test;
        double[] x =
        {
            1.0, 2.0, 3.0,
            -2.0, 4.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("MONO_VALUE_TEST");
        Console.WriteLine("  MONO_VALUE evaluates a monomial.");

        const int n = 6;

        Console.WriteLine("");
        Console.WriteLine("  Let M = " + m + "");
        Console.WriteLine("      N = " + n + "");

        int seed = 123456789;
        const int test_num = 5;

        for (test = 1; test <= test_num; test++)
        {
            int[] f = Monomial.mono_upto_random(m, n, ref seed, ref rank);
            Console.WriteLine("");
            Monomial.mono_print(m, f, "  M(X) = ");
            double[] v = Monomial.mono_value(m, nx, f, x);
            int j;
            for (j = 0; j < nx; j++)
            {
                Console.WriteLine("  M(" + x[0 + j * m]
                                         + "," + x[1 + j * m]
                                         + "," + x[2 + j * m]
                                         + ") = " + v[j] + "");
            }
        }
    }
}