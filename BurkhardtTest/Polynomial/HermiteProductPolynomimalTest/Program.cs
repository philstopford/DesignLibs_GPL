using System;
using Burkardt;
using Burkardt.PolynomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HermiteProductPolynomimalTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for HERMITE_PRODUCT_POLYNOMIAL_TEST.
            //
            //  Discussion:
            //
            //    HERMITE_PRODUCT_POLYNOMIAL_TEST tests the HERMITE_PRODUCT_POLYNOMIAL library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("HERMITE_PRODUCT_POLYNOMIAL_TEST:");
            Console.WriteLine("  Test the HERMITE_PRODUCT_POLYNOMIAL library.");

            hpp_test01();
            hpp_test015();
            hpp_test02();
            hpp_test03();
            hpp_test04();

            Console.WriteLine("");
            Console.WriteLine("HERMITE_PRODUCT_POLYNOMIAL_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void hpp_test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HPP_TEST01 tests routines for the GRLEX ordering of compositions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int k = 2;
            int rank;
            int rank1;
            int rank2;
            int seed;
            int test;
            int[] x;
            int x_sum;
            int x_sum_old;

            string cout = "";

            x = new int[k];

            Console.WriteLine("");
            Console.WriteLine("HPP_TEST01:");
            Console.WriteLine("  COMP_NEXT_GRLEX is given a composition, and computes the ");
            Console.WriteLine("  next composition in grlex order.");

            Console.WriteLine("");
            Console.WriteLine("  Rank   Sum   Components");

            for (i = 0; i < k; i++)
            {
                x[i] = 0;
            }

            x_sum_old = -1;
            rank = 1;

            for (;;)
            {
                x_sum = typeMethods.i4vec_sum(k, x);

                if (x_sum_old < x_sum)
                {
                    x_sum_old = x_sum;
                    Console.WriteLine("");
                }

                cout = rank.ToString().PadLeft(6) + "  "
                    + x_sum.ToString().PadLeft(6);
                for (i = 0; i < k; i++)
                {
                    cout += x[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);

                if (20 <= rank)
                {
                    break;
                }

                Comp.comp_next_grlex(k, ref x);
                rank = rank + 1;
            }

            Console.WriteLine("");
            Console.WriteLine("  COMP_UNRANK_GRLEX is given a rank and returns the");
            Console.WriteLine("  corresponding set of multinomial exponents.");
            Console.WriteLine("");
            Console.WriteLine("  Rank   Sum   Components");
            Console.WriteLine("");

            seed = 123456789;

            for (test = 1; test <= 5; test++)
            {
                rank = UniformRNG.i4_uniform_ab(1, 20, ref seed);
                x = Comp.comp_unrank_grlex(k, rank);
                x_sum = typeMethods.i4vec_sum(k, x);
                cout = rank.ToString().PadLeft(6) + "  "
                    + x_sum.ToString().PadLeft(6);
                for (i = 0; i < k; i++)
                {
                    cout += x[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  COMP_RANDOM_GRLEX randomly selects a composition");
            Console.WriteLine("  between given lower and upper ranks.");
            Console.WriteLine("");
            Console.WriteLine("  Rank   Sum   Components");
            Console.WriteLine("");

            seed = 123456789;
            rank1 = 5;
            rank2 = 20;

            for (test = 1; test <= 5; test++)
            {
                x = Comp.comp_random_grlex(k, rank1, rank2, ref seed, ref rank);
                x_sum = typeMethods.i4vec_sum(k, x);
                cout = rank.ToString().PadLeft(6) + "  "
                    + x_sum.ToString().PadLeft(6);
                for (i = 0; i < k; i++)
                {
                    cout += x[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  COMP_RANK_GRLEX returns the rank of a given composition.");
            Console.WriteLine("");
            Console.WriteLine("  Rank   Sum   Components");
            Console.WriteLine("");

            x = new int[k];

            x[0] = 4;
            x[1] = 0;
            rank = Comp.comp_rank_grlex(k, x);
            x_sum = typeMethods.i4vec_sum(k, x);
            cout = rank.ToString().PadLeft(6) + "  "
                + x_sum.ToString().PadLeft(6);
            for (i = 0; i < k; i++)
            {
                cout += x[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            ;

            x[0] = 11;
            x[1] = 5;
            rank = Comp.comp_rank_grlex(k, x);
            x_sum = typeMethods.i4vec_sum(k, x);
            cout = rank.ToString().PadLeft(6) + "  "
                + x_sum.ToString().PadLeft(6);
            for (i = 0; i < k; i++)
            {
                cout += x[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }

        static void hpp_test015()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HPP_TEST015 tests HEP_COEFFICIENTS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c;
            int[] e;
            int[] f;
            int[] l = new int[1];
            int m;
            int n;
            int o;
            int o_max;
            string title;

            m = 1;

            Console.WriteLine("");
            Console.WriteLine("HPP_TEST015:");
            Console.WriteLine("  HEP_COEFFICIENTS computes the coefficients and");
            Console.WriteLine("  exponents of the Hermite polynomial He(n,x).");

            for (n = 1; n <= 5; n++)
            {
                o = (n + 2) / 2;
                c = new double[o];
                e = new int[o];
                f = new int[o];

                Hermite.hep_coefficients(n, ref o, ref c, ref f);

                l[0] = n;
                o_max = o;

                Hermite.hepp_to_polynomial(m, l, o_max, o, ref c, ref e);

                Console.WriteLine("");
                title = "  He(" + n + ",x) =";
                Polynomial.polynomial_print(m, o, c, e, title);
            }
        }

        static void hpp_test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HPP_TEST02 tests HEP_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double e;
            int n;
            int n_data;
            int o = 0;
            double x = 0;
            double[] xvec = new double[1];
            double fx1 = 0;
            double[] fx2;

            n = 1;

            Console.WriteLine("");
            Console.WriteLine("HPP_TEST02:");
            Console.WriteLine("  HEP_VALUES stores values of");
            Console.WriteLine("  the Hermite polynomial He(o,x).");
            Console.WriteLine("  HEP_VALUE evaluates a Hermite polynomial.");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated                 Computed");
            Console.WriteLine("     O        X          He(O,X)                   He(O,X)" +
                              "                   Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.TestValues.Hermite.hep_values(ref n_data, ref o, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                xvec[0] = x;

                fx2 = Hermite.hep_value(n, o, ref xvec);

                e = fx1 - fx2[0];

                Console.WriteLine(o.ToString().PadLeft(6) + "  "
                    + x.ToString().PadLeft(12) + "  "
                    + fx1.ToString().PadLeft(24) + "  "
                    + fx2[0].ToString().PadLeft(24) + "  "
                    + e.ToString().PadLeft(8) + "");
            }
        }

        static void hpp_test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HPP_TEST03 tests HEPP_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c;
            int[] e;
            int i;
            int[] l;
            int m = 3;
            int n = 1;
            int o = 0;
            int o_max;
            int rank;
            int seed;
            double[] v1;
            double[] v2;
            double[] x;
            double xhi;
            double xlo;
            string cout = "";

            Console.WriteLine("");
            Console.WriteLine("HPP_TEST03:");
            Console.WriteLine("  HEPP_VALUE evaluates a Hermite product polynomial.");
            Console.WriteLine("  POLYNOMIAL_VALUE evaluates a polynomial.");

            xlo = -1.0;
            xhi = +1.0;
            seed = 123456789;
            x = UniformRNG.r8vec_uniform_ab_new(m, xlo, xhi, ref seed);

            Console.WriteLine("");
            cout = "  Evaluate at X = ";
            for (i = 0; i < m; i++)
            {
                cout += "  " + x[i + 0 * m];
            }

            Console.WriteLine(cout);
            Console.WriteLine("");
            Console.WriteLine("  Rank  I1  I2  I3:  He(I1,X1)*He(I2,X2)*He(I3,X3)    P(X1,X2,X3)");
            Console.WriteLine("");

            for (rank = 1; rank <= 20; rank++)
            {
                l = Comp.comp_unrank_grlex(m, rank);
                //
                //  Evaluate the HePP directly.
                //
                v1 = Hermite.hepp_value(m, n, l, x);
                //
                //  Convert the HePP to a polynomial.
                //
                o_max = 1;
                for (i = 0; i < m; i++)
                {
                    o_max = o_max * (l[i] + 2) / 2;
                }

                c = new double[o_max];
                e = new int[o_max];

                Hermite.hepp_to_polynomial(m, l, o_max, o, ref c, ref e);
                //
                //  Evaluate the polynomial.
                //
                v2 = Polynomial.polynomial_value(m, o, c, e, n, x);
                //
                //  Compare results.
                //
                Console.WriteLine(rank.ToString().PadLeft(6) + "  "
                    + l[0].ToString().PadLeft(2) + "  "
                    + l[1].ToString().PadLeft(2) + "  "
                    + l[2].ToString().PadLeft(2) + "  "
                    + v1[0].ToString().PadLeft(14) + "  "
                    + v2[0].ToString().PadLeft(14) + "");
            }
        }

        static void hpp_test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HPP_TEST04 tests HEPP_TO_POLYNOMIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c;
            int[] e;
            int i;
            int[] l;
            string label;
            int m = 2;
            int o = 0;
            int o_max;
            int rank;

            Console.WriteLine("");
            Console.WriteLine("HPP_TEST04:");
            Console.WriteLine("  HEPP_TO_POLYNOMIAL is given a Hermite product polynomial");
            Console.WriteLine("  and determines its polynomial representation.");

            Console.WriteLine("");
            Console.WriteLine("  Using spatial dimension M = " + m + "");

            for (rank = 1; rank <= 11; rank++)
            {
                l = Comp.comp_unrank_grlex(m, rank);

                o_max = 1;
                for (i = 0; i < m; i++)
                {
                    o_max = o_max * (l[i] + 2) / 2;
                }

                c = new double[o_max];
                e = new int[o_max];

                Hermite.hepp_to_polynomial(m, l, o_max, o, ref c, ref e);

                label = "  HePP #" + rank
                    + " = He(" + l[0]
                    + ",X)*He(" + l[1]
                    + ",Y) =";

                Console.WriteLine("");
                Polynomial.polynomial_print(m, o, c, e, label);
            }
        }
    }
}