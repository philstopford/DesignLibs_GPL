using System;
using Burkardt.PolynomialNS;

namespace LegendreProductPolynomialTest
{
    public static class lpTest
    {
        public static void lp_coefficients_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LP_COEFFICIENTS_TEST tests LP_COEFFICIENTS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c;
            int[] e;
            int[] f;
            int i;
            string label;
            int m = 1;
            int n;
            int n_max = 10;
            int o = 0;

            Console.WriteLine("");
            Console.WriteLine("LP_COEFFICIENTS_TEST");
            Console.WriteLine("  LP_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).");
            Console.WriteLine("");

            for (n = 0; n <= n_max; n++)
            {
                c = new double[n + 1];
                f = new int[n + 1];

                Legendre.lp_coefficients(n, ref o, ref c, ref f);

                e = new int[o];
                for (i = 0; i < o; i++)
                {
                    e[i] = f[i] + 1;
                }

                label = "  P(" + n + ",x) = ";
                Polynomial.polynomial_print(m, o, c, e, label);
            }
        }

        public static void lp_value_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LP_VALUE_TEST tests LP_VALUE.
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
            Console.WriteLine("LP_VALUE_TEST:");
            Console.WriteLine("  LP_VALUE evaluates a Legendre polynomial.");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated                 Computed");
            Console.WriteLine("     O        X           L(O,X)                    L(O,X)" +
                              "                   Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Legendre.lp_values(ref n_data, ref o, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                xvec[0] = x;

                fx2 = Legendre.lp_value(n, o, xvec);

                e = fx1 - fx2[0];

                Console.WriteLine(o.ToString().PadLeft(6) + "  "
                                                          + x.ToString().PadLeft(12) + "  "
                                                          + fx1.ToString().PadLeft(24) + "  "
                                                          + fx2[0].ToString().PadLeft(24) + "  "
                                                          + e.ToString().PadLeft(8) + "");
            }

        }

        public static void lp_values_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LP_VALUES_TEST tests LP_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data;
            int o = 0;
            double x = 0;
            double fx = 0;

            Console.WriteLine("");
            Console.WriteLine("LP_VALUES_TEST:");
            Console.WriteLine("  LP_VALUES stores values of");
            Console.WriteLine("  the Legendre polynomial P(o,x).");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated");
            Console.WriteLine("     O        X           L(O,X)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Legendre.lp_values(ref n_data, ref o, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine(o.ToString().PadLeft(6) + "  "
                                                          + x.ToString().PadLeft(12) + "  "
                                                          + fx.ToString().PadLeft(24) + "");
            }
        }
    }
}