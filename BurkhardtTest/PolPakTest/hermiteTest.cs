using System;
using Burkardt.PolynomialNS;

namespace PolPakTest
{
    public static class hermiteTest
    {
        public static void gen_hermite_poly_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEN_HERMITE_POLY_TEST tests GEN_HERMITE_POLY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;
            int N_TEST = 6;

            double[] c = new double[N + 1];
            int i;
            int j;
            double mu;
            double[] mu_test = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
            double x;
            double[] x_test = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

            Console.WriteLine("");
            Console.WriteLine("GEN_HERMITE_POLY_TEST");
            Console.WriteLine("  GEN_HERMITE_POLY evaluates the generalized Hermite");
            Console.WriteLine("  polynomial.");

            for (i = 0; i < N_TEST; i++)
            {

                x = x_test[i];
                mu = mu_test[i];

                Console.WriteLine("");
                Console.WriteLine("  Table of H(N,MU)(X) for");
                Console.WriteLine("");
                Console.WriteLine("    N(max) = " + N + "");
                Console.WriteLine("    MU =     " + mu + "");
                Console.WriteLine("    X =      " + x + "");
                Console.WriteLine("");

                Hermite.gen_hermite_poly(N, x, mu, ref c);

                for (j = 0; j <= N; j++)
                {
                    Console.WriteLine("  "
                                      + j.ToString().PadLeft(6) + "  "
                                      + c[j].ToString().PadLeft(14) + "");
                }
            }

        }

        public static void hermite_poly_phys_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLY_PHYS_TEST tests HERMITE_POLY_PHYS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 12;

            double fx = 0;
            double[] fx2 = new double[N_MAX + 1];
            int n = 0;
            int n_data = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("HERMITE_POLY_PHYS_TEST:");
            Console.WriteLine("  HERMITE_POLY_PHYS evaluates the physicist's Hermite polynomial.");
            Console.WriteLine("");
            Console.WriteLine("     N      X        Exact F       H(N)(X)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Hermite.hermite_poly_phys_values(ref n_data, ref n, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                Hermite.hermite_poly_phys(n, x, ref fx2);

                Console.WriteLine("  "
                     + n.ToString().PadLeft(8)+ "  "
                     + x.ToString().PadLeft(8) + "  "
                     + fx.ToString().PadLeft(14) + "  "
                     + fx2[n].ToString().PadLeft(14) + "");
            }

        }

        public static void hermite_poly_phys_coef_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLY_PHYS_COEF_TEST tests HERMITE_POLY_PHYS_COEF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            double[] c = new double[(N + 1) * (N + 1)];
            int i;
            int j;

            Console.WriteLine("");
            Console.WriteLine("HERMITE_POLY_PHYS_COEF_TEST");
            Console.WriteLine("  HERMITE_POLY_PHYS_COEF: physicist's Hermite polynomial coefficients.");

            Hermite.hermite_poly_phys_coef(N, ref c);

            for (i = 0; i <= N; i++)
            {
                Console.WriteLine("");
                Console.WriteLine("  H(" + i + ")");
                Console.WriteLine("");
                for (j = i; 0 <= j; j--)
                {
                    if (j == 0)
                    {
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + "");
                        ;
                    }
                    else if (j == 1)
                    {
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x");
                    }
                    else
                    {
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                    }
                }
            }

        }

    }
}