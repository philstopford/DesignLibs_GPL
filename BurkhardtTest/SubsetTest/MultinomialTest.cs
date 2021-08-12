using System;
using Burkardt.PolynomialNS;

namespace SubsetTest
{
    public static class MultinomialTest
    {
        public static void multinomial_coef1_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTINOMIAL_COEF1_TEST tests MULTINOMIAL_COEF1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int MAXFACTOR = 5;

            int[] factor = new int[MAXFACTOR];
            int i;
            int j;
            int n;
            int ncomb;
            int nfactor;

            Console.WriteLine("");
            Console.WriteLine("MULTINOMIAL_COEF1_TEST");
            Console.WriteLine("  MULTINOMIAL_COEF1 computes multinomial");
            Console.WriteLine("  coefficients using the Gamma function;");

            Console.WriteLine("");
            Console.WriteLine("  Line 10 of the BINOMIAL table:");
            Console.WriteLine("");

            n = 10;
            nfactor = 2;

            for (i = 0; i <= n; i++)
            {
                factor[0] = i;
                factor[1] = n - i;

                ncomb = Multinomial.multinomial_coef1(nfactor, factor);

                Console.WriteLine(factor[0].ToString().PadLeft(4) + "  "
                                  + factor[1].ToString().PadLeft(4) + "  "
                                  + "    "
                                  + ncomb.ToString().PadLeft(5) + "");

            }

            Console.WriteLine("");
            Console.WriteLine("  Level 5 of the TRINOMIAL coefficients:");

            n = 5;
            nfactor = 3;

            for (i = 0; i <= n; i++)
            {
                factor[0] = i;

                Console.WriteLine("");

                for (j = 0; j <= n - factor[0]; j++)
                {
                    factor[1] = j;
                    factor[2] = n - factor[0] - factor[1];

                    ncomb = Multinomial.multinomial_coef1(nfactor, factor);

                    Console.WriteLine(factor[0].ToString().PadLeft(4) + "  "
                                      + factor[1].ToString().PadLeft(4) + "  "
                                      + factor[2].ToString().PadLeft(4) + "  "
                                      + "    "
                                      + ncomb.ToString().PadLeft(5) + "");
                }
            }
        }

        public static void multinomial_coef2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTINOMIAL_COEF2_TEST tests MULTINOMIAL_COEF2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int MAXFACTOR = 5;

            int[] factor = new int[MAXFACTOR];
            int i;
            int j;
            int n;
            int ncomb;
            int nfactor;

            Console.WriteLine("");
            Console.WriteLine("MULTINOMIAL_COEF2_TEST");
            Console.WriteLine("  MULTINOMIAL_COEF2 computes multinomial");
            Console.WriteLine("  coefficients directly.");

            Console.WriteLine("");
            Console.WriteLine("  Line 10 of the BINOMIAL table:");
            Console.WriteLine("");

            n = 10;
            nfactor = 2;

            for (i = 0; i <= n; i++)
            {
                factor[0] = i;
                factor[1] = n - i;

                ncomb = Multinomial.multinomial_coef2(nfactor, factor);

                Console.WriteLine(factor[0].ToString().PadLeft(4) + "  "
                                  + factor[1].ToString().PadLeft(4) + "  "
                                  + "    "
                                  + ncomb.ToString().PadLeft(5) + "");

            }

            Console.WriteLine("");
            Console.WriteLine("  Level 5 of the TRINOMIAL coefficients:");

            n = 5;
            nfactor = 3;

            for (i = 0; i <= n; i++)
            {
                factor[0] = i;

                Console.WriteLine("");

                for (j = 0; j <= n - factor[0]; j++)
                {
                    factor[1] = j;
                    factor[2] = n - factor[0] - factor[1];

                    ncomb = Multinomial.multinomial_coef2(nfactor, factor);

                    Console.WriteLine(factor[0].ToString().PadLeft(4) + "  "
                                      + factor[1].ToString().PadLeft(4) + "  "
                                      + factor[2].ToString().PadLeft(4) + "  "
                                      + "    "
                                      + ncomb.ToString().PadLeft(5) + "");

                }

            }
        }

    }
}