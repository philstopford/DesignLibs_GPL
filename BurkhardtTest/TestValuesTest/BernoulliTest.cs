using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public static class BernoulliTest
    {
        public static void bernoulli_number_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BERNOULLI_NUMBER_VALUES_TEST tests BERNOULLI_NUMBER_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double c = 0;
            int n = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("BERNOULLI_NUMBER_VALUES_TEST:");
            Console.WriteLine("  BERNOULLI_NUMBER_VALUES returns values of ");
            Console.WriteLine("  the Bernoulli numbers.");
            Console.WriteLine("");
            Console.WriteLine("     N              B(N)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Bernoulli.bernoulli_number_values ( ref n_data, ref n, ref c );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  " + n.ToString().PadLeft(6)
                                       + "  " + c.ToString("0.################").PadLeft(24) + "");
            }
        }
        public static void bernoulli_poly_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BERNOULLI_POLY_VALUES_TEST tests BERNOULLI_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double b = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BERNOULLI_POLY_VALUES_TEST:");
            Console.WriteLine("  BERNOULLI_POLY_VALUES returns values of ");
            Console.WriteLine("  the Bernoulli Polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N     X      BERNOULLI(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Bernoulli.bernoulli_poly_values ( ref n_data, ref n, ref x, ref b );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + b.ToString().PadLeft(12) + "");
            }
        }

    }
}