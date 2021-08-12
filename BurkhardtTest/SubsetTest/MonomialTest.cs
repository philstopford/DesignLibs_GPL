using System;
using Burkardt;

namespace SubsetTest
{
    public static class MonomialTest
    {
        public static void monomial_count_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_COUNT_TEST tests MONOMIAL_COUNT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree_max = 9;
            int dim;
            int total;

            Console.WriteLine("");
            Console.WriteLine("MONOMIAL_COUNT_TEST");
            Console.WriteLine("  MONOMIAL_COUNT counts the number of monomials of");
            Console.WriteLine("  degrees 0 through DEGREE_MAX in a space of dimension DIM.");

            Console.WriteLine("");
            Console.WriteLine("  DIM  Count");
            Console.WriteLine("");

            for (dim = 1; dim <= 6; dim++)
            {
                total = Monomial.monomial_count(degree_max, dim);
                Console.WriteLine("  " + dim.ToString().PadLeft(2)
                                  + "  " + total.ToString().PadLeft(8) + "");
            }
        }

        public static void monomial_counts_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_COUNTS_TEST tests MONOMIAL_COUNTS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] counts;
            int degree;
            int degree_max = 9;
            int dim;
            int total;

            Console.WriteLine("");
            Console.WriteLine("MONOMIAL_COUNTS");
            Console.WriteLine("  MONOMIAL_COUNTS counts the number of monomials of");
            Console.WriteLine("  degrees 0 through DEGREE_MAX in a space of dimension DIM.");

            for (dim = 1; dim <= 6; dim++)
            {
                counts = Monomial.monomial_counts(degree_max, dim);

                Console.WriteLine("");
                Console.WriteLine("  DIM = " + dim + "");
                Console.WriteLine("");

                for (degree = 0; degree <= degree_max; degree++)
                {
                    Console.WriteLine("  " + degree.ToString().PadLeft(8)
                                      + "  " + counts[degree].ToString().PadLeft(8) + "");
                }

                total = 0;
                for (degree = 0; degree <= degree_max; degree++)
                {
                    total = total + counts[degree];
                }

                Console.WriteLine("");
                Console.WriteLine("     Total"
                                  + "  " + total.ToString().PadLeft(8) + "");
            }
        }

    }
}