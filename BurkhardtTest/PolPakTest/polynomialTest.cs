using System;
using Burkardt.PolynomialNS;

namespace PolPakTest;

public static class polynomialTest
{
    public static void poly_coef_count_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_COEF_COUNT_TEST tests POLY_COEF_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int degree;
        int dim;

        Console.WriteLine("");
        Console.WriteLine("POLY_COEF_COUNT_TEST");
        Console.WriteLine("  POLY_COEF_COUNT counts the number of coefficients");
        Console.WriteLine("  in a polynomial of degree DEGREE and dimension DIM.");
        Console.WriteLine("");
        Console.WriteLine(" Dimension    Degree     Count");

        for ( dim = 1; dim <= 10; dim += 3 )
        {
            Console.WriteLine("");
            for ( degree = 0; degree <= 5; degree++ )
            {
                Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                       + "  " + degree.ToString().PadLeft(8)
                                       + "  " + Polynomial.poly_coef_count ( dim, degree ).ToString().PadLeft(8) + "");
            }
        }

    }
        
    public static void trinomial_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRINOMIAL_TEST tests TRINOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        int t;

        Console.WriteLine("");
        Console.WriteLine("TRINOMIAL_TEST");
        Console.WriteLine("  TRINOMIAL evaluates the trinomial coefficient:");
        Console.WriteLine("");
        Console.WriteLine("  T(I,J,K) = (I+J+K)! / I! / J! / K!");
        Console.WriteLine("");
        Console.WriteLine("     I     J     K    T(I,J,K)");
        Console.WriteLine("");

        for ( k = 0; k <= 4; k++ )
        {
            for ( j = 0; j <= 4; j++ )
            {
                for ( i = 0; i <= 4; i++ )
                {
                    t = Trinomial.trinomial ( i, j, k );
                    Console.WriteLine("  " + i.ToString().PadLeft(4)
                                           + "  " + j.ToString().PadLeft(4)
                                           + "  " + k.ToString().PadLeft(4)
                                           + "  " + t.ToString().PadLeft(8) + "");
                }
            }
        }

    }

}