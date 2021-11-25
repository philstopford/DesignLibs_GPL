using System;
using System.Globalization;
using Burkardt.PolynomialNS;

namespace PolPakTest;

public static class charlierTest
{
    public static void charlier_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHARLIER_TEST tests CHARLIER.
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
        const int TEST_NUM = 5;
        const int N = 5;

        double[] a_test = { 0.25, 0.5, 1.0, 2.0, 10.0 };
        int test;
        double[] value = new double[N+1];

        Console.WriteLine("");
        Console.WriteLine("CHARLIER_TEST:");
        Console.WriteLine("  CHARLIER evaluates Charlier polynomials.");
        Console.WriteLine("");
        Console.WriteLine("       N      A         X        P(N,A,X)");
        Console.WriteLine("");

        for ( test = 0; test < TEST_NUM; test++ )
        {
            int n = N;
            double a = a_test[test];

            Console.WriteLine("");

            int j;
            for ( j = 0; j <= 5; j++ )
            {
                double x = j / 2.0;

                Charlier.charlier ( n, a, x, ref value );

                Console.WriteLine("");
                int i;
                for ( i = 0; i <= 5; i++ )
                {

                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                           + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + value[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }

    }

}