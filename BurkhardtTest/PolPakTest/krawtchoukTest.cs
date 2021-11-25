using System;
using System.Globalization;
using Burkardt.PolynomialNS;

namespace PolPakTest;

public static class krawtchoukTest
{
    public static void krawtchouk_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KRAWTCHOUK_TEST tests KRAWTCHOUK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 2;
        const int N = 5;

        double[] p_test = { 0.25, 0.5 };
        int test;
        double[] value = new double[N+1];

        Console.WriteLine("");
        Console.WriteLine("KRAWTCHOUK_TEST:");
        Console.WriteLine("  KRAWTCHOUK evaluates Krawtchouk polynomials.");
        Console.WriteLine("");
        Console.WriteLine("        N         P         X          M      K(N,P,X,M)");
        Console.WriteLine("");

        const int m = 5;

        for ( test = 0; test < TEST_NUM; test++ )
        {
            double p = p_test[test];

            Console.WriteLine("");

            int j;
            for ( j = 0; j <= 5; j++ )
            {
                double x = j / 2.0;

                Krawtchouk.krawtchouk ( N, p, x, m, ref value );

                Console.WriteLine("");
                int i;
                for ( i = 0; i <= 5; i++ )
                {

                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + p.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + value[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

}