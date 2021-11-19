using System;
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
        int  TEST_NUM = 2;
        int  N = 5;

        int i = 0;
        int j = 0;
        int m = 0;
        int n = 0;
        double p = 0;
        double[] p_test = { 0.25, 0.5 };
        int test = 0;
        double x = 0;
        double[] value = new double[N+1];

        Console.WriteLine("");
        Console.WriteLine("KRAWTCHOUK_TEST:");
        Console.WriteLine("  KRAWTCHOUK evaluates Krawtchouk polynomials.");
        Console.WriteLine("");
        Console.WriteLine("        N         P         X          M      K(N,P,X,M)");
        Console.WriteLine("");

        m = 5;
        n = N;

        for ( test = 0; test < TEST_NUM; test++ )
        {
            p = p_test[test];

            Console.WriteLine("");

            for ( j = 0; j <= 5; j++ )
            {
                x = j / 2.0;

                Krawtchouk.krawtchouk ( n, p, x, m, ref value );

                Console.WriteLine("");
                for ( i = 0; i <= 5; i++ )
                {

                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + p.ToString().PadLeft(8)
                                           + "  " + x.ToString().PadLeft(8)
                                           + "  " + m.ToString().PadLeft(8)
                                           + "  " + value[i].ToString().PadLeft(14) + "");
                }
            }
        }
    }

}