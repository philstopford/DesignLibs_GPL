using System;
using Burkardt.PolynomialNS;

namespace PolPakTest;

public static class meixnerTest
{
    public static void meixner_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MEIXNER_TEST tests MEIXNER.
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
        int  N = 5;
        int TEST_NUM = 3;

        double beta;
        double[] beta_test = { 0.5, 1.0, 2.0 };
        double c;
        double[] c_test = { 0.125, 0.25, 0.5 };
        int i;
        int j;
        int n;
        int test;
        double[] v = new double[N+1];
        double x;

        Console.WriteLine("");
        Console.WriteLine("MEIXNER_TEST:");
        Console.WriteLine("  MEIXNER evaluates Meixner polynomials.");
        Console.WriteLine("");
        Console.WriteLine("       N      BETA         C         X        M(N,BETA,C,X)");

        for ( test = 0; test < TEST_NUM; test++ )
        {
            n = N;
            beta = beta_test[test];
            c = c_test[test];

            for ( j = 0; j <= 5; j++ )
            {
                x = j / 2.0;

                Meixner.meixner ( n, beta, c, x, ref v );

                Console.WriteLine("");

                for ( i = 0; i <= n; i++ )
                {
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + beta.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + c.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }
}