using System;
using System.Globalization;
using Burkardt.CircleNS;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace CircleMonteCarloTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CIRCLE_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    CIRCLE_MONTE_CARLO_TEST tests the CIRCLE_MONTE_CARLO library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 June 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CIRCLE_MONTE_CARLO_TEST");

        Console.WriteLine("  Test the CIRCLE_MONTE_CARLO library.");

        circle01_sample_random_test();
        circle01_sample_ergodic_test();

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void circle01_sample_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE01_SAMPLE_RANDOM_TEST uses CIRCLE01_SAMPLE_RANDOM with an increasing number of points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[2];
        int[] e_test =  {
                0, 0,
                2, 0,
                0, 2,
                4, 0,
                2, 2,
                0, 4,
                6, 0
            }
            ;
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE0_SAMPLE_RANDOM_TEST");
        Console.WriteLine("  CIRCLE01_SAMPLE_RANDOM randomly samples the unit circle.");
        Console.WriteLine("  Use it to estimate integrals.");

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N        1              X^2             Y^2" + 
                          "             X^4           X^2Y^2          Y^4          X^6");
        Console.WriteLine("");

        int n = 1;

        while (n <= 65536)
        {
            double[] x = MonteCarlo.circle01_sample_random(n, ref seed);
            string cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 7; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    e[i] = e_test[i + j * 2];
                }

                double[] value = Monomial.monomial_value(2, n, e, x);

                double result = Integrals.circle01_length() * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
                
            n = 2 * n;
        }

        Console.WriteLine("");
        string cout2 = "     Exact";
        for (j = 0; j < 7; j++)
        {
            for (i = 0; i < 2; i++)
            {
                e[i] = e_test[i + j * 2];
            }

            double exact = Integrals.circle01_monomial_integral(e);
            cout2 += "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout2);
    }

    private static void circle01_sample_ergodic_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE01_SAMPLE_ERGODIC_TEST uses CIRCLE01_SAMPLE_ERGODIC with an increasing number of points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 June 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[2];
        int[] e_test =  {
                0, 0,
                2, 0,
                0, 2,
                4, 0,
                2, 2,
                0, 4,
                6, 0
            }
            ;
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE0_SAMPLE_ERGODIC_TEST");
        Console.WriteLine("  CIRCLE01_SAMPLE_ERGODIC ergodically samples the unit circle.");
        Console.WriteLine("  Use it to estimate integrals.");

        double angle = 0.0;

        Console.WriteLine("");
        Console.WriteLine("         N        1              X^2             Y^2" + 
                          "             X^4           X^2Y^2          Y^4          X^6");
        Console.WriteLine("");

        int n = 1;

        while (n <= 65536)
        {
            double[] x = MonteCarlo.circle01_sample_ergodic(n, ref angle);
            string cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 7; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    e[i] = e_test[i + j * 2];
                }

                double[] value = Monomial.monomial_value(2, n, e, x);

                double result = Integrals.circle01_length() * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine("");
        string cout2 = "     Exact";
        for (j = 0; j < 7; j++)
        {
            for (i = 0; i < 2; i++)
            {
                e[i] = e_test[i + j * 2];
            }

            double exact = Integrals.circle01_monomial_integral(e);
            cout2 += "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout2);
    }
}