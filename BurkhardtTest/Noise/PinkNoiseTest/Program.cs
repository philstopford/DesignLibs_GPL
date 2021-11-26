using System;
using System.Globalization;
using Burkardt.Noise;

namespace PinkNoiseTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PINK_NOISE_TEST.
        //
        //  Discussion:
        //
        //    PINK_NOISE_TEST tests the PINK_NOISE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("PINK_NOISE_TEST:");
        Console.WriteLine("  Test the PINK_NOISE library.");

        test01();
        test02();
        test03();
        test04();

        Console.WriteLine("");
        Console.WriteLine("PINK_NOISE_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests WRAP2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  WRAP2 performs a circular wrap.");
        Console.WriteLine("  Q is expected to range between 0 and M.");
        Console.WriteLine("  WRAP2 takes an input value of Q, and either");
        Console.WriteLine("  increments it by M+1 until in the range, or");
        Console.WriteLine("  decrements it by M+1 until in the range,");
        Console.WriteLine("  and returns the result as the function value.");

        for (m = 2; m <= 4; m++)
        {
            Console.WriteLine("");
            Console.WriteLine("   M  Qin  Qout");
            Console.WriteLine("");
            int i;
            for (i = -5; i < 3 * m; i++)
            {
                int q = i;
                int q_in = q;
                Pink.wrap2(m, ref q);
                Console.WriteLine("  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q_in.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "");
            }
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests CDELAY2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  CDELAY2 is a circular buffer implementation");
        Console.WriteLine("  of an M-fold delay.  Q is a counter");
        Console.WriteLine("  which is decremented by CDELAY2, but reset to M");
        Console.WriteLine("  after it reaches 0.");

        for (m = 2; m <= 4; m++)
        {
            Console.WriteLine("");
            Console.WriteLine("   I   M  Qin  Qout");
            Console.WriteLine("");
            int q = m;
            int i;
            for (i = 1; i <= 3 * (m + 1); i++)
            {
                int q_in = q;
                Pink.cdelay2(m, ref q);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q_in.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "");
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests RANH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  RANH is a random hold function.");
        Console.WriteLine("  Given a value U and a delay D, it returns the value");
        Console.WriteLine("  U for D calls, then resets U.");

        for (d = 5; 1 <= d; d--)
        {
            Console.WriteLine("");
            Console.WriteLine("   I   D   Q      U           Y");
            Console.WriteLine("");
            double u = 0.5;
            int q = 3;
            int i;
            for (i = 1; i <= 20; i++)
            {
                double y = Pink.ranh(d, ref u, ref q);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + d.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests RAN1F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int b;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  RAN1F generates random values with an approximate");
        Console.WriteLine("  1/F distribution.");

        for (b = 1; b < 32; b *= 2)
        {
            double[] u = new double[b];
            int[] q = new int[b];
            int rep;
            for (rep = 1; rep <= 4; rep++)
            {
                int i;
                for (i = 0; i < b; i++)
                {
                    u[i] = entropyRNG.RNG.nextdouble() - 0.5;
                }

                for (i = 0; i < b; i++)
                {
                    q[i] = 0;
                }

                Console.WriteLine("");
                Console.WriteLine("   B   I      Y");
                Console.WriteLine("");

                for (i = 1; i <= 20; i++)
                {
                    double y = Pink.ran1f(b, ref u, ref q);
                    Console.WriteLine("  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                           + "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                           + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                }
            }
        }
    }
}