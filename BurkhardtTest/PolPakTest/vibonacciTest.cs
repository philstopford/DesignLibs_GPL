using System;
using Burkardt.Sequence;

namespace PolPakTest;

public static class vibonacciTest
{
    public static void vibonacci_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VIBONACCI_TEST tests VIBONACCI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;

        int i;
        int seed;
        int[] v1 = new int[N];
        int[] v2 = new int[N];
        int[] v3 = new int[N];

        Console.WriteLine("");
        Console.WriteLine("VIBONACCI_TEST");
        Console.WriteLine("  VIBONACCI computes a Vibonacci sequence.");
        Console.WriteLine("");
        Console.WriteLine("  We compute the series 3 times.");
        Console.WriteLine("");
        Console.WriteLine("     I      V1      V2      V3");
        Console.WriteLine("");

        seed = 123456789;

        Vibonacci.vibonacci(N, ref seed, ref v1);
        Vibonacci.vibonacci(N, ref seed, ref v2);
        Vibonacci.vibonacci(N, ref seed, ref v3);

        for (i = 0; i < N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + v1[i].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + v2[i].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + v3[i].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

    }

}