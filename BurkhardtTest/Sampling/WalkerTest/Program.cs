using System;
using Burkardt.Sampling;

namespace WalkerTest;

internal static class Program
{
    private static void Main()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for WALKER_SAMPLE_TEST.
        //
        //  Discussion:
        //
        //    WALKER_SAMPLE_TEST tests WALKER_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WALKER_SAMPLE_TEST:");
        Console.WriteLine("  Test the WALKER_SAMPLE library.");

        Walker.walker_sampler_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("WALKER_SAMPLE_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}