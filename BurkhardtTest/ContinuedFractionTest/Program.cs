using System;
using Burkardt;

namespace ContinuedFractionTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONTINUED_FRACTION_TEST tests the CONTINUED_FRACTION library.
        //
        //  Licensing:
        //
        //    I don't care what you do with this code.
        //
        //  Modified:
        //
        //    07 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CONTINUED_FRACTION_TEST:");
        Console.WriteLine("  CONTINUED_FRACTION is a library for dealing with");
        Console.WriteLine("  expresssions representing a continued fraction.");

        ContinuedFraction.i4cf_evaluate_test ( );
        ContinuedFraction.i4scf_evaluate_test ( );
        ContinuedFraction.i8cf_evaluate_test ( );
        ContinuedFraction.i8scf_evaluate_test ( );
        ContinuedFraction.r8_to_i4scf_test ( );
        ContinuedFraction.r8_to_i8scf_test ( );
        ContinuedFraction.r8cf_evaluate_test ( );
        ContinuedFraction.r8scf_evaluate_test ( );

        Console.WriteLine("");
        Console.WriteLine("CONTINUED_FRACTION_TEST:");
        Console.WriteLine("  Normal end of execution.");
    }
}