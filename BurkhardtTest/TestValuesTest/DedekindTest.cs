using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class DedekindTest
{
    public static void dedekind_sum_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEDEKIND_SUM_VALUES_TEST tests DEDEKIND_SUM_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d = 0;
        int n = 0;
        int n_data;
        int p = 0;
        int q = 0;
        Console.WriteLine("");
        Console.WriteLine("DEDEKIND_SUM_VALUES_TEST:");
        Console.WriteLine("  DEDEKIND_SUM_VALUES stores values of the Dedekind sum");
        Console.WriteLine("  (N/D) = Dedekind_Sum(P,Q).");
        Console.WriteLine("");
        Console.WriteLine("       P       Q       N       D");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            Dedekind.dedekind_sum_values ( ref n_data, ref p, ref q, ref n, ref d );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + p.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + q.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + d.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }
}