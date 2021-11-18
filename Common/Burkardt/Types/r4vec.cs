﻿using System;
using System.Globalization;
using System.Linq;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static float r4vec_max(int n, float[] dvec)
    {
        switch (dvec.Length)
        {
            case <= 0:
                return 0;
        }

        // Limit to the number of items in the array as a maximum
        n = Math.Min(n, dvec.Length);

        return n == dvec.Length ? dvec.Max() : dvec.Take(n).Max();
    }

    public static float r4vec_mean ( int n, float[] x )
    {
        switch (x.Length)
        {
            case <= 0:
                return 0;
        }

        // Limit to the number of items in the array as a maximum
        n = Math.Min(n, x.Length);

        return n == x.Length ? x.Average() : x.Take(n).Average();
    }

    public static float r4vec_min ( int n, float[] dvec )
    {
        switch (dvec.Length)
        {
            case <= 0:
                return 0;
        }

        // Limit to the number of items in the array as a maximum
        n = Math.Min(n, dvec.Length);

        return n == dvec.Length ? dvec.Min() : dvec.Take(n).Min();
    }

    public static float r4vec_variance ( int n, float[] x )
    {
        float mean = r4vec_mean ( n, x );

        float variance = 0.0f;
        for (int i = 0; i < n; i++ )
        {
            variance += ( x[i] - mean ) * ( x[i] - mean );
        }

        switch (n)
        {
            case > 1:
                variance /= n - 1;
                break;
            default:
                variance = 0.0f;
                break;
        }

        return variance;
    }

        
    public static void r4vec_print_part ( int n, float[] a, int max_print, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_PRINT_PART prints "part" of an R4VEC.
        //
        //  Discussion:
        //
        //    The user specifies MAX_PRINT, the maximum number of lines to print.
        //
        //    If N, the size of the vector, is no more than MAX_PRINT, then
        //    the entire vector is printed, one entry per line.
        //
        //    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
        //    followed by a line of periods suggesting an omission,
        //    and the last entry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of the vector.
        //
        //    Input, float A[N], the vector to be printed.
        //
        //    Input, int MAX_PRINT, the maximum number of lines
        //    to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        switch (max_print)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 0:
                return;
        }

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        if ( n <= max_print )
        {
            for ( i = 0; i < n; i++ )
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
        else
        {
            switch (max_print)
            {
                case >= 3:
                {
                    for ( i = 0; i < max_print - 2; i++ )
                    {
                        Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + "  " + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    }
                    Console.WriteLine("  ........  ..............");
                    i = n - 1;
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
                }
                default:
                {
                    for ( i= 0; i < max_print - 1; i++ )
                    {
                        Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + "  " + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    }
                    i = max_print - 1;
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + ""
                                           + "  " + "...more entries...");
                    break;
                }
            }
        }
    }
        
}