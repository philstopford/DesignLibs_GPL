using System;
using System.Linq;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static float r4vec_max(int n, float[] dvec)
        {
            if (dvec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, dvec.Length);

            if (n == dvec.Length)
            {
                return dvec.Max();
            }
            
            return dvec.Take(n).Max();
        }

        public static float r4vec_mean ( int n, float[] x )
        {
            if (x.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, x.Length);

            if (n == x.Length)
            {
                return x.Average();
            }
                    
            return x.Take(n).Average();
        }

        public static float r4vec_min ( int n, float[] dvec )
        {
            if (dvec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, dvec.Length);

            if (n == dvec.Length)
            {
                return dvec.Min();
            }
                    
            return dvec.Take(n).Min();
        }

        public static float r4vec_variance ( int n, float[] x )
        {
            float mean = r4vec_mean ( n, x );

            float variance = 0.0f;
            for (int i = 0; i < n; i++ )
            {
                variance = variance + ( x[i] - mean ) * ( x[i] - mean );
            }

            if ( 1 < n )
            {
                variance = variance / ( float ) ( n - 1 );
            }
            else
            {
                variance = 0.0f;
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

            if ( max_print <= 0 )
            {
                return;
            }

            if ( n <= 0 )
            {
                return;
            }

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");

            if ( n <= max_print )
            {
                for ( i = 0; i < n; i++ )
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                        + "  " + a[i].ToString().PadLeft(14) + "");
                }
            }
            else if ( 3 <= max_print )
            {
                for ( i = 0; i < max_print - 2; i++ )
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + a[i].ToString().PadLeft(14) + "");
                }
                Console.WriteLine("  ........  ..............");
                i = n - 1;
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + a[i].ToString().PadLeft(14) + "");
            }
            else
            {
                for ( i= 0; i < max_print - 1; i++ )
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + a[i].ToString().PadLeft(14) + "");
                }
                i = max_print - 1;
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + a[i].ToString().PadLeft(14) + ""
                    + "  " + "...more entries...");
            }

            return;
        }
        
    }
}