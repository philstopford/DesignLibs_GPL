using System;
using System.Linq;

namespace Burkardt.Types
{
    public class i8
    {
        public bool error { get; set; }
        public long val { get; set; }
        public int lchar { get; set; }
    }

    public class i8vec
    {
        public bool error { get; set; }
        public long[] ivec { get; set; }
    }

    public static partial class typeMethods
    {
        public static i8vec s_to_i8vec ( string s, int n )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_TO_I4VEC reads an I4VEC from a string.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string S, the string to be read.
            //
            //    Input, int N, the number of values expected.
            //
            //    Output, long IVEC[N], the values read from the string.
            //
            //    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
            //
        {
            i8vec ret = new i8vec() {ivec = new long[n]};

            string[] tokens = Helpers.splitStringByWhitespace(s);

            for (int i = 0; i < n; i++)
            {
                ret.ivec[i] = s_to_i8(tokens[i]).val;
            }

            return ret;
        }
        
        public static long i8vec_max(int n, long[] ivec)
        {
            if (ivec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, ivec.Length);

            if (n == ivec.Length)
            {
                return ivec.Max();
            }
            
            return ivec.Take(n).Max();
        }
        
        public static long i8vec_min ( int n, long[] ivec )
        {
            if (ivec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, ivec.Length);

            if (n == ivec.Length)
            {
                return ivec.Min();
            }
                    
            return ivec.Take(n).Min();
        }
        
        public static void i8vec_print ( int n, long[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8VEC_PRINT prints an I8VEC.
        //
        //  Discussion:
        //
        //    An I8VEC is a vector of I8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 August 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, long long int A[N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
        {
            int i;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");
            for ( i = 0; i < n; i++ )
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                    + ": " + a[i].ToString().PadLeft(24)  + "");
            }
        }

        public static double i8vec_variance(int n, long[] x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_VARIANCE returns the variance of an I4VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 May 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, int X[N], the vector whose variance is desired.
            //
            //    Output, double I4VEC_VARIANCE, the variance of the vector entries.
            //
        {
            double mean = i8vec_mean(n, x);

            double variance = 0.0;
            for (int i = 0; i < n; i++)
            {
                variance = variance + ((double) x[i] - mean) * ((double) x[i] - mean);
            }

            if (1 < n)
            {
                variance = variance / (double) (n - 1);
            }
            else
            {
                variance = 0.0;
            }

            return variance;
        }
        
        public static double i8vec_mean ( int n, long[] x )
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

    }
    
}
