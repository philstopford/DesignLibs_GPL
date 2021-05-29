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
        public static i8 s_to_i8 ( string s )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_I8 reads an I8 from a string.
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
        //    Input, string S, a string to be examined.
        //
        //    Output, int *LAST, the last character of S used to make IVAL.
        //
        //    Output, bool *ERROR is TRUE if an error occurred.
        //
        //    Output, int *S_TO_I4, the integer value read from the string.
        //    If the string is blank, then IVAL will be returned 0.
        //
        {
            i8 ival = new i8() {val = 0};

            int istate = 0;
            long isgn = 1;
            int i = 0;

            while ( i < s.Length )
            {
                char c = s[i];
                i = i + 1;
                //
                //  Haven't read anything.
                //
                if ( istate == 0 )
                {
                    if ( c == ' ' )
                    {
                    }
                    else if ( c == '-' )
                    {
                        istate = 1;
                        isgn = -1;
                    }
                    else if ( c == '+' )
                    {
                        istate = 1;
                        isgn = + 1;
                    }
                    else if ( '0' <= c && c <= '9' )
                    {
                        istate = 2;
                        ival.val = c - '0';
                    }
                    else
                    {
                        ival.error = true;
                        return ival;
                    }
                }
                //
                //  Have read the sign, expecting digits.
                //
                else if ( istate == 1 )
                {
                    if ( c == ' ' )
                    {
                    }
                    else if ( '0' <= c && c <= '9' )
                    {
                        istate = 2;
                        ival.val = c - '0';
                    }
                    else
                    {
                        ival.error = true;
                        return ival;
                    }
                }
                //
                //  Have read at least one digit, expecting more.
                //
                else if ( istate == 2 )
                {
                    if ( '0' <= c && c <= '9' )
                    {
                        ival.val = 10 * (ival.val) + c - '0';
                    }
                    else
                    {
                        ival.val = isgn * ival.val;
                        ival.lchar = i - 1;
                        return ival;
                    }
                }
            }
            //
            //  If we read all the characters in the string, see if we're OK.
            //
            if ( istate == 2 )
            {
                ival.val = isgn * ival.val;
                ival.lchar = typeMethods.s_len_trim ( s );
            }
            else
            {
                ival.error = true;
                ival.lchar = 0;
            }

            return ival;
        }
     
        
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

            int begin = 0;
            int length = s.Length;

            for ( int i = 0; i < n; i++ )
            {
                i8 res = s_to_i8 ( s.Substring(begin,length) );

                ret.ivec[i] = res.val;
                int lchar = res.lchar;
                
                if ( res.error )
                {
                    return ret;
                }
                begin = begin + lchar;
                length = length - lchar;
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
