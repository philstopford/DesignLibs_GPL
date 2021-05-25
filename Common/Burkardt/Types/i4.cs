using System;

namespace Burkardt.Types
{
    public class i4
    {
        public bool error { get; set; }
        public int val { get; set; }
        public int lchar { get; set; }
    }

    public class i4vec
    {
        public bool error { get; set; }
        public int[] ivec { get; set; }
    }

    public static partial class typeMethods
    {
        public static i4 s_to_i4 ( string s )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_I4 reads an I4 from a string.
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
            i4 ival = new i4() {val = 0};
            char c;
            int i;
            int isgn;
            int istate;

            istate = 0;
            isgn = 1;
            i = 0;

            for ( ; ; )
            {
                c = s[i];
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
     
        
        public static i4vec s_to_i4vec ( string s, int n )
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
            //    Output, int IVEC[N], the values read from the string.
            //
            //    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
            //
        {
            i4vec ret = new i4vec() {ivec = new int[n]};

            int begin = 0;
            int length = s.Length;

            for ( int i = 0; i < n; i++ )
            {
                i4 res = s_to_i4 ( s.Substring(begin,length) );

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
        

        public static int i4_wrap ( int ival, int ilo, int ihi )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_WRAP forces an I4 to lie between given limits by wrapping.
        //
        //  Example:
        //
        //    ILO = 4, IHI = 8
        //
        //    I   Value
        //
        //    -2     8
        //    -1     4
        //     0     5
        //     1     6
        //     2     7
        //     3     8
        //     4     4
        //     5     5
        //     6     6
        //     7     7
        //     8     8
        //     9     4
        //    10     5
        //    11     6
        //    12     7
        //    13     8
        //    14     4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int IVAL, an integer value.
        //
        //    Input, int ILO, IHI, the desired bounds for the integer value.
        //
        //    Output, int I4_WRAP, a "wrapped" version of IVAL.
        //
        {
            int value;

            int jlo = Math.Min ( ilo, ihi );
            int jhi = Math.Max ( ilo, ihi );

            int wide = jhi + 1 - jlo;

            if ( wide == 1 )
            {
                value = jlo;
            }
            else
            {
                value = jlo + Math.Abs( (ival - jlo) % wide );
            }

            return value;
        }
        
        public static void i4block_print ( int l, int m, int n, int[] a, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4BLOCK_PRINT prints an I4BLOCK.
        //
        //  Discussion:
        //
        //    An I4BLOCK is a 3D array of I4 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int L, M, N, the dimensions of the block.
        //
        //    Input, int A[L*M*N], the matrix to be printed.
        //
        //    Input, string TITLE, a title.
        //
        {
            Console.WriteLine();
            Console.WriteLine(title);

            for (int k = 0; k < n; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  K = " + k);
                for (int jlo = 0; jlo < m; jlo = jlo + 10 )
                {
                    int jhi = Math.Min ( jlo + 10, m );
                    Console.WriteLine();
                    string cout = "        J:";
                    for (int j = jlo; j < jhi; j++ )
                    {
                        cout += "  " + j.ToString().PadLeft(6);
                    }
                    Console.WriteLine(cout);
                    Console.WriteLine("       I:");
                    for (int i = 0; i < l; i++ )
                    {
                        cout = "  " + i.ToString().PadLeft(6) + ": ";
                        for ( int j = jlo; j < jhi; j++ )
                        {
                            cout += "  " + a[i+j*l+k*l*m].ToString().PadLeft(6);
                        }
                        Console.WriteLine(cout);
                    }
                }
            }
        }
        
        
        
    }
}