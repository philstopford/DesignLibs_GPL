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
    }
}