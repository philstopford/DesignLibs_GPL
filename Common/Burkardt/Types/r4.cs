using System;
using System.Linq;

namespace Burkardt.Types
{
    public class r4
    {
        public bool error { get; set; }
        public float val { get; set; }
        public int lchar { get; set; }
    }

    public class r4vec
    {
        public bool error { get; set; }
        public float[] rvec { get; set; }
    }

    public static partial class typeMethods
    {
        public static r4 s_to_r4 ( string s )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R4 reads an R4 from a string.
        //
        //  Discussion:
        //
        //    This routine will read as many characters as possible until it reaches
        //    the end of the string, or encounters a character which cannot be
        //    part of the real number.
        //
        //    Legal input is:
        //
        //       1 blanks,
        //       2 '+' or '-' sign,
        //       2.5 spaces
        //       3 integer part,
        //       4 decimal point,
        //       5 fraction part,
        //       6 'E' or 'e' or 'D' or 'd', exponent marker,
        //       7 exponent sign,
        //       8 exponent integer part,
        //       9 exponent decimal point,
        //      10 exponent fraction part,
        //      11 blanks,
        //      12 final comma or semicolon.
        //
        //    with most quantities optional.
        //
        //  Example:
        //
        //    S                 R
        //
        //    '1'               1.0
        //    '     1   '       1.0
        //    '1A'              1.0
        //    '12,34,56'        12.0
        //    '  34 7'          34.0
        //    '-1E2ABCD'        -100.0
        //    '-1X2ABCD'        -1.0
        //    ' 2E-1'           0.2
        //    '23.45'           23.45
        //    '-4.2E+2'         -420.0
        //    '17d2'            1700.0
        //    '-14e-2'         -0.14
        //    'e2'              100.0
        //    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
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
        //    Input, string S, the string containing the
        //    data to be read.  Reading will begin at position 1 and
        //    terminate at the end of the string, or when no more
        //    characters can be read to form a legal real.  Blanks,
        //    commas, or other nonnumeric data will, in particular,
        //    cause the conversion to halt.
        //
        //    Output, int *LCHAR, the number of characters read from
        //    the string to form the number, including any terminating
        //    characters such as a trailing comma or blanks.
        //
        //    Output, bool *ERROR, is true if an error occurred.
        //
        //    Output, double S_TO_R8, the real value that was read from the string.
        //
        {
            r4 ret = new r4 {lchar = -1};
            float rexp;
            char TAB = (char) 9;

            int nchar = typeMethods.s_len_trim ( s );
            int isgn = 1;
            float rtop = 0.0f;
            float rbot = 1.0f;
            int jsgn = 1;
            int jtop = 0;
            int jbot = 1;
            int ihave = 1;
            int iterm = 0;

            for ( ; ; )
            {
                char c = s[ret.lchar+1];
                ret.lchar = ret.lchar + 1;
                //
                //  Blank or TAB character.
                //
                if ( c == ' ' || c == TAB )
                {
                    if ( ihave == 2 )
                    {
                    }
                    else if ( ihave == 6 || ihave == 7 )
                    {
                        iterm = 1;
                    }
                    else if ( 1 < ihave )
                    {
                        ihave = 11;
                    }
                }
                //
                //  Comma.
                //
                else if ( c == ',' || c == ';' )
                {
                    if ( ihave != 1 )
                    {
                        iterm = 1;
                        ihave = 12;
                        ret.lchar = ret.lchar + 1;
                    }
                }
                //
                //  Minus sign.
                //
                else if ( c == '-' )
                {
                    if ( ihave == 1 )
                    {
                        ihave = 2;
                        isgn = -1;
                    }
                    else if ( ihave == 6 )
                    {
                        ihave = 7;
                        jsgn = -1;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Plus sign.
                //
                else if ( c == '+' )
                {
                    if ( ihave == 1 )
                    {
                        ihave = 2;
                    }
                    else if ( ihave == 6 )
                    {
                        ihave = 7;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Decimal point.
                //
                else if ( c == '.' )
                {
                    if ( ihave < 4 )
                    {
                        ihave = 4;
                    }
                    else if ( 6 <= ihave && ihave <= 8 )
                    {
                        ihave = 9;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Exponent marker.
                //
                else if ( ( Char.ToUpper(c) == 'E' ) || ( Char.ToUpper(c) == 'D' ) )
                {
                    if ( ihave < 6 )
                    {
                        ihave = 6;
                    }
                    else
                    {
                        iterm = 1;
                    }
                }
                //
                //  Digit.
                //
                else if ( ihave < 11 && '0' <= c && c <= '9' )
                {
                    if ( ihave <= 2 )
                    {
                        ihave = 3;
                    }
                    else if ( ihave == 4 )
                    {
                        ihave = 5;
                    }
                    else if ( ihave == 6 || ihave == 7 )
                    {
                        ihave = 8;
                    }
                    else if ( ihave == 9 )
                    {
                        ihave = 10;
                    }

                    int ndig = typeMethods.ch_to_digit ( c );

                    if ( ihave == 3 )
                    {
                        rtop = (float) 10.0 * rtop + ndig;
                    }
                    else if ( ihave == 5 )
                    {
                        rtop = (float) 10.0 * rtop + ndig;
                        rbot = (float) 10.0 * rbot;
                    }
                    else if ( ihave == 8 )
                    {
                        jtop = 10 * jtop + ndig;
                    }
                    else if ( ihave == 10 )
                    {
                        jtop = 10 * jtop + ndig;
                        jbot = 10 * jbot;
                    }
                }
                //
                //  Anything else is regarded as a terminator.
                //
                else
                {
                    iterm = 1;
                }
                //
                //  If we haven't seen a terminator, and we haven't examined the
                //  entire string, go get the next character.
                //
                if ( iterm == 1 || nchar <= ret.lchar + 1 )
                {
                    break;
                }

            }
            //
            //  If we haven't seen a terminator, and we have examined the
            //  entire string, then we're done, and LCHAR is equal to NCHAR.
            //
            if ( iterm != 1 && (ret.lchar) + 1 == nchar )
            {
                ret.lchar = nchar;
            }
            //
            //  Number seems to have terminated.  Have we got a legal number?
            //  Not if we terminated in states 1, 2, 6 or 7!
            //
            if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
            {
                ret.error = true;
                return ret;
            }
            //
            //  Number seems OK.  Form it.
            //
            if ( jtop == 0 )
            {
                rexp = 1.0f;
            }
            else
            {
                if ( jbot == 1 )
                {
                    rexp = (float) Math.Pow ( 10.0, jsgn * jtop );
                }
                else
                {
                    rexp = jsgn * jtop;
                    rexp = rexp / jbot;
                    rexp = (float) Math.Pow ( 10.0, rexp );
                }
            }

            ret.val = isgn * rexp * rtop / rbot;

            return ret;
        }

        
        public static r4vec s_to_r4vec ( string s, int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R4VEC reads an R4VEC from a string.
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
        //    Output, float RVEC[N], the values read from the string.
        //
        //    Output, bool S_TO_R4VEC, is true if an error occurred.
        //
        {
            int begin = 0;
            int length = s.Length;

            r4vec ret = new r4vec() {rvec = new float[n]} ;

            for (int i = 0; i < n; i++ )
            {
                r4 res = s_to_r4 ( s.Substring(begin,length));

                ret.rvec[i] = res.val;
                int lchar = res.lchar;
                
                if ( ret.error )
                {
                    return ret;
                }
                begin = begin + lchar;
                length = length - lchar;
            }

            return ret;
        }
        
        public static int r4_nint ( float x )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R4_NINT returns the nearest integer to an R4.
            //
            //  Example:
            //
            //        X         R4_NINT
            //
            //      1.3         1
            //      1.4         1
            //      1.5         1 or 2
            //      1.6         2
            //      0.0         0
            //     -0.7        -1
            //     -1.1        -1
            //     -1.6        -2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 November 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, float X, the value.
            //
            //    Output, int R4_NINT, the nearest integer to X.
            //
        {
            int value = ( int ) ( Math.Abs( x ) + 0.5 );

            if ( x < 0.0 )
            {
                value = -value;
            }

            return value;
        }
 
        public static float r4_huge ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_HUGE returns a "huge" R4.
        //
        //  Discussion:
        //
        //    The value returned by this function is NOT required to be the
        //    maximum representable R4.  This value varies from machine to machine,
        //    from compiler to compiler, and may cause problems when being printed.
        //    We simply want a "very large" but non-infinite number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, float R4_HUGE, a "huge" R4 value.
        //
        {
            return 1.0E+30f;
        }
        
        public static float r4poly_value(int n, float[] a, float x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4POLY_VALUE evaluates a real polynomial.
        //
        //  Discussion:
        //
        //    For sanity's sake, the value of N indicates the NUMBER of 
        //    coefficients, or more precisely, the ORDER of the polynomial,
        //    rather than the DEGREE of the polynomial.  The two quantities
        //    differ by 1, but cause a great deal of confusion.
        //
        //    Given N and A, the form of the polynomial is:
        //
        //      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the polynomial.
        //
        //    Input, float A[N], the coefficients of the polynomial.
        //    A[0] is the constant term.
        //
        //    Input, float X, the point at which the polynomial is to be evaluated.
        //
        //    Output, float R4POLY_VALUE, the value of the polynomial at X.
        //
        {
            float value = 0.0f;

            for (int i = n - 1; 0 <= i; i--)
            {
                value = value * x + a[i];
            }

            return value;
        }
        
        public static float r4_sign ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN returns the sign of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose sign is desired.
//
//    Output, float R4_SIGN, the sign of X.
//
        {
            float value;

            if ( x < 0.0 )
            {
                value = -1.0f;
            }
            else
            {
                value = 1.0f;
            }
            return value;
        }

    }
}