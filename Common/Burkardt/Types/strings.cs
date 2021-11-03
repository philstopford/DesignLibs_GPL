using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static char[] s_substring ( char[] s, int a, int b )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_SUBSTRING returns a substring of a given string.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char *S, a pointer to a string.
            //
            //    Input, int A, B, the indices of the first and last character of S to copy.
            //    These are 1-based indices!  B should be 
            //
            //    Output, char *S_SUBSTRING, a pointer to the substring.
            //
        {
            int i;
            int j;
            char[] t;

            t = new char[b+2-a];

            j = 0;
            for ( i = a; i <= b; i++ )
            {
                t[j] = s[i-1];
                j = j + 1;
            }
            t[j] = '\0';

            return t;
        }

        public static void s_to_format(char[] s, ref int r, ref char code, ref int w, ref int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_TO_FORMAT reads a FORTRAN format from a string.
            //
            //  Discussion:
            //
            //    This routine will read as many characters as possible until it reaches
            //    the end of the string, or encounters a character which cannot be
            //    part of the format.  This routine is limited in its ability to
            //    recognize FORTRAN formats.  In particular, we are only expecting
            //    a single format specification, and cannot handle extra features
            //    such as 'ES' and 'EN' codes, '5X' spacing, and so on.
            //
            //    Legal input is:
            //
            //       0 nothing
            //       1 blanks
            //       2 optional '('
            //       3 blanks
            //       4 optional repeat factor R
            //       5 blanks
            //       6 CODE ( 'A', 'B', 'E', 'F', 'G', 'I', 'L', 'O', 'Z', '*' )
            //       7 blanks
            //       8 width W
            //       9 optional decimal point
            //      10 optional mantissa M
            //      11 blanks
            //      12 optional ')'
            //      13 blanks
            //
            //  Example:
            //
            //    S                 R   CODE   W    M
            //
            //    'I12              1   I      12   0
            //    'E8.0'            1   E       8   0
            //    'F10.5'           1   F      10   5
            //    '2G14.6'          2   G      14   6
            //    '*'               1   *      -1  -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char *S, the string containing the
            //    data to be read.  Reading will begin at position 1 and
            //    terminate at the end of the string, or when no more
            //    characters can be read.
            //
            //    Output, int *R, the repetition factor, which defaults to 1.
            //
            //    Output, char *CODE, the format code.
            //
            //    Output, int *W, the field width.
            //
            //    Output, int *M, the mantissa width.
            //
        {
            char c;
            int d;
            bool debug = true;
            int LEFT = 1;
            int paren_sum;
            int pos;
            int RIGHT = -1;
            int s_length;
            int state;

            state = 0;
            paren_sum = 0;
            pos = 0;
            s_length = s_len_trim(s);

            r = 0;
            w = 0;
            code = '?';
            m = 0;

            while (pos < s_length)
            {
                c = s[pos];
                pos = pos + 1;
                //
                //  BLANK character:
                //
                if (c == ' ')
                {
                    if (state == 4)
                    {
                        state = 5;
                    }
                    else if (state == 6)
                    {
                        state = 7;
                    }
                    else if (state == 10)
                    {
                        state = 11;
                    }
                    else if (state == 12)
                    {
                        state = 13;
                    }
                }
                //
                //  LEFT PAREN
                //
                else if (c == '(')
                {
                    if (state < 2)
                    {
                        paren_sum = paren_sum + LEFT;
                    }
                    else
                    {
                        if (debug)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                        }

                        state = -1;
                        break;
                    }
                }
                //
                //  DIGIT (R, F, or W)
                //
                else if (ch_is_digit(c))
                {
                    if (state <= 3)
                    {
                        state = 4;
                        r = ch_to_digit(c);
                    }
                    else if (state == 4)
                    {
                        d = ch_to_digit(c);
                        r = 10 * (r) + d;
                    }
                    else if (state == 6 || state == 7)
                    {
                        if (code == '*')
                        {
                            if (debug)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("S_TO_FORMAT - Fatal error!");
                                Console.WriteLine("  Current state = " + state + "");
                                Console.WriteLine("  Current code = '" + code + "'.");
                                Console.WriteLine("  Input character = '" + c + "'.");
                            }

                            state = -1;
                            break;
                        }

                        state = 8;
                        w = ch_to_digit(c);
                    }
                    else if (state == 8)
                    {
                        d = ch_to_digit(c);
                        w = 10 * (w) + d;
                    }
                    else if (state == 9)
                    {
                        state = 10;
                        m = ch_to_digit(c);
                    }
                    else if (state == 10)
                    {
                        d = ch_to_digit(c);
                        m = 10 * (m) + d;
                    }
                    else
                    {
                        if (debug)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                        }

                        state = -1;
                        break;
                    }
                }
                //
                //  DECIMAL POINT
                //
                else if (c == '.')
                {
                    if (state == 8)
                    {
                        state = 9;
                    }
                    else
                    {
                        if (debug)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                        }

                        state = -1;
                        break;
                    }
                }
                //
                //  RIGHT PAREN
                //
                else if (c == ')')
                {
                    paren_sum = paren_sum + RIGHT;

                    if (paren_sum != 0)
                    {
                        if (debug)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current paren sum = " + paren_sum + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                        }

                        state = -1;
                        break;
                    }

                    if (state == 6 && code == '*')
                    {
                        state = 12;
                    }
                    else if (6 <= state)
                    {
                        state = 12;
                    }
                    else
                    {
                        if (debug)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                        }

                        state = -1;
                        break;
                    }
                }
                //
                //  Code
                //
                else if (ch_is_format_code(c))
                {
                    if (state < 6)
                    {
                        state = 6;
                        code = c;
                    }
                    else
                    {
                        if (debug)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                        }

                        state = -1;
                        break;
                    }
                }
                //
                //  Unexpected character
                //
                else
                {
                    if (debug)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("S_TO_FORMAT - Fatal error!");
                        Console.WriteLine("  Current state = " + state + "");
                        Console.WriteLine("  Input character = '" + c + "'.");
                    }

                    state = -1;
                    break;
                }
            }

            if (paren_sum != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("S_TO_FORMAT - Fatal error!");
                Console.WriteLine("  Parentheses mismatch.");
                return;
            }

            if (state < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("S_TO_FORMAT - Fatal error!");
                Console.WriteLine("  Parsing error.");
                return;
            }

            if (r == 0)
            {
                r = 1;
            }

            return;
        }

        public static void s_trim ( ref char[] s )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_TRIM promotes the final null forward through trailing blanks.
            //
            //  Discussion:
            //
            //    What we're trying to say is that we reposition the null character
            //    so that trailing blanks are no longer visible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, char *S, the string to be trimmed.
            //
        {
            s = String.Join("",s).TrimEnd().ToCharArray();
        }
        
        public static int s_word_count ( string s )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_WORD_COUNT counts the number of "words" in a string.
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
            //    Input, string S, the string to be examined.
            //
            //    Output, int S_WORD_COUNT, the number of "words" in the string.
            //    Words are presumed to be separated by one or more blanks.
            //
        {
            bool blank;
            int char_count;
            int i;
            int word_count;

            word_count = 0;
            blank = true;

            char_count = s.Length;

            for ( i = 0; i < char_count; i++ )
            {
                if ( ( s[i] == ' ' ) )
                {
                    blank = true;
                }
                else if ( blank )
                {
                    word_count = word_count + 1;
                    blank = false;
                }
            }

            return word_count;
        }
        
        public static void s_word_extract_first ( string s, ref string s1, ref string s2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_WORD_EXTRACT_FIRST extracts the first word from a string.
        //
        //  Discussion:
        //
        //    A "word" is a string of characters terminated by a blank or
        //    the end of the string.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string S, the string.
        //
        //    Output, string &S1, the first word (initial blanks removed).
        //
        //    Output, string &S2, the remainder of the string, after removing
        //    the first word (initial blanks removed).
        //
        {
            int i;
            int mode;
            int s_len;

            s_len = s.Length;
            s1 = "";
            s2 = "";
            mode = 1;

            for ( i = 0; i < s_len; i++ )
            {
                if ( mode == 1 )
                {
                    if ( s[i] != ' ' )
                    {
                        mode = 2;
                    }
                }
                else if ( mode == 2 )
                {
                    if ( s[i] == ' ' )
                    {
                        mode = 3;
                    }
                }
                else if ( mode == 3 )
                {
                    if ( s[i] != ' ' )
                    {
                        mode = 4;
                    }
                }
                if ( mode == 2 )
                {
                    s1 = s1 + s[i];
                }
                else if ( mode == 4 )
                {
                    s2 = s2 + s[i];
                }
            }
        }
        public static bool s_eqi(string a, string b)
        {
            return a.ToLower() == b.ToLower();
        }

        public static int s_len_trim ( char[] s )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_LEN_TRIM returns the length of a string to the last nonblank.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char *S, a pointer to a string.
            //
            //    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
            //    If S_LEN_TRIM is 0, then the string is entirely blank.
            //
        {
            int n;
            int t;

            n =  ( s.Length );
            t = n  - 1;

            while ( 0 < n ) 
            {
                if ( s[t] != ' ' )
                {
                    return n;
                }
                t--;
                n--;
            }

            return n;
        }
        public static int s_len_trim(string line)
        {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_LEN_TRIM returns the length of a string to the last nonblank.
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
            //    Input, string S, a string.
            //
            //    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
            //    If S_LEN_TRIM is 0, then the string is entirely blank.
            //
            if (line == "")
            {
                return 0;
            }

            string tmp = line.TrimEnd();
            tmp = line.TrimStart();

            int index = tmp.Length;

            /*
            int index = line.Length - 1;
            while ((index >= 0) && (line[index] == ' '))
            {
            index--;
            }
            */
            return index;
        }

        public static bool ch_is_digit ( char c )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CH_IS_DIGIT returns TRUE if a character is a decimal digit.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 December 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char C, the character to be analyzed.
            //
            //    Output, bool CH_IS_DIGIT, is TRUE if C is a digit.
            //
        {
            if ( '0' <= c && c <= '9' )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        public static int ch_to_digit(char ch)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CH_TO_DIGIT returns the integer value of a base 10 digit.
            //
            //  Example:
            //
            //     CH  DIGIT
            //    ---  -----
            //    '0'    0
            //    '1'    1
            //    ...  ...
            //    '9'    9
            //    ' '    0
            //    'X'   -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
            //
            //    Output, int CH_TO_DIGIT, the corresponding integer value.  If the
            //    character was 'illegal', then DIGIT is -1.
            //
        {
            int digit;

            if ('0' <= ch && ch <= '9')
            {
                digit = ch - '0';
            }
            else if (ch == ' ')
            {
                digit = 0;
            }
            else
            {
                digit = -1;
            }

            return digit;
        }
        
        public static bool ch_is_format_code ( char c )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CH_IS_FORMAT_CODE returns TRUE if a character is a FORTRAN format code.
            //
            //  Discussion:
            //
            //    The format codes accepted here are not the only legal format
            //    codes in FORTRAN90.  However, they are more than sufficient
            //    for my needs!
            //
            //  Table:
            //
            //    A  Character
            //    B  Binary digits
            //    D  Real number, exponential representation
            //    E  Real number, exponential representation
            //    F  Real number, fixed point
            //    G  General format
            //    I  Integer
            //    L  Logical variable
            //    O  Octal digits
            //    Z  Hexadecimal digits
            //    *  Free format
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char C, the character to be analyzed.
            //
            //    Output, bool CH_IS_FORMAT_CODE, is TRUE if C is a FORTRAN format code.
            //
        {
            if ( ch_eqi ( c, 'A' ) ) 
            {
                return true;
            }
            else if ( ch_eqi ( c, 'B' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'D' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'E' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'F' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'G' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'I' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'L' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'O' ) )
            {
                return true;
            }
            else if ( ch_eqi ( c, 'Z' ) )
            {
                return true;
            }
            else if ( c == '*' )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
        public static bool ch_eqi ( char c1, char c2 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CH_EQI is true if two characters are equal, disregarding case.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 June 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, char C1, C2, the characters to compare.
            //
            //    Output, bool CH_EQI, is true if the two characters are equal,
            //    disregarding case.
            //
        {
            return ( c1.ToString().ToLower() == c2.ToString().ToLower() );
        }

        public static int s_to_i4(string st, ref int last, ref bool error )

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
        //    Output, int &LAST, the last character of S used to make IVAL.
        //
        //    Output, bool &ERROR is TRUE if an error occurred.
        //
        //    Output, int *S_TO_I4, the integer value read from the string.
        //    If the string is blank, then IVAL will be returned 0.
        //
        {
            char c;
            int i;
            int isgn;
            int istate;
            int ival;

            error = false;
            istate = 0;
            isgn = 1;
            i = 0;
            ival = 0;

            char[] s = st.ToCharArray();

            for (;;)
            {
                try
                {
                    c = s[i];
                }
                catch (Exception e)
                {
                    break;
                }
                i = i + 1;
                //
                //  Haven't read anything.
                //
                if (istate == 0)
                {
                    if (c == ' ')
                    {
                    }
                    else if (c == '-')
                    {
                        istate = 1;
                        isgn = -1;
                    }
                    else if (c == '+')
                    {
                        istate = 1;
                        isgn = +1;
                    }
                    else if ('0' <= c && c <= '9')
                    {
                        istate = 2;
                        ival = c - '0';
                    }
                    else
                    {
                        error = true;
                        return ival;
                    }
                }
                //
                //  Have read the sign, expecting digits.
                //
                else if (istate == 1)
                {
                    if (c == ' ')
                    {
                    }
                    else if ('0' <= c && c <= '9')
                    {
                        istate = 2;
                        ival = c - '0';
                    }
                    else
                    {
                        error = true;
                        return ival;
                    }
                }
                //
                //  Have read at least one digit, expecting more.
                //
                else if (istate == 2)
                {
                    if ('0' <= c && c <= '9')
                    {
                        ival = 10 * ival + c - '0';
                    }
                    else
                    {
                        ival = isgn * ival;
                        last = i - 1;
                        return ival;
                    }

                }
            }

            //
            //  If we read all the characters in the string, see if we're OK.
            //
            if (istate == 2)
            {
                ival = isgn * ival;
                last = s_len_trim(st);
            }
            else
            {
                error = true;
                last = 0;
            }

            return ival;
        }

        public static double s_to_r8(string s, ref int lchar, ref bool error )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    S_TO_R8 reads an R8 from a string.
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
        //    Output, int &LCHAR, the number of characters read from
        //    the string to form the number, including any terminating
        //    characters such as a trailing comma or blanks.
        //
        //    Output, bool &ERROR, is true if an error occurred.
        //
        //    Output, double S_TO_R8, the real value that was read from the string.
        //
        {
            char c;
            int ihave;
            int isgn;
            int iterm;
            int jbot;
            int jsgn;
            int jtop;
            int nchar;
            int ndig;
            double r;
            double rbot;
            double rexp;
            double rtop;
            int TAB = 9;

            nchar = s_len_trim(s);
            error = false;
            r = 0.0;
            lchar = -1;
            isgn = 1;
            rtop = 0.0;
            rbot = 1.0;
            jsgn = 1;
            jtop = 0;
            jbot = 1;
            ihave = 1;
            iterm = 0;

            for (;;)
            {
                c = s[lchar + 1];
                lchar = lchar + 1;
                //
                //  Blank or TAB character.
                //
                if (c == ' ' || c == TAB)
                {
                    if (ihave == 2)
                    {
                    }
                    else if (ihave == 6 || ihave == 7)
                    {
                        iterm = 1;
                    }
                    else if (1 < ihave)
                    {
                        ihave = 11;
                    }
                }
                //
                //  Comma.
                //
                else if (c == ',' || c == ';')
                {
                    if (ihave != 1)
                    {
                        iterm = 1;
                        ihave = 12;
                        lchar = lchar + 1;
                    }
                }
                //
                //  Minus sign.
                //
                else if (c == '-')
                {
                    if (ihave == 1)
                    {
                        ihave = 2;
                        isgn = -1;
                    }
                    else if (ihave == 6)
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
                else if (c == '+')
                {
                    if (ihave == 1)
                    {
                        ihave = 2;
                    }
                    else if (ihave == 6)
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
                else if (c == '.')
                {
                    if (ihave < 4)
                    {
                        ihave = 4;
                    }
                    else if (6 <= ihave && ihave <= 8)
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
                else if ((c == 'E') || (c == 'D'))
                {
                    if (ihave < 6)
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
                else if (ihave < 11 && '0' <= c && c <= '9')
                {
                    if (ihave <= 2)
                    {
                        ihave = 3;
                    }
                    else if (ihave == 4)
                    {
                        ihave = 5;
                    }
                    else if (ihave == 6 || ihave == 7)
                    {
                        ihave = 8;
                    }
                    else if (ihave == 9)
                    {
                        ihave = 10;
                    }

                    ndig = ch_to_digit(c);

                    if (ihave == 3)
                    {
                        rtop = 10.0 * rtop + (double) ndig;
                    }
                    else if (ihave == 5)
                    {
                        rtop = 10.0 * rtop + (double) ndig;
                        rbot = 10.0 * rbot;
                    }
                    else if (ihave == 8)
                    {
                        jtop = 10 * jtop + ndig;
                    }
                    else if (ihave == 10)
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
                if (iterm == 1 || nchar <= lchar + 1)
                {
                    break;
                }

            }

            //
            //  If we haven't seen a terminator, and we have examined the
            //  entire string, then we're done, and LCHAR is equal to NCHAR.
            //
            if (iterm != 1 && lchar + 1 == nchar)
            {
                lchar = nchar;
            }

            //
            //  Number seems to have terminated.  Have we got a legal number?
            //  Not if we terminated in states 1, 2, 6 or 7!
            //
            if (ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7)
            {
                error = true;
                return r;
            }

            //
            //  Number seems OK.  Form it.
            //
            if (jtop == 0)
            {
                rexp = 1.0;
            }
            else
            {
                if (jbot == 1)
                {
                    rexp = Math.Pow(10.0, jsgn * jtop);
                }
                else
                {
                    rexp = jsgn * jtop;
                    rexp = rexp / jbot;
                    rexp = Math.Pow(10.0, rexp);
                }

            }

            r = isgn * rexp * rtop / rbot;

            return r;
        }
        
        public static string s_escape_tex ( string s1 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    S_ESCAPE_TEX de-escapes TeX escape sequences.
            //
            //  Discussion:
            //
            //    In particular, every occurrence of the characters '\', '_',
            //    '^', '{' and '}' will be replaced by '\\', '\_', '\^',
            //    '\{' and '\}'.  A TeX interpreter, on seeing these character
            //    strings, is then likely to return the original characters.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string S1, the string to be de-escaped.
            //
            //    Output, string S_ESCAPE_TEX, a copy of the string,
            //    modified to avoid TeX escapes.
            //
        {
            char ch;
            int s1_length;
            int s1_pos;
            char[] s2;
            int s2_pos;
            string s3;
            int slash_count;

            s1_length = s1.Length;

            slash_count = 0;
            for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
            {
                ch = s1[s1_pos];

                if ( ch == '\\' ||
                     ch == '_' ||
                     ch == '^' ||
                     ch == '{' ||
                     ch == '}' )
                {
                    slash_count = slash_count + 1;
                }
            }
            s2 = new char[s1_length + slash_count + 1];
            //
            //  Now copy S1 into S2.
            //
            s1_pos = 0;
            s2_pos = 0;

            for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
            {
                ch = s1[s1_pos];

                if ( ch == '\\' ||
                     ch == '_' ||
                     ch == '^' ||
                     ch == '{' ||
                     ch == '}' )
                {
                    s2[s2_pos] = '\\';
                    s2_pos = s2_pos + 1;
                }

                s2[s2_pos] = ch;
                s2_pos = s2_pos + 1;
            }

            s2[s2_pos] = '\0';
            s2_pos = s2_pos + 1;

            s3 = s2.ToString();

            return s3;
        }

    }
}