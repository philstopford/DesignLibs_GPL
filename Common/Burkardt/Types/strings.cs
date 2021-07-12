using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static bool s_eqi(string a, string b)
        {
            return a.ToLower() == b.ToLower();
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

        public static int s_to_i4(string s, ref int last, ref bool error )

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

            for (;;)
            {
                c = s[i];
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
                last = s_len_trim(s);
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