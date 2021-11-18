using System;
using System.Globalization;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public class WordData
    {
        public int lenc;
        public int next;

    }
    public static string word_next_read(ref WordData data, string s, ref bool done)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WORD_NEXT_READ "reads" words from a string, one at a time.
        //
        //  Discussion:
        //
        //    This routine was written to process tokens in a file.
        //    A token is considered to be an alphanumeric string delimited
        //    by whitespace, or any of various "brackets".
        //
        //    The following characters are considered to be a single word,
        //    whether surrounded by spaces or not:
        //
        //      " ( ) { } [ ]
        //
        //    Also, if there is a trailing comma on the word, it is stripped off.
        //    This is to facilitate the reading of lists.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string S, a string, presumably containing words
        //    separated by spaces.
        //
        //    Input/output, bool *DONE.
        //    On input with a fresh string, set DONE to TRUE.
        //    On output, the routine sets DONE:
        //      FALSE if another word was read,
        //      TRUE if no more words could be read.
        //
        //    Output, string WORD_NEXT_READ.
        //    If DONE is FALSE, then WORD contains the "next" word read.
        //    If DONE is TRUE, then WORD is NULL, because there was no more to read.
        //
    {
        int i;
        int j;
        const char TAB = '9';
        string word;
        char[] word_chstar;
        switch (done)
        {
            //
            //  We "remember" LENC and NEXT from the previous call.
            //
            //  An input value of DONE = TRUE signals a new line of text to examine.
            //
            case true:
            {
                data.next = 0;
                done = false;
                data.lenc = s.Length;
                switch (data.lenc)
                {
                    case <= 0:
                        done = true;
                        word = "\n";
                        return word;
                }

                break;
            }
        }

        //
        //  Beginning at index NEXT, search the string for the next nonblank,
        //  which signals the beginning of a word.
        //
        int ilo = data.next;
        //
        //  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
        //
        for (;;)
        {
            if (data.lenc < ilo)
            {
                word = "\n";
                done = true;
                data.next = data.lenc + 1;
                return word;
            }

            //
            //  If the current character is blank, skip to the next one.
            //
            if (s[ilo] != ' ' && s[ilo] != TAB)
            {
                break;
            }

            ilo += 1;
        }

        switch (s[ilo])
        {
            //
            //  ILO is the index of the next nonblank character in the string.
            //
            //  If this initial nonblank is a special character,
            //  then that's the whole word as far as we're concerned,
            //  so return immediately.
            //
            case '"':
                word = "\"\"";
                data.next = ilo + 1;
                return word;
            case '(':
                word = "(";
                data.next = ilo + 1;
                return word;
            case ')':
                word = ")";
                data.next = ilo + 1;
                return word;
            case '{':
                word = "{";
                data.next = ilo + 1;
                return word;
            case '}':
                word = "}";
                data.next = ilo + 1;
                return word;
            case '[':
                word = "[";
                data.next = ilo + 1;
                return word;
            case ']':
                word = "]";
                data.next = ilo + 1;
                return word;
        }

        //
        //  Now search for the last contiguous character that is not a
        //  blank, TAB, or special character.
        //
        data.next = ilo + 1;

        while (data.next <= data.lenc)
        {
            if (s[data.next] == ' ')
            {
                break;
            }

            if (s[data.next] == TAB)
            {
                break;
            }
            if (s[data.next] == '"')
            {
                break;
            }
            if (s[data.next] == '(')
            {
                break;
            }
            if (s[data.next] == ')')
            {
                break;
            }
            if (s[data.next] == '{')
            {
                break;
            }
            if (s[data.next] == '}')
            {
                break;
            }
            if (s[data.next] == '[')
            {
                break;
            }
            if (s[data.next] == ']')
            {
                break;
            }

            data.next += 1;
        }

        switch (s[data.next - 1])
        {
            //
            //  Allocate WORD, copy characters, and return.
            //
            case ',':
            {
                word_chstar = new char[data.next - ilo];
                i = 0;
                for (j = ilo; j <= data.next - 2; j++)
                {
                    word_chstar[i] = s[j];
                    i += 1;
                }

                word_chstar[i] = '\0';
                word = new string( word_chstar );
                break;
            }
            default:
            {
                word_chstar = new char[data.next + 1 - ilo];
                i = 0;
                for (j = ilo; j <= data.next - 1; j++)
                {
                    word_chstar[i] = s[j];
                    i += 1;
                }

                word_chstar[i] = '\0';
                word = new string(word_chstar);
                break;
            }
        }

        return word;
    }

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

        char[] t = new char[b+2-a];

        int j = 0;
        for ( i = a; i <= b; i++ )
        {
            t[j] = s[i-1];
            j += 1;
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
        const bool debug = true;
        const int LEFT = 1;
        const int RIGHT = -1;

        int state = 0;
        int paren_sum = 0;
        int pos = 0;
        int s_length = s_len_trim(s);

        r = 0;
        w = 0;
        code = '?';
        m = 0;

        while (pos < s_length)
        {
            char c = s[pos];
            pos += 1;
            //
            //  BLANK character:
            //
            if (c == ' ')
            {
                state = state switch
                {
                    4 => 5,
                    6 => 7,
                    10 => 11,
                    12 => 13,
                    _ => state
                };
            }
            //
            //  LEFT PAREN
            //
            else if (c == '(')
            {
                if (state < 2)
                {
                    paren_sum += LEFT;
                }
                else
                {
                    switch (debug)
                    {
                        case true:
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                            break;
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
                else
                {
                    int d;
                    if (state == 4)
                    {
                        d = ch_to_digit(c);
                        r = 10 * r + d;
                    }
                    else if (state is 6 or 7)
                    {
                        if (code == '*')
                        {
                            switch (debug)
                            {
                                case true:
                                    Console.WriteLine("");
                                    Console.WriteLine("S_TO_FORMAT - Fatal error!");
                                    Console.WriteLine("  Current state = " + state + "");
                                    Console.WriteLine("  Current code = '" + code + "'.");
                                    Console.WriteLine("  Input character = '" + c + "'.");
                                    break;
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
                        w = 10 * w + d;
                    }
                    else if (state == 9)
                    {
                        state = 10;
                        m = ch_to_digit(c);
                    }
                    else if (state == 10)
                    {
                        d = ch_to_digit(c);
                        m = 10 * m + d;
                    }
                    else
                    {
                        switch (debug)
                        {
                            case true:
                                Console.WriteLine("");
                                Console.WriteLine("S_TO_FORMAT - Fatal error!");
                                Console.WriteLine("  Current state = " + state + "");
                                Console.WriteLine("  Input character = '" + c + "'.");
                                break;
                        }

                        state = -1;
                        break;
                    }
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
                    switch (debug)
                    {
                        case true:
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                            break;
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
                paren_sum += RIGHT;

                if (paren_sum != 0)
                {
                    switch (debug)
                    {
                        case true:
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current paren sum = " + paren_sum + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                            break;
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
                    switch (debug)
                    {
                        case true:
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                            break;
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
                    switch (debug)
                    {
                        case true:
                            Console.WriteLine("");
                            Console.WriteLine("S_TO_FORMAT - Fatal error!");
                            Console.WriteLine("  Current state = " + state + "");
                            Console.WriteLine("  Input character = '" + c + "'.");
                            break;
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
                switch (debug)
                {
                    case true:
                        Console.WriteLine("");
                        Console.WriteLine("S_TO_FORMAT - Fatal error!");
                        Console.WriteLine("  Current state = " + state + "");
                        Console.WriteLine("  Input character = '" + c + "'.");
                        break;
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

        switch (state)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("S_TO_FORMAT - Fatal error!");
                Console.WriteLine("  Parsing error.");
                return;
        }

        r = r switch
        {
            0 => 1,
            _ => r
        };
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
        s = string.Join("",s).TrimEnd().ToCharArray();
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
        int i;

        int word_count = 0;
        bool blank = true;

        int char_count = s.Length;

        for ( i = 0; i < char_count; i++ )
        {
            switch (s[i])
            {
                case ' ':
                    blank = true;
                    break;
                default:
                {
                    switch (blank)
                    {
                        case true:
                            word_count += 1;
                            blank = false;
                            break;
                    }

                    break;
                }
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

        int s_len = s.Length;
        s1 = "";
        s2 = "";
        int mode = 1;

        for ( i = 0; i < s_len; i++ )
        {
            switch (mode)
            {
                case 1:
                {
                    if ( s[i] != ' ' )
                    {
                        mode = 2;
                    }

                    break;
                }
                case 2:
                {
                    mode = s[i] switch
                    {
                        ' ' => 3,
                        _ => mode
                    };

                    break;
                }
                case 3:
                {
                    if ( s[i] != ' ' )
                    {
                        mode = 4;
                    }

                    break;
                }
            }

            switch (mode)
            {
                case 2:
                    s1 += s[i];
                    break;
                case 4:
                    s2 += s[i];
                    break;
            }
        }
    }
    public static bool s_eqi(string a, string b)
    {
        return string.Equals(a, b, StringComparison.CurrentCultureIgnoreCase);
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
        int n = s.Length;
        int t = n  - 1;

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
        switch (line)
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
            case "":
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
        return c switch
        {
            >= '0' and <= '9' => true,
            _ => false
        };
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
        int digit = ch switch
        {
            >= '0' and <= '9' => ch - '0',
            ' ' => 0,
            _ => -1
        };

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

        if ( ch_eqi ( c, 'B' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'D' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'E' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'F' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'G' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'I' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'L' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'O' ) )
        {
            return true;
        }
        if ( ch_eqi ( c, 'Z' ) )
        {
            return true;
        }

        return c switch
        {
            '*' => true,
            _ => false
        };
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
        return string.Equals(c1.ToString(CultureInfo.InvariantCulture), c2.ToString(CultureInfo.InvariantCulture), StringComparison.CurrentCultureIgnoreCase);
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
        error = false;
        int istate = 0;
        int isgn = 1;
        int i = 0;
        int ival = 0;

        char[] s = st.ToCharArray();

        for (;;)
        {
            char c;
            try
            {
                c = s[i];
            }
            catch (Exception)
            {
                break;
            }
            i += 1;
            switch (istate)
            {
                //
                //  Haven't read anything.
                //
                case 0 when c == ' ':
                    break;
                case 0 when c == '-':
                    istate = 1;
                    isgn = -1;
                    break;
                case 0 when c == '+':
                    istate = 1;
                    isgn = +1;
                    break;
                case 0 when c is >= '0' and <= '9':
                    istate = 2;
                    ival = c - '0';
                    break;
                case 0:
                    error = true;
                    return ival;
                //
                //  Have read the sign, expecting digits.
                //
                case 1 when c == ' ':
                    break;
                case 1 when c is >= '0' and <= '9':
                    istate = 2;
                    ival = c - '0';
                    break;
                case 1:
                    error = true;
                    return ival;
                //
                //  Have read at least one digit, expecting more.
                //
                case 2 when c is >= '0' and <= '9':
                    ival = 10 * ival + c - '0';
                    break;
                case 2:
                    ival = isgn * ival;
                    last = i - 1;
                    return ival;
            }
        }

        switch (istate)
        {
            //
            //  If we read all the characters in the string, see if we're OK.
            //
            case 2:
                ival = isgn * ival;
                last = s_len_trim(st);
                break;
            default:
                error = true;
                last = 0;
                break;
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
        double rexp;
        const int TAB = 9;

        int nchar = s_len_trim(s);
        error = false;
        double r = 0.0;
        lchar = -1;
        int isgn = 1;
        double rtop = 0.0;
        double rbot = 1.0;
        int jsgn = 1;
        int jtop = 0;
        int jbot = 1;
        int ihave = 1;
        int iterm = 0;

        for (;;)
        {
            char c = s[lchar + 1];
            lchar += 1;
            //
            //  Blank or TAB character.
            //
            if (c == ' ' || c == TAB)
            {
                switch (ihave)
                {
                    case 2:
                        break;
                    case 6:
                    case 7:
                        iterm = 1;
                        break;
                    case > 1:
                        ihave = 11;
                        break;
                }
            }
            else
            {
                switch (c)
                {
                    //
                    //  Comma.
                    //
                    case ',':
                    case ';':
                    {
                        if (ihave != 1)
                        {
                            iterm = 1;
                            ihave = 12;
                            lchar += 1;
                        }

                        break;
                    }
                    //
                    //  Minus sign.
                    //
                    case '-' when ihave == 1:
                        ihave = 2;
                        isgn = -1;
                        break;
                    case '-' when ihave == 6:
                        ihave = 7;
                        jsgn = -1;
                        break;
                    case '-':
                        iterm = 1;
                        break;
                    //
                    //  Plus sign.
                    //
                    case '+' when ihave == 1:
                        ihave = 2;
                        break;
                    case '+' when ihave == 6:
                        ihave = 7;
                        break;
                    case '+':
                        iterm = 1;
                        break;
                    //
                    //  Decimal point.
                    //
                    case '.' when ihave < 4:
                        ihave = 4;
                        break;
                    case '.' when ihave is >= 6 and <= 8:
                        ihave = 9;
                        break;
                    case '.':
                        iterm = 1;
                        break;
                    //
                    //  Exponent marker.
                    //
                    case 'E':
                    case 'D':
                    {
                        switch (ihave)
                        {
                            case < 6:
                                ihave = 6;
                                break;
                            default:
                                iterm = 1;
                                break;
                        }

                        break;
                    }
                    //
                    default:
                    {
                        switch (ihave)
                        {
                            case < 11 when c is >= '0' and <= '9':
                                switch (ihave)
                                {
                                    case <= 2:
                                        ihave = 3;
                                        break;
                                    case 4:
                                        ihave = 5;
                                        break;
                                    case 6:
                                    case 7:
                                        ihave = 8;
                                        break;
                                    case 9:
                                        ihave = 10;
                                        break;
                                }

                                int ndig = ch_to_digit(c);

                                switch (ihave)
                                {
                                    case 3:
                                        rtop = 10.0 * rtop + ndig;
                                        break;
                                    case 5:
                                        rtop = 10.0 * rtop + ndig;
                                        rbot = 10.0 * rbot;
                                        break;
                                    case 8:
                                        jtop = 10 * jtop + ndig;
                                        break;
                                    case 10:
                                        jtop = 10 * jtop + ndig;
                                        jbot = 10 * jbot;
                                        break;
                                }

                                break;
                            //
                            default:
                                iterm = 1;
                                break;
                        }

                        break;
                    }
                }
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

        switch (ihave)
        {
            //
            //  Number seems to have terminated.  Have we got a legal number?
            //  Not if we terminated in states 1, 2, 6 or 7!
            //
            case 1:
            case 2:
            case 6:
            case 7:
                error = true;
                return r;
        }

        switch (jtop)
        {
            //
            //  Number seems OK.  Form it.
            //
            case 0:
                rexp = 1.0;
                break;
            default:
            {
                switch (jbot)
                {
                    case 1:
                        rexp = Math.Pow(10.0, jsgn * jtop);
                        break;
                    default:
                        rexp = jsgn * jtop;
                        rexp /= jbot;
                        rexp = Math.Pow(10.0, rexp);
                        break;
                }

                break;
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
        int s1_pos;

        int s1_length = s1.Length;

        int slash_count = 0;
        for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
        {
            ch = s1[s1_pos];

            switch (ch)
            {
                case '\\':
                case '_':
                case '^':
                case '{':
                case '}':
                    slash_count += 1;
                    break;
            }
        }
        char[] s2 = new char[s1_length + slash_count + 1];
        //
        //  Now copy S1 into S2.
        //
        s1_pos = 0;
        int s2_pos = 0;

        for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
        {
            ch = s1[s1_pos];

            switch (ch)
            {
                case '\\':
                case '_':
                case '^':
                case '{':
                case '}':
                    s2[s2_pos] = '\\';
                    s2_pos += 1;
                    break;
            }

            s2[s2_pos] = ch;
            s2_pos += 1;
        }

        s2[s2_pos] = '\0';

        string s3 = s2.ToString();

        return s3;
    }

}