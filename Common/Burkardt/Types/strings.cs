namespace Burkardt.Types
{
    public static partial class typeMethods
    {
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
        
        public static int ch_to_digit ( char ch )
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

            if ( '0' <= ch && ch <= '9' )
            {
                digit = ch - '0';
            }
            else if ( ch == ' ' )
            {
                digit = 0;
            }
            else
            {
                digit = -1;
            }

            return digit;
        }        
        
    }
}