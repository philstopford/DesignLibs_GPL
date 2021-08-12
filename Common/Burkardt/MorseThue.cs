using Burkardt.Types;

namespace Burkardt
{
    public static class MorseThue
    {
        public static int morse_thue ( int i )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MORSE_THUE generates a Morse_Thue number.
            //
            //  Discussion:
            //
            //    The Morse_Thue sequence can be defined in a number of ways.
            //
            //    A) Start with the string containing the single letter '0'; then
            //       repeatedly apply the replacement rules '0' -> '01' and
            //       '1' -> '10' to the letters of the string.  The Morse_Thue sequence
            //       is the resulting letter sequence.
            //
            //    B) Starting with the string containing the single letter '0',
            //       repeatedly append the binary complement of the string to itself.
            //       Thus, '0' becomes '0' + '1' or '01', then '01' becomes
            //       '01' + '10', which becomes '0110' + '1001', and so on.
            //
            //    C) Starting with I = 0, the I-th Morse-Thue number is determined
            //       by taking the binary representation of I, adding the digits,
            //       and computing the remainder modulo 2.
            //
            //  Example:
            //
            //     I  binary   S
            //    --  ------  --
            //     0       0   0
            //     1       1   1
            //     2      10   1
            //     3      11   0
            //     4     100   1
            //     5     101   0
            //     6     110   0
            //     7     111   1
            //     8    1000   1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int I, the index of the Morse-Thue number.
            //
            //    Output, int MORSE_THUE, the Morse-Thue number of index I.
            //
        {
            int NBITS = 32;

            uint[] b = new uint[NBITS];
            uint s;
            uint ui;
            int value;
            //
            //  Expand I into binary form.
            //
            ui = (uint)i;

            typeMethods.ui4_to_ubvec ( ui, NBITS, ref b );
            //
            //  Sum the 1's in the binary representation.
            //
            s = 0;
            for ( i = 0; i < NBITS; i++ )
            {
                s = s + b[i];
            }
            //
            //  Take the value modulo 2.
            //
            s = s % 2;

            value = ( int ) s;

            return value;
        }
    }
}