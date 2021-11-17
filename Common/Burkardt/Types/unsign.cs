using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static uint nim_sum ( uint ui, uint uj )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NIM_SUM computes the Nim sum of two integers.
        //
        //  Discussion:
        //
        //    If K is the Nim sum of I and J, then each bit of K is the exclusive
        //    OR of the corresponding bits of I and J.
        //
        //  Example:
        //
        //     I     J     K       I_2        J_2         K_2
        //   ----  ----  ----  ----------  ----------  ----------
        //      0     0     0           0           0           0
        //      1     0     1           1           0           1
        //      1     1     0           1           1           0
        //      2     7     5          10         111         101
        //     11    28    23        1011       11100       10111
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, unsigned int UI, UJ, the integers to be Nim-summed.
        //
        //    Output, unsigned int NIM_SUM, the Nim sum of I and J.
        //
    {
        const int NBITS = 32;

        uint[] bvec1 = new uint[NBITS];
        uint[] bvec2 = new uint[NBITS];
        uint[] bvec3 = new uint[NBITS];

        ui4_to_ubvec ( ui, NBITS, ref bvec1 );

        ui4_to_ubvec ( uj, NBITS, ref bvec2 );

        ubvec_xor ( NBITS, bvec1, bvec2, bvec3 );

        uint value = ubvec_to_ui4 ( NBITS, bvec3 );

        return value;
    }
        
    public static bool ubvec_add(int n, uint[] bvec1, uint[] bvec2,
            ref uint[] bvec3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_ADD adds two unsigned binary vectors.
        //
        //  Discussion:
        //
        //    A UBVEC is a vector of binary digits representing an unsigned integer.  
        //
        //    UBVEC[N-1] contains the units digit, UBVEC[N-2]
        //    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
        //    so that printing the digits in order gives the binary form of the number.
        //
        //  Example:
        //
        //    N = 4
        //
        //     UBVEC1       +  UBVEC2       =  UBVEC3
        //
        //    ( 0 0 0 1 )   + ( 0 0 1 1 )   = ( 0 1 0 0 )
        //
        //      1           +   3           =   4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the length of the vectors.
        //
        //    Input, unsigned int BVEC1[N], BVEC2[N], the vectors to be added.
        //
        //    Output, unsigned int BVEC3[N], the sum of the two input vectors.
        //
        //    Output, bool UBVEC_ADD, is TRUE if an error occurred.
        //
    {
        const uint base_ = 2;
        int i;

        bool overflow = false;

        for (i = 0; i < n; i++)
        {
            bvec3[i] = bvec1[i] + bvec2[i];
        }

        for (i = n - 1; 0 <= i; i--)
        {
            while (base_ <= bvec3[i])
            {
                bvec3[i] -= base_;
                switch (i)
                {
                    case > 0:
                        bvec3[i - 1] += 1;
                        break;
                    default:
                        overflow = true;
                        break;
                }
            }
        }

        return overflow;
    }

    public static void ubvec_print(int n, uint[] bvec, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_PRINT prints a UBVEC, with an optional title.
        //
        //  Discussion:
        //
        //    A UBVEC is a vector of binary digits representing an unsigned integer.  
        //
        //    UBVEC[N-1] contains the units digit, UBVEC[N-2]
        //    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
        //    so that printing the digits in order gives the binary form of the number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, unsigned int BVEC[N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int ilo;

        switch (title.Length)
        {
            case > 0:
                Console.WriteLine("");
                Console.WriteLine(title + "");
                break;
        }

        for (ilo = 0; ilo < n; ilo += 70)
        {
            int ihi = Math.Min(ilo + 70 - 1, n - 1);
            string cout = "  ";

            int i;
            for (i = ilo; i <= ihi; i++)
            {
                cout += bvec[i];
            }

            Console.WriteLine(cout);
        }

    }

    public static uint ubvec_to_ui4(int n, uint[] bvec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_TO_UI4 makes an unsigned integer from an unsigned binary vector.
        //
        //  Discussion:
        //
        //    A UBVEC is a vector of binary digits representing an unsigned integer.  
        //
        //    UBVEC[N-1] contains the units digit, UBVEC[N-2]
        //    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
        //    so that printing the digits in order gives the binary form of the number.
        //
        //  Example:
        //
        //    N = 4
        //
        //         BVEC   binary  I
        //    ----------  -----  --
        //    1  2  3  4
        //    ----------
        //    0  0  0  1       1  1
        //    0  0  1  0      10  2
        //    0  0  1  1      11  3
        //    0  1  0  0     100  4
        //    1  0  0  1    1001  9
        //    1  1  1  1    1111 15
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //
        //    Input, unsigned int BVEC[N], the binary representation.
        //
        //    Input, unsigned int UBVEC_TO_UI4, the integer value.
        //
    {
        const uint base_ = 2;
        int i;

        uint ui4 = 0;
        for (i = 0; i < n; i++)
        {
            ui4 = base_ * ui4 + bvec[i];
        }

        return ui4;
    }

    public static void ubvec_xor(int n, uint[] bvec1, uint[] bvec2,
            uint[] bvec3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UBVEC_XOR computes the exclusive OR of two UBVECs.
        //
        //  Discussion:
        //
        //    A UBVEC is a vector of binary digits representing an unsigned integer.  
        //
        //    UBVEC[N-1] contains the units digit, UBVEC[N-2]
        //    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
        //    so that printing the digits in order gives the binary form of the number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the length of the vectors.
        //
        //    Input, unsigned int BVEC1[N], BVEC2[N], the binary vectors to be XOR'ed.
        //
        //    Input, unsigned int BVEC3[N], the exclusive OR of the two vectors.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            bvec3[i] = (bvec1[i] + bvec2[i]) % 2;
        }

    }

    public static void ui4_to_ubvec(uint ui4, int n, ref uint[] bvec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UI4_TO_UBVEC makes a unsigned binary vector from an integer.
        //
        //  Discussion:
        //
        //    A UBVEC is a vector of binary digits representing an unsigned integer.  
        //
        //    UBVEC[N-1] contains the units digit, UBVEC[N-2]
        //    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
        //    so that printing the digits in order gives the binary form of the number.
        //
        //    To guarantee that there will be enough space for any
        //    value of I, it would be necessary to set N = 32.
        //
        //  Example:
        //
        //    UI4      BVEC         binary
        //        0  1  2  3  4  5
        //        1  2  4  8 16 32
        //    --  ----------------  ------
        //     1  1  0  0  0  0  1       1
        //     2  0  1  0  0  1  0      10
        //     3  1  1  0  0  1  1      11
        //     4  0  0  0  1  0  0     100
        //     9  0  0  1  0  0  1    1001
        //    57  1  1  1  0  1  1  110111
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, unsigned int UI4, an integer to be represented.
        //
        //    Input, int N, the dimension of the vector.
        //
        //    Output, unsigned int BVEC[N], the binary representation.
        //
    {
        const uint base_ = 2;
        int i;

        for (i = n - 1; 0 <= i; i--)
        {
            bvec[i] = ui4 % base_;
            ui4 /= base_;
        }

    }
}