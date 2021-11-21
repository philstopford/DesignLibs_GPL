using System;

namespace Burkardt.FullertonFnLib;

public static partial class FullertonLib
{
    public static int i4_abs(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_ABS returns the absolute value of an I4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, an integer.
        //
        //    Output, int I4_ABS, the absolute value of the integer.
        //
    {
        int value = i switch
        {
            >= 0 => i,
            _ => -i
        };

        return value;
    }
        
    public static int i4_binom ( int n, int k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOM computes the binomial coefficient.
        //
        //  Discussion:
        //
        //    This is ACM algorithm 160 translated to Fortran.
        //
        //    It calculates the number of combinations of N things taken K at a time.
        //
        //  Modified:
        //
        //    01 April 2016
        //
        //  Author:
        //
        //    Bill Buckles, Matthew Lybanon
        //
        //  Reference:
        //
        //    Bill Buckles, Matthew Lybanon,
        //    Algorithm 515: Generation of a Vector from the Lexicographical Index,
        //    ACM Transactions on Mathematical Software,
        //    Volume 3, Number 2, June 1977, pages 180-182.
        //
        //  Parameters:
        //
        //    Input, int N, K, the parameters for the binomial 
        //    coefficient.
        //
        //    Output, int BINOM, the binomial coefficient.
        //
    {
        int i;

        int k1 = k;
        int p = n - k1;

        if ( k1 < p )
        {
            p = k1;
            k1 = n - p;
        }

        int r = p switch
        {
            0 => 1,
            _ => k1 + 1
        };

        for ( i = 2; i <= p; i++ )
        {
            r = r * ( k1 + i ) / i;
        }

        return r;
    }

    public static int i4_mach(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MACH returns integer machine constants.
        //
        //  Discussion:
        //
        //    Input/output unit numbers.
        //
        //      I4_MACH(1) = the standard input unit.
        //      I4_MACH(2) = the standard output unit.
        //      I4_MACH(3) = the standard punch unit.
        //      I4_MACH(4) = the standard error message unit.
        //
        //    Words.
        //
        //      I4_MACH(5) = the number of bits per integer storage unit.
        //      I4_MACH(6) = the number of characters per integer storage unit.
        //
        //    Integers.
        //
        //    Assume integers are represented in the S digit base A form:
        //
        //      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
        //
        //    where 0 <= X(1:S-1) < A.
        //
        //      I4_MACH(7) = A, the base.
        //      I4_MACH(8) = S, the number of base A digits.
        //      I4_MACH(9) = A^S-1, the largest integer.
        //
        //    Floating point numbers
        //
        //    Assume floating point numbers are represented in the T digit 
        //    base B form:
        //
        //      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
        //
        //    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
        //
        //      I4_MACH(10) = B, the base.
        //
        //    Single precision
        //
        //      I4_MACH(11) = T, the number of base B digits.
        //      I4_MACH(12) = EMIN, the smallest exponent E.
        //      I4_MACH(13) = EMAX, the largest exponent E.
        //
        //    Double precision
        //
        //      I4_MACH(14) = T, the number of base B digits.
        //      I4_MACH(15) = EMIN, the smallest exponent E.
        //      I4_MACH(16) = EMAX, the largest exponent E.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 April 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Phyllis Fox, Andrew Hall, Norman Schryer,
        //    Algorithm 528,
        //    Framework for a Portable Library,
        //    ACM Transactions on Mathematical Software,
        //    Volume 4, Number 2, June 1978, page 176-188.
        //
        //  Parameters:
        //
        //    Input, int I, chooses the parameter to be returned.
        //    1 <= I <= 16.
        //
        //    Output, int I4_MACH, the value of the chosen parameter.
        //
    {
        int value;

        switch (i)
        {
            case 1:
                value = 5;
                break;
            case 2:
                value = 6;
                break;
            case 3:
                value = 7;
                break;
            case 4:
                value = 6;
                break;
            case 5:
                value = 32;
                break;
            case 6:
                value = 4;
                break;
            case 7:
                value = 2;
                break;
            case 8:
                value = 31;
                break;
            case 9:
                value = 2147483647;
                break;
            case 10:
                value = 2;
                break;
            case 11:
                value = 24;
                break;
            case 12:
                value = -125;
                break;
            case 13:
                value = 128;
                break;
            case 14:
                value = 53;
                break;
            case 15:
                value = -1021;
                break;
            case 16:
                value = 1024;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("I4_MACH - Fatal error!");
                Console.WriteLine("  The input argument I is out of bounds.");
                Console.WriteLine("  Legal values satisfy 1 <= I <= 16.");
                Console.WriteLine("  I = " + i + "");
                value = 0;
                return 1;
        }

        return value;
    }

    public static int i4_max(int i1, int i2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MAX returns the maximum of two I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1, I2, are two integers to be compared.
        //
        //    Output, int I4_MAX, the larger of I1 and I2.
        //
    {
        int value = i2 < i1 ? i1 : i2;

        return value;
    }

    public static int i4_min(int i1, int i2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MIN returns the minimum of two I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I1, I2, two integers to be compared.
        //
        //    Output, int I4_MIN, the smaller of I1 and I2.
        //
    {
        int value = i1 < i2 ? i1 : i2;

        return value;
    }

    public static int i4_pow(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_POW returns the value of I^J.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the base and the power.  J should be nonnegative.
        //
        //    Output, int I4_POW, the value of I^J.
        //
    {
        int value;

        switch (j)
        {
            case < 0 when i == 1:
                value = 1;
                break;
            case < 0 when i == 0:
                Console.WriteLine("");
                Console.WriteLine("I4_POW - Fatal error!");
                Console.WriteLine("  I^J requested, with I = 0 and J negative.");
                return 1;
            case < 0:
                value = 0;
                break;
            case 0 when i == 0:
                Console.WriteLine("");
                Console.WriteLine("I4_POW - Fatal error!");
                Console.WriteLine("  I^J requested, with I = 0 and J = 0.");
                return 1;
            case 0:
                value = 1;
                break;
            case 1:
                value = i;
                break;
            default:
            {
                value = 1;
                int k;
                for (k = 1; k <= j; k++)
                {
                    value *= i;
                }

                break;
            }
        }

        return value;
    }
}