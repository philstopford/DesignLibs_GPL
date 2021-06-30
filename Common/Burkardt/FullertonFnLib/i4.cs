using System;

namespace Burkardt.FullertonFnLib
{
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
            int value;

            if (0 <= i)
            {
                value = i;
            }
            else
            {
                value = -i;
            }

            return value;
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

            if (i == 1)
            {
                value = 5;
            }
            else if (i == 2)
            {
                value = 6;
            }
            else if (i == 3)
            {
                value = 7;
            }
            else if (i == 4)
            {
                value = 6;
            }
            else if (i == 5)
            {
                value = 32;
            }
            else if (i == 6)
            {
                value = 4;
            }
            else if (i == 7)
            {
                value = 2;
            }
            else if (i == 8)
            {
                value = 31;
            }
            else if (i == 9)
            {
                value = 2147483647;
            }
            else if (i == 10)
            {
                value = 2;
            }
            else if (i == 11)
            {
                value = 24;
            }
            else if (i == 12)
            {
                value = -125;
            }
            else if (i == 13)
            {
                value = 128;
            }
            else if (i == 14)
            {
                value = 53;
            }
            else if (i == 15)
            {
                value = -1021;
            }
            else if (i == 16)
            {
                value = 1024;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("I4_MACH - Fatal error!");
                Console.WriteLine("  The input argument I is out of bounds.");
                Console.WriteLine("  Legal values satisfy 1 <= I <= 16.");
                Console.WriteLine("  I = " + i + "");
                value = 0;
                return (1);
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
            int value;

            if (i2 < i1)
            {
                value = i1;
            }
            else
            {
                value = i2;
            }

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
            int value;

            if (i1 < i2)
            {
                value = i1;
            }
            else
            {
                value = i2;
            }

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
            int k;
            int value;

            if (j < 0)
            {
                if (i == 1)
                {
                    value = 1;
                }
                else if (i == 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("I4_POW - Fatal error!");
                    Console.WriteLine("  I^J requested, with I = 0 and J negative.");
                    return (1);
                }
                else
                {
                    value = 0;
                }
            }
            else if (j == 0)
            {
                if (i == 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("I4_POW - Fatal error!");
                    Console.WriteLine("  I^J requested, with I = 0 and J = 0.");
                    return (1);
                }
                else
                {
                    value = 1;
                }
            }
            else if (j == 1)
            {
                value = i;
            }
            else
            {
                value = 1;
                for (k = 1; k <= j; k++)
                {
                    value = value * i;
                }
            }

            return value;
        }
    }
}