using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void dec_add(int mantissa1, int exponent1, int mantissa2, int exponent2,
            int dec_digit, ref int mantissa, ref int exponent )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEC_ADD adds two decimal quantities.
        //
        //  Discussion:
        //
        //    A decimal value is represented as MANTISSA * 10^EXPONENT.
        //
        //    The routine computes
        //
        //      MANTISSA * 10^EXPONENT = MANTISSA1 * 10^EXPONENT1 + MANTISSA2 * 10^EXPONENT2
        //
        //    using DEC_DIGIT arithmetic.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MANTISSA1, EXPONENT1, the first number to be added.
        //
        //    Input, int MANTISSA2, EXPONENT2, the second number to be added.
        //
        //    Input, int DEC_DIGIT, the number of decimal digits.
        //
        //    Output, int &MANTISSA, &EXPONENT, the sum.
        //
        {
            if (mantissa1 == 0)
            {
                mantissa = mantissa2;
                exponent = exponent2;
                dec_round(mantissa, exponent, dec_digit, ref mantissa, ref exponent);
                return;
            }
            else if (mantissa2 == 0)
            {
                mantissa = mantissa1;
                exponent = exponent1;
                dec_round(mantissa, exponent, dec_digit, ref mantissa, ref exponent);
                return;
            }

            //
            //  Line up the exponents.
            //
            if (exponent1 < exponent2)
            {
                mantissa2 = mantissa2 * (int)Math.Pow(10, exponent2 - exponent1);
                exponent2 = exponent1;
                mantissa = mantissa1 + mantissa2;
                exponent = exponent1;
            }
            else if (exponent1 == exponent2)
            {
                mantissa = mantissa1 + mantissa2;
                exponent = exponent1;
            }
            else if (exponent2 < exponent1)
            {
                mantissa1 = mantissa1 * (int)Math.Pow(10, exponent1 - exponent2);
                exponent1 = exponent2;
                mantissa = mantissa1 + mantissa2;
                exponent = exponent2;
            }

            //
            //  Clean up the result.
            //
            dec_round(mantissa, exponent, dec_digit, ref mantissa, ref exponent);
        }

        public static void dec_div(int mantissa1, int exponent1, int mantissa2, int exponent2,
            int dec_digit, ref int mantissa, ref int exponent, ref bool error )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEC_DIV divides two decimal values.
        //
        //  Discussion:
        //
        //    A decimal value is represented as MANTISSA * 10^EXPONENT.
        //
        //    The routine computes
        //
        //      MANTISSA * 10^EXPONENT 
        //      = (MANTISSA1 * 10^EXPONENT1) / (MANTISSA2 * 10^EXPONENT2)
        //      = (MANTISSA1/MANTISSA2) * 10^(EXPONENT1-EXPONENT2)
        //
        //    while avoiding integer overflow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MANTISSA1, EXPONENT1, the numerator.
        //
        //    Input, int MANTISSA2, EXPONENT2, the denominator.
        //
        //    Input, int DEC_DIGIT, the number of decimal digits.
        //
        //    Output, int &MANTISSA, &EXPONENT, the result.
        //
        //    Output, bool &ERROR is true if an error occurred.
        //
        {
            double dval;
            int exponent3 = 0;
            int mantissa3 = 0;

            error = false;
            //
            //  First special case, top fraction is 0.
            //
            if (mantissa1 == 0)
            {
                mantissa = 0;
                exponent = 0;
                return;
            }

            //
            //  First error, bottom of fraction is 0.
            //
            if (mantissa2 == 0)
            {
                error = true;
                mantissa = 0;
                exponent = 0;
                return;
            }

            //
            //  Second special case, result is 1.
            //
            if (mantissa1 == mantissa2 && exponent1 == exponent2)
            {
                mantissa = 1;
                exponent = 0;
                return;
            }

            //
            //  Third special case, result is power of 10.
            //
            if (mantissa1 == mantissa2)
            {
                mantissa = 1;
                exponent = exponent1 - exponent2;
                return;
            }

            //
            //  Fourth special case: MANTISSA1/MANTISSA2 is exact.
            //
            if ((mantissa1 / mantissa2) * mantissa2 == mantissa1)
            {
                mantissa = mantissa1 / mantissa2;
                exponent = exponent1 - exponent2;
                return;
            }

            //
            //  General case.
            //
            dval = (double)mantissa1 / (double)mantissa2;

            r8_to_dec(dval, dec_digit, ref mantissa3, ref exponent3);

            mantissa = mantissa3;
            exponent = exponent3 + exponent1 - exponent2;
        }

        public static void dec_mul(int mantissa1, int exponent1, int mantissa2, int exponent2,
            int dec_digit, ref int mantissa, ref int exponent )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEC_MUL multiplies two decimals.
        //
        //  Discussion:
        //
        //    A decimal value is represented as MANTISSA * 10^EXPONENT.
        //
        //    The routine computes
        //
        //      MANTISSA * 10^EXPONENT = (MANTISSA1 * 10^EXPONENT1) 
        //                             * (MANTISSA2 * 10^EXPONENT2)
        //                      = (MANTISSA1*MANTISSA2) * 10^(EXPONENT1+EXPONENT2)
        //
        //    while avoiding integer overflow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MANTISSA1, EXPONENT1, the first multiplier.
        //
        //    Input, int MANTISSA2, EXPONENT2, the second multiplier.
        //
        //    Input, int DEC_DIGIT, the number of decimal digits.
        //
        //    Output, int &MANTISSA, &EXPONENT, the product.
        //
        {
            double dval;
            int exponent3 = 0;
            int mantissa3 = 0;
            double temp;
            //
            //  The result is zero if either MANTISSA1 or MANTISSA2 is zero.
            //
            if (mantissa1 == 0 || mantissa2 == 0)
            {
                mantissa = 0;
                exponent = 0;
                return;
            }

            //
            //  The result is simple if either MANTISSA1 or MANTISSA2 is one.
            //
            if (Math.Abs(mantissa1) == 1 || Math.Abs(mantissa2) == 1)
            {
                mantissa = mantissa1 * mantissa2;
                exponent = exponent1 + exponent2;
                return;
            }

            temp = Math.Log((double)Math.Abs(mantissa1))
                   + Math.Log((double)Math.Abs(mantissa2));

            if (temp < Math.Log(typeMethods.i4_huge()))
            {
                mantissa = mantissa1 * mantissa2;
                exponent = exponent1 + exponent2;
            }
            else
            {
                dval = (double)mantissa1 * (double)mantissa2;

                r8_to_dec(dval, dec_digit, ref mantissa3, ref exponent3);

                mantissa = mantissa3;
                exponent = exponent3 + (exponent1 + exponent2);
            }

            dec_round(mantissa, exponent, dec_digit, ref mantissa, ref exponent);
        }

        public static void dec_round(int mantissa1, int exponent1, int dec_digit,
            ref int mantissa2, ref int exponent2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEC_ROUND rounds a decimal fraction to a given number of digits.
        //
        //  Discussion:
        //
        //    A decimal value is represented as MANTISSA * 10^EXPONENT.
        //
        //    The routine takes an arbitrary decimal value and makes sure that MANTISSA
        //    has no more than DEC_DIGIT digits.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MANTISSA1, EXPONENT1, the coefficient and exponent
        //    of a decimal fraction to be rounded.
        //
        //    Input, int DEC_DIGIT, the number of decimal digits.
        //
        //    Output, int &MANTISSA2, &EXPONENT2, the rounded coefficient and exponent
        //    of a decimal fraction.  MANTISSA2 has no more than
        //    DEC_DIGIT decimal digits.
        //
        {
            int i;
            int limit;
            int sgn;

            mantissa2 = mantissa1;
            exponent2 = exponent1;
            //
            //  Watch out for the special case of 0.
            //
            if (mantissa2 == 0)
            {
                exponent2 = 0;
                return;
            }

            //
            //  Record the sign of MANTISSA.
            //
            sgn = 1;
            if (mantissa2 < 0)
            {
                mantissa2 = -(mantissa2);
                sgn = -sgn;
            }

            //
            //  If MANTISSA is too big, knock it down.
            //
            limit = 1;
            for (i = 1; i <= dec_digit; i++)
            {
                limit = limit * 10;
            }

            while (limit <= Math.Abs(mantissa2))
            {
                mantissa2 = (mantissa2 + 5) / 10;
                exponent2 = exponent2 + 1;
            }

            //
            //  Absorb trailing 0's into the exponent.
            //
            if (0 < mantissa2)
            {
                while ((mantissa2 / 10) * 10 == mantissa2)
                {
                    mantissa2 = mantissa2 / 10;
                    exponent2 = exponent2 + 1;
                }
            }

            mantissa2 = sgn * (mantissa2);
        }

        public static double dec_to_r8(int mantissa, int exponent)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEC_TO_R8 converts a decimal value to an R8.
            //
            //  Discussion:
            //
            //    A decimal value is represented as MANTISSA * 10^EXPONENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int MANTISSA, EXPONENT, the coefficient and exponent
            //    of the decimal value.
            //
            //    Output, double DEC_TO_R8, the real value of the decimal.
            //
        {
            double value;

            value = mantissa * Math.Pow(10.0, exponent);

            return value;
        }

        public static void dec_to_rat(int mantissa, int exponent, ref int rat_top, ref int rat_bot )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEC_TO_RAT converts a decimal to a rational representation.
        //
        //  Discussion:
        //
        //    A decimal value is represented as MANTISSA * 10^EXPONENT.
        //
        //    A rational value is represented by RAT_TOP / RAT_BOT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MANTISSA, EXPONENT, the decimal number.
        //
        //    Output, int &RAT_TOP, &RAT_BOT, the rational value.
        //
        {
            int gcd;
            int i;

            if (exponent == 0)
            {
                rat_top = mantissa;
                rat_bot = 1;
            }
            else if (0 < exponent)
            {
                rat_top = mantissa;
                for (i = 1; i <= exponent; i++)
                {
                    rat_top = rat_top * 10;
                }

                rat_bot = 1;
            }
            else
            {
                rat_top = mantissa;
                rat_bot = 1;
                for (i = 1; i <= -exponent; i++)
                {
                    rat_bot = rat_bot * 10;
                }

                gcd = i4_gcd(rat_top, rat_bot);
                rat_top = rat_top / gcd;
                rat_bot = rat_bot / gcd;
            }
        }

        public static string dec_to_s(int mantissa, int exponent)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEC_TO_S converts a decimal number to a string.
            //
            //  Discussion:
            //
            //    This is a draft version that is NOT WORKING YET.
            //
            //  Example:
            //
            //    Mantissa  Exponent Representation:
            //
            //         523      -1              5.23
            //         134       2          13400
            //           0      10              0
            //
            //      123456       3      123456000
            //      123456       2       12345600
            //      123456       1        1234560
            //      123456       0         123456
            //      123456      -1          12345.6
            //      123456      -2           1234.56
            //      123456      -3            123.456
            //      123456      -4             12.3456
            //      123456      -5             1.23456
            //      123456      -6             0.123456
            //      123456      -7             0.0123456
            //      123456      -8             0.00123456
            //      123456      -9             0.000123456
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 July 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int MANTISSA, EXPONENT, integers which represent the decimal.
            //
            //    Output, char *S, the representation of the value.
            //
        {
            int digit;
            int i;
            int mantissa_exponent;
            int mantissa_exponent_copy;
            int mantissa_10;
            int pos;
            char[] s;
            int s_length;

            s_length = dec_width(mantissa, exponent) + 1;

            s = new char[s_length];

            for (i = 0; i < s_length - 1; i++)
            {
                s[i] = '0';
            }

            s[s_length - 1] = '\0';

            if (mantissa == 0)
            {
                return String.Join("",s);
            }

            pos = 0;

            if (mantissa < 0)
            {
                s[pos] = '-';
                pos = pos + 1;
                mantissa = -mantissa;
            }

            mantissa_exponent = (int)Math.Log10(mantissa) + 1;
            mantissa_10 = (int)Math.Pow(10, mantissa_exponent - 1);
            //
            //  Are the next characters "0."?
            //
            if (mantissa_exponent + exponent <= 0)
            {
                s[pos] = '0';
                pos = pos + 1;
                s[pos] = '.';
                pos = pos + 1;

                for (i = mantissa_exponent + exponent; i < 0; i++)
                {
                    s[pos] = '0';
                    pos = pos + 1;
                }
            }

            //
            //  Print the digits of the mantissa.
            //
            mantissa_exponent_copy = mantissa_exponent;

            for (i = 0; i < mantissa_exponent; i++)
            {
                digit = mantissa / mantissa_10;
                mantissa = mantissa % mantissa_10;
                s[pos] = (digit).ToString()[0];
                pos = pos + 1;
                mantissa_10 = mantissa_10 / 10;
                mantissa_exponent_copy = mantissa_exponent_copy - 1;
                if (exponent < 0)
                {
                    if (mantissa_exponent_copy + exponent == 0)
                    {
                        s[pos] = '.';
                        pos = pos + 1;
                    }
                }
            }
            //
            //  Print any trailing zeros.
            //

            if (0 < exponent)
            {
                for (i = exponent; 0 < i; i--)
                {
                    s[pos] = '0';
                    pos = pos + 1;
                }
            }

            return String.Join("",s);
        }

        public static int dec_width(int mantissa, int exponent)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEC_WIDTH returns the "width" of a decimal number.
            //
            //  Discussion:
            //
            //    A decimal value is represented as MANTISSA * 10^EXPONENT.
            //
            //    The "width" of a decimal number is the number of characters
            //    required to print it.
            //
            //  Example:
            //
            //    Mantissa  Exponent Width  Representation:
            //
            //         523      -1       4           5.23
            //         134       2       5       13400
            //           0      10       1           0
            //
            //      123456       3       9    123456000
            //      123456       2       8     12345600
            //      123456       1       7      1234560
            //      123456       0       6       123456
            //      123456      -1       7        12345.6
            //      123456      -2       7         1234.56
            //      123456      -3       7          123.456
            //      123456      -4       7           12.3456
            //      123456      -5       7           1.23456
            //      123456      -6       8           0.123456
            //      123456      -7       9           0.0123456
            //      123456      -8      10           0.00123456
            //      123456      -9      11           0.000123456
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 July 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int MANTISSA, EXPONENT, the decimal number.
            //
            //    Output, int DEC_WIDTH, the "width" of the decimal number.
            //
        {
            int mantissa_abs;
            int ten_pow;
            int value;
            //
            //  Special case of 0.
            //
            value = 1;

            if (mantissa == 0)
            {
                return value;
            }

            //
            //  Determine a power of 10 that is strictly bigger than MANTISSA.
            //  The exponent of that power of 10 is our first estimate for 
            //  the number of places.
            //
            ten_pow = 10;
            mantissa_abs = Math.Abs(mantissa);

            while (ten_pow <= mantissa_abs)
            {
                value = value + 1;
                ten_pow = ten_pow * 10;
            }

            //
            //  If the exponent is nonnegative, that just adds more places.
            //
            if (0 <= exponent)
            {
                value = value + exponent;
            }
            //
            //  If the exponent is a little negative, then we are essentially
            //  just inserting a decimal point, and moving it around.
            //
            else if (-value < exponent)
            {
                value = value + 1;
            }
            //
            //  A very negative value of B means we have a leading 0 and decimal,
            //  and B trailing places.
            //
            else if (exponent <= -value)
            {
                value = 2 - exponent;
            }

            //
            //  Take care of sign of MANTISSA.
            //
            if (mantissa < 0)
            {
                value = value + 1;
            }

            return value;
        }
    }
}