﻿using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_round(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ROUND rounds an R8 to the nearest integral value.
        //
        //  Example:
        //
        //        X         Value
        //
        //      1.3         1.0
        //      1.4         1.0
        //      1.5         1.0 or 2.0
        //      1.6         2.0
        //      0.0         0.0
        //     -0.7        -1.0
        //     -1.1        -1.0
        //     -1.6        -2.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the value.
        //
        //    Output, double R8_ROUND, the rounded value.
        //
    {
        double value = x switch
        {
            < 0.0 => -Math.Floor(-x + 0.5),
            _ => Math.Floor(x + 0.5)
        };

        return value;
    }

    public static int r8_round_i4(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ROUND_I4 rounds an R8, returning an I4.
        //
        //  Example:
        //
        //        X         Value
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
        //    25 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the value.
        //
        //    Output, int R8_ROUND_I4, the rounded value.
        //
    {
        int value = x switch
        {
            < 0.0 => -(int) Math.Floor(-x + 0.5),
            _ => (int) Math.Floor(x + 0.5)
        };

        return value;
    }

    public static double r8_round2(int nplace, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ROUND2 rounds an R8 in base 2.
        //
        //  Discussion:
        //
        //    Assume that the input quantity X has the form
        //
        //      X = S * J * 2^L
        //
        //    where S is plus or minus 1, L is an integer, and J is a binary
        //    mantissa which is either exactly zero, or greater than or equal
        //    to 0.5 and less than 1.0.
        //
        //    Then on return, XROUND = R8_ROUND2 ( NPLACE, X ) will satisfy
        //
        //      XROUND = S * K * 2^L
        //
        //    where S and L are unchanged, and K is a binary mantissa which
        //    agrees with J in the first NPLACE binary digits and is zero
        //    thereafter.
        //
        //    If NPLACE is 0, XROUND will always be zero.
        //
        //    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
        //
        //    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
        //    or 0.75.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPLACE, the number of binary digits to
        //    preserve.  NPLACE should be 0 or positive.
        //
        //    Input, double X, the real number to be decomposed.
        //
        //    Output, double R8_ROUND2, the rounded value of X.
        //
    {
        int s;
        double xtemp;

        double value = 0.0;
        switch (x)
        {
            //
            //  1: Handle the special case of 0.
            //
            case 0.0:
                return value;
        }

        switch (nplace)
        {
            case <= 0:
                return value;
        }

        switch (x)
        {
            //
            //  2: Determine the sign S.
            //
            case > 0.0:
                s = 1;
                xtemp = x;
                break;
            default:
                s = -1;
                xtemp = -x;
                break;
        }

        //
        //  3: Force XTEMP to lie between 1 and 2, and compute the
        //  logarithm L.
        //
        int l = 0;

        while (2.0 <= xtemp)
        {
            xtemp /= 2.0;
            l += 1;
        }

        while (xtemp < 1.0)
        {
            xtemp *= 2.0;
            l -= 1;
        }

        //
        //  4: Strip out the digits of the mantissa as XMANT, and decrease L.
        //
        double xmant = 0.0;
        int iplace = 0;

        for (;;)
        {
            xmant = 2.0 * xmant;

            switch (xtemp)
            {
                case >= 1.0:
                    xmant += 1.0;
                    xtemp -= 1.0;
                    break;
            }

            iplace += 1;

            if (xtemp == 0.0 || nplace <= iplace)
            {
                value = s * xmant * Math.Pow(2.0, l);
                break;
            }

            l -= 1;
            xtemp *= 2.0;
        }

        return value;
    }

    public static double r8_roundb(int base_, int nplace, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ROUNDB rounds an R8 in a given base.
        //
        //  Discussion:
        //
        //    The code does not seem to do a good job of rounding when
        //    the base is negative.
        //
        //    Assume that the input quantity X has the form
        //
        //      X = S * J * BASE^L
        //
        //    where S is plus or minus 1, L is an integer, and J is a
        //    mantissa base BASE which is either exactly zero, or greater
        //    than or equal to (1/BASE) and less than 1.0.
        //
        //    Then on return, XROUND will satisfy
        //
        //      XROUND = S * K * BASE^L
        //
        //    where S and L are unchanged, and K is a mantissa base BASE
        //    which agrees with J in the first NPLACE digits and is zero
        //    thereafter.
        //
        //    Note that because of rounding, for most bases, most numbers
        //    with a fractional quantities cannot be stored exactly in the
        //    computer, and hence will have trailing "bogus" digits.
        //
        //    If NPLACE is 0, XROUND will always be zero.
        //
        //    If NPLACE is 1, the mantissa of XROUND will be 0,
        //    1/BASE, 2/BASE, ..., (BASE-1)/BASE.
        //
        //    If NPLACE is 2, the mantissa of XROUND will be 0,
        //    BASE/BASE^2, (BASE+1)/BASE^2, ...,
        //    BASE^2-2/BASE^2, BASE^2-1/BASE^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int BASE, the base of the arithmetic.
        //    BASE must not be zero.  Theoretically, BASE may be negative.
        //
        //    Input, int NPLACE, the number of digits base BASE to
        //    preserve.  NPLACE should be 0 or positive.
        //
        //    Input, double X, the number to be decomposed.
        //
        //    Output, double R8_ROUNDB, the rounded value of X.
        //
    {
        int is_;
        double xtemp;

        double value = 0.0;
        double r8_base = base_;
        switch (base_)
        {
            //
            //  0: Error checks.
            //
            case 0:
                Console.WriteLine("");
                Console.WriteLine("R8_ROUNDB - Fatal error!");
                Console.WriteLine("  The base BASE cannot be zero.");
                return 1;
        }

        switch (x)
        {
            //
            //  1: Handle the special case of 0.
            //
            case 0.0:
                return value;
        }

        switch (nplace)
        {
            case <= 0:
                return value;
        }

        switch (x)
        {
            //
            //  2: Determine the sign IS.
            //
            case > 0.0:
                is_ = 1;
                xtemp = x;
                break;
            default:
                is_ = -1;
                xtemp = -x;
                break;
        }

        //
        //  3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
        //  logarithm L.
        //
        int l = 0;

        while (Math.Abs(r8_base) <= Math.Abs(xtemp))
        {
            xtemp /= r8_base;

            switch (xtemp)
            {
                case < 0.0:
                    is_ = -is_;
                    xtemp = -xtemp;
                    break;
            }

            l += 1;
        }

        while (Math.Abs(xtemp) < 1.0)
        {
            xtemp *= r8_base;

            switch (xtemp)
            {
                case < 0.0:
                    is_ = -is_;
                    xtemp = -xtemp;
                    break;
            }

            l -= 1;
        }

        //
        //  4: Now strip out the digits of the mantissa as XMANT, and
        //  decrease L.
        //
        double xmant = 0.0;
        int iplace = 0;
        int js = is_;

        for (;;)
        {
            xmant = r8_base * xmant;

            switch (xmant)
            {
                case < 0.0:
                    js = -js;
                    xmant = -xmant;
                    break;
            }

            switch (xtemp)
            {
                case >= 1.0:
                    xmant += (int) xtemp;
                    xtemp -= (int) xtemp;
                    break;
            }

            iplace += 1;

            if (xtemp == 0.0 || nplace <= iplace)
            {
                value = js * xmant * Math.Pow(r8_base, l);
                break;
            }

            l -= 1;
            xtemp *= r8_base;

            switch (xtemp)
            {
                case < 0.0:
                    is_ = -is_;
                    xtemp = -xtemp;
                    break;
            }
        }

        return value;
    }

    public static double r8_roundx(int nplace, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ROUNDX rounds an R8 in base 10.
        //
        //  Discussion:
        //
        //    Assume that the input quantity X has the form
        //
        //      X = S * J * 10^L
        //
        //    where S is plus or minus 1, L is an integer, and J is a decimal
        //    mantissa which is either exactly zero, or greater than or equal
        //    to 0.1 and less than 1.0.
        //
        //    Then on return, XROUND will satisfy
        //
        //      XROUND = S * K * 10^L
        //
        //    where S and L are unchanged, and K is a decimal mantissa which
        //    agrees with J in the first NPLACE decimal digits and is zero
        //    thereafter.
        //
        //    Note that because of rounding, most decimal fraction quantities
        //    cannot be stored exactly in the computer, and hence will have
        //    trailing "bogus" digits.
        //
        //    If NPLACE is 0, XROUND will always be zero.
        //
        //    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
        //    0.2, ..., or 0.9.
        //
        //    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
        //    0.03, ..., 0.98, 0.99.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPLACE, the number of decimal digits to
        //    preserve.  NPLACE should be 0 or positive.
        //
        //    Input, double X, the number to be decomposed.
        //
        //    Output, double R8_ROUNDX, the rounded value of X.
        //
    {
        int is_;
        double xtemp;

        double xround = 0.0;
        switch (x)
        {
            //
            //  1: Handle the special case of 0.
            //
            case 0.0:
                return xround;
        }

        switch (nplace)
        {
            case <= 0:
                return xround;
        }

        switch (x)
        {
            //
            //  2: Determine the sign IS.
            //
            case > 0.0:
                is_ = 1;
                xtemp = x;
                break;
            default:
                is_ = -1;
                xtemp = -x;
                break;
        }

        //
        //  3: Force XTEMP to lie between 1 and 10, and compute the
        //  logarithm L.
        //
        int l = 0;

        while (10.0 <= x)
        {
            xtemp /= 10.0;
            l += 1;
        }

        while (xtemp < 1.0)
        {
            xtemp *= 10.0;
            l -= 1;
        }

        //
        //  4: Now strip out the digits of the mantissa as XMANT, and
        //  decrease L.
        //
        double xmant = 0.0;
        int iplace = 0;

        for (;;)
        {
            xmant = 10.0 * xmant;

            switch (xtemp)
            {
                case >= 1.0:
                    xmant += (int) xtemp;
                    xtemp -= (int) xtemp;
                    break;
            }

            iplace += 1;

            if (xtemp == 0.0 || nplace <= iplace)
            {
                xround =  is_ * xmant * Math.Pow(10.0, l);
                break;
            }

            l -= 1;
            xtemp *= 10.0;
        }

        return xround;
    }

}