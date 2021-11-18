using System;

namespace Burkardt.FullertonFnLib;

public static partial class FullertonLib
{
    public static double r8_mach(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MACH returns double precision real machine constants.
        //
        //  Discussion:
        //
        //    Assuming that the internal representation of a double precision real
        //    number is in base B, with T the number of base-B digits in the mantissa,
        //    and EMIN the smallest possible exponent and EMAX the largest possible 
        //    exponent, then
        //
        //      R8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
        //      R8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
        //      R8_MACH(3) = B^(-T), the smallest relative spacing.
        //      R8_MACH(4) = B^(1-T), the largest relative spacing.
        //      R8_MACH(5) = log10(B).
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
        //    Algorithm 528:
        //    Framework for a Portable Library,
        //    ACM Transactions on Mathematical Software,
        //    Volume 4, Number 2, June 1978, page 176-188.
        //
        //  Parameters:
        //
        //    Input, int I, chooses the parameter to be returned.
        //    1 <= I <= 5.
        //
        //    Output, double R8_MACH, the value of the chosen parameter.
        //
    {
        double value = 0;

        switch (i)
        {
            case 1:
                value = 4.450147717014403E-308;
                break;
            case 2:
                value = 8.988465674311579E+307;
                break;
            case 3:
                value = 1.110223024625157E-016;
                break;
            case 4:
                value = 2.220446049250313E-016;
                break;
            case 5:
                value = 0.301029995663981E+000;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("R8_MACH - Fatal error!");
                Console.WriteLine("  The input argument I is out of bounds.");
                Console.WriteLine("  Legal values satisfy 1 <= I <= 5.");
                Console.WriteLine("  I = " + i + "");
                value = 0.0;
                return 1;
        }

        return value;
    }

    public static void r8_machar(ref long ibeta, ref long it, ref long irnd, ref long ngrd,
            ref long machep, ref long negep, ref long iexp, ref long minexp,
            ref long maxexp, ref double eps, ref double epsneg, ref double xmin, ref double xmax )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MACHAR computes machine constants for R8 arithmetic.
        //
        //  Discussion:
        //
        //    This routine determines the parameters of the floating-point 
        //    arithmetic system specified below.  The determination of the first 
        //    three uses an extension of an algorithm due to Malcolm, 
        //    incorporating some of the improvements suggested by Gentleman and 
        //    Marovich.  
        //
        //    A FORTRAN version of this routine appeared as ACM algorithm 665.
        //
        //    This routine is a C translation of the FORTRAN code, and appeared
        //    as part of ACM algorithm 722.
        //
        //    An earlier version of this program was published in Cody and Waite.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 April 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody,
        //    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
        //      machine parameters,
        //    ACM Transactions on Mathematical Software,
        //    Volume 14, Number 4, pages 303-311, 1988.
        //
        //    William Cody and W Waite,
        //    Software Manual for the Elementary Functions,
        //    Prentice Hall, 1980.
        //
        //    M Gentleman and S Marovich,
        //    Communications of the ACM,
        //    Volume 17, pages 276-277, 1974.
        //
        //    M. Malcolm,
        //    Communications of the ACM,
        //    Volume 15, pages 949-951, 1972.
        //
        //  Parameters:
        //
        //    Output, long int* ibeta, the radix for the floating-point representation.
        //
        //    Output, long int* it, the number of base IBETA digits in the floating-point
        //    significand.
        //
        //    Output, long int* irnd:
        //    0, if floating-point addition chops.
        //    1, if floating-point addition rounds, but not in the IEEE style.
        //    2, if floating-point addition rounds in the IEEE style.
        //    3, if floating-point addition chops, and there is partial underflow.
        //    4, if floating-point addition rounds, but not in the IEEE style, and 
        //      there is partial underflow.
        //    5, if floating-point addition rounds in the IEEE style, and there is 
        //      partial underflow.
        //
        //    Output, long int* ngrd, the number of guard digits for multiplication with
        //    truncating arithmetic.  It is
        //    0, if floating-point arithmetic rounds, or if it truncates and only 
        //      IT base IBETA digits participate in the post-normalization shift of the
        //      floating-point significand in multiplication;
        //   1, if floating-point arithmetic truncates and more than IT base IBETA
        //      digits participate in the post-normalization shift of the floating-point
        //      significand in multiplication.
        //
        //    Output, long int* MACHEP, the largest negative integer such that
        //      1.0 + ( double ) IBETA ^ MACHEP != 1.0, 
        //    except that MACHEP is bounded below by - ( IT + 3 ).
        //
        //    Output, long int* NEGEPS, the largest negative integer such that
        //      1.0 - ( double ) IBETA ) ^ NEGEPS != 1.0, 
        //    except that NEGEPS is bounded below by - ( IT + 3 ).
        //
        //    Output, long int* IEXP, the number of bits (decimal places if IBETA = 10)
        //    reserved for the representation of the exponent (including the bias or
        //    sign) of a floating-point number.
        //
        //    Output, long int* MINEXP, the largest in magnitude negative integer such 
        //    that
        //      ( double ) IBETA ^ MINEXP 
        //    is positive and normalized.
        //
        //    Output, long int* MAXEXP, the smallest positive power of BETA that overflows.
        // 
        //    Output, double* EPS, the smallest positive floating-point number such
        //    that  
        //      1.0 + EPS != 1.0. 
        //    in particular, if either IBETA = 2  or IRND = 0, 
        //      EPS = ( double ) IBETA ^ MACHEP.
        //    Otherwise,  
        //      EPS = ( ( double ) IBETA ^ MACHEP ) / 2.
        //
        //    Output, double* EPSNEG, a small positive floating-point number such that
        //      1.0 - EPSNEG != 1.0. 
        //    In particular, if IBETA = 2 or IRND = 0, 
        //      EPSNEG = ( double ) IBETA ^ NEGEPS.
        //    Otherwise,  
        //      EPSNEG = ( double ) IBETA ^ NEGEPS ) / 2.  
        //    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
        //    smallest number that can alter 1.0 by subtraction.
        //
        //    Output, double* XMIN, the smallest non-vanishing normalized floating-point
        //    power of the radix:
        //      XMIN = ( double ) IBETA ^ MINEXP
        //
        //    Output, float* XMAX, the largest finite floating-point number.  In
        //    particular,
        //      XMAX = ( 1.0 - EPSNEG ) * ( double ) IBETA ^ MAXEXP
        //    On some machines, the computed value of XMAX will be only the second, 
        //    or perhaps third, largest number, being too small by 1 or 2 units in 
        //    the last digit of the significand.
        //
    {
        double a;
        double b;
        double beta;
        double betah;
        double betain;
        int i;
        int itmp;
        int iz;
        int j;
        int k;
        int mx;
        int nxres;
        double one;
        double t;
        double tmp;
        double tmp1;
        double tmpa;
        double two;
        double y;
        double z;
        double zero;

        irnd = 1;
        one = irnd;
        two = one + one;
        a = two;
        b = a;
        zero = 0.0e0;
        //
        //  Determine IBETA and BETA ala Malcolm.
        //
        tmp = a + one - a - one;

        while (Math.Abs(tmp - zero) <= double.Epsilon)
        {
            a += a;
            tmp = a + one;
            tmp1 = tmp - a;
            tmp = tmp1 - one;
        }

        tmp = a + b;
        itmp = (int) (tmp - a);

        while (itmp == 0)
        {
            b += b;
            tmp = a + b;
            itmp = (int) (tmp - a);
        }

        ibeta = itmp;
        beta = ibeta;
        //
        //  Determine IRND, IT.
        //
        it = 0;
        b = one;
        tmp = b + one - b - one;

        while (Math.Abs(tmp - zero) <= double.Epsilon)
        {
            it += 1;
            b *= beta;
            tmp = b + one;
            tmp1 = tmp - b;
            tmp = tmp1 - one;
        }

        irnd = 0;
        betah = beta / two;
        tmp = a + betah;
        tmp1 = tmp - a;

        if (Math.Abs(tmp1 - zero) > double.Epsilon)
        {
            irnd = 1;
        }

        tmpa = a + beta;
        tmp = tmpa + betah;

        irnd = irnd switch
        {
            0 when Math.Abs(tmp - tmpa - zero) > double.Epsilon => 2,
            _ => irnd
        };

        //
        //  Determine NEGEP, EPSNEG.
        //
        negep = it + 3;
        betain = one / beta;
        a = one;

        for (i = 1; i <= negep; i++)
        {
            a *= betain;
        }

        b = a;
        tmp = one - a;
        tmp -= one;

        while (Math.Abs(tmp - zero) <= double.Epsilon)
        {
            a *= beta;
            negep -= 1;
            tmp1 = one - a;
            tmp = tmp1 - one;
        }

        negep = -negep;
        epsneg = a;
        //
        //  Determine MACHEP, EPS.
        //

        machep = -it - 3;
        a = b;
        tmp = one + a;

        while (Math.Abs(tmp - one - zero) <= double.Epsilon)
        {
            a *= beta;
            machep += 1;
            tmp = one + a;
        }

        eps = a;
        //
        //  Determine NGRD.
        //
        ngrd = 0;
        tmp = one + eps;
        tmp *= one;

        ngrd = irnd switch
        {
            0 when Math.Abs(tmp - one - zero) > double.Epsilon => 1,
            _ => ngrd
        };
        //
        //  Determine IEXP, MINEXP and XMIN.
        //
        //  Loop to determine largest I such that (1/BETA) ** (2**(I))
        //  does not underflow.  Exit from loop is signaled by an underflow.
        //

        i = 0;
        k = 1;
        z = betain;
        t = one + eps;
        nxres = 0;

        for (;;)
        {
            y = z;
            z = y * y;
            //
            //  Check for underflow
            //

            a = z * one;
            tmp = z * t;

            if (Math.Abs(a + a - zero) <= double.Epsilon || Math.Abs(z) > y)
            {
                break;
            }

            tmp1 = tmp * betain;

            if (Math.Abs(tmp1 * beta - z) <= double.Epsilon)
            {
                break;
            }

            i += 1;
            k += k;
        }

        //
        //  Determine K such that (1/BETA)**K does not underflow.
        //  First set  K = 2 ^ I.
        //
        iexp = i + 1;
        mx = k + k;
        switch (ibeta)
        {
            //
            //  For decimal machines only
            //
            case 10:
            {
                iexp = 2;
                iz = (int)ibeta;
                while (iz <= k)
                {
                    iz *= (int)ibeta;
                    iexp += 1;
                }

                mx = iz + iz - 1;
                break;
            }
        }

        //
        //  Loop to determine MINEXP, XMIN.
        //  Exit from loop is signaled by an underflow.
        //
        for (;;)
        {
            xmin = y;
            y *= betain;
            a = y * one;
            tmp = y * t;
            tmp1 = a + a;

            if (Math.Abs(tmp1 - zero) <= double.Epsilon || Math.Abs(y) >= xmin)
            {
                break;
            }

            k += 1;
            tmp1 = tmp * betain;
            tmp1 *= beta;

            if (Math.Abs(tmp1 - y) <= double.Epsilon && Math.Abs(tmp - y) > double.Epsilon)
            {
                nxres = 3;
                xmin = y;
                break;
            }

        }

        minexp = -k;
        //
        //  Determine MAXEXP, XMAX.
        //
        if (mx <= k + k - 3 && ibeta != 10)
        {
            mx += mx;
            iexp += 1;
        }

        maxexp = mx + minexp;
        //
        //  Adjust IRND to reflect partial underflow.
        //
        irnd += nxres;
        switch (irnd)
        {
            //
            //  Adjust for IEEE style machines.
            //
            case >= 2:
                maxexp -= 2;
                break;
        }

        //
        //  Adjust for machines with implicit leading bit in binary
        //  significand and machines with radix point at extreme
        //  right of significand.
        //
        i = (int)(maxexp + minexp);

        switch (ibeta)
        {
            case 2 when i == 0:
                maxexp -= 1;
                break;
        }

        switch (i)
        {
            case > 20:
                maxexp -= 1;
                break;
        }

        if (Math.Abs(a - y) > double.Epsilon)
        {
            maxexp -= 2;
        }

        xmax = one - epsneg;
        tmp = xmax * one;

        if (Math.Abs(tmp - xmax) > double.Epsilon)
        {
            xmax = one - beta * epsneg;
        }

        xmax /= (beta * beta * beta * xmin);
        i = (int)(maxexp + minexp + 3);

        switch (i)
        {
            case > 0:
            {
                for (j = 1; j <= i; j++)
                {
                    switch (ibeta)
                    {
                        case 2:
                            xmax += xmax;
                            break;
                    }

                    if (ibeta != 2)
                    {
                        xmax *= beta;
                    }
                }

                break;
            }
        }
    }

    public static double r8_max(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MAX returns the maximum of two R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the quantities to compare.
        //
        //    Output, double R8_MAX, the maximum of X and Y.
        //
    {
        double value = 0;

        if (y < x)
        {
            value = x;
        }
        else
        {
            value = y;
        }

        return value;
    }

    public static double r8_min(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MIN returns the minimum of two R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the quantities to compare.
        //
        //    Output, double R8_MIN, the minimum of X and Y.
        //
    {
        double value = 0;

        if (y < x)
        {
            value = y;
        }
        else
        {
            value = x;
        }

        return value;
    }

    public static double r8_mod(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MOD returns the remainder of R8 division.
        //
        //  Discussion:
        //
        //    If
        //      REM = R8_MOD ( X, Y )
        //      RMULT = ( X - REM ) / Y
        //    then
        //      X = Y * RMULT + REM
        //    where REM has the same sign as X, and abs ( REM ) < Y.
        //
        //  Example:
        //
        //        X         Y     R8_MOD   R8_MOD  Factorization
        //
        //      107        50       7     107 =  2 *  50 + 7
        //      107       -50       7     107 = -2 * -50 + 7
        //     -107        50      -7    -107 = -2 *  50 - 7
        //     -107       -50      -7    -107 =  2 * -50 - 7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number to be divided.
        //
        //    Input, double Y, the number that divides X.
        //
        //    Output, double R8_MOD, the remainder when X is divided by Y.
        //
    {
        double value = 0;

        switch (y)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_MOD - Fatal error!");
                Console.WriteLine("  R8_MOD ( X, Y ) called with Y = " + y + "");
                return 1;
        }

        value = x - (int) (x / y) * y;

        switch (x)
        {
            case < 0.0 when 0.0 < value:
                value -= Math.Abs(y);
                break;
            case > 0.0 when value < 0.0:
                value += Math.Abs(y);
                break;
        }

        return value;
    }

    public static double r8_mop(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MOP returns the I-th power of -1 as an R8 value.
        //
        //  Discussion:
        //
        //    An R8 is an double value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the power of -1.
        //
        //    Output, double R8_MOP, the I-th power of -1.
        //
    {
        double value = (i % 2) switch
        {
            0 => 1.0,
            _ => -1.0
        };

        return value;
    }
}