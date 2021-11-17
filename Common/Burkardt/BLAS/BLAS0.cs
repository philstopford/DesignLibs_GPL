using System;
using System.Numerics;

namespace Burkardt.BLAS;

public static partial class BLAS0
{
    public static Complex c4_uniform_01(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4_UNIFORM_01 returns a unit complex pseudorandom number.
        //
        //  Discussion:
        //
        //    The angle should be uniformly distributed between 0 and 2 * PI,
        //    the square root of the radius uniformly distributed between 0 and 1.
        //
        //    This results in a uniform distribution of values in the unit circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, complex <float> C4_UNIFORM_01, a pseudorandom complex value.
        //
    {
        int k = seed / 127773;

        seed = 16807 * (seed - k * 127773) - k * 2836;

        switch (seed)
        {
            case < 0:
                seed += 2147483647;
                break;
        }

        float r = (float) Math.Sqrt((float) ((double) seed * 4.656612875E-10f));

        k = seed / 127773;

        seed = 16807 * (seed - k * 127773) - k * 2836;

        switch (seed)
        {
            case < 0:
                seed += 2147483647;
                break;
        }

        float theta = (float) (2.0f * Math.PI * (float)
            ((double) seed * 4.656612875E-10f));

        Complex value = new(r * Math.Cos(theta), r * Math.Sin(theta));

        return value;
    }

    public static double zabs1(Complex z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZABS1 returns the L1 norm of a complex <double>.
        //
        //  Discussion:
        //
        //    This routine uses double precision complex arithmetic.
        //
        //    The L1 norm of a complex number is the sum of the absolute values
        //    of the real and imaginary components.
        //
        //    ZABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979.
        //
        //    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
        //    Basic Linear Algebra Subprograms for FORTRAN usage,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, pages 308-323, 1979.
        //
        //  Parameters:
        //
        //    Input, complex <double> Z, the number whose norm is desired.
        //
        //    Output, double ZABS1, the L1 norm of Z.
        //
    {
        double value = Math.Abs(z.Real) + Math.Abs(z.Imaginary);

        return value;
    }

    public static double zabs2(Complex z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZABS2 returns the L2 norm of a complex <double>.
        //
        //  Discussion:
        //
        //    This routine uses double precision complex arithmetic.
        //
        //    The L2 norm of a complex number is the square root of the sum
        //    of the squares of the real and imaginary components.
        //
        //    ZABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979.
        //
        //    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
        //    Basic Linear Algebra Subprograms for FORTRAN usage,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, pages 308-323, 1979.
        //
        //  Parameters:
        //
        //    Input, complex <double> Z, the number whose norm is desired.
        //
        //    Output, float ZABS2, the L2 norm of Z.
        //
    {
        double zr = z.Real;
        double zi = z.Imaginary;
        double value = Math.Sqrt(zr * zr + zi * zi);

        return value;
    }

    public static double zmach(int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZMACH computes machine parameters for complex <double> arithmetic.
        //
        //  Discussion:
        //
        //    Assume the computer has
        //
        //      B = base of arithmetic;
        //      T = number of base B digits;
        //      L = smallest possible exponent;
        //      U = largest possible exponent;
        //
        //    then
        //
        //      EPS = B^(1-T)
        //      TINY = 100.0 * B^(-L+T)
        //      HUGE = 0.01 * B^(U-T)
        //
        //    If complex division is done by
        //
        //      1 / (X+i*Y) = (X-i*Y) / (X^2+Y^2)
        //
        //    then
        //
        //      TINY = sqrt ( TINY )
        //      HUGE = sqrt ( HUGE )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979.
        //
        //    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
        //    Basic Linear Algebra Subprograms for FORTRAN usage,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 3, pages 308-323, 1979.
        //
        //  Parameters:
        //
        //    Input, int JOB:
        //    1, EPS is desired;
        //    2, TINY is desired;
        //    3, HUGE is desired.
        //
        //    Output, double ZMACH, the requested value.
        //
    {
        double eps;
        double huge;
        double s;
        Complex temp1;
        Complex temp2;
        Complex temp3;
        double tiny;
        double value = 0;

        eps = 1.0;

        for (;;)
        {
            eps /= 2.0;
            s = 1.0 + eps;
            if (s <= 1.0)
            {
                break;
            }
        }

        eps = 2.0 * eps;

        s = 1.0;

        for (;;)
        {
            tiny = s;
            s /= 16.0;

            if (s * 1.0 == 0.0)
            {
                break;
            }
        }

        tiny = tiny / eps * 100.0;
        //
        //  Had to insert this manually!
        //
        tiny = Math.Sqrt(tiny);

        switch (false)
        {
            case true:
            {
                temp1 = new Complex(1.0, 0.0);
                temp2 = new Complex(tiny, 0.0);
                temp3 = temp1 / temp2;

                s = temp3.Real;

                tiny = Math.Abs(s - 1.0 / tiny) switch
                {
                    > double.Epsilon => Math.Sqrt(tiny),
                    _ => tiny
                };

                break;
            }
        }

        huge = 1.0 / tiny;

        value = job switch
        {
            1 => eps,
            2 => tiny,
            3 => huge,
            _ => 0.0
        };

        return value;
    }
    //****************************************************************************80

    public static Complex zsign1(Complex z1, Complex z2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZSIGN1 is a complex <double> transfer-of-sign function.
        //
        //  Discussion:
        //
        //    This routine uses double precision complex arithmetic.
        //
        //    The L1 norm is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, complex <double> Z1, Z2, the arguments.
        //
        //    Output, complex <double> ZSIGN1,  a complex value, with the magnitude of
        //    Z1, and the argument of Z2.
        //
    {
        Complex value;

        if (zabs1(z2) == 0.0)
        {
            value = new Complex(0.0, 0.0);
        }
        else
        {
            value = zabs1(z1) * (z2 / zabs1(z2));
        }

        return value;
    }

    public static Complex zsign2(Complex z1, Complex z2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZSIGN2 is a complex <double> transfer-of-sign function.
        //
        //  Discussion:
        //
        //    This routine uses double precision complex arithmetic.
        //
        //    The L2 norm is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, complex <double> Z1, Z2, the arguments.
        //
        //    Output, complex <double> ZSIGN2,  a complex value, with the magnitude of
        //    Z1, and the argument of Z2.
        //
    {
        Complex value;

        if (zabs2(z2) == 0.0)
        {
            value = new Complex(0.0, 0.0);
        }
        else
        {
            value = zabs2(z1) * (z2 / zabs2(z2));
        }

        return value;
    }

}