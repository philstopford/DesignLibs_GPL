using System;
using System.Numerics;

namespace Burkardt.FullertonFnLib;

public static partial class FullertonLib
{
    public static Complex cos(Complex z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4_COS evaluates the cosine of a C4 argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 August 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, complex <float> Z, the argument.
        //
        //    Output, complex <float> C4_COS, the cosine of Z.
        //
    {
        double x = z.Real;
        double y = z.Imaginary;

        double cs = Math.Cos(x);

        Complex value = new(cs * Math.Cosh(y), -Math.Sin(x) * Math.Sinh(y));

        return value;
    }

    public static Complex c4_sin(Complex z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4_SIN evaluates the sine of a C4 argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 August 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, complex <float> Z, the argument.
        //
        //    Output, complex <float> C4_SIN, the sine of Z.
        //
    {
        double x = z.Real;
        double y = z.Imaginary;

        double sn = Math.Sin(x);

        Complex value = new(sn * Math.Cosh(y), Math.Cos(x) * Math.Sinh(y));

        return value;
    }
}