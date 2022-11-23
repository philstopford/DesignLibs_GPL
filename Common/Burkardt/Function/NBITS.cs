﻿using System;
using Burkardt.Types;

namespace Burkardt.Function;

public static class NBITS
{
    public static int nbits_compute ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NBITS_COMPUTE computes the mantissa length minus one.
        //
        //  Discussion:
        //
        //    NBITS is the number of bits (less 1) in the mantissa of the
        //    floating point number number representation of your machine.
        //    It is used to determine the level of accuracy to which the W
        //    function should be calculated.
        //
        //    Most machines use a 24-bit matissa for single precision and
        //    53-56 bits for double. The IEEE standard is 53
        //    bits. The Fujitsu VP2200 uses 56 bits. Long word length
        //    machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
        //    single precision.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
        //    Patricia Culligan-Hensley.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
        //    Algorithm 743: WAPR - A Fortran routine for calculating real 
        //    values of the W-function,
        //    ACM Transactions on Mathematical Software,
        //    Volume 21, Number 2, June 1995, pages 172-181.
        //
        //  Parameters:
        //
        //    Output, int NBITS_COMPUTE, the mantissa length, in bits, 
        //    minus one.
        //
    {
        int nbits = 0;

        double b = 1.0;

        for ( ; ; )
        {
            b /= 2.0;
            double v = b + 1.0;

            if ( Math.Abs(v - 1.0) <= typeMethods.r8_epsilon())
            {
                break;
            }
            nbits += 1;
        }

        return nbits;
    }
}