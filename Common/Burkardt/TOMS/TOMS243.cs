﻿using System;
using System.Numerics;
using Burkardt.Types;

namespace Burkardt.TOMSNS;

public static partial class TOMS
{
    public static Complex toms243 ( Complex z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOMS243 computes the natural logarithm for complex values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Collens,
        //    Algorithm 243: Logarithm of a Complex Number,
        //    Communications of the Association for Computing Machinery,
        //    Volume 7, Number 11, November 1964, page 660.
        //
        //  Parameters:
        //
        //    Input, complex <double> Z, the argument of the function.
        //
        //    Output, complex <double> TOMS243, the value of the function.
        //
    {
        double c;
        double d;

        double a = z.Real;
        double b = z.Imaginary;
        switch (a)
        {
            //
            //  Ugly hack to get NaN values.
            //
            case 0.0 when b == 0.0:
                c = 1.0 / a;
                d = 1.0 / b;
                break;
            default:
            {
                double e = a / 2.0;
                double f = b / 2.0;
                switch (Math.Abs ( e ))
                {
                    case < 0.5 when Math.Abs ( f ) < 0.5:
                        c = Math.Abs ( 2.0 * a ) + Math.Abs ( 2.0 * b );
                        d = 8.0 * ( a / c ) * a + 8.0 * ( b / c ) * b;
                        c = 0.5 * ( Math.Log ( c ) + Math.Log ( d ) ) - Math.Log ( Math.Sqrt ( 8.0 ) );
                        break;
                    default:
                        c = Math.Abs ( e / 2.0 ) + Math.Abs ( f / 2.0 );
                        d = 0.5 * ( e / c ) * e + 0.5 * ( f / c ) * f;
                        c = 0.5 * ( Math.Log ( c ) + Math.Log ( d ) ) + Math.Log ( Math.Sqrt ( 8.0 ) );
                        break;
                }

                if ( a != 0.0 && Math.Abs ( f ) <= Math.Abs ( e ) )
                {
                    if ( Math.Abs(Math.CopySign ( 1.0, a ) - -1.0) > typeMethods.r8_epsilon() )
                    {
                        d = Math.Atan ( b / a );
                    }
                    else if ( Math.Abs(Math.CopySign ( 1.0, b ) - -1.0) > typeMethods.r8_epsilon() )
                    {
                        d = Math.Atan ( b / a ) + Math.PI;
                    }
                    else
                    {
                        d = Math.Atan ( b / a ) - Math.PI;
                    }
                }
                else
                {
                    d = - Math.Atan ( a / b ) + Math.PI / 2.0 * Math.CopySign ( 1.0, b );
                }

                break;
            }
        }

        Complex value = new( c, d );

        return value;
    }
}