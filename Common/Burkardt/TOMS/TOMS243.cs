﻿using System;
using System.Numerics;

namespace Burkardt.TOMSNS
{
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
            double a;
            double b;
            double c;
            double d;
            double e;
            double f;
            double r8_pi = 3.141592653589793;
            Complex value;

            a = ( z.Real );
            b = ( z.Imaginary );
            //
            //  Ugly hack to get NaN values.
            //
            if ( a == 0.0 && b == 0.0 )
            {
                c = 1.0 / a;
                d = 1.0 / b;
            }
            else
            {
                e = a / 2.0;
                f = b / 2.0;
                if ( Math.Abs ( e ) < 0.5 && Math.Abs ( f ) < 0.5 )
                {
                    c = Math.Abs ( 2.0 * a ) + Math.Abs ( 2.0 * b );
                    d = 8.0 * ( a / c ) * a + 8.0 * ( b / c ) * b;
                    c = 0.5 * ( Math.Log ( c ) + Math.Log ( d ) ) - Math.Log ( Math.Sqrt ( 8.0 ) );
                }
                else
                {
                    c = Math.Abs ( e / 2.0 ) + Math.Abs ( f / 2.0 );
                    d = 0.5 * ( e / c ) * e + 0.5 * ( f / c ) * f;
                    c = 0.5 * ( Math.Log ( c ) + Math.Log ( d ) ) + Math.Log ( Math.Sqrt ( 8.0 ) );
                }

                if ( ( a != 0.0 ) && Math.Abs ( f ) <= Math.Abs ( e ) )
                {
                    if ( Math.CopySign ( 1.0, a ) != -1.0 )
                    {
                        d = Math.Atan ( b / a );
                    }
                    else if ( Math.CopySign ( 1.0, b ) != -1.0 )
                    {
                        d = Math.Atan ( b / a ) + r8_pi;
                    }
                    else
                    {
                        d = Math.Atan ( b / a ) - r8_pi;
                    }
                }
                else
                {
                    d = - Math.Atan ( a / b ) + r8_pi / 2.0 * Math.CopySign ( 1.0, b );
                }

            }

            value = new Complex ( c, d );

            return value;
        }
    }
}