using System;
using System.Numerics;

namespace Burkardt.BLAS
{
    public static partial class BLAS0
    {
        public static Complex c4_uniform_01 ( ref int seed )

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

            seed = 16807 * ( seed - k * 127773 ) - k * 2836;

            if ( seed < 0 )
            {
                seed = seed + 2147483647;
            }

            float r = (float) Math.Sqrt ( ( float ) ( ( double ) ( seed ) * 4.656612875E-10f ) );

            k = seed / 127773;

            seed = 16807 * ( seed - k * 127773 ) - k * 2836;

            if ( seed < 0 )
            {
                seed = seed + 2147483647;
            }

            float theta = (float) ( 2.0f * Math.PI * ( float )
                ( ( double ) ( seed ) * 4.656612875E-10f ));

            Complex value = new Complex ( r * Math.Cos ( theta ), r * Math.Sin ( theta ) );

            return value;
        }
    }
}
