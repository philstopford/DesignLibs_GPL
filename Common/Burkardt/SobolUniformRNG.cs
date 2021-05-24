using System;

namespace Burkardt.Sobol
{
    public static partial class SobolUniformRNG
    {
        static int i4_uniform(int a, int b, int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_UNIFORM returns a scaled pseudorandom I4.
        //
        //  Discussion:
        //
        //    The pseudorandom number should be uniformly distributed
        //    between A and B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Linus Schrage,
        //    A Guide to Simulation,
        //    Springer Verlag, pages 201-202, 1983.
        //
        //    Pierre L'Ecuyer,
        //    Random Number Generation,
        //    in Handbook of Simulation,
        //    edited by Jerry Banks,
        //    Wiley Interscience, page 95, 1998.
        //
        //    Bennett Fox,
        //    Algorithm 647:
        //    Implementation and Relative Efficiency of Quasirandom
        //    Sequence Generators,
        //    ACM Transactions on Mathematical Software,
        //    Volume 12, Number 4, pages 362-376, 1986.
        //
        //    Peter Lewis, Allen Goodman, James Miller
        //    A Pseudo-Random Number Generator for the System/360,
        //    IBM Systems Journal,
        //    Volume 8, pages 136-143, 1969.
        //
        //  Parameters:
        //
        //    Input, int A, B, the limits of the interval.
        //
        //    Input/output, int *SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, int I4_UNIFORM, a number between A and B.
        //
        {
            if (seed == 0)
            {
                Console.WriteLine();
                Console.WriteLine("I4_UNIFORM - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return 1;
            }

            int k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + 2147483647;
            }

            float r = seed * (float)4.656612875E-10;
            //
            //  Scale R to lie between A-0.5 and B+0.5.
            //
            r = (1.0f - r) * (Math.Min(a, b) - 0.5f)
                + r * (Math.Max(a, b) + 0.5f);
            //
            //  Use rounding to convert R to an integer between A and B.
            //
            int value = r4_nint(r);

            value = Math.Max(value, Math.Min(a, b));
            value = Math.Min(value, Math.Max(a, b));

            return value;
        }


        static int r4_nint ( float x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_NINT returns the nearest integer to an R4.
        //
        //  Example:
        //
        //        X         R4_NINT
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
        //    14 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, float X, the value.
        //
        //    Output, int R4_NINT, the nearest integer to X.
        //
        {
            int value = ( int ) ( Math.Abs( x ) + 0.5 );

            if ( x < 0.0 )
            {
                value = -value;
            }

            return value;
        }
    }
}