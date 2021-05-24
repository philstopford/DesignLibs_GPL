using System;

namespace Burkardt.Sobol
{
    public static partial class SobolUniformRNG
    {
        public static long i8_uniform ( long a, long b, long seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_UNIFORM returns a scaled pseudorandom I8.
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
        //    12 May 2007
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
        //    Input, long long int A, B, the limits of the interval.
        //
        //    Input/output, int *SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, long long int I8_UNIFORM, a number between A and B.
        //
        {
            if ( seed == 0 )
            {
                Console.WriteLine();
                Console.WriteLine("I8_UNIFORM - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return 1;
            }

            long k = seed / 127773;

            seed = 16807 * ( seed - k * 127773 ) - k * 2836;

            if ( seed < 0 )
            {
                seed = seed + 2147483647;
            }

            double r = seed * 4.656612875E-10;
            //
            //  Scale R to lie between A-0.5 and B+0.5.
            //
            r = ( 1.0 - r ) * ( Math.Min ( a, b ) - 0.5 ) 
            +         r   * ( Math.Max ( a, b ) + 0.5 );
            //
            //  Use rounding to convert R to an integer between A and B.
            //
            long value = r8_nint( r );

            value = Math.Max ( value, Math.Min ( a, b ) );
            value = Math.Min ( value, Math.Max ( a, b ) );

            return value;
        }
        
        
        static long r8_nint ( double x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_NINT returns the nearest integer to an R8.
        //
        //  Examples:
        //
        //        X         R8_NINT
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
        //    Input, double X, the value.
        //
        //    Output, int R8_NINT, the nearest integer to X.
        //
        {
            long value =  (long) ( Math.Abs ( x ) + 0.5 );

            if ( x < 0.0 )
            {
                value = -value;
            }

            return value;
        }
    }
}