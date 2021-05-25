using System;

namespace Burkardt.Uniform
{
    public static partial class UniformRNG
    {
        public static double[] r8mat_uniform_01 ( int m, int n, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8's.
        //
        //    This routine implements the recursion
        //
        //      seed = 16807 * seed mod ( 2^31 - 1 )
        //      unif = seed / ( 2^31 - 1 )
        //
        //    The integer arithmetic never requires more than 32 bits,
        //    including a sign bit.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2005
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
        //    Bennett Fox,
        //    Algorithm 647:
        //    Implementation and Relative Efficiency of Quasirandom
        //    Sequence Generators,
        //    ACM Transactions on Mathematical Software,
        //    Volume 12, Number 4, pages 362-376, 1986.
        //
        //    Peter Lewis, Allen Goodman, James Miller,
        //    A Pseudo-Random Number Generator for the System/360,
        //    IBM Systems Journal,
        //    Volume 8, pages 136-143, 1969.
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, int &SEED, the "seed" value.  Normally, this
        //    value should not be 0.  On output, SEED has
        //    been updated.
        //
        //    Output, double R8MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
        //
        {
            double[] r = new double[m*n];

            if ( seed == 0 )
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_UNIFORM_01 - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return r;
            }
            
            for (int j = 0; j < n; j++ )
            {
                for (int i = 0; i < m; i++ )
                {
                    int k = seed / 127773;

                    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

                    if ( seed < 0 )
                    {
                        seed = seed + 2147483647;
                    }

                    r[i+j*m] = seed * 4.656612875E-10;
                }
            }

            return r;
        }
        
    }
}