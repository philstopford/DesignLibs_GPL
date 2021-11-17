using System;
using Burkardt.Types;

namespace Burkardt.Uniform;

public static partial class UniformRNG
{
    public static float r4_uniform_01(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_UNIFORM_01 returns a unit pseudorandom R4.
        //
        //  Discussion:
        //
        //    This routine implements the recursion
        //
        //      seed = 16807 * seed mod ( 2^31 - 1 )
        //      r4_uniform_01 = seed / ( 2^31 - 1 )
        //
        //    The integer arithmetic never requires more than 32 bits,
        //    including a sign bit.
        //
        //    If the initial seed is 12345, then the first three computations are
        //
        //      Input     Output      R4_UNIFORM_01
        //      SEED      SEED
        //
        //         12345   207482415  0.096616
        //     207482415  1790989824  0.833995
        //    1790989824  2035175616  0.947702
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2004
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
        //    in Handbook of Simulation
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
        //    Peter Lewis, Allen Goodman, James Miller,
        //    A Pseudo-Random Number Generator for the System/360,
        //    IBM Systems Journal,
        //    Volume 8, pages 136-143, 1969.
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, the "seed" value.  Normally, this
        //    value should not be 0.  On output, SEED has been updated.
        //
        //    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
        //    0 and 1.
        //
    {
        switch (seed)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("R4_UNIFORM_01 - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return 1;
        }

        int k = seed / 127773;

        seed = 16807 * (seed - k * 127773) - k * 2836;

        switch (seed)
        {
            case < 0:
                seed += 2147483647;
                break;
        }

        //
        //  Although SEED can be represented exactly as a 32 bit integer,
        //  it generally cannot be represented exactly as a 32 bit real number.
        //
        float value = seed * 4.656612875E-10f;

        return value;
    }

    public static float r4_uniform_ab(float a, float b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4_UNIFORM_AB returns a scaled pseudorandom R4.
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
        //    21 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, float A, B, the limits of the interval.
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, float R4_UNIFORM_AB, a number strictly between A and B.
        //
    {
        switch (seed)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("R4_UNIFORM_AB - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return 1;
        }

        int k = seed / 127773;

        seed = 16807 * (seed - k * 127773) - k * 2836;

        switch (seed)
        {
            case < 0:
                seed += typeMethods.i4_huge();
                break;
        }

        float value = seed * 4.656612875E-10f;

        value = a + (b - a) * value;

        return value;
    }

    public static float[] r4vec_uniform_ab_new(int n, float b, float c, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R4VEC.
        //
        //  Discussion:
        //
        //    An R4VEC is a vector of R4's.
        //
        //    This routine implements the recursion
        //
        //      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
        //      u = seed / ( 2^31 - 1 )
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
        //    23 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Linus Schrage,
        //    A Guide to Simulation,
        //    Second Edition,
        //    Springer, 1987,
        //    ISBN: 0387964673,
        //    LC: QA76.9.C65.B73.
        //
        //    Bennett Fox,
        //    Algorithm 647:
        //    Implementation and Relative Efficiency of Quasirandom
        //    Sequence Generators,
        //    ACM Transactions on Mathematical Software,
        //    Volume 12, Number 4, December 1986, pages 362-376.
        //
        //    Pierre L'Ecuyer,
        //    Random Number Generation,
        //    in Handbook of Simulation,
        //    edited by Jerry Banks,
        //    Wiley, 1998,
        //    ISBN: 0471134031,
        //    LC: T57.62.H37.
        //
        //    Peter Lewis, Allen Goodman, James Miller,
        //    A Pseudo-Random Number Generator for the System/360,
        //    IBM Systems Journal,
        //    Volume 8, Number 2, 1969, pages 136-143.
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, float B, C, the lower and upper limits of the pseudorandom values.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, float R4VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
        //
    {
        int i;

        switch (seed)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("R4VEC_UNIFORM_AB_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return null;
        }

        float[] r = new float[n];

        for (i = 0; i < n; i++)
        {
            int k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            switch (seed)
            {
                case < 0:
                    seed += typeMethods.i4_huge();
                    break;
            }

            r[i] = (float)(b + (c - b) * seed * 4.656612875E-10);
        }

        return r;
    }

}