using System;
using Burkardt.Types;

namespace Burkardt.Uniform
{
    public static partial class UniformRNG
    {
        public static int i4_uniform_ab(int a, int b, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
            //    02 October 2012
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
            //    Input, int A, B, the limits of the interval.
            //
            //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
            //    On output, SEED has been updated.
            //
            //    Output, int I4_UNIFORM, a number between A and B.
            //
        {
            int c;
            const int i4_huge = 2147483647;
            int k;
            float r;
            int value;

            if (seed == 0)
            {
                Console.WriteLine();
                Console.WriteLine("R4_UNIFORM_AB - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return 1;
            }

            //
            //  Guarantee A <= B.
            //
            if (b < a)
            {
                c = a;
                a = b;
                b = c;
            }

            k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + i4_huge;
            }

            r = seed * (float) 4.656612875E-10;
            //
            //  Scale R to lie between A-0.5 and B+0.5.
            //
            r = (1.0f - r) * (a - 0.5f)
                + r * (b + 0.5f);
            //
            //  Use rounding to convert R to an integer between A and B.
            //
            value = typeMethods.r4_nint(r);
            //
            //  Guarantee A <= VALUE <= B.
            //
            if (value < a)
            {
                value = a;
            }

            if (b < value)
            {
                value = b;
            }

            return value;
        }

        public static int i4_uniform(int a, int b, ref int seed)
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

            float r = seed * (float) 4.656612875E-10;
            //
            //  Scale R to lie between A-0.5 and B+0.5.
            //
            r = (1.0f - r) * (Math.Min(a, b) - 0.5f)
                + r * (Math.Max(a, b) + 0.5f);
            //
            //  Use rounding to convert R to an integer between A and B.
            //
            int value = typeMethods.r4_nint(r);

            value = Math.Max(value, Math.Min(a, b));
            value = Math.Min(value, Math.Max(a, b));

            return value;
        }

        public static int[] i4vec_uniform_ab_new(int n, int a, int b, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    The pseudorandom numbers should be uniformly distributed
            //    between A and B.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 May 2012
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
            //    Input, int N, the dimension of the vector.
            //
            //    Input, int A, B, the limits of the interval.
            //
            //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
            //    On output, SEED has been updated.
            //
            //    Output, int IVEC_UNIFORM_AB_NEW[N], a vector of random values 
            //    between A and B.
            //
        {
            int c;
            const int i4_huge = 2147483647;
            int k;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_UNIFORM_AB_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return new int[0];
            }

            //
            //  Guarantee A <= B.
            //
            if (b < a)
            {
                c = a;
                a = b;
                b = c;
            }

            int[] x = new int[n];

            for (int i = 0; i < n; i++)
            {
                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                float r = (float) (seed) * 4.656612875E-10f;
                //
                //  Scale R to lie between A-0.5 and B+0.5.
                //
                r = (1.0f - r) * ((float) a - 0.5f)
                    + r * ((float) b + 0.5f);
                //
                //  Use rounding to convert R to an integer between A and B.
                //
                int value = (int) Math.Round(r);
                //
                //  Guarantee A <= VALUE <= B.
                //
                if (value < a)
                {
                    value = a;
                }

                if (b < value)
                {
                    value = b;
                }

                x[i] = value;
            }

            return x;
        }

        public static int[] i4vec_uniform(int n, int a, int b, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
            //
            //  Discussion:
            //
            //    The pseudorandom numbers should be uniformly distributed
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
            //    Input, integer N, the dimension of the vector.
            //
            //    Input, int A, B, the limits of the interval.
            //
            //    Input/output, int *SEED, the "seed" value, which should NOT be 0.
            //    On output, SEED has been updated.
            //
            //    Output, int IVEC_UNIFORM[N], a vector of random values between A and B.
            //
        {
            int i;
            int k;
            float r;
            int value;
            int[] x;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_UNIFORM - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return null;
            }

            x = new int[n];

            for (i = 0; i < n; i++)
            {
                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + 2147483647;
                }

                r = (float) (seed) * 4.656612875E-10f;
                //
                //  Scale R to lie between A-0.5 and B+0.5.
                //
                r = (float)((1.0 - r) * ((float) (Math.Min(a, b)) - 0.5)
                    + r * ((float) (Math.Max(a, b)) + 0.5));
                //
                //  Use rounding to convert R to an integer between A and B.
                //
                value = (int)(r);

                value = Math.Max(value, Math.Min(a, b));
                value = Math.Min(value, Math.Max(a, b));

                x[i] = value;
            }

            return x;
        }

    }
}