using System;
using System.Numerics;

namespace Burkardt.Uniform
{
    public static partial class UniformRNG
    {
        public static Complex[] c4mat_uniform_01_new(int m, int n, ref int seed)

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01_NEW returns a new unit pseudorandom C4MAT.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 June 2015
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
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4MAT_UNIFORM_01[M*N], the pseudorandom complex matrix.
//
        {
            const int i4_huge = 2147483647;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("C4MAT_UNIFORM_01_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return new Complex[1];
            }

            Complex[] c = new Complex[m * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    int k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + i4_huge;
                    }

                    float r = (float)Math.Sqrt((float) (seed) * 4.656612875E-10);

                    k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + i4_huge;
                    }

                    float theta = (float)(2.0 * Math.PI * ((float) (seed) * 4.656612875E-10));

                    c[i + j * m] = r * new Complex(Math.Cos(theta), Math.Sin(theta));
                }
            }

            return c;
        }

        public static Complex[] c4vec_uniform_01_new(int n, ref int seed)

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01_NEW returns a unit pseudorandom C4VEC.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 June 2015
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
//    Input, int N, the number of values to compute.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4VEC_UNIFORM_01_NEW[N], the pseudorandom 
//    complex vector.
//
        {
            const int i4_huge = 2147483647;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("C4VEC_UNIFORM_01_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return new Complex[1];
            }

            Complex[] c = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                int k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                float r = (float)Math.Sqrt((float) (seed) * 4.656612875E-10);

                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                float theta = (float)(2.0 * Math.PI * ((float) (seed) * 4.656612875E-10));

                c[i] = r * new Complex(Math.Cos(theta), Math.Sin(theta));
            }

            return c;
        }

        public static Complex c8_uniform_01(ref int seed)

//****************************************************************************80
//
//  Purpose:
//
//    C8_UNIFORM_01 returns a unit double complex pseudorandom number.
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
//    17 April 2006
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
//    Output, complex <double> C8_UNIFORM_01, a pseudorandom complex value.
//
        {
            int k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + 2147483647;
            }

            double r = Math.Sqrt(((double) (seed) * 4.656612875E-10));

            k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + 2147483647;
            }

            double theta = 2.0 * Math.PI * ((double) (seed) * 4.656612875E-10);

            Complex value = new Complex(r * Math.Cos(theta), r * Math.Sin(theta));

            return value;
        }
    }
}