using System;
using Burkardt.Types;

namespace Burkardt.Uniform
{
    public static partial class UniformRNG
    {
        public static double[] r8mat_uniform_01(int m, int n, ref int seed)
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
            double[] r = new double[m * n];

            if (seed == 0)
            {
                Console.WriteLine();
                Console.WriteLine("R8MAT_UNIFORM_01 - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return r;
            }

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    int k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + 2147483647;
                    }

                    r[i + j * m] = seed * 4.656612875E-10;
                }
            }

            return r;
        }

        public static double r8_uniform_ab(double a, double b, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_UNIFORM_AB returns a scaled pseudorandom R8.
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
            //    09 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, the limits of the interval.
            //
            //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
            //    On output, SEED has been updated.
            //
            //    Output, double R8_UNIFORM_AB, a number strictly between A and B.
            //
        {
            int i4_huge = 2147483647;
            int k;
            double value;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_UNIFORM_AB - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
            }

            k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + i4_huge;
            }

            value = (double) (seed) * 4.656612875E-10;

            value = a + (b - a) * value;

            return value;
        }

        public static void r8vec_uniform_01 ( int n, ref int seed, ref double[] r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
        //
        //  Discussion:
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
        //    19 August 2004
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R[N], the vector of pseudorandom values.
        //
        {
            int i;
            int i4_huge = 2147483647;
            int k;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_UNIFORM_01 - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                r[i] = (double) (seed) * 4.656612875E-10;
            }

            return;
        }

        public static double[] r8vec_uniform_01_new(int n, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC
            //
            //  Discussion:
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
            //    19 August 2004
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
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
            //
        {
            double[] r = new double[n];

            for (int i = 0; i < n; i++)
            {
                int k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + 2147483647;
                }

                r[i] = (double) (seed) * 4.656612875E-10;
            }

            return r;
        }

        public static double r8_uniform(double a, double b, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_UNIFORM returns a scaled pseudorandom R8.
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
            //    Input, double A, B, the limits of the interval.
            //
            //    Input/output, int *SEED, the "seed" value, which should NOT be 0.
            //    On output, SEED has been updated.
            //
            //    Output, double R8_UNIFORM, a number strictly between A and B.
            //
        {
            int i4_huge = 2147483647;
            int k;
            double value;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_UNIFORM - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return (1);
            }

            k = seed / 127773;

            seed = 16807 * (seed - k * 127773) - k * 2836;

            if (seed < 0)
            {
                seed = seed + i4_huge;
            }

            value = (double) (seed) * 4.656612875E-10;

            value = a + (b - a) * value;

            return value;
        }

        public static double[] r8vec_uniform_ab_new(int n, double a, double b, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.
            //
            //  Discussion:
            //
            //    Each dimension ranges from A to B.
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
            //    09 April 2012
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
            //    Input, double A, B, the lower and upper limits of the pseudorandom values.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
            //
        {
            const int i4_huge = 2147483647;

            if (seed == 0)
            {
                Console.WriteLine();
                Console.WriteLine("R8VEC_UNIFORM_AB_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return new double[1];
            }

            double[] r = new double[n];

            for (int i = 0; i < n; i++)
            {
                int k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                r[i] = a + (b - a) * (double) (seed) * 4.656612875E-10;
            }

            return r;
        }

        public static double[] r8mat_uniform_01_new(int m, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
            //    in column-major order.
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
            //    Philip Lewis, Allen Goodman, James Miller,
            //    A Pseudo-Random Number Generator for the System/360,
            //    IBM Systems Journal,
            //    Volume 8, pages 136-143, 1969.
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input/output, int &SEED, the "seed" value.  Normally, this
            //    value should not be 0, otherwise the output value of SEED
            //    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
            //    been updated.
            //
            //    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
            //
        {
            int i;
            const int i4_huge = 2147483647;
            int j;
            int k;
            double[] r;

            r = new double[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + i4_huge;
                    }

                    r[i + j * m] = (double) (seed) * 4.656612875E-10;
                }
            }

            return r;
        }

        public static void r8mat_uniform_ab(int m, int n, double a, double b, ref int seed, ref double[] r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8's.
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
        //    09 April 2012
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
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A, B, the limits of the pseudorandom values.
        //
        //    Input/output, int &SEED, the "seed" value.  Normally, this
        //    value should not be 0.  On output, SEED has 
        //    been updated.
        //
        //    Output, double R[M*N], a matrix of pseudorandom values.
        //
        {
            int i;
            const int i4_huge = 2147483647;
            int j;
            int k;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8MAT_UNIFORM_AB - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + i4_huge;
                    }

                    r[i + j * m] = a + (b - a) * (double) (seed) * 4.656612875E-10;
                }
            }
        }

        public static double[] r8mat_uniform_ab_new(int m, int n, double a, double b, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_UNIFORM_AB_NEW returns a new scaled pseudorandom R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is an array of R8's.
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
            //    09 April 2012
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
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double A, B, the limits of the pseudorandom values.
            //
            //    Input/output, int &SEED, the "seed" value.  Normally, this
            //    value should not be 0.  On output, SEED has 
            //    been updated.
            //
            //    Output, double R8MAT_UNIFORM_AB_NEW[M*N], a matrix of pseudorandom values.
            //
        {
            int i;
            const int i4_huge = 2147483647;
            int j;
            int k;
            double[] r;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8MAT_UNIFORM_AB_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return null;
            }

            r = new double[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + i4_huge;
                    }

                    r[i + j * m] = a + (b - a) * (double) (seed) * 4.656612875E-10;
                }
            }

            return r;
        }

        public static void r8mat_uniform_abvec(int m, int n, double[] a, double[] b, ref int seed,
        ref double[] r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UNIFORM_ABVEC returns a scaled pseudorandom R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8's.
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
        //    09 April 2012
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
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M], B[M], the limits of the pseudorandom values.
        //
        //    Input/output, int &SEED, the "seed" value.  Normally, this
        //    value should not be 0.  On output, SEED has 
        //    been updated.
        //
        //    Output, double R[M*N], a matrix of pseudorandom values.
        //
        {
            int i;
            const int i4_huge = 2147483647;
            int j;
            int k;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8MAT_UNIFORM_ABVEC - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + i4_huge;
                    }

                    r[i + j * m] = a[i] + (b[i] - a[i]) * (double) (seed) * 4.656612875E-10;
                }
            }
        }

        public static double[] r8mat_uniform_abvec_new(int m, int n, double[] a, double[] b,
        ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UNIFORM_ABVEC_NEW returns a new scaled pseudorandom R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8's.
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
        //    09 April 2012
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
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M], B[M], the limits of the pseudorandom values.
        //
        //    Input/output, int &SEED, the "seed" value.  Normally, this
        //    value should not be 0.  On output, SEED has 
        //    been updated.
        //
        //    Output, double R8MAT_UNIFORM_ABVEC_NEW[M*N], a matrix of
        //    pseudorandom values.
        //
        {
            int i;
            const int i4_huge = 2147483647;
            int j;
            int k;
            double[] r;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8MAT_UNIFORM_ABVEC_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return null;
            }

            r = new double[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    k = seed / 127773;

                    seed = 16807 * (seed - k * 127773) - k * 2836;

                    if (seed < 0)
                    {
                        seed = seed + i4_huge;
                    }

                    r[i + j * m] = a[i] + (b[i] - a[i]) * (double) (seed) * 4.656612875E-10;
                }
            }

            return r;
        }
        
        public static double[] r8vec_uniform_unit_new ( int m, ref int seed )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIFORM_UNIT_NEW generates a random unit vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the space.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8VEC_UNIFORM_UNIT_NEW[M], a random direction vector, 
            //    with unit norm.
            //
        {
            double[] a;
            int i;
            double norm;
            //
            //  Take M random samples from the normal distribution.
            //
            a = typeMethods.r8vec_normal_01_new ( m, ref seed );
            //
            //  Compute the norm.
            //
            norm = 0.0;
            for ( i = 0; i < m; i++ )
            {
                norm = norm + a[i] * a[i];
            }
            norm = Math.Sqrt ( norm );
            //
            //  Normalize.
            //
            for ( i = 0; i < m; i++ )
            {
                a[i] = a[i] / norm;
            }

            return a;
        }
        
        public static void r8vec_uniform_abvec ( int n, double[] a, double[] b, ref int seed, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_UNIFORM_ABVEC returns a scaled pseudorandom R8VEC.
        //
        //  Discussion:
        //
        //    Dimension I ranges from A[I] to B[I].
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
        //    09 April 2012
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
        //    Input, double A[N], B[N], the lower and upper limits of the 
        //    pseudorandom values.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double X[N], the vector of pseudorandom values.
        //
        {
            int i;
            const int i4_huge = 2147483647;
            int k;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_UNIFORM_ABVEC - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                x[i] = a[i] + (b[i] - a[i]) * (double) (seed) * 4.656612875E-10;
            }
        }

        public static void r8vec_uniform_ab(int n, double a, double b, ref int seed, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
            //
            //  Discussion:
            //
            //    Each dimension ranges from A to B.
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
            //    09 April 2012
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
            //    Input, double A, B, the lower and upper limits of the pseudorandom values.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double X[N], the vector of pseudorandom values.
            //
        {
            int i;
            const int i4_huge = 2147483647;
            int k;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_UNIFORM_AB - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                x[i] = a + (b - a) * (double) (seed) * 4.656612875E-10;
            }
        }

        public static double[] r8vec_uniform_abvec_new(int n, double[] a, double[] b, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_UNIFORM_ABVEC_NEW returns a scaled pseudorandom R8VEC.
            //
            //  Discussion:
            //
            //    Dimension I ranges from A[I] to B[I].
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
            //    09 April 2012
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
            //    Input, double A[N], B[N], the lower and upper limits of the 
            //    pseudorandom values.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8VEC_UNIFORM_ABVEC_NEW[N], the vector of 
            //    pseudorandom values.
            //
        {
            int i;
            const int i4_huge = 2147483647;
            int k;
            double[] r;

            if (seed == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_UNIFORM_ABVEC_NEW - Fatal error!");
                Console.WriteLine("  Input value of SEED = 0.");
                return null;
            }

            r = new double[n];

            for (i = 0; i < n; i++)
            {
                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                if (seed < 0)
                {
                    seed = seed + i4_huge;
                }

                r[i] = a[i] + (b[i] - a[i]) * (double) (seed) * 4.656612875E-10;
            }

            return r;
        }

    }
}