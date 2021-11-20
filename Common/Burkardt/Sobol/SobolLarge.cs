﻿using System;

namespace Burkardt.Sobol;

public static partial class SobolSampler
{
    public static double[] i8_sobol_generate ( int m, int n, int skip )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SOBOL_GENERATE generates a Sobol dataset.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 December 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points to generate.
        //
        //    Input, int SKIP, the number of initial points to skip.
        //
        //    Output, float I4_SOBOL_GENERATE[M*N], the points.
        //
    {
        double[] r = new double[m*n];

        SobolConfigLarge config = new(m) {seed = skip};

        for (int j = 0; j < n; j++ )
        {
            i8_sobol ( m, ref config );
            for (int i = 0; i < config.quasi.Length; i++)
            {
                r[j * m + i] = config.quasi[i];
            }
        }

        return r;
    }

    public static int i8_sobol ( int dim_num, ref SobolConfigLarge config )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_SOBOL generates a new quasirandom Sobol vector with each call.
        //
        //  Discussion:
        //
        //    The routine adapts the ideas of Antonov and Saleev.
        //
        //    This routine uses LONG LONG INT for integers and DOUBLE for real values.
        //
        //    Thanks to Steffan Berridge for supplying (twice) the properly
        //    formatted V data needed to extend the original routine's dimension
        //    limit from 40 to 1111, 05 June 2007.
        //
        //    Thanks to Francis Dalaudier for pointing out that the range of allowed
        //    values of DIM_NUM should start at 1, not 2!  17 February 2009.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2009
        //
        //  Author:
        //
        //    FORTRAN77 original version by Bennett Fox.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    IA Antonov, VM Saleev,
        //    An Economic Method of Computing LP Tau-Sequences,
        //    USSR Computational Mathematics and Mathematical Physics,
        //    Volume 19, 1980, pages 252 - 256.
        //
        //    Paul Bratley, Bennett Fox,
        //    Algorithm 659:
        //    Implementing Sobol's Quasirandom Sequence Generator,
        //    ACM Transactions on Mathematical Software,
        //    Volume 14, Number 1, pages 88-100, 1988.
        //
        //    Bennett Fox,
        //    Algorithm 647:
        //    Implementation and Relative Efficiency of Quasirandom 
        //    Sequence Generators,
        //    ACM Transactions on Mathematical Software,
        //    Volume 12, Number 4, pages 362-376, 1986.
        //
        //    Stephen Joe, Frances Kuo
        //    Remark on Algorithm 659:
        //    Implementing Sobol's Quasirandom Sequence Generator,
        //    ACM Transactions on Mathematical Software,
        //    Volume 29, Number 1, pages 49-57, March 2003.
        //
        //    Ilya Sobol,
        //    USSR Computational Mathematics and Mathematical Physics,
        //    Volume 16, pages 236-242, 1977.
        //
        //    Ilya Sobol, YL Levitan, 
        //    The Production of Points Uniformly Distributed in a Multidimensional 
        //    Cube (in Russian),
        //    Preprint IPM Akad. Nauk SSSR, 
        //    Number 40, Moscow 1976.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //    DIM_NUM must satisfy 1 <= DIM_NUM <= 1111.
        //
        //    Input/output, long long int *SEED, the "seed" for the sequence.
        //    This is essentially the index in the sequence of the quasirandom
        //    value to be generated.  On output, SEED has been set to the
        //    appropriate next value, usually simply SEED+1.
        //    If SEED is less than 0 on input, it is treated as though it were 0.
        //    An input value of 0 requests the first (0-th) element of the sequence.
        //
        //    Output, double QUASI[DIM_NUM], the next quasirandom vector.
        //
    {
        //
        //  Here, we have commented out the definition of ATMOST, because
        //  in some cases, a compiler was complaining that the value of ATMOST could not
        //  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
        //  so as long as we set MAXCOL (below) to what we expect it should be, we
        //  may be able to get around this difficulty.
        //  JVB, 24 January 2006.
        //
        //static long long int atmost = 4611686018427387903;
        //
        long i;
        long l = 1;
        bool[] includ = new bool[SobolConfigLarge.LOG_MAX];

        config.recipd = 0.5 ;

        if ( !config.initialized || dim_num != config.dim_num_save )
        {
            config.initialized = true;

            config.initv();

            //
            //  Check parameters.
            //
            if ( dim_num < 1 || SobolConfigLarge.DIM_MAX2 < dim_num )
            {
                Console.WriteLine();
                Console.WriteLine("I8_SOBOL - Fatal error!");
                Console.WriteLine("  The spatial dimension DIM_NUM should satisfy:");
                Console.WriteLine("    1 <= DIM_NUM <= " + SobolConfigLarge.DIM_MAX2);
                Console.WriteLine("  But this input value is DIM_NUM = " + dim_num);
                return 1;
            }

            config.dim_num_save = dim_num;
            //
            //  Find the number of bits in ATMOST.
            //
            //  Here, we have short-circuited the computation of MAXCOL from ATMOST, because
            //  in some cases, a compiler was complaining that the value of ATMOST could not
            //  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
            //  so if we know what the answer should be we can try to get it this way!
            //  JVB, 24 January 2006.
            //
            //  maxcol = i8_bit_hi1 ( atmost );
            //
            config.maxcol = 62;
            //
            //  Initialize row 1 of V.
            //
            long j;
            for ( j = 0; j < config.maxcol; j++ )
            {
                config.v[0,j] = 1;
            }
            //
            //  Initialize the remaining rows of V.
            //
            for ( i = 1; i < dim_num; i++ )
            {
                //
                //  The bit pattern of the integer POLY(I) gives the form
                //  of polynomial I.
                //
                //  Find the degree of polynomial I from binary encoding.
                //
                j = config.poly[i];
                long  m = 0;

                while ( true )
                {
                    j /= 2;
                    if ( j <= 0 )
                    {
                        break;
                    }
                    m += 1;
                }
                //
                //  We expand this bit pattern to separate components
                //  of the logical array INCLUD.
                //
                j = config.poly[i];
                long k;
                for ( k = m-1; 0 <= k; k-- )
                {
                    long j2 = j / 2;
                    includ[k] = j != 2 * j2;
                    j = j2;
                }
                //
                //  Calculate the remaining elements of row I as explained
                //  in Bratley and Fox, section 2.
                //
                for ( j = m; j < config.maxcol; j++ )
                {
                    long newv = config.v[i,j-m];
                    l = 1;

                    for ( k = 0; k < m; k++ )
                    {
                        l = 2 * l;

                        switch (includ[k])
                        {
                            case true:
                                newv ^= l * config.v[i,j-k-1];
                                break;
                        }
                    }
                    config.v[i,j] = newv;
                }
            }
            //
            //  Multiply columns of V by appropriate power of 2.
            //
            l = 1;
            for ( j = config.maxcol - 2; 0 <= j; j-- )
            {
                l = 2 * l;
                for ( i = 0; i < dim_num; i++ )
                {
                    config.v[i,j] *= l;
                }
            }
            //
            //  RECIPD is 1/(common denominator of the elements in V).
            //
            config.recipd = (double)1 / ( 2 * l );
        }

        config.seed = config.seed switch
        {
            < 0 => 0,
            _ => config.seed
        };

        switch (config.seed)
        {
            case 0:
            {
                l = 1;
                for ( i = 0; i < dim_num; i++ )
                {
                    config.lastq[i] = 0;
                }

                break;
            }
            default:
            {
                if ( config.seed == config.seed_save + 1 )
                {
                    l = i8_bit_lo0 ( config.seed );
                }
                else
                {
                    long seed_temp;
                    if ( config.seed <= config.seed_save )
                    {
                        config.seed_save = 0;
                        l = 1;
                        for ( i = 0; i < dim_num; i++ )
                        {
                            config.lastq[i] = 0;
                        }

                        for ( seed_temp = config.seed_save; seed_temp <= config.seed-1; seed_temp++ )
                        {

                            l = i8_bit_lo0 ( seed_temp );

                            for ( i = 0; i < dim_num; i++ )
                            {
                                config.lastq[i] ^= config.v[i,l-1];
                            }
                        }
                        l = i8_bit_lo0 ( config.seed );
                    }
                    else if ( config.seed_save+1 < config.seed )
                    {
                        for ( seed_temp = config.seed_save+1; seed_temp <= config.seed-1; seed_temp++ )
                        {
                            l = i8_bit_lo0 ( seed_temp );

                            for ( i = 0; i < dim_num; i++ )
                            {
                                config.lastq[i] ^= config.v[i,l-1];
                            }
                        }
                        l = i8_bit_lo0 ( config.seed );
                    }
                }

                break;
            }
        }
        //
        //  Check that the user is not calling too many times!
        //
        if ( config.maxcol < l )
        {
            Console.WriteLine();
            Console.WriteLine("I8_SOBOL - Fatal error!");
            Console.WriteLine("  The value of SEED seems to be too large!");
            Console.WriteLine("  SEED =   " + config.seed);
            Console.WriteLine("  MAXCOL = " + config.maxcol);
            Console.WriteLine("  L =      " + l);
            return 2;
        }
        //
        //  Calculate the new components of QUASI.
        //  The caret indicates the bitwise exclusive OR.
        //
        for ( i = 0; i < dim_num; i++ )
        {
            config.quasi[i] = config.lastq[i] * config.recipd;
            config.lastq[i] ^= config.v[i,l-1];
        }

        config.seed_save = config.seed;
        config.seed += 1;

        return 0;
    }

        
    public static int i8_bit_hi1 ( long n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_BIT_HI1 returns the position of the high 1 bit base 2 in an integer.
        //
        //  Example:
        //
        //       N    Binary    Hi 1
        //    ----    --------  ----
        //       0           0     0
        //       1           1     1
        //       2          10     2
        //       3          11     2 
        //       4         100     3
        //       5         101     3
        //       6         110     3
        //       7         111     3
        //       8        1000     4
        //       9        1001     4
        //      10        1010     4
        //      11        1011     4
        //      12        1100     4
        //      13        1101     4
        //      14        1110     4
        //      15        1111     4
        //      16       10000     5
        //      17       10001     5
        //    1023  1111111111    10
        //    1024 10000000000    11
        //    1025 10000000001    11
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
        //  Parameters:
        //
        //    Input, long long int N, the integer to be measured.
        //    N should be nonnegative.  If N is nonpositive, I8_BIT_HI1
        //    will always be 0.
        //
        //    Output, int I8_BIT_HI1, the number of bits base 2.
        //
    {
        int bit = 0;

        while ( 0 < n )
        {
            bit += 1;
            n /= 2;
        }

        return bit;
    }

        
    public static int i8_bit_lo0 ( long n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I8_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
        //
        //  Example:
        //
        //       N    Binary    Lo 0
        //    ----    --------  ----
        //       0           0     1
        //       1           1     2
        //       2          10     1
        //       3          11     3 
        //       4         100     1
        //       5         101     2
        //       6         110     1
        //       7         111     4
        //       8        1000     1
        //       9        1001     2
        //      10        1010     1
        //      11        1011     3
        //      12        1100     1
        //      13        1101     2
        //      14        1110     1
        //      15        1111     5
        //      16       10000     1
        //      17       10001     2
        //    1023  1111111111    11
        //    1024 10000000000     1
        //    1025 10000000001     2
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
        //  Parameters:
        //
        //    Input, long long int N, the integer to be measured.
        //    N should be nonnegative.
        //
        //    Output, int I8_BIT_LO0, the position of the low 1 bit.
        //
    {
        int bit = 0;

        while ( true )
        {
            bit += 1;
            long n2 = n / 2;

            if ( n == 2 * n2 )
            {
                break;
            }
            n = n2;
        }

        return bit;
    }

}