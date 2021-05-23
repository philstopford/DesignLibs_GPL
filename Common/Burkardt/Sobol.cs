using System;

namespace Sobol
{
    public static partial class SobolSampler
    {
        public class SobolVector
        {
            public int seed { get; set; }
            public float[] quasi { get; set; }
            
            public int resultCode { get; set; }
        }

        public class SobolVectorLarge
        {
            public long seed { get; set; }
            public double[] quasi { get; set; }
            
            public int resultCode { get; set; }
        }

        private const int DIM_MAX2 = 1111;
        private const int LOG_MAX = 30;
        private static bool[] includ = new bool[LOG_MAX];
        private static float recipd;
        private static int seed_save = - 1;
        private static int[,] v = new int[DIM_MAX2,LOG_MAX];

        public static SobolVector i4_sobol ( int dim_num, int seed )
        
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SOBOL generates a new quasirandom Sobol vector with each call.
        //
        //  Discussion:
        //
        //    The routine adapts the ideas of Antonov and Saleev.
        //
        //    This routine uses INT for integers and FLOAT for real values.
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
        //    Original FORTRAN77 version by Bennett Fox.
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
        //    Input/output, int *SEED, the "seed" for the sequence.
        //    This is essentially the index in the sequence of the quasirandom
        //    value to be generated.  On output, SEED has been set to the
        //    appropriate next value, usually simply SEED+1.
        //    If SEED is less than 0 on input, it is treated as though it were 0.
        //    An input value of 0 requests the first (0-th) element of the sequence.
        //
        //    Output, float QUASI[DIM_NUM], the next quasirandom vector.
        //
        {
            SobolVector ret = new SobolVector() {seed = seed, quasi = new float[dim_num]};

            int dim_num_save = 0;
            bool initialized = false;
            int i;
            int l = 1;
            int[] lastq = new int[DIM_MAX2];
            int maxcol = 1;
            int seed_temp;
        
            if ( !initialized || dim_num != dim_num_save )
            {
                initialized = true;

                initv();

                //
                //  Check parameters.
                //
                if ( dim_num < 1 || DIM_MAX2 < dim_num )
                {
                    Console.WriteLine();
                    Console.WriteLine("I4_SOBOL - Fatal error!");
                    Console.WriteLine("  The spatial dimension DIM_NUM should satisfy:");
                    Console.WriteLine("    1 <= DIM_NUM <= " + DIM_MAX2);
                    Console.WriteLine("  But this input value is DIM_NUM = " + dim_num);
                    ret.resultCode = 1;
                    return ret;
                }
        
                dim_num_save = dim_num;
                //
                //  Set ATMOST = 2^LOG_MAX - 1.
                //
                int atmost = 0;
                for (i = 1; i <= LOG_MAX; i++ )
                {
                    atmost = 2 * atmost + 1;
                }
                //
                //  Find the highest 1 bit in ATMOST (should be LOG_MAX).
                //
                maxcol = i4_bit_hi1 ( atmost );
                //
                //  Initialize row 1 of V.
                //
                int j;
                for (j = 0; j < maxcol; j++ )
                {
                    v[0,j] = 1;
                }
                //
                //  Initialize the remaining rows of V.
                //
                for (i = 1; i < dim_num; i++ )
                {
                    //
                    //  The bit pattern of the integer POLY(I) gives the form
                    //  of polynomial I.
                    //
                    //  Find the degree of polynomial I from binary encoding.
                    //
                    j = poly[i];
                    int m = 0;
        
                    while ( true )
                    {
                        j = j / 2;
                        if ( j <= 0 )
                        {
                            break;
                        }
                        m = m + 1;
                    }
                    //
                    //  We expand this bit pattern to separate components
                    //  of the logical array INCLUD.
                    //
                    j = poly[i];
                    int k;
                    for (k = m-1; 0 <= k; k-- )
                    {
                        int j2 = j / 2;
                        includ[k] = ( j != ( 2 * j2 ) );
                        j = j2;
                    }
                    //
                    //  Calculate the remaining elements of row I as explained
                    //  in Bratley and Fox, section 2.
                    //
                    for (j = m; j < maxcol; j++ )
                    {
                        int newv = v[i,j-m];
                        l = 1;

                        for (k = 0; k < m; k++ )
                        {
                            l = 2 * l;

                            if ( includ[k] )
                            {
                                newv = ( newv ^ ( l * v[i,j-k-1] ) );
                            }
                        }
                        v[i,j] = newv;
                    }
                }
                //
                //  Multiply columns of V by appropriate power of 2.
                //
                l = 1;
                for (j = maxcol-2; 0 <= j; j-- )
                {
                    l = 2 * l;
                    for (i = 0; i < dim_num; i++ )
                    {
                        v[i,j] = v[i,j] * l;
                    }
                }
                //
                //  RECIPD is 1/(common denominator of the elements in V).
                //
                recipd = 1.0f / ( 2 * l );
            }
        
            if ( seed < 0 )
            {
                seed = 0;
            }
        
            if ( seed == 0 )
            {
                l = 1;
                for (i = 0; i < dim_num; i++ )
                {
                    lastq[i] = 0;
                }
            }
            else if ( seed == seed_save + 1 )
            {
                l = i4_bit_lo0 ( seed );
            }
            else if ( seed <= seed_save )
            {
                seed_save = 0;
                l = 1;
                for (i = 0; i < dim_num; i++ )
                {
                    lastq[i] = 0;
                }
        
                for (seed_temp = seed_save; seed_temp <= (seed)-1; seed_temp++ )
                {
                    l = i4_bit_lo0 ( seed_temp );

                    for (i = 0; i < dim_num; i++ )
                    {
                        lastq[i] = ( lastq[i] ^ v[i,l-1] );
                    }
                }
                l = i4_bit_lo0 ( seed );
            }
            else if ( seed_save+1 < seed )
            {
                for (seed_temp = seed_save+1; seed_temp <= (seed)-1; seed_temp++ )
                {
                    l = i4_bit_lo0 ( seed_temp );
                    for (i = 0; i < dim_num; i++ )
                    {
                        lastq[i] = ( lastq[i] ^ v[i,l-1] );
                    }
                }
                l = i4_bit_lo0 ( seed );
            }
            //
            //  Check that the user is not calling too many times!
            //
            if ( maxcol < l )
            {
                Console.WriteLine();
                Console.WriteLine("I4_SOBOL - Fatal error!");
                Console.WriteLine("  The value of SEED seems to be too large!");
                Console.WriteLine("  SEED =   " + seed);
                Console.WriteLine("  MAXCOL = " + maxcol);
                Console.WriteLine("  L =      " + l);
                ret.resultCode = 2;
                return ret;
            }
            //
            //  Calculate the new components of QUASI.
            //  The caret indicates the bitwise exclusive OR.
            //
            for (i = 0; i < dim_num; i++ )
            {
                ret.quasi[i] = ( ( float ) lastq[i] ) * recipd;
                lastq[i] = ( lastq[i] ^ v[i,l-1] );
            }

            seed_save = seed;
            seed = seed + 1;

            ret.seed = seed;

            return ret;
        }

        public static int i4_bit_hi1 ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_HI1 returns the position of the high 1 bit base 2 in an integer.
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
        //    13 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be measured.
        //    N should be nonnegative.  If N is nonpositive, I4_BIT_HI1
        //    will always be 0.
        //
        //    Output, int I4_BIT_HI1, the location of the high order bit.
        //
        {
            int bit;

            bit = 0;

            while ( 0 < n )
            {
                bit = bit + 1;
                n = n / 2;
            }

            return bit;
        }

        public static int i4_bit_lo0 ( int n )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
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
        //    08 February 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer to be measured.
        //    N should be nonnegative.
        //
        //    Output, int I4_BIT_LO0, the position of the low 1 bit.
        //
        {
            int bit;
            int n2;
            bit = 0;
            while ( true )
            {
                bit = bit + 1;
                n2 = n / 2;

                if ( n == 2 * n2 )
                {
                    break;
                }
                n = n2;
            }
            return bit;
        }
        
        public static SobolVectorLarge i8_sobol ( int dim_num, long seed )
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
            SobolVectorLarge ret = new SobolVectorLarge() {seed = seed, quasi = new double[dim_num]};

            int DIM_MAX = 40;
            int DIM_MAX2 = 1111;
            int LOG_MAX = 62;
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
            int dim_num_save = 0;
            long i;
            bool[] includ = new bool[LOG_MAX];
            bool initialized = false;
            long j;
            long j2;
            long k;
            long l = 1;
            long[] lastq = new long[DIM_MAX2];
            long  m;
            long maxcol = 62;
            long newv;

            double recipd = 1.0E+00 / ( ( double ) ( 2 * l ) );
            long seed_save = - 1;
            long seed_temp;
            long[,] v = new long[DIM_MAX2,LOG_MAX];

            if ( !initialized || dim_num != dim_num_save )
            {
                initialized = true;
                for ( i = 0; i < DIM_MAX2; i++ )
                {
                    for ( j = 0; j < LOG_MAX; j++ )
                    {
                        v[i,j] = 0;
                    }
                }
                //
                //  Initialize (part of) V.
                //
                initv();

                //
                //  Check parameters.
                //
                if ( dim_num < 1 || DIM_MAX2 < dim_num )
                {
                    Console.WriteLine();;
                    Console.WriteLine("I8_SOBOL - Fatal error!");
                    Console.WriteLine("  The spatial dimension DIM_NUM should satisfy:");
                    Console.WriteLine("    1 <= DIM_NUM <= " + DIM_MAX2);
                    Console.WriteLine("  But this input value is DIM_NUM = " + dim_num);
                    ret.resultCode = 1;
                    return ret;
                }

                dim_num_save = dim_num;
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
                // maxcol = 62;
                //
                //  Initialize row 1 of V.
                //
                for ( j = 0; j < maxcol; j++ )
                {
                    v[0,j] = 1;
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
                    j = poly[i];
                    m = 0;

                    while ( true )
                    {
                        j = j / 2;
                        if ( j <= 0 )
                        {
                            break;
                        }
                        m = m + 1;
                    }
                    //
                    //  We expand this bit pattern to separate components
                    //  of the logical array INCLUD.
                    //
                    j = poly[i];
                    for ( k = m-1; 0 <= k; k-- )
                    {
                        j2 = j / 2;
                        includ[k] = ( j != ( 2 * j2 ) );
                        j = j2;
                    }
                    //
                    //  Calculate the remaining elements of row I as explained
                    //  in Bratley and Fox, section 2.
                    //
                    for ( j = m; j < maxcol; j++ )
                    {
                        newv = v[i,j-m];
                        l = 1;

                        for ( k = 0; k < m; k++ )
                        {
                            l = 2 * l;

                            if ( includ[k] )
                            {
                                newv = ( newv ^ ( l * v[i,j-k-1] ) );
                            }
                        }
                        v[i,j] = newv;
                    }
                }
                //
                //  Multiply columns of V by appropriate power of 2.
                //
                l = 1;
                for ( j = maxcol - 2; 0 <= j; j-- )
                {
                    l = 2 * l;
                    for ( i = 0; i < dim_num; i++ )
                    {
                        v[i,j] = v[i,j] * l;
                    }
                }
                //
                //  RECIPD is 1/(common denominator of the elements in V).
                //
                recipd = 1.0E+00 / ( ( double ) ( 2 * l ) );
            }

            if ( seed < 0 )
            {
                seed = 0;
            }

            if ( seed == 0 )
            {
                l = 1;
                for ( i = 0; i < dim_num; i++ )
                {
                    lastq[i] = 0;
                }
            }
            else if ( seed == seed_save + 1 )
            {
                l = i8_bit_lo0 ( seed );
            }
            else if ( seed <= seed_save )
            {
                seed_save = 0;
                l = 1;
                for ( i = 0; i < dim_num; i++ )
                {
                    lastq[i] = 0;
                }

                for ( seed_temp = seed_save; seed_temp <= (seed)-1; seed_temp++ )
                {

                    l = i8_bit_lo0 ( seed_temp );

                    for ( i = 0; i < dim_num; i++ )
                    {
                        lastq[i] = ( lastq[i] ^ v[i,l-1] );
                    }
                }
                l = i8_bit_lo0 ( seed );
            }
            else if ( seed_save+1 < seed )
            {
                for ( seed_temp = seed_save+1; seed_temp <= (seed)-1; seed_temp++ )
                {

                    l = i8_bit_lo0 ( seed_temp );

                    for ( i = 0; i < dim_num; i++ )
                    {
                        lastq[i] = ( lastq[i] ^ v[i,l-1] );
                    }
                }
                l = i8_bit_lo0 ( seed );
            }
            //
            //  Check that the user is not calling too many times!
            //
            if ( maxcol < l )
            {
                Console.WriteLine();
                Console.WriteLine("I8_SOBOL - Fatal error!");
                Console.WriteLine("  The value of SEED seems to be too large!");
                Console.WriteLine("  SEED =   " + seed);
                Console.WriteLine("  MAXCOL = " + maxcol);
                Console.WriteLine("  L =      " + l);
                ret.resultCode = 2 ;
                return ret;
            }
            //
            //  Calculate the new components of QUASI.
            //  The caret indicates the bitwise exclusive OR.
            //
            for ( i = 0; i < dim_num; i++ )
            {
                ret.quasi[i] = ( ( double ) lastq[i] ) * recipd;

                lastq[i] = ( lastq[i] ^ v[i,l-1] );
            }

            seed_save = seed;
            seed = seed + 1;

            ret.seed = seed;
            return ret;
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
                bit = bit + 1;
                n = n / 2;
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
                bit = bit + 1;
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
}