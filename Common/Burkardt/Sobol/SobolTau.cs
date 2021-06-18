namespace Burkardt.Sobol
{
    public static class Tau
    {
        public static int tau_sobol ( int dim_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TAU_SOBOL defines favorable starting seeds for Sobol sequences.
        //
        //  Discussion:
        //
        //    For spatial dimensions 1 through 13, this routine returns 
        //    a "favorable" value TAU by which an appropriate starting point
        //    in the Sobol sequence can be determined.
        //
        //    These starting points have the form N = 2**K, where
        //    for integration problems, it is desirable that
        //      TAU + DIM_NUM - 1 <= K
        //    while for optimization problems, it is desirable that
        //      TAU < K.
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
        //    Original FORTRAN77 version by Bennett Fox.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    IA Antonov, VM Saleev,
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
        //    Input, int DIM_NUM, the spatial dimension.  Only values
        //    of 1 through 13 will result in useful responses.
        //
        //    Output, int TAU_SOBOL, the value TAU.
        //
        {
            int DIM_MAX = 13;

            int tau;
            int[] tau_table = {
                                0,  0,  1,  3,  5, 
                                8, 11, 15, 19, 23, 
                                27, 31, 35
                              };

            if ( 1 <= dim_num && dim_num <= DIM_MAX )
            {
                tau = tau_table[dim_num-1];
            }
            else
            {
                tau = - 1;
            }

            return tau;
        }
    }
}