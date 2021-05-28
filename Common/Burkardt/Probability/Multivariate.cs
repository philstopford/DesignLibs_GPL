namespace Burkardt.Probability
{
    public static class Multivariate
    {
        public static double[] multivariate_normal_sample ( int n, double[] mean, 
        double[] covar_factor, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTIVARIATE_NORMAL_SAMPLE samples the Multivariate Normal PDF.
        //
        //  Discussion:
        //
        //    PDF ( Mean(1:N), S(1:N,1:N); X(1:N) ) = 
        //      1 / ( 2 * pi ) ^ ( N / 2 ) * 1 / det ( S )
        //      * exp ( - ( X - Mean )' * inverse ( S ) * ( X - Mean ) / 2 )
        //
        //    Here,
        //
        //      X is the argument vector of length N,
        //      Mean is the mean vector of length N,
        //      S is an N by N positive definite symmetric covariance matrix.
        //
        //    The properties of S guarantee that it has a lower triangular
        //    matrix L, the Cholesky factor, such that S = L * L'.  It is the
        //    matrix L, rather than S, that is required by this routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jerry Banks, editor,
        //    Handbook of Simulation,
        //    Engineering and Management Press Books, 1998, page 167.
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double MEAN[N], the mean vector.
        //
        //    Input, double COVAR_FACTOR[N*N], the lower triangular Cholesky
        //    factor L of the covariance matrix S.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double MULTIVARIATE_NORMAL_SAMPLE[N], a sample point
        //    of the distribution.
        //
        {
            double[] x = new double[n];

            for (int i = 0; i < n; i++ )
            {
                double z = Normal.normal_01_sample ( ref seed );

                x[i] = mean[i];

                for (int j = 0; j <= i; j++ )
                {
                    x[i] = x[i] + covar_factor[i+j*n] * z;
                }
            }

            return x;
        }
    }
}