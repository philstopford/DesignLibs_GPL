using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.RandomMatrix;

public static class PositiveDefiniteSymmetric
{
    public static double[] pds_random ( int n, ref typeMethods.r8NormalData data, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PDS_RANDOM returns the PDS_RANDOM matrix.
        //
        //  Discussion:
        //
        //    The matrix is a "random" positive definite symmetric matrix.
        //
        //    The matrix returned will have eigenvalues in the range [0,1].
        //
        //  Properties:
        //
        //    A is symmetric: A' = A.
        //
        //    A is positive definite: 0 < x'*A*x for nonzero x.
        //
        //    The eigenvalues of A will be real.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double PDS_RANDOM[N*N], the matrix.
        //
    {
        int j;

        double[] a = new double[n*n];
        //
        //  Get a random set of eigenvalues.
        //
        double[] lambda = UniformRNG.r8vec_uniform_01_new ( n, ref seed );
        //
        //  Get a random orthogonal matrix Q.
        //
        double[] q = Orthogonal.orth_random ( n, ref data, ref seed );
        //
        //  Set A = Q * Lambda * Q'.
        //
        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < n; i++ )
            {
                a[i+j*n] = 0.0;
                int k;
                for ( k = 0; k < n; k++ )
                {
                    a[i+j*n] += q[i+k*n] * lambda[k] * q[j+k*n];
                }
            }
        }
        return a;
    }
}