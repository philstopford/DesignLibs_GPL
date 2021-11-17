namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_polarize ( int n, double[] a, double[] p, ref double[] a_normal,
            ref double[] a_parallel, int aIndex = 0, int pIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The (nonzero) vector P defines a direction.
        //
        //    The vector A can be written as the sum
        //
        //      A = A_normal + A_parallel
        //
        //    where A_parallel is a linear multiple of P, and A_normal
        //    is perpendicular to P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], the vector to be polarized.
        //
        //    Input, double P[N], the polarizing direction.
        //
        //    Output, double A_NORMAL[N], A_PARALLEL[N], the normal
        //    and parallel components of A.
        //
    {
        int i;

        double p_norm = r8vec_norm ( n, p, pIndex );

        switch (p_norm)
        {
            case 0.0:
            {
                for ( i = 0; i < n; i++ )
                {
                    a_normal[i] = a[aIndex + i];
                }
                for ( i = 0; i < n; i++ )
                {
                    a_parallel[i] = 0.0;
                }
                return;
            }
        }
        double a_dot_p = r8vec_dot_product ( n, a, p, aIndex, pIndex ) / p_norm;

        for ( i = 0; i < n; i++ )
        {
            a_parallel[i] = a_dot_p * p[pIndex + i] / p_norm;
        }

        for ( i = 0; i < n; i++ )
        {
            a_normal[i] = a[aIndex + i] - a_parallel[i];
        }
    }

}