using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static double[] basis_map_3d ( double[] u, double[] v )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MAP_3D computes the matrix which maps one basis to another.
        //
        //  Discussion:
        //
        //    As long as the vectors U1, U2 and U3 are linearly independent,
        //    a matrix A will be computed that maps U1 to V1, U2 to V2, and
        //    U3 to V3.
        //
        //    Depending on the values of the vectors, A may represent a
        //    rotation, reflection, dilation, project, or a combination of these
        //    basic linear transformations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double U[3*3], the matrix whose columns are three "domain" or "preimage"
        //    vectors, which should be linearly independent.
        //
        //    Input, double V[3*3], the matrix whose columns are the three "range" or
        //    "image" vectors.
        //
        //    Output, double BASIS_MAP_3D[3*3], a matrix with the property that A * U1 = V1,
        //    A * U2 = V2 and A * U3 = V3.
        //
    {
        double[] a;
        double[] c;
        int i;
        int j;
        int k;
        //
        //  Compute C = the inverse of U.
        //
        c = typeMethods.r8mat_inverse_3d ( u );

        switch (c)
        {
            case null:
                return null;
        }
        //
        //  A = V * inverse ( U ).
        //
        a = new double[3*3];

        for ( j = 0; j < 3; j++ )
        {
            for ( i = 0; i < 3; i++ )
            {
                a[i+j*3] = 0.0;
                for ( k = 0; k < 3; k++ )
                {
                    a[i+j*3] += v[i+k*3] * c[k+j*3];
                }
            }
        }
            
        return a;
    }
}