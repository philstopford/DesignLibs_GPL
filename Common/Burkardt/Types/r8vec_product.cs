namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8vec_product(int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRODUCT returns the product of the entries of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A[N], the vector.
        //
        //    Output, double R8VEC_PRODUCT, the product of the vector.
        //
    {
        int i;

        double product = 1.0;
        for (i = 0; i < n; i++)
        {
            product *= a[i];
        }

        return product;
    }

    public static double[] r8vec_vector_triple_product ( double[] v1, double[] v2, double[] v3 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_VECTOR_TRIPLE_PRODUCT computes the vector triple product.
        //
        //  Discussion:
        //
        //    VTRIPLE = V1 x (V2 x V3)
        //
        //    VTRIPLE is a vector perpendicular to V1, lying in the plane
        //    spanned by V2 and V3.  The norm of VTRIPLE is the product
        //    of the norms of V1, V2 and V3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double V1[3], V2[3], V3[3], the coordinates
        //    of the three vectors.
        //
        //    Output, double R8VEC_VECTOR_TRIPLE_PRODUCT[3], the vector triple product.
        //
    {
        double[] v23 = r8vec_cross_product_3d ( v2, v3 );

        double[] v123 = r8vec_cross_product_3d ( v1, v23 );

        return v123;
    }
        
    public static double r8vec_scalar_triple_product ( double[] v1, double[] v2, double[] v3 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SCALAR_TRIPLE_PRODUCT computes the scalar triple product.
        //
        //  Discussion:
        //
        //    STRIPLE = V1 dot ( V2 x V3 ).
        //
        //    STRIPLE is the volume of the parallelogram whose sides are
        //    formed by V1, V2 and V3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double V1[3], V2[3], V3[3], the three vectors.
        //
        //    Output, double R8VEC_SCALAR_TRIPLE_PRODUCT, the scalar
        //    triple product.
        //
    {
        double value = v1[0] * ( v2[1] * v3[2] - v2[2] * v3[1] )
                       + v1[1] * ( v2[2] * v3[0] - v2[0] * v3[2] )
                       + v1[2] * ( v2[0] * v3[1] - v2[1] * v3[0] );

        return value;
    }


}