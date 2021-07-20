namespace Burkardt.Types
{
    public static partial class typeMethods
    {
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
            double[] v123;
            double[] v23;

            v23 = r8vec_cross_product_3d ( v2, v3 );

            v123 = r8vec_cross_product_3d ( v1, v23 );

            return v123;
        }

    }
}