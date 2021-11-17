namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8vec_cross_product ( double[] v1, double[] v2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CROSS_PRODUCT computes the cross product of two R8VEC's in 3D.
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
        //    07 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double V1[3], V2[3], the coordinates of the vectors.
        //
        //    Output, double R8VEC_CROSS_PRODUCT[3], the cross product vector.
        //
    {
        double[] v3 = new double[3];

        v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

        return v3;
    }
        
    public static double r8vec_cross_product_2d(double[] v1, double[] v2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CROSS_PRODUCT_2D finds the cross product of a pair of R8VEC's in 2D.
        //
        //  Discussion:
        //
        //    Strictly speaking, the vectors lie in the (X,Y) plane, and
        //    the cross product here is a vector in the Z direction.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double V1[2], V2[2], the vectors.
        //
        //    Output, double R8VEC_CROSS_PRODUCT_2D, the Z component of the cross product
        //    of V1 and V2.
        //
    {
        double value = v1[0] * v2[1] - v1[1] * v2[0];

        return value;
    }

    public static double r8vec_cross_product_affine_2d(double[] v0, double[] v1,
            double[] v2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.
        //
        //  Discussion:
        //
        //    Strictly speaking, the vectors lie in the (X,Y) plane, and
        //    the cross product here is a vector in the Z direction.
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
        //    Input, double V0[2], the base vector.
        //
        //    Input, double V1[2], V2[2], the vectors.
        //
        //    Output, double R8VEC_CROSS_PRODUCT_AFFINE_2D, the Z component of the
        //    cross product of V1 and V2.
        //
    {
        double value = (v1[0] - v0[0]) * (v2[1] - v0[1])
                       - (v2[0] - v0[0]) * (v1[1] - v0[1]);

        return value;
    }

    public static double[] r8vec_cross_product_3d(double[] v1, double[] v2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.
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
        //    07 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double V1[3], V2[3], the coordinates of the vectors.
        //
        //    Output, double R8VEC_CROSS_PRODUCT_3D[3], the cross product vector.
        //
    {
        double[] v3 = new double[3];

        v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

        return v3;
    }

    public static double[] r8vec_cross_product_affine_3d(double[] v0, double[] v1,
            double[] v2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.
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
        //    27 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double V0[3], the base vector.
        //
        //    Input, double V1[3], V2[3], the coordinates of the vectors.
        //
        //    Output, double R8VEC_CROSS_PRODUCT_AFFINE_3D[3], the cross product vector.
        //
    {
        double[] v3 = new double[3];

        v3[0] =
            (v1[1] - v0[1]) * (v2[2] - v0[2])
            - (v2[1] - v0[1]) * (v1[2] - v0[2]);

        v3[1] =
            (v1[2] - v0[2]) * (v2[0] - v0[0])
            - (v2[2] - v0[2]) * (v1[0] - v0[0]);

        v3[2] =
            (v1[0] - v0[0]) * (v2[1] - v0[1])
            - (v2[0] - v0[0]) * (v1[1] - v0[1]);

        return v3;
    }

}