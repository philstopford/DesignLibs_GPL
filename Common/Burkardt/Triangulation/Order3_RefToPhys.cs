namespace Burkardt.TriangulationNS;

public static partial class Triangulation
{
    public static void triangle_order3_reference_to_physical ( double[] t, int n,
            double[] ref_, ref double[] phy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL maps reference points to physical points.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 3 physical triangle and a point
        //    (XSI,ETA) in the reference triangle, the routine computes the value
        //    of the corresponding image point (X,Y) in physical space.
        //
        //    Note that this routine may also be appropriate for an order 6
        //    triangle, if the mapping between reference and physical space
        //    is linear.  This implies, in particular, that the sides of the
        //    image triangle are straight and that the "midside" nodes in the
        //    physical triangle are halfway along the sides of
        //    the physical triangle.
        //
        //  Reference Element T3:
        //
        //    |
        //    1  3
        //    |  |.
        //    |  | .
        //    S  |  .
        //    |  |   .
        //    |  |    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the coordinates of the vertices.
        //    The vertices are assumed to be the images of (0,0), (1,0) and
        //    (0,1) respectively.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double REF[2*N], points in the reference triangle.
        //
        //    Output, double PHY[2*N], corresponding points in the
        //    physical triangle.
        //
    {
        int i;

        for ( i = 0; i < 2; i++ )
        {
            int j;
            for ( j = 0; j < n; j++ )
            {
                phy[i+j*2] = t[i+0*2] * ( 1.0 - ref_[0+j*2] - ref_[1+j*2] )
                             + t[i+1*2] *       + ref_[0+j*2]
                             + t[i+2*2] *                    + ref_[1+j*2];
            }
        }
    }

}