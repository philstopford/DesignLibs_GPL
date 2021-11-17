namespace Burkardt.TriangulationNS;

public static partial class Triangulation
{
    public static void triangle_order3_physical_to_reference(double[] t, int n,
            double[] phy, ref double[] ref_ )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE maps physical points to reference points.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 3 physical triangle and a point
        //    (X,Y) in the physical triangle, the routine computes the value
        //    of the corresponding image point (XSI,ETA) in reference space.
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
        //    Input, double T[2*3], the X and Y coordinates
        //    of the vertices.  The vertices are assumed to be the images of
        //    (0,0), (1,0) and (0,1) respectively.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double PHY[2*N], the coordinates of physical points
        //    to be transformed.
        //
        //    Output, double REF[2*N], the coordinates of the corresponding
        //    points in the reference space.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {

            ref_[
                    0 + j * 2] = ((t[1 + 2 * 2] - t[1 + 0 * 2]) * (phy[0 + j * 2] - t[0 + 0 * 2])
                                  - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (phy[1 + j * 2] - t[1 + 0 * 2]))
                                 / ((t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2])
                                    - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2]));

            ref_[
                    1 + j * 2] = ((t[0 + 1 * 2] - t[0 + 0 * 2]) * (phy[1 + j * 2] - t[1 + 0 * 2])
                                  - (t[1 + 1 * 2] - t[1 + 0 * 2]) * (phy[0 + j * 2] - t[0 + 0 * 2]))
                                 / ((t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2])
                                    - (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2]));
        }
    }
}