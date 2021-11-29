using Burkardt.Types;

namespace Burkardt.FEM;

public static class PhysicalToRef
{
    public static void physical_to_reference_t3(double[] t, int n, double[] phy, ref double[] ref_ )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PHYSICAL_TO_REFERENCE_T3 maps physical points to reference points.
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
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  .   .
        //    |  .    .
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
        //    Input, double[] t, the X and Y coordinates
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
        
    public static double[] physical_to_reference_tet4 ( double[] t, int n, double[] phy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PHYSICAL_TO_REFERENCE_TET4 maps physical points to reference points.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 4 physical tetrahedron and a point 
        //    (X,Y,Z) in the physical tetrahedron, the routine computes the value 
        //    of the corresponding point (R,S,T) in the reference tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*4], the coordinates of the vertices of the
        //    physical tetrahedron.  The vertices are assumed to be the images of
        //    (1,0,0), (0,1,0), (0,0,1) and (0,0,0) respectively.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double PHY[3*N], the coordinates of physical points
        //    to be transformed.
        //
        //    Output, double PHYSICAL_TO_REFERENCE[3*N], the coordinates of the 
        //    corresponding points in the reference tetrahedron.
        //
    {
        double[] a = new double[3*3];
        int i;
        int j;

        for ( j = 0; j < 3; j++ )
        {
            for ( i = 0; i < 3; i++ )
            {
                a[i+j*3] = t[i+j*3] - t[i+3*3];
            }
        }

        double[] ref_ = new double[3*n];

        for ( j = 0; j < n; j++ )
        {
            for ( i = 0; i < 3; i++ )
            {
                ref_[i+j*3] = phy[i+j*3] - t[i+3*3];
            }
        }

        typeMethods.r8ge_fss ( 3, ref a, n, ref ref_ );

        return ref_;
    }
}