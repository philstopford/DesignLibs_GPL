using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TetrahedronNS;

public static class Geometry
{
    public static double[] tetrahedron_barycentric_3d(double[] tetra, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_BARYCENTRIC_3D returns the barycentric coordinates of a point in 3D.
        //
        //  Discussion:
        //
        //    The barycentric coordinates of a point P with respect to
        //    a tetrahedron are a set of four values C(1:4), each associated
        //    with a vertex of the tetrahedron.  The values must sum to 1.
        //    If all the values are between 0 and 1, the point is contained
        //    within the tetrahedron.
        //
        //    The barycentric coordinate of point X related to vertex A can be
        //    interpreted as the ratio of the volume of the tetrahedron with
        //    vertex A replaced by vertex X to the volume of the original
        //    tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Input, double P[3], the point to be checked.
        //
        //    Output, double C[4], the barycentric coordinates of the point with
        //    respect to the tetrahedron.
        //
    {
        int N = 3;
        int RHS_NUM = 1;

        double[] a = new double[N * (N + RHS_NUM)];
        double[] c;
        int info;
        //
        //  Set up the linear system
        //
        //    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
        //    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
        //    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
        //
        //  which is satisfied by the barycentric coordinates.
        //

        a[0 + 0 * N] = tetra[0 + 1 * 3] - tetra[0 + 0 * 3];
        a[1 + 0 * N] = tetra[1 + 1 * 3] - tetra[1 + 0 * 3];
        a[2 + 0 * N] = tetra[2 + 1 * 3] - tetra[2 + 0 * 3];

        a[0 + 1 * N] = tetra[0 + 2 * 3] - tetra[0 + 0 * 3];
        a[1 + 1 * N] = tetra[1 + 2 * 3] - tetra[1 + 0 * 3];
        a[2 + 1 * N] = tetra[2 + 2 * 3] - tetra[2 + 0 * 3];

        a[0 + 2 * N] = tetra[0 + 3 * 3] - tetra[0 + 0 * 3];
        a[1 + 2 * N] = tetra[1 + 3 * 3] - tetra[1 + 0 * 3];
        a[2 + 2 * N] = tetra[2 + 3 * 3] - tetra[2 + 0 * 3];

        a[0 + 3 * N] = p[0] - tetra[0 + 0 * 3];
        a[1 + 3 * N] = p[1] - tetra[1 + 0 * 3];
        a[2 + 3 * N] = p[2] - tetra[2 + 0 * 3];
        //
        //  Solve the linear system.
        //
        info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_BARYCENTRIC_3D - Fatal error!");
            Console.WriteLine("  The linear system is singular.");
            Console.WriteLine("  The input data does not form a proper tetrahedron.");
            return null;
        }

        c = new double[4];

        c[1] = a[0 + 3 * N];
        c[2] = a[1 + 3 * N];
        c[3] = a[2 + 3 * N];

        c[0] = 1.0 - c[1] - c[2] - c[3];

        return c;
    }

    public static double[] tetrahedron_centroid_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_CENTROID_3D computes the centroid of a tetrahedron in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Output, double TETRAHEDRON_CENTROID_3D[3], the coordinates of the centroid.
        //
    {
        int DIM_NUM = 3;

        double[] centroid;

        centroid = new double[3];

        centroid[0] = 0.25 * (tetra[0 + 0 * DIM_NUM] + tetra[0 + 1 * DIM_NUM]
                                                     + tetra[0 + 2 * DIM_NUM] + tetra[0 + 3 * DIM_NUM]);
        centroid[1] = 0.25 * (tetra[1 + 0 * DIM_NUM] + tetra[1 + 1 * DIM_NUM]
                                                     + tetra[1 + 2 * DIM_NUM] + tetra[1 + 3 * DIM_NUM]);
        centroid[2] = 0.25 * (tetra[2 + 0 * DIM_NUM] + tetra[2 + 1 * DIM_NUM]
                                                     + tetra[2 + 2 * DIM_NUM] + tetra[2 + 3 * DIM_NUM]);

        return centroid;
    }

    public static void tetrahedron_circumsphere_3d(double[] tetra, ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
        //
        //  Discussion:
        //
        //    The circumsphere, or circumscribed sphere, of a tetrahedron is the sphere that
        //    passes through the four vertices.  The circumsphere is not necessarily
        //    the smallest sphere that contains the tetrahedron.
        //
        //    Surprisingly, the diameter of the sphere can be found by solving
        //    a 3 by 3 linear system.  This is because the vectors P2 - P1,
        //    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
        //    right triangle with the diameter through P1.  Hence, the dot product of
        //    P2 - P1 with that diameter is equal to the square of the length
        //    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
        //    the diameter vector originating at P1, and hence the radius and
        //    center.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Output, double *R, PC[3], the coordinates of the center of the
        //    circumscribed sphere, and its radius.  If the linear system is
        //    singular, then R = -1, PC[] = 0.
        //
    {
        int DIM_NUM = 3;
        int RHS_NUM = 1;

        double[] a = new double[DIM_NUM * (DIM_NUM + RHS_NUM)];
        int info;
        //
        //  Set up the linear system.
        //
        a[0 + 0 * 3] = tetra[0 + 1 * 3] - tetra[0 + 0 * 3];
        a[0 + 1 * 3] = tetra[1 + 1 * 3] - tetra[1 + 0 * 3];
        a[0 + 2 * 3] = tetra[2 + 1 * 3] - tetra[2 + 0 * 3];
        a[0 + 3 * 3] = Math.Pow(tetra[0 + 1 * 3] - tetra[0 + 0 * 3], 2)
                       + Math.Pow(tetra[1 + 1 * 3] - tetra[1 + 0 * 3], 2)
                       + Math.Pow(tetra[2 + 1 * 3] - tetra[2 + 0 * 3], 2);

        a[1 + 0 * 3] = tetra[0 + 2 * 3] - tetra[0 + 0 * 3];
        a[1 + 1 * 3] = tetra[1 + 2 * 3] - tetra[1 + 0 * 3];
        a[1 + 2 * 3] = tetra[2 + 2 * 3] - tetra[2 + 0 * 3];
        a[1 + 3 * 3] = Math.Pow(tetra[0 + 2 * 3] - tetra[0 + 0 * 3], 2)
                       + Math.Pow(tetra[1 + 2 * 3] - tetra[1 + 0 * 3], 2)
                       + Math.Pow(tetra[2 + 2 * 3] - tetra[2 + 0 * 3], 2);

        a[2 + 0 * 3] = tetra[0 + 3 * 3] - tetra[0 + 0 * 3];
        a[2 + 1 * 3] = tetra[1 + 3 * 3] - tetra[1 + 0 * 3];
        a[2 + 2 * 3] = tetra[2 + 3 * 3] - tetra[2 + 0 * 3];
        a[2 + 3 * 3] = Math.Pow(tetra[0 + 3 * 3] - tetra[0 + 0 * 3], 2)
                       + Math.Pow(tetra[1 + 3 * 3] - tetra[1 + 0 * 3], 2)
                       + Math.Pow(tetra[2 + 3 * 3] - tetra[2 + 0 * 3], 2);
        //
        //  Solve the linear system.
        //
        info = typeMethods.r8mat_solve(DIM_NUM, RHS_NUM, ref a);
        //
        //  If the system was singular, return a consolation prize.
        //
        if (info != 0)
        {
            r = -1.0;
            typeMethods.r8vec_zero(DIM_NUM, ref pc);
            return;
        }

        //
        //  Compute the radius and center.
        //
        r = 0.5 * Math.Sqrt
        (a[0 + 3 * 3] * a[0 + 3 * 3]
         + a[1 + 3 * 3] * a[1 + 3 * 3]
         + a[2 + 3 * 3] * a[2 + 3 * 3]);

        pc[0] = tetra[0 + 0 * 3] + 0.5 * a[0 + 3 * 3];
        pc[1] = tetra[1 + 0 * 3] + 0.5 * a[1 + 3 * 3];
        pc[2] = tetra[2 + 0 * 3] + 0.5 * a[2 + 3 * 3];
    }

    public static bool tetrahedron_contains_point_3d(double[] tetra, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_CONTAINS_POINT_3D: a tetrahedron contains a point in 3D.
        //
        //  Discussion:
        //
        //    Thanks to Saiful Akbar for pointing out that the array of barycentric
        //    coordinated was not being deleted!  29 January 2006
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Input, double P[3], the point to be checked.
        //
        //    Output, bool TETRAHEDRON_CONTAINS_POINT_3D, is TRUE if the point is inside
        //    the tetrahedron or on its boundary, and FALSE otherwise.
        //
    {
        double[] c;
        bool value;

        c = tetrahedron_barycentric_3d(tetra, p);
        //
        //  If the point is in the tetrahedron, its barycentric coordinates
        //  must be nonnegative.
        //
        value = 0.0 <= c[0] &&
                0.0 <= c[1] &&
                0.0 <= c[2] &&
                0.0 <= c[3];

        return value;
    }

    public static double[] tetrahedron_dihedral_angles_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_DIHEDRAL_ANGLES_3D computes dihedral angles of a tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron,
        //    which can be labeled as A, B, C and D.
        //
        //    Output, double TETRAHEDRON_DIHEDRAL_ANGLES_3D[6], the dihedral angles
        //    along the axes AB, AC, AD, BC, BD and CD, respectively.
        //
    {
        double[] ab = new double[3];
        double[] abc_normal;
        double[] abd_normal;
        double[] ac = new double[3];
        double[] acd_normal;
        double[] ad = new double[3];
        double[] angle;
        double[] bc = new double[3];
        double[] bcd_normal;
        double[] bd = new double[3];
        int i;

        for (i = 0; i < 3; i++)
        {
            ab[i] = tetra[i + 1 * 3] - tetra[i + 0 * 3];
            ac[i] = tetra[i + 2 * 3] - tetra[i + 0 * 3];
            ad[i] = tetra[i + 3 * 3] - tetra[i + 0 * 3];
            bc[i] = tetra[i + 2 * 3] - tetra[i + 1 * 3];
            bd[i] = tetra[i + 3 * 3] - tetra[i + 1 * 3];
        }

        abc_normal = typeMethods.r8vec_cross_product_3d(ac, ab);
        abd_normal = typeMethods.r8vec_cross_product_3d(ab, ad);
        acd_normal = typeMethods.r8vec_cross_product_3d(ad, ac);
        bcd_normal = typeMethods.r8vec_cross_product_3d(bc, bd);

        angle = new double[6];

        angle[0] = typeMethods.r8vec_angle_3d(abc_normal, abd_normal);
        angle[1] = typeMethods.r8vec_angle_3d(abc_normal, acd_normal);
        angle[2] = typeMethods.r8vec_angle_3d(abd_normal, acd_normal);
        angle[3] = typeMethods.r8vec_angle_3d(abc_normal, bcd_normal);
        angle[4] = typeMethods.r8vec_angle_3d(abd_normal, bcd_normal);
        angle[5] = typeMethods.r8vec_angle_3d(acd_normal, bcd_normal);

        for (i = 0; i < 6; i++)
        {
            angle[i] = Math.PI - angle[i];
        }

        return angle;
    }

    public static double[] tetrahedron_edge_length_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the tetrahedron vertices.
        //
        //    Output, double EDGE_LENGTH[6], the length of the edges.
        //
    {
        int DIM_NUM = 3;

        double[] edge_length;
        int i;
        int j1;
        int j2;
        int k;
        double[] v = new double[DIM_NUM];

        edge_length = new double[6];

        k = 0;
        for (j1 = 0; j1 < 3; j1++)
        {
            for (j2 = j1 + 1; j2 < 4; j2++)
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    v[i] = tetra[i + j2 * DIM_NUM] - tetra[i + j1 * DIM_NUM];
                }

                edge_length[k] = typeMethods.r8vec_norm(DIM_NUM, v);
                k += 1;
            }
        }

        return edge_length;
    }

    public static void tetrahedron_face_angles_3d(double[] tetra, ref double[] angles)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_FACE_ANGLES_3D returns the 12 face angles of a tetrahedron 3D.
        //
        //  Discussion:
        //
        //    The tetrahedron has 4 triangular faces.  This routine computes the
        //    3 planar angles associated with each face.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4] the tetrahedron vertices.
        //
        //    Output, double ANGLES[3*4], the face angles.
        //
    {
        double[] tri;

        tri = new double[3 * 3];
        //
        //  Face 123
        //
        tri[0 + 0 * 3] = tetra[0 + 0 * 3];
        tri[1 + 0 * 3] = tetra[1 + 0 * 3];
        tri[2 + 0 * 3] = tetra[2 + 0 * 3];
        tri[0 + 1 * 3] = tetra[0 + 1 * 3];
        tri[1 + 1 * 3] = tetra[1 + 1 * 3];
        tri[2 + 1 * 3] = tetra[2 + 1 * 3];
        tri[0 + 2 * 3] = tetra[0 + 2 * 3];
        tri[1 + 2 * 3] = tetra[1 + 2 * 3];
        tri[2 + 2 * 3] = tetra[2 + 2 * 3];

        TriangleNS.Geometry.triangle_angles_3d(tri, ref angles);
        //
        //  Face 124
        //
        tri[0 + 0 * 3] = tetra[0 + 0 * 3];
        tri[1 + 0 * 3] = tetra[1 + 0 * 3];
        tri[2 + 0 * 3] = tetra[2 + 0 * 3];
        tri[0 + 1 * 3] = tetra[0 + 1 * 3];
        tri[1 + 1 * 3] = tetra[1 + 1 * 3];
        tri[2 + 1 * 3] = tetra[2 + 1 * 3];
        tri[0 + 2 * 3] = tetra[0 + 3 * 3];
        tri[1 + 2 * 3] = tetra[1 + 3 * 3];
        tri[2 + 2 * 3] = tetra[2 + 3 * 3];

        TriangleNS.Geometry.triangle_angles_3d(tri, ref angles, +3);
        //
        //  Face 134
        //
        tri[0 + 0 * 3] = tetra[0 + 0 * 3];
        tri[1 + 0 * 3] = tetra[1 + 0 * 3];
        tri[2 + 0 * 3] = tetra[2 + 0 * 3];
        tri[0 + 1 * 3] = tetra[0 + 2 * 3];
        tri[1 + 1 * 3] = tetra[1 + 2 * 3];
        tri[2 + 1 * 3] = tetra[2 + 2 * 3];
        tri[0 + 2 * 3] = tetra[0 + 3 * 3];
        tri[1 + 2 * 3] = tetra[1 + 3 * 3];
        tri[2 + 2 * 3] = tetra[2 + 3 * 3];

        TriangleNS.Geometry.triangle_angles_3d(tri, ref angles, +6);
        //
        //  Face 234
        //
        tri[0 + 0 * 3] = tetra[0 + 1 * 3];
        tri[1 + 0 * 3] = tetra[1 + 1 * 3];
        tri[2 + 0 * 3] = tetra[2 + 1 * 3];
        tri[0 + 1 * 3] = tetra[0 + 2 * 3];
        tri[1 + 1 * 3] = tetra[1 + 2 * 3];
        tri[2 + 1 * 3] = tetra[2 + 2 * 3];
        tri[0 + 2 * 3] = tetra[0 + 3 * 3];
        tri[1 + 2 * 3] = tetra[1 + 3 * 3];
        tri[2 + 2 * 3] = tetra[2 + 3 * 3];

        TriangleNS.Geometry.triangle_angles_3d(tri, ref angles, +9);
    }

    public static void tetrahedron_face_areas_3d(double[] tetra, ref double[] areas)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_FACE_AREAS_3D returns the 4 face areas of a tetrahedron 3D.
        //
        //  Discussion:
        //
        //    The tetrahedron has 4 triangular faces.  This routine computes the
        //    areas associated with each face.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4] the tetrahedron vertices.
        //
        //    Output, double AREAS[4], the face areas.
        //
    {
        double[] tri;

        tri = new double[3 * 3];
        //
        //  Face 123
        //
        tri[0 + 0 * 3] = tetra[0 + 0 * 3];
        tri[1 + 0 * 3] = tetra[1 + 0 * 3];
        tri[2 + 0 * 3] = tetra[2 + 0 * 3];
        tri[0 + 1 * 3] = tetra[0 + 1 * 3];
        tri[1 + 1 * 3] = tetra[1 + 1 * 3];
        tri[2 + 1 * 3] = tetra[2 + 1 * 3];
        tri[0 + 2 * 3] = tetra[0 + 2 * 3];
        tri[1 + 2 * 3] = tetra[1 + 2 * 3];
        tri[2 + 2 * 3] = tetra[2 + 2 * 3];

        areas[0] = TriangleNS.Geometry.triangle_area_3d(tri);
        //
        //  Face 124
        //
        tri[0 + 0 * 3] = tetra[0 + 0 * 3];
        tri[1 + 0 * 3] = tetra[1 + 0 * 3];
        tri[2 + 0 * 3] = tetra[2 + 0 * 3];
        tri[0 + 1 * 3] = tetra[0 + 1 * 3];
        tri[1 + 1 * 3] = tetra[1 + 1 * 3];
        tri[2 + 1 * 3] = tetra[2 + 1 * 3];
        tri[0 + 2 * 3] = tetra[0 + 3 * 3];
        tri[1 + 2 * 3] = tetra[1 + 3 * 3];
        tri[2 + 2 * 3] = tetra[2 + 3 * 3];

        areas[1] = TriangleNS.Geometry.triangle_area_3d(tri);
        //
        //  Face 134
        //
        tri[0 + 0 * 3] = tetra[0 + 0 * 3];
        tri[1 + 0 * 3] = tetra[1 + 0 * 3];
        tri[2 + 0 * 3] = tetra[2 + 0 * 3];
        tri[0 + 1 * 3] = tetra[0 + 2 * 3];
        tri[1 + 1 * 3] = tetra[1 + 2 * 3];
        tri[2 + 1 * 3] = tetra[2 + 2 * 3];
        tri[0 + 2 * 3] = tetra[0 + 3 * 3];
        tri[1 + 2 * 3] = tetra[1 + 3 * 3];
        tri[2 + 2 * 3] = tetra[2 + 3 * 3];

        areas[2] = TriangleNS.Geometry.triangle_area_3d(tri);
        //
        //  Face 234
        //
        tri[0 + 0 * 3] = tetra[0 + 1 * 3];
        tri[1 + 0 * 3] = tetra[1 + 1 * 3];
        tri[2 + 0 * 3] = tetra[2 + 1 * 3];
        tri[0 + 1 * 3] = tetra[0 + 2 * 3];
        tri[1 + 1 * 3] = tetra[1 + 2 * 3];
        tri[2 + 1 * 3] = tetra[2 + 2 * 3];
        tri[0 + 2 * 3] = tetra[0 + 3 * 3];
        tri[1 + 2 * 3] = tetra[1 + 3 * 3];
        tri[2 + 2 * 3] = tetra[2 + 3 * 3];

        areas[3] = TriangleNS.Geometry.triangle_area_3d(tri);

    }

    public static void tetrahedron_insphere_3d(double[] tetra, ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
        //
        //  Discussion:
        //
        //    The insphere of a tetrahedron is the inscribed sphere, which touches
        //    each face of the tetrahedron at a single point.
        //
        //    The points of contact are the centroids of the triangular faces
        //    of the tetrahedron.  Therefore, the point of contact for a face
        //    can be computed as the average of the vertices of that face.
        //
        //    The sphere can then be determined as the unique sphere through
        //    the four given centroids.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Schneider, David Eberly,
        //    Geometric Tools for Computer Graphics,
        //    Elsevier, 2002,
        //    ISBN: 1558605940,
        //    LC: T385.G6974.
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Output, double *R, PC[3], the radius and the center
        //    of the sphere.
        //
    {
        int DIM_NUM = 3;

        double[] b = new double[4 * 4];
        double gamma;
        int i;
        int j;
        double l123;
        double l124;
        double l134;
        double l234;
        double[] n123;
        double[] n124;
        double[] n134;
        double[] n234;
        double[] v21 = new double[DIM_NUM];
        double[] v31 = new double[DIM_NUM];
        double[] v41 = new double[DIM_NUM];
        double[] v32 = new double[DIM_NUM];
        double[] v42 = new double[DIM_NUM];
        //double v43[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            v21[i] = tetra[i + 1 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v31[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v41[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v32[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v42[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
        }
        //for ( i = 0; i < DIM_NUM; i++ )
        //{
        //  v43[i] = tetra[i+3*DIM_NUM] - tetra[i+2*DIM_NUM];
        //}

        n123 = typeMethods.r8vec_cross_product_3d(v21, v31);
        n124 = typeMethods.r8vec_cross_product_3d(v41, v21);
        n134 = typeMethods.r8vec_cross_product_3d(v31, v41);
        n234 = typeMethods.r8vec_cross_product_3d(v42, v32);

        l123 = typeMethods.r8vec_norm(DIM_NUM, n123);
        l124 = typeMethods.r8vec_norm(DIM_NUM, n124);
        l134 = typeMethods.r8vec_norm(DIM_NUM, n134);
        l234 = typeMethods.r8vec_norm(DIM_NUM, n234);

        for (i = 0; i < DIM_NUM; i++)
        {
            pc[i] = (l234 * tetra[i + 0 * DIM_NUM]
                     + l134 * tetra[i + 1 * DIM_NUM]
                     + l124 * tetra[i + 2 * DIM_NUM]
                     + l123 * tetra[i + 3 * DIM_NUM])
                    / (l234 + l134 + l124 + l123);
        }

        for (j = 0; j < 4; j++)
        {
            for (i = 0; i < DIM_NUM; i++)
            {
                b[i + j * 4] = tetra[i + j * DIM_NUM];
            }

            b[3 + j * 4] = 1.0;
        }

        gamma = Math.Abs(typeMethods.r8mat_det_4d(b));

        r = gamma / (l234 + l134 + l124 + l123);
    }

    public static void tetrahedron_lattice_layer_point_next(int[] c, ref int[] v, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_LATTICE_LAYER_POINT_NEXT: next tetrahedron lattice layer point.
        //
        //  Discussion:
        //
        //    The tetrahedron lattice layer L is bounded by the lines
        //
        //      0 <= X,
        //      0 <= Y,
        //      0 <= Z,
        //      L - 1 < X / C[0] + Y / C[1] + Z/C[2] <= L.
        //
        //    In particular, layer L = 0 always contains the single point (0,0).
        //
        //    This function returns, one at a time, the points that lie within
        //    a given tetrahedron lattice layer.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int C[4], coefficients defining the
        //    lattice layer in the first 3 entries, and the laver index in C[3].
        //    The coefficients should be positive, and C[3] must be nonnegative.
        //
        //    Input/output, int V[3].  On first call for a given layer,
        //    the input value of V is not important.  On a repeated call for the same
        //    layer, the input value of V should be the output value from the previous
        //    call.  On output, V contains the next lattice layer point.
        //
        //    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
        //    that this is the first call for a given layer.  Thereafter, the input
        //    value should be the output value from the previous call.  On output,
        //    MORE is TRUE if the returned value V is a new point.
        //    If the output value is FALSE, then no more points were found,
        //    and V was reset to 0, and the lattice layer has been exhausted.
        //
    {
        int c1n;
        int lhs;
        int n = 3;
        int rhs1;
        int rhs2;
        switch (c[3])
        {
            //
            //  Treat layer C[3] = 0 specially.
            //
            case 0:
            {
                switch (more)
                {
                    case false:
                        v[0] = 0;
                        v[1] = 0;
                        v[2] = 0;
                        more = true;
                        break;
                    default:
                        more = false;
                        break;
                }

                return;
            }
        }

        switch (more)
        {
            //
            //  Compute the first point.
            //
            case false:
                v[0] = (c[n] - 1) * c[0] + 1;
                v[1] = 0;
                v[2] = 0;
                more = true;
                break;
            default:
            {
                c1n = typeMethods.i4vec_lcm(n, c);

                rhs1 = c1n * (c[n] - 1);
                rhs2 = c1n * c[n];
                //
                //  Can we simply increase X?
                //
                v[0] += 1;

                lhs = c1n / c[0] * v[0]
                      + c1n / c[1] * v[1]
                      + c1n / c[2] * v[2];

                if (lhs <= rhs2)
                {
                }
                //
                //  No.  Increase Y, and set X so we just exceed RHS1, if possible.
                //
                else
                {
                    v[1] += 1;

                    v[0] = c[0] * (rhs1 - c1n / c[1] * v[1]
                                        - c1n / c[2] * v[2]) / c1n;
                    v[0] = Math.Max(v[0], 0);

                    lhs = c1n / c[0] * v[0]
                          + c1n / c[1] * v[1]
                          + c1n / c[2] * v[2];

                    if (lhs <= rhs1)
                    {
                        v[0] += 1;
                        lhs += c1n / c[0];
                    }

                    //
                    //  We have increased Y by 1.  Have we stayed below the upper bound?
                    //
                    if (lhs <= rhs2)
                    {
                    }
                    //
                    //  No.  Increase Z, and set X so we just exceed RHS1, if possible.
                    //
                    else
                    {
                        v[2] += 1;
                        v[1] = 0;
                        v[0] = c[0] * (rhs1 - c1n / c[1] * v[1]
                                            - c1n / c[2] * v[2]) / c1n;
                        v[0] = Math.Max(v[0], 0);

                        lhs = c1n / c[0] * v[0]
                              + c1n / c[1] * v[1]
                              + c1n / c[2] * v[2];

                        if (lhs <= rhs1)
                        {
                            v[0] += 1;
                            lhs += c1n / c[0];
                        }

                        if (lhs <= rhs2)
                        {
                        }
                        else
                        {
                            more = false;
                            v[0] = 0;
                            v[1] = 0;
                            v[2] = 0;
                        }
                    }
                }

                break;
            }
        }
    }

    public static void tetrahedron_lattice_point_next(int[] c, ref int[] v, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_LATTICE_POINT_NEXT returns the next tetrahedron lattice point.
        //
        //  Discussion:
        //
        //    The lattice tetrahedron is defined by the vertices:
        //
        //      (0,0,0), (C[3]/C[0],0,0), (0,C[3]/C[1],0) and (0,0,C[3]/C[2])
        //
        //    The lattice tetrahedron is bounded by the lines
        //
        //      0 <= X,
        //      0 <= Y
        //      0 <= Z,
        //      X / C[0] + Y / C[1] + Z / C[2] <= C[3]
        //
        //    Lattice points are listed one at a time, starting at the origin,
        //    with X increasing first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int C[4], coefficients defining the
        //    lattice tetrahedron.  These should be positive.
        //
        //    Input/output, int V[3].  On first call, the input
        //    value is not important.  On a repeated call, the input value should
        //    be the output value from the previous call.  On output, V contains
        //    the next lattice point.
        //
        //    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
        //    that this is the first call for a given tetrahedron.  Thereafter, the input
        //    value should be the output value from the previous call.  On output,
        //    MORE is TRUE if not only is the returned value V a lattice point,
        //    but the routine can be called again for another lattice point.
        //    If the output value is FALSE, then no more lattice points were found,
        //    and V was reset to 0, and the routine should not be called further
        //    for this tetrahedron.
        //
    {
        int c1n;
        int lhs;
        int n = 3;
        int rhs;

        switch (more)
        {
            case false:
                v[0] = 0;
                v[1] = 0;
                v[2] = 0;
                more = true;
                break;
            default:
            {
                c1n = typeMethods.i4vec_lcm(n, c);

                rhs = c1n * c[n];

                lhs = c[1] * c[2] * v[0]
                      + c[0] * c[2] * v[1]
                      + c[0] * c[1] * v[2];

                if (lhs + c1n / c[0] <= rhs)
                {
                    v[0] += 1;
                }
                else
                {
                    lhs -= c1n * v[0] / c[0];
                    v[0] = 0;
                    if (lhs + c1n / c[1] <= rhs)
                    {
                        v[1] += 1;
                    }
                    else
                    {
                        lhs -= c1n * v[1] / c[1];
                        v[1] = 0;
                        if (lhs + c1n / c[2] <= rhs)
                        {
                            v[2] += 1;
                        }
                        else
                        {
                            v[2] = 0;
                            more = false;
                        }
                    }
                }

                break;
            }
        }
    }

    public static double tetrahedron_quality1_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
        //
        //  Discussion:
        //
        //    The quality of a tetrahedron is 3.0 times the ratio of the radius of
        //    the inscribed sphere divided by that of the circumscribed sphere.
        //
        //    An equilateral tetrahredron achieves the maximum possible quality of 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the tetrahedron vertices.
        //
        //    Output, double TETRAHEDRON_QUALITY1_3D, the quality of the tetrahedron.
        //
    {
        int DIM_NUM = 3;

        double[] pc = new double[DIM_NUM];
        double quality;
        double r_in = 0;
        double r_out = 0;

        tetrahedron_circumsphere_3d(tetra, ref r_out, ref pc);

        tetrahedron_insphere_3d(tetra, ref r_in, ref pc);

        quality = 3.0 * r_in / r_out;

        return quality;
    }

    public static double tetrahedron_quality2_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_QUALITY2_3D: "quality" of a tetrahedron in 3D.
        //
        //  Discussion:
        //
        //    The quality measure #2 of a tetrahedron is:
        //
        //      QUALITY2 = 2 * Math.Sqrt ( 6 ) * RIN / LMAX
        //
        //    where
        //
        //      RIN = radius of the inscribed sphere;
        //      LMAX = length of longest side of the tetrahedron.
        //
        //    An equilateral tetrahredron achieves the maximum possible quality of 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Qiang Du, Desheng Wang,
        //    The Optimal Centroidal Voronoi Tesselations and the Gersho's
        //      Conjecture in the Three-Dimensional Space,
        //    Computers and Mathematics with Applications,
        //    Volume 49, 2005, pages 1355-1373.
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the tetrahedron vertices.
        //
        //    Output, double TETRAHEDRON_QUALITY2_3D, the quality of the tetrahedron.
        //
    {
        int DIM_NUM = 3;

        double[] edge_length;
        double l_max;
        double[] pc = new double[DIM_NUM];
        double quality2;
        double r_in = 0;

        edge_length = tetrahedron_edge_length_3d(tetra);

        l_max = typeMethods.r8vec_max(6, edge_length);

        tetrahedron_insphere_3d(tetra, ref r_in, ref pc);

        quality2 = 2.0 * Math.Sqrt(6.0) * r_in / l_max;

        return quality2;
    }

    public static double tetrahedron_quality3_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_QUALITY3_3D computes the mean ratio of a tetrahedron.
        //
        //  Discussion:
        //
        //    This routine computes QUALITY3, the eigenvalue or mean ratio of
        //    a tetrahedron.
        //
        //      QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
        //
        //    This value may be used as a shape quality measure for the tetrahedron.
        //
        //    For an equilateral tetrahedron, the value of this quality measure
        //    will be 1.  For any other tetrahedron, the value will be between
        //    0 and 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Barry Joe.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Barry Joe,
        //    GEOMPACK - a software package for the generation of meshes
        //    using geometric algorithms,
        //    Advances in Engineering Software,
        //    Volume 13, pages 325-331, 1991.
        //
        //  Parameters:
        //
        //    Input, double TETRA(3,4), the vertices of the tetrahedron.
        //
        //    Output, double TETRAHEDRON_QUALITY3_3D, the mean ratio of the tetrahedron.
        //
    {
        int DIM_NUM = 3;

        double[] ab = new double[DIM_NUM];
        double[] ac = new double[DIM_NUM];
        double[] ad = new double[DIM_NUM];
        double[] bc = new double[DIM_NUM];
        double[] bd = new double[DIM_NUM];
        double[] cd = new double[DIM_NUM];
        double denom;
        int i;
        double lab;
        double lac;
        double lad;
        double lbc;
        double lbd;
        double lcd;
        double quality3;
        double volume;
        //
        //  Compute the vectors representing the sides of the tetrahedron.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            ab[i] = tetra[i + 1 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
            ac[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
            ad[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
            bc[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
            bd[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
            cd[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 2 * DIM_NUM];
        }

        //
        //  Compute the squares of the lengths of the sides.
        //
        lab = Math.Pow(ab[0], 2) + Math.Pow(ab[1], 2) + Math.Pow(ab[2], 2);
        lac = Math.Pow(ac[0], 2) + Math.Pow(ac[1], 2) + Math.Pow(ac[2], 2);
        lad = Math.Pow(ad[0], 2) + Math.Pow(ad[1], 2) + Math.Pow(ad[2], 2);
        lbc = Math.Pow(bc[0], 2) + Math.Pow(bc[1], 2) + Math.Pow(bc[2], 2);
        lbd = Math.Pow(bd[0], 2) + Math.Pow(bd[1], 2) + Math.Pow(bd[2], 2);
        lcd = Math.Pow(cd[0], 2) + Math.Pow(cd[1], 2) + Math.Pow(cd[2], 2);
        //
        //  Compute the volume.
        //
        volume = Math.Abs(
            ab[0] * (ac[1] * ad[2] - ac[2] * ad[1])
            + ab[1] * (ac[2] * ad[0] - ac[0] * ad[2])
            + ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])) / 6.0;

        denom = lab + lac + lad + lbc + lbd + lcd;

        quality3 = denom switch
        {
            0.0 => 0.0,
            _ => 12.0 * Math.Pow(3.0 * volume, 2.0 / 3.0) / denom
        };

        return quality3;
    }

    public static double tetrahedron_quality4_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_QUALITY4_3D computes the minimum solid angle of a tetrahedron.
        //
        //  Discussion:
        //
        //    This routine computes a quality measure for a tetrahedron, based
        //    on the sine of half the minimum of the four solid angles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Barry Joe.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Barry Joe,
        //    GEOMPACK - a software package for the generation of meshes
        //    using geometric algorithms,
        //    Advances in Engineering Software,
        //    Volume 13, pages 325-331, 1991.
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Output, double QUALITY4, the value of the quality measure.
        //
    {
        int DIM_NUM = 3;

        double[] ab = new double[DIM_NUM];
        double[] ac = new double[DIM_NUM];
        double[] ad = new double[DIM_NUM];
        double[] bc = new double[DIM_NUM];
        double[] bd = new double[DIM_NUM];
        double[] cd = new double[DIM_NUM];
        double denom;
        int i;
        double l1;
        double l2;
        double l3;
        double lab;
        double lac;
        double lad;
        double lbc;
        double lbd;
        double lcd;
        double quality4;
        double volume;
        //
        //  Compute the vectors that represent the sides.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            ab[i] = tetra[i + 1 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
            ac[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
            ad[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
            bc[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
            bd[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
            cd[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 2 * DIM_NUM];
        }

        //
        //  Compute the lengths of the sides.
        //
        lab = typeMethods.r8vec_norm(DIM_NUM, ab);
        lac = typeMethods.r8vec_norm(DIM_NUM, ac);
        lad = typeMethods.r8vec_norm(DIM_NUM, ad);
        lbc = typeMethods.r8vec_norm(DIM_NUM, bc);
        lbd = typeMethods.r8vec_norm(DIM_NUM, bd);
        lcd = typeMethods.r8vec_norm(DIM_NUM, cd);
        //
        //  Compute the volume.
        //
        volume = Math.Abs(
            ab[0] * (ac[1] * ad[2] - ac[2] * ad[1])
            + ab[1] * (ac[2] * ad[0] - ac[0] * ad[2])
            + ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])) / 6.0;

        quality4 = 1.0;

        l1 = lab + lac;
        l2 = lab + lad;
        l3 = lac + lad;

        denom = (l1 + lbc) * (l1 - lbc)
                           * (l2 + lbd) * (l2 - lbd)
                           * (l3 + lcd) * (l3 - lcd);

        quality4 = denom switch
        {
            <= 0.0 => 0.0,
            _ => Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom))
        };

        l1 = lab + lbc;
        l2 = lab + lbd;
        l3 = lbc + lbd;

        denom = (l1 + lac) * (l1 - lac)
                           * (l2 + lad) * (l2 - lad)
                           * (l3 + lcd) * (l3 - lcd);

        quality4 = denom switch
        {
            <= 0.0 => 0.0,
            _ => Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom))
        };

        l1 = lac + lbc;
        l2 = lac + lcd;
        l3 = lbc + lcd;

        denom = (l1 + lab) * (l1 - lab)
                           * (l2 + lad) * (l2 - lad)
                           * (l3 + lbd) * (l3 - lbd);

        quality4 = denom switch
        {
            <= 0.0 => 0.0,
            _ => Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom))
        };

        l1 = lad + lbd;
        l2 = lad + lcd;
        l3 = lbd + lcd;

        denom = (l1 + lab) * (l1 - lab)
                           * (l2 + lac) * (l2 - lac)
                           * (l3 + lbc) * (l3 - lbc);

        quality4 = denom switch
        {
            <= 0.0 => 0.0,
            _ => Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom))
        };

        quality4 = quality4 * 1.5 * Math.Sqrt(6.0);

        return quality4;
    }

    public static void tetrahedron_rhombic_shape_3d(int point_num, int face_num,
            int face_order_max, ref double[] point_coord, ref int[] face_order,
            ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_RHOMBIC_SHAPE_3D describes a rhombic tetrahedron in 3D.
        //
        //  Discussion:
        //
        //    Call TETRAHEDRON_RHOMBIC_SIZE_3D first, to get dimension information.
        //
        //    The tetrahedron is described using 10 nodes.  If we label the vertices
        //    P0, P1, P2 and P3, then the extra nodes lie halfway between vertices,
        //    and have the labels P01, P02, P03, P12, P13 and P23.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Anwei Liu, Barry Joe,
        //    Quality Local Refinement of Tetrahedral Meshes Based
        //    on 8-Subtetrahedron Subdivision,
        //    Mathematics of Computation,
        //    Volume 65, Number 215, July 1996, pages 1183-1200.
        //
        //  Parameters:
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int FACE_NUM, the number of faces.
        //
        //    Input, int FACE_ORDER_MAX, the maximum number of vertices per face.
        //
        //    Output, double POINT_COORD[3*POINT_NUM], the vertices.
        //
        //    Output, int FACE_ORDER[FACE_NUM], the number of vertices
        //    for each face.
        //
        //    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
        //    contains the index of the I-th point in the J-th face.  The
        //    points are listed in the counter clockwise direction defined
        //    by the outward normal at the face.
        //
    {
        double a;
        double b;
        double c;
        double d;
        int face;
        double z = 0.0;

        a = 1.0 / Math.Sqrt(3.0);
        b = Math.Sqrt(2.0) / Math.Sqrt(3.0);
        c = Math.Sqrt(3.0) / 6.0;
        d = 1.0 / Math.Sqrt(6.0);
        //
        //  Set the point coordinates.
        //
        point_coord[0 + 0 * 3] = -b;
        point_coord[1 + 0 * 3] = z;
        point_coord[2 + 0 * 3] = z;

        point_coord[0 + 1 * 3] = z;
        point_coord[1 + 1 * 3] = -a;
        point_coord[2 + 1 * 3] = z;

        point_coord[0 + 2 * 3] = z;
        point_coord[1 + 2 * 3] = a;
        point_coord[2 + 2 * 3] = z;

        point_coord[0 + 3 * 3] = z;
        point_coord[1 + 3 * 3] = z;
        point_coord[2 + 3 * 3] = b;

        point_coord[0 + 4 * 3] = -d;
        point_coord[1 + 4 * 3] = -c;
        point_coord[2 + 4 * 3] = z;

        point_coord[0 + 5 * 3] = -d;
        point_coord[1 + 5 * 3] = c;
        point_coord[2 + 5 * 3] = z;

        point_coord[0 + 6 * 3] = -d;
        point_coord[1 + 6 * 3] = z;
        point_coord[2 + 6 * 3] = d;

        point_coord[0 + 7 * 3] = z;
        point_coord[1 + 7 * 3] = z;
        point_coord[2 + 7 * 3] = z;

        point_coord[0 + 8 * 3] = z;
        point_coord[1 + 8 * 3] = -c;
        point_coord[2 + 8 * 3] = d;

        point_coord[0 + 9 * 3] = z;
        point_coord[1 + 9 * 3] = c;
        point_coord[2 + 9 * 3] = d;
        //
        //  Set the face orders.
        //
        for (face = 0; face < face_num; face++)
        {
            face_order[face] = 6;
        }

        //
        //  Set faces.
        //
        face_point[0 + 0 * face_order_max] = 1;
        face_point[1 + 0 * face_order_max] = 5;
        face_point[2 + 0 * face_order_max] = 2;
        face_point[3 + 0 * face_order_max] = 9;
        face_point[4 + 0 * face_order_max] = 4;
        face_point[5 + 0 * face_order_max] = 7;

        face_point[0 + 1 * face_order_max] = 2;
        face_point[1 + 1 * face_order_max] = 8;
        face_point[2 + 1 * face_order_max] = 3;
        face_point[3 + 1 * face_order_max] = 10;
        face_point[4 + 1 * face_order_max] = 4;
        face_point[5 + 1 * face_order_max] = 9;

        face_point[0 + 2 * face_order_max] = 3;
        face_point[1 + 2 * face_order_max] = 6;
        face_point[2 + 2 * face_order_max] = 1;
        face_point[3 + 2 * face_order_max] = 7;
        face_point[4 + 2 * face_order_max] = 4;
        face_point[5 + 2 * face_order_max] = 10;

        face_point[0 + 3 * face_order_max] = 1;
        face_point[1 + 3 * face_order_max] = 6;
        face_point[2 + 3 * face_order_max] = 3;
        face_point[3 + 3 * face_order_max] = 8;
        face_point[4 + 3 * face_order_max] = 2;
        face_point[5 + 3 * face_order_max] = 5;

    }

    public static void tetrahedron_rhombic_size_3d(ref int point_num, ref int edge_num,
            ref int face_num, ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_RHOMBIC_SIZE_3D gives "sizes" for a rhombic tetrahedron in 3D.
        //
        //  Discussion:
        //
        //    Call this routine first, in order to learn the required dimensions
        //    of arrays to be set up by TETRAHEDRON_RHOMBIC_SHAPE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int *POINT_NUM, the number of vertices.
        //
        //    Output, int *EDGE_NUM, the number of edges.
        //
        //    Output, int *FACE_NUM, the number of faces.
        //
        //    Output, int *FACE_ORDER_MAX, the maximum order of any face.
        //
    {
        point_num = 10;
        edge_num = 6;
        face_num = 4;
        face_order_max = 6;
    }

    public static void tetrahedron_sample_3d(double[] tetra, int n, ref int seed, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_SAMPLE_3D returns random points in a tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the tetrahedron vertices.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double P[3*N], random points in the tetrahedron.
        //
    {
        int DIM_NUM = 3;

        double alpha;
        double beta;
        double gamma;
        int i;
        int j;
        int k;
        double[] p12;
        double[] p13;
        double r;
        double[] t;

        p12 = new double[DIM_NUM];
        p13 = new double[DIM_NUM];
        t = new double[DIM_NUM * 3];

        for (k = 0; k < n; k++)
        {
            r = UniformRNG.r8_uniform_01(ref seed);
            //
            //  Interpret R as a percentage of the tetrahedron's volume.
            //
            //  Imagine a plane, parallel to face 1, so that the volume between
            //  vertex 1 and the plane is R percent of the full tetrahedron volume.
            //
            //  The plane will intersect sides 12, 13, and 14 at a fraction
            //  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
            //
            alpha = Math.Pow(r, 1.0 / 3.0);
            //
            //  Determine the coordinates of the points on sides 12, 13 and 14 intersected
            //  by the plane, which form a triangle TR.
            //
            for (i = 0; i < DIM_NUM; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    t[i + j * 3] = (1.0 - alpha) * tetra[i + 0 * 3]
                                   + alpha * tetra[i + (j + 1) * 3];
                }
            }

            //
            //  Now choose, uniformly at random, a point in this triangle.
            //
            r = UniformRNG.r8_uniform_01(ref seed);
            //
            //  Interpret R as a percentage of the triangle's area.
            //
            //  Imagine a line L, parallel to side 1, so that the area between
            //  vertex 1 and line L is R percent of the full triangle's area.
            //
            //  The line L will intersect sides 2 and 3 at a fraction
            //  ALPHA = Math.Sqrt ( R ) of the distance from vertex 1 to vertices 2 and 3.
            //
            beta = Math.Sqrt(r);
            //
            //  Determine the coordinates of the points on sides 2 and 3 intersected
            //  by line L.
            //
            for (i = 0; i < DIM_NUM; i++)
            {
                p12[i] = (1.0 - beta) * t[i + 0 * 3]
                         + beta * t[i + 1 * 3];

                p13[i] = (1.0 - beta) * t[i + 0 * 3]
                         + beta * t[i + 2 * 3];
            }

            //
            //  Now choose, uniformly at random, a point on the line L.
            //
            gamma = UniformRNG.r8_uniform_01(ref seed);

            for (i = 0; i < DIM_NUM; i++)
            {
                p[i + k * 3] = gamma * p12[i] + (1.0 - gamma) * p13[i];
            }
        }
    }

    public static void tetrahedron_shape_3d(int point_num, int face_num, int face_order_max,
            ref double[] point_coord, ref int[] face_order, ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_SHAPE_3D describes a tetrahedron in 3D.
        //
        //  Discussion:
        //
        //    The vertices lie on the unit sphere.
        //
        //    The dual of the tetrahedron is the tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int FACE_NUM, the number of faces.
        //
        //    Input, int FACE_ORDER_MAX, the maximum number of vertices per face.
        //
        //    Output, double POINT_COORD[3*POINT_NUM], the point coordinates.
        //
        //    Output, int FACE_ORDER[FACE_NUM], the number of vertices
        //    for each face.
        //
        //    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
        //    contains the index of the I-th point in the J-th face.  The
        //    points are listed in the counter clockwise direction defined
        //    by the outward normal at the face.
        //
    {
        int DIM_NUM = 3;

        int[] face_order_save =
        {
            3, 3, 3, 3
        };
        int[] face_point_save =
        {
            1, 3, 2,
            1, 2, 4,
            1, 4, 3,
            2, 3, 4
        };
        double[] point_coord_save =
        {
            0.942809, 0.000000, -0.333333,
            -0.471405, 0.816497, -0.333333,
            -0.471405, -0.816497, -0.333333,
            0.000000, 0.000000, 1.000000
        };

        typeMethods.i4vec_copy(face_num, face_order_save, ref face_order);
        typeMethods.i4vec_copy(face_order_max * face_num, face_point_save, ref face_point);
        typeMethods.r8vec_copy(DIM_NUM * point_num, point_coord_save, ref point_coord);

    }

    public static void tetrahedron_size_3d(ref int point_num, ref int edge_num, ref int face_num,
            ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_SIZE_3D gives "sizes" for a tetrahedron in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int *POINT_NUM, the number of points.
        //
        //    Output, int *EDGE_NUM, the number of edges.
        //
        //    Output, int *FACE_NUM, the number of faces.
        //
        //    Output, int *FACE_ORDER_MAX, the maximum order of any face.
        //
    {
        point_num = 4;
        edge_num = 12;
        face_num = 4;
        face_order_max = 3;
    }

    public static double[] tetrahedron_solid_angles_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_SOLID_ANGLES_3D computes solid angles of a tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Output, double TETRAHEDRON_SOLID_ANGLES_3D[4], the solid angles.
        //
    {
        double[] angle;
        double[] dihedral_angle;

        dihedral_angle = tetrahedron_dihedral_angles_3d(tetra);

        angle = new double[4];

        angle[0] = dihedral_angle[0] + dihedral_angle[1] + dihedral_angle[2] - Math.PI;
        angle[1] = dihedral_angle[0] + dihedral_angle[3] + dihedral_angle[4] - Math.PI;
        angle[2] = dihedral_angle[1] + dihedral_angle[3] + dihedral_angle[5] - Math.PI;
        angle[3] = dihedral_angle[2] + dihedral_angle[4] + dihedral_angle[5] - Math.PI;

        return angle;
    }

    public static int tetrahedron_unit_lattice_point_num_3d(int s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_UNIT_LATTICE_POINT_NUM_3D: count lattice points.
        //
        //  Discussion:
        //
        //    The tetrahedron is assumed to be the unit tetrahedron:
        //
        //    ( (0,0,0), (1,0,0), (0,1,0), (0,0,1) )
        //
        //    or a copy of this tetrahedron scaled by an integer S:
        //
        //    ( (0,0,0), (S,0,0), (0,S,0), (0,0,S) ).
        //
        //    The routine returns the number of integer lattice points that appear
        //    inside the tetrahedron, or on its faces, edges or vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Matthias Beck, Sinai Robins,
        //    Computing the Continuous Discretely,
        //    Springer, 2006,
        //    ISBN13: 978-0387291390,
        //    LC: QA640.7.B43.
        //
        //  Parameters:
        //
        //    Input, int S, the scale factor.
        //
        //    Output, int TETRAHEDRON_UNIT_LATTICE_POINT_NUM_3D, the number of lattice points.
        //
    {
        int n;

        n = (s + 3) * (s + 2) * (s + 1) / 6;

        return n;
    }

    public static double tetrahedron_volume_3d(double[] tetra)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TETRA[3*4], the vertices of the tetrahedron.
        //
        //    Output, double TETRAHEDRON_VOLUME_3D, the volume of the tetrahedron.
        //
    {
        double[] a = new double[4 * 4];
        int i;
        int j;
        double volume;

        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 4; j++)
            {
                a[i + j * 4] = tetra[i + j * 3];
            }
        }

        i = 3;
        for (j = 0; j < 4; j++)
        {
            a[i + j * 4] = 1.0;
        }

        volume = Math.Abs(typeMethods.r8mat_det_4d(a)) / 6.0;

        return volume;
    }

}