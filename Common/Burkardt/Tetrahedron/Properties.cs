using System;
using Burkardt.Types;

namespace Burkardt.TetrahedronNS
{
    public static class Properties
    {
        public static double[] tetrahedron_centroid(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_CENTROID computes the centroid of a tetrahedron.
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
            //    Output, double TETRAHEDRON_CENTROID[3], the coordinates of the centroid.
            //
        {
            double[] centroid;
            int i;

            centroid = new double[3];

            centroid[0] = (tetra[0 + 0 * 3] + tetra[0 + 1 * 3]
                                            + tetra[0 + 2 * 3] + tetra[0 + 3 * 3]);
            centroid[1] = (tetra[1 + 0 * 3] + tetra[1 + 1 * 3]
                                            + tetra[1 + 2 * 3] + tetra[1 + 3 * 3]);
            centroid[2] = (tetra[2 + 0 * 3] + tetra[2 + 1 * 3]
                                            + tetra[2 + 2 * 3] + tetra[2 + 3 * 3]);

            for (i = 0; i < 3; i++)
            {
                centroid[i] = centroid[i] / 4.0;
            }

            return centroid;
        }

        public static void tetrahedron_circumsphere(double[] tetra, ref double r, ref double[] pc)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_CIRCUMSPHERE computes the circumsphere of a tetrahedron.
            //
            //  Discussion:
            //
            //    The circumsphere, or circumscribed sphere, of a tetrahedron is the 
            //    sphere that passes through the four vertices.  The circumsphere is not
            //    necessarily the smallest sphere that contains the tetrahedron.
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
            //    Output, double &R, PC[3], the coordinates of the center of the
            //    circumscribed sphere, and its radius.  If the linear system is
            //    singular, then R = -1, PC[] = 0.
            //
        {
            double[] a = new double[3 * 4];
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
            info = typeMethods.r8mat_solve(3, 1, ref a);
            //
            //  If the system was singular, return a consolation prize.
            //
            if (info != 0)
            {
                r = -1.0;
                typeMethods.r8vec_zero(3, ref pc);
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

            return;
        }
        
        public static void tetrahedron_circumsphere_3d ( double[] tetra, ref double r, ref double[] pc )

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
        //    Input, double TETRA[3*4], the coordinates of the vertices.
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

        public static double[] tetrahedron_dihedral_angles(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_DIHEDRAL_ANGLES computes dihedral angles of a tetrahedron.
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
            //    Output, double TETRAHEDRON_DIHEDRAL_ANGLES[6], the dihedral angles 
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
            double[] cd = new double[3];
            int i;
            const double r8_pi = 3.141592653589793;

            tetrahedron_edges(tetra, ref ab, ref ac, ref ad, ref bc, ref bd, ref cd);

            abc_normal = typeMethods.r8vec_cross_3d(ac, ab);
            abd_normal = typeMethods.r8vec_cross_3d(ab, ad);
            acd_normal = typeMethods.r8vec_cross_3d(ad, ac);
            bcd_normal = typeMethods.r8vec_cross_3d(bc, bd);

            angle = new double[6];

            angle[0] = typeMethods.r8vec_angle_3d(abc_normal, abd_normal);
            angle[1] = typeMethods.r8vec_angle_3d(abc_normal, acd_normal);
            angle[2] = typeMethods.r8vec_angle_3d(abd_normal, acd_normal);
            angle[3] = typeMethods.r8vec_angle_3d(abc_normal, bcd_normal);
            angle[4] = typeMethods.r8vec_angle_3d(abd_normal, bcd_normal);
            angle[5] = typeMethods.r8vec_angle_3d(acd_normal, bcd_normal);

            for (i = 0; i < 6; i++)
            {
                angle[i] = r8_pi - angle[i];
            }

            return angle;
        }

        public static double[] tetrahedron_edge_length(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_EDGE_LENGTH returns edge lengths of a tetrahedron.
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
            double[] edge_length;
            int i;
            int j1;
            int j2;
            int k;
            double[] v = new double[3];

            edge_length = new double[6];

            k = 0;
            for (j1 = 0; j1 < 3; j1++)
            {
                for (j2 = j1 + 1; j2 < 4; j2++)
                {
                    for (i = 0; i < 3; i++)
                    {
                        v[i] = tetra[i + j2 * 3] - tetra[i + j1 * 3];
                    }

                    edge_length[k] = typeMethods.r8vec_length(3, v);
                    k = k + 1;
                }
            }

            return edge_length;
        }

        public static double[] tetrahedron_edge_length_3d ( double[] tetra )

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
        //    Input, double TETRA[3*4], the coordinates of the vertices.
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

                    edge_length[k] = typeMethods.r8vec_length(DIM_NUM, v);
                    k = k + 1;
                }
            }

            return edge_length;
        }

        public static void tetrahedron_insphere_3d ( double[] tetra, ref double r, ref double[] pc )

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
        //    Input, double TETRA[3*4], the coordinates of the vertices.
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

            n123 = typeMethods.r8vec_cross_3d(v21, v31);
            n124 = typeMethods.r8vec_cross_3d(v41, v21);
            n134 = typeMethods.r8vec_cross_3d(v31, v41);
            n234 = typeMethods.r8vec_cross_3d(v42, v32);

            l123 = typeMethods.r8vec_length(DIM_NUM, n123);
            l124 = typeMethods.r8vec_length(DIM_NUM, n124);
            l134 = typeMethods.r8vec_length(DIM_NUM, n134);
            l234 = typeMethods.r8vec_length(DIM_NUM, n234);

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

        public static void tetrahedron_edges(double[] tetra, ref double[] ab, ref double[] ac,
                ref double[] ad, ref double[] bc, ref double[] bd, ref double[] cd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_EDGES returns the edges of a tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double TETRA[3*4], the tetrahedron vertices.
            //
            //    Output, double AB[3], AC[3], AD[3], BC[3], BD[3], CD[3], the edges.
            //
        {
            int i;
            //
            //  Compute the vectors that represent the sides.
            //
            for (i = 0; i < 3; i++)
            {
                ab[i] = tetra[i + 1 * 3] - tetra[i + 0 * 3];
                ac[i] = tetra[i + 2 * 3] - tetra[i + 0 * 3];
                ad[i] = tetra[i + 3 * 3] - tetra[i + 0 * 3];
                bc[i] = tetra[i + 2 * 3] - tetra[i + 1 * 3];
                bd[i] = tetra[i + 3 * 3] - tetra[i + 1 * 3];
                cd[i] = tetra[i + 3 * 3] - tetra[i + 2 * 3];
            }

        }

        public static void tetrahedron_face_angles(double[] tetra, ref double[] angles)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_FACE_ANGLES returns the 12 face angles of a tetrahedron.
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

            typeMethods.triangle_angles_3d(tri, ref angles);
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

            typeMethods.triangle_angles_3d(tri, ref angles, angleIndex: +3);
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

            typeMethods.triangle_angles_3d(tri, ref angles, angleIndex: +6);
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

            typeMethods.triangle_angles_3d(tri, ref angles, angleIndex: +9);

        }

        public static void tetrahedron_face_areas(double[] tetra, ref double[] areas)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_FACE_AREAS returns the 4 face areas of a tetrahedron.
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

            areas[0] = typeMethods.triangle_area_3d(tri);
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

            areas[1] = typeMethods.triangle_area_3d(tri);
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

            areas[2] = typeMethods.triangle_area_3d(tri);
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

            areas[3] = typeMethods.triangle_area_3d(tri);

        }

        public static void tetrahedron_insphere(double[] tetra, ref double r, ref double[] pc)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_INSPHERE finds the insphere of a tetrahedron.
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
            //    Output, double &R, PC[3], the radius and the center
            //    of the sphere.
            //
        {
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
            double[] v21 = new double[3];
            double[] v31 = new double[3];
            double[] v41 = new double[3];
            double[] v32 = new double[3];
            double[] v42 = new double[3];
            double[] v43 = new double[3];

            tetrahedron_edges(tetra, ref v21, ref v31, ref v41, ref v32, ref v42, ref v43);

            n123 = typeMethods.r8vec_cross_3d(v21, v31);
            n124 = typeMethods.r8vec_cross_3d(v41, v21);
            n134 = typeMethods.r8vec_cross_3d(v31, v41);
            n234 = typeMethods.r8vec_cross_3d(v42, v32);

            l123 = typeMethods.r8vec_length(3, n123);
            l124 = typeMethods.r8vec_length(3, n124);
            l134 = typeMethods.r8vec_length(3, n134);
            l234 = typeMethods.r8vec_length(3, n234);

            for (i = 0; i < 3; i++)
            {
                pc[i] = (l234 * tetra[i + 0 * 3]
                         + l134 * tetra[i + 1 * 3]
                         + l124 * tetra[i + 2 * 3]
                         + l123 * tetra[i + 3 * 3])
                        / (l234 + l134 + l124 + l123);
            }

            for (j = 0; j < 4; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    b[i + j * 4] = tetra[i + j * 3];
                }

                b[3 + j * 4] = 1.0;
            }

            gamma = Math.Abs(typeMethods.r8mat_det_4d(b));

            r = gamma / (l234 + l134 + l124 + l123);

        }

        public static double tetrahedron_quality1(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_QUALITY1: "quality" of a tetrahedron.
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
            //    Output, double TETRAHEDRON_QUALITY1, the quality of the tetrahedron.
            //
        {
            double[] pc = new double[3];
            double quality;
            double r_in = 0;
            double r_out = 0;

            tetrahedron_circumsphere(tetra, ref r_out, ref pc);

            tetrahedron_insphere(tetra, ref r_in, ref pc);

            quality = 3.0 * r_in / r_out;

            return quality;
        }

        public static double tetrahedron_quality2(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_QUALITY2: "quality" of a tetrahedron.
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
            //    Conjecture in the Three-Dimensional Space,
            //    Computers and Mathematics with Applications,
            //    Volume 49, 2005, pages 1355-1373.
            //
            //  Parameters:
            //
            //    Input, double TETRA[3*4], the tetrahedron vertices.
            //
            //    Output, double TETRAHEDRON_QUALITY2, the quality of the tetrahedron.
            //
        {
            double[] edge_length;
            double l_max;
            double[] pc = new double[3];
            double quality2;
            double r_in = 0;

            edge_length = tetrahedron_edge_length(tetra);

            l_max = typeMethods.r8vec_max(6, edge_length);

            tetrahedron_insphere(tetra, ref r_in, ref pc);

            quality2 = 2.0 * Math.Sqrt(6.0) * r_in / l_max;

            return quality2;
        }

        public static double tetrahedron_quality3(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_QUALITY3 computes the mean ratio of a tetrahedron.
            //
            //  Discussion:
            //
            //    This routine computes QUALITY3, the eigenvalue or mean ratio of
            //    a tetrahedron.
            //
            //      QUALITY3 = 12 * ( 3 * volume )^(2/3) / (sum of square of edge lengths).
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
            //    Input, double TETRA[3*4], the vertices of the tetrahedron.
            //
            //    Output, double TETRAHEDRON_QUALITY3, the mean ratio of the tetrahedron.
            //
        {
            double[] ab = new double[3];
            double[] ac = new double[3];
            double[] ad = new double[3];
            double[] bc = new double[3];
            double[] bd = new double[3];
            double[] cd = new double[3];
            double denom;
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
            tetrahedron_edges(tetra, ref ab, ref ac, ref ad, ref bc, ref bd, ref cd);
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

            if (denom == 0.0)
            {
                quality3 = 0.0;
            }
            else
            {
                quality3 = 12.0 * Math.Pow(3.0 * volume, 2.0 / 3.0) / denom;
            }

            return quality3;
        }

        public static double tetrahedron_quality4(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_QUALITY4 computes the minimum solid angle of a tetrahedron.
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
            double[] ab = new double[3];
            double[] ac = new double[3];
            double[] ad = new double[3];
            double[] bc = new double[3];
            double[] bd = new double[3];
            double[] cd = new double[3];
            double denom;
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
            //  Compute the edges.
            //
            tetrahedron_edges(tetra, ref ab, ref ac, ref ad, ref bc, ref bd, ref cd);
            //
            //  Compute the lengths of the sides.
            //
            lab = typeMethods.r8vec_length(3, ab);
            lac = typeMethods.r8vec_length(3, ac);
            lad = typeMethods.r8vec_length(3, ad);
            lbc = typeMethods.r8vec_length(3, bc);
            lbd = typeMethods.r8vec_length(3, bd);
            lcd = typeMethods.r8vec_length(3, cd);
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

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            l1 = lab + lbc;
            l2 = lab + lbd;
            l3 = lbc + lbd;

            denom = (l1 + lac) * (l1 - lac)
                               * (l2 + lad) * (l2 - lad)
                               * (l3 + lcd) * (l3 - lcd);

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            l1 = lac + lbc;
            l2 = lac + lcd;
            l3 = lbc + lcd;

            denom = (l1 + lab) * (l1 - lab)
                               * (l2 + lad) * (l2 - lad)
                               * (l3 + lbd) * (l3 - lbd);

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            l1 = lad + lbd;
            l2 = lad + lcd;
            l3 = lbd + lcd;

            denom = (l1 + lab) * (l1 - lab)
                               * (l2 + lac) * (l2 - lac)
                               * (l3 + lbc) * (l3 - lbc);

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            quality4 = quality4 * 1.5 * Math.Sqrt(6.0);

            return quality4;
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
            //    Input, double TETRA[3*4], the coordinates of the vertices.
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
            //      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
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
            //    Conjecture in the Three-Dimensional Space,
            //    Computers and Mathematics with Applications,
            //    Volume 49, 2005, pages 1355-1373.
            //
            //  Parameters:
            //
            //    Input, double TETRA[3*4], the coordinates of the vertices.
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
            //    Input, double TETRA(3,4), the coordinates of the vertices.
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

            if (denom == 0.0)
            {
                quality3 = 0.0;
            }
            else
            {
                quality3 = 12.0 * Math.Pow(3.0 * volume, 2.0 / 3.0) / denom;
            }

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
            //    Input, double TETRA[3*4], the coordinates of the vertices.
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
            lab = typeMethods.r8vec_length(DIM_NUM, ab);
            lac = typeMethods.r8vec_length(DIM_NUM, ac);
            lad = typeMethods.r8vec_length(DIM_NUM, ad);
            lbc = typeMethods.r8vec_length(DIM_NUM, bc);
            lbd = typeMethods.r8vec_length(DIM_NUM, bd);
            lcd = typeMethods.r8vec_length(DIM_NUM, cd);
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

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            l1 = lab + lbc;
            l2 = lab + lbd;
            l3 = lbc + lbd;

            denom = (l1 + lac) * (l1 - lac)
                               * (l2 + lad) * (l2 - lad)
                               * (l3 + lcd) * (l3 - lcd);

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            l1 = lac + lbc;
            l2 = lac + lcd;
            l3 = lbc + lcd;

            denom = (l1 + lab) * (l1 - lab)
                               * (l2 + lad) * (l2 - lad)
                               * (l3 + lbd) * (l3 - lbd);

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            l1 = lad + lbd;
            l2 = lad + lcd;
            l3 = lbd + lcd;

            denom = (l1 + lab) * (l1 - lab)
                               * (l2 + lac) * (l2 - lac)
                               * (l3 + lbc) * (l3 - lbc);

            if (denom <= 0.0)
            {
                quality4 = 0.0;
            }
            else
            {
                quality4 = Math.Min(quality4, 12.0 * volume / Math.Sqrt(denom));
            }

            quality4 = quality4 * 1.5 * Math.Sqrt(6.0);

            return quality4;
        }

        public static double[] tetrahedron_solid_angles(double[] tetra)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_SOLID_ANGLES computes solid angles of a tetrahedron.
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
            //    Output, double TETRAHEDRON_SOLID_ANGLES[4], the solid angles.
            //
        {
            double[] angle;
            double[] dihedral_angles;
            const double r8_pi = 3.141592653589793;

            dihedral_angles = tetrahedron_dihedral_angles(tetra);

            angle = new double[4];

            angle[0] = dihedral_angles[0]
                + dihedral_angles[1]
                + dihedral_angles[2] - r8_pi;

            angle[1] = dihedral_angles[0]
                + dihedral_angles[3]
                + dihedral_angles[4] - r8_pi;

            angle[2] = dihedral_angles[1]
                + dihedral_angles[3]
                + dihedral_angles[5] - r8_pi;

            angle[3] = dihedral_angles[2]
                + dihedral_angles[4]
                + dihedral_angles[5] - r8_pi;

            return angle;
        }
    }
}