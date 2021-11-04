using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SphereNS
{
    public static class Triangle
    {
        public static double sphere01_triangle_angles_to_area(double a, double b, double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_ANGLES_TO_AREA: area of a spherical triangle on the unit sphere.
            //
            //  Discussion:
            //
            //    A unit sphere centered at 0 in 3D satisfies the equation:
            //
            //      X^2 + Y^2 + Z^2 = 1
            //
            //    A spherical triangle is specified by three points on the surface
            //    of the sphere.
            //
            //    The area formula is known as Girard's formula.
            //
            //    The area of a spherical triangle on the unit sphere is:
            //
            //      AREA = ( A + B + C - PI )
            //
            //    where A, B and C are the (surface) angles of the triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, C, the angles of the triangle.
            //
            //    Output, double STRI_ANGLES_TO_AREA, the area of the spherical triangle.
            //
        {
            double area;
            

            area = a + b + c - Math.PI;

            return area;
        }

        public static double[] sphere01_triangle_project(double[] a_xyz, double[] b_xyz, double[] c_xyz,
                int f1, int f2, int f3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_PROJECT projects from a plane triangle to a spherical triangle.
            //
            //  Discussion:
            //
            //    We assume that points A, B and C lie on the unit sphere, and they
            //    thus define a spherical triangle.
            //
            //    They also, of course, define a planar triangle.
            //
            //    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
            //    planar triangle.
            //
            //    This function determines the coordinates of the point in the planar
            //    triangle identified by the barycentric coordinates, and returns the
            //    coordinates of the projection of that point onto the unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the coordinates
            //    of the points A, B, and C.
            //
            //    Input, int F1, F2, F3, the barycentric coordinates
            //    of a point in the triangle ABC.  Normally, these coordinates would
            //    be real numbers, and would sum to 1.  For convenience, we allow these
            //    to be integers which must be divided by F1+F2+F3.
            //
            //    Output, double NODE_XYZ[3], the coordinates of the 
            //    point on the unit sphere which is the projection of the point on the plane
            //    whose barycentric coordinates with respect to A, B, and C is
            //    (F1,F2,F3)/(F1+F2+F3).
            //
        {
            int i;
            double[] node_xyz;
            double norm;

            node_xyz = new double[3];

            for (i = 0; i < 3; i++)
            {
                node_xyz[i] =
                    ((double)(f1) * a_xyz[i]
                     + (double)(f2) * b_xyz[i]
                     + (double)(f3) * c_xyz[i])
                    / (double)(f1 + f2 + f3);
            }

            norm = typeMethods.r8vec_norm(3, node_xyz);

            for (i = 0; i < 3; i++)
            {
                node_xyz[i] = node_xyz[i] / norm;
            }

            return node_xyz;
        }

        public static double[] sphere01_triangle_project2(double[] a_xyz, double[] b_xyz, double[] c_xyz,
                int f1, int f2, int f3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_PROJECT2 projects from a plane triangle to a spherical triangle.
            //
            //  Discussion:
            //
            //    We assume that points A, B and C lie on the unit sphere, and they
            //    thus define a spherical triangle.
            //
            //    They also, of course, define a planar triangle.
            //
            //    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
            //    planar triangle.
            //
            //    This function determines the coordinates of the point in the planar
            //    triangle identified by the barycentric coordinates, and returns the
            //    coordinates of the projection of that point onto the unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
            //    of the points A, B, and C.
            //
            //    Input, int F1, F2, F3, the barycentric coordinates
            //    of a point in the triangle ABC.  Normally, these coordinates would
            //    be real numbers, and would sum to 1.  For convenience, we allow these
            //    to be integers which must be divided by F1+F2+F3.
            //
            //    Output, double SPHERE01_TRIANGLE_PROJECT2[3], the coordinates of the 
            //    point on the unit sphere which is the projection of the point on the 
            //    plane whose barycentric coordinates with respect to A, B, and C is
            //    (F1,F2,F3)/(F1+F2+F3).
            //
        {
            double[] ab = new double[3];
            double[] ac = new double[3];
            double[] acn = new double[3];
            double[] acp = new double[3];
            double angle;
            double[] bn = new double[3];
            double[] bp = new double[3];
            double[] cn = new double[3];
            double[] cp = new double[3];
            int i;
            double[] node_xyz;
            double norm;
            double theta_ab;
            double theta_ac;
            double theta_bc;

            node_xyz = new double[3];
            //
            //  This check avoids 0/0 calculations later.
            //
            if (f2 == 0 && f3 == 0)
            {
                for (i = 0; i < 3; i++)
                {
                    node_xyz[i] = a_xyz[i];
                }

                return node_xyz;
            }
            else if (f1 == 0 && f3 == 0)
            {
                for (i = 0; i < 3; i++)
                {
                    node_xyz[i] = b_xyz[i];
                }

                return node_xyz;
            }
            else if (f1 == 0 && f2 == 0)
            {
                for (i = 0; i < 3; i++)
                {
                    node_xyz[i] = c_xyz[i];
                }

                return node_xyz;
            }

            //
            //  Determine the angular distances (A,B) and (A,C).
            //
            theta_ab = Distance.sphere01_distance_xyz(a_xyz, b_xyz);

            theta_ac = Distance.sphere01_distance_xyz(a_xyz, c_xyz);
            //
            //  Polarize B = BP + BN
            //  Normalize BN, 
            //  Same for C.
            //
            typeMethods.r8vec_polarize(3, b_xyz, a_xyz, ref bn, ref bp);
            norm = typeMethods.r8vec_norm(3, bn);
            for (i = 0; i < 3; i++)
            {
                bn[i] = bn[i] / norm;
            }

            typeMethods.r8vec_polarize(3, c_xyz, a_xyz, ref cn, ref cp);
            norm = typeMethods.r8vec_norm(3, cn);
            for (i = 0; i < 3; i++)
            {
                cn[i] = cn[i] / norm;
            }

            //
            //  Determine AB and AC that use Math.Cos ( ( F2 + F3 ) / ( F1 + F2 + F3 ) ) of A
            //  and Math.Cos ( F1 / ( F1 + F2 + F3 ) ) of B or C.
            //
            angle = ((double)(f2 + f3) * theta_ab) / (double)(f1 + f2 + f3);
            for (i = 0; i < 3; i++)
            {
                ab[i] = Math.Cos(angle) * a_xyz[i] + Math.Sin(angle) * bn[i];
            }

            angle = ((double)(f2 + f3) * theta_ac) / (double)(f1 + f2 + f3);
            for (i = 0; i < 3; i++)
            {
                ac[i] = Math.Cos(angle) * a_xyz[i] + Math.Sin(angle) * cn[i];
            }

            //
            //  Determine the angular distance between AB and AC.
            //
            theta_bc = Distance.sphere01_distance_xyz(ab, ac);
            //
            //  Polarize AC = ACP + ACN, normalize ACN.
            //
            typeMethods.r8vec_polarize(3, ac, ab, ref acn, ref acp);
            norm = typeMethods.r8vec_norm(3, acn);
            for (i = 0; i < 3; i++)
            {
                acn[i] = acn[i] / norm;
            }

            //
            //  The interval between AB and AC is marked by F2+F3+1 vertices 0 through F2+F3.
            //
            angle = ((double)(f3) * theta_bc) / (double)(f2 + f3);

            for (i = 0; i < 3; i++)
            {
                node_xyz[i] = Math.Cos(angle) * ab[i] + Math.Sin(angle) * acn[i];
            }

            return node_xyz;
        }

        public static double[] sphere01_triangle_sample(int n, double[] v1, double[] v2, double[] v3,
                ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_SAMPLE: sample points from triangle on unit sphere.
            //
            //  Discussion:
            //
            //    The sphere has center 0 and radius 1.
            //
            //    A spherical triangle on the surface of the unit sphere contains those 
            //    points with radius R = 1, bounded by the vertices V1, V2, V3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    James Arvo,
            //    Stratified sampling of spherical triangles,
            //    Computer Graphics Proceedings, Annual Conference Series, 
            //    ACM SIGGRAPH '95, pages 437-438, 1995.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double V1[3], V2[3], V3[3], the XYZ coordinates of
            //    the vertices of the spherical triangle.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, double SPHERE01_TRIANGLE_SAMPLE[3*N], the XYZ coordinates of the 
            //    sample points.
            //
        {
            double a = 0;
            double alpha = 0;
            double area = 0;
            double area_hat = 0;
            double b = 0;
            double beta = 0;
            double c = 0;
            double gamma = 0;
            int i;
            int j;
            double norm;
            double q;
            double s;
            double t;
            double u;
            double v;
            double[] v31 = new double[3];
            double[] v4 = new double[3];
            double[] v42 = new double[3];
            double w;
            double[] x;
            double xsi1;
            double xsi2;
            double z;

            sphere01_triangle_vertices_to_sides(v1, v2, v3, ref a, ref b, ref c);

            sphere01_triangle_sides_to_angles(a, b, c, ref alpha, ref beta, ref gamma);

            area = sphere01_triangle_angles_to_area(alpha, beta, gamma);

            x = new double[3 * n];

            for (j = 0; j < n; j++)
            {
                //
                //  Select the new area.
                //
                xsi1 = UniformRNG.r8_uniform_01(ref seed);
                area_hat = xsi1 * area;
                //
                //  Compute the sine and cosine of the angle phi.
                //
                s = Math.Sin(area_hat - alpha);
                t = Math.Cos(area_hat - alpha);
                //
                //  Compute the pair that determines beta_hat.
                //
                u = t - Math.Cos(alpha);
                v = s + Math.Sin(alpha) * Math.Cos(c);
                //
                //  Q is the cosine of the new edge length b_hat.
                //
                q = ((v * t - u * s) * Math.Cos(alpha) - v)
                    / ((v * s + u * t) * Math.Sin(alpha));
                //
                //  We very occasionally get a Q value out of bounds.
                //
                q = Math.Max(q, -1.0);
                q = Math.Min(q, +1.0);
                //
                //  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
                //
                w = typeMethods.r8vec_dot_product(3, v3, v1);
                for (i = 0; i < 3; i++)
                {
                    v31[i] = v3[i] - w * v1[i];
                }

                norm = typeMethods.r8vec_norm(3, v31);
                for (i = 0; i < 3; i++)
                {
                    v31[i] = v31[i] / norm;
                }

                //
                //  V4 is the third vertex of the subtriangle V1, V2, V4.
                //
                for (i = 0; i < 3; i++)
                {
                    v4[i] = q * v1[i] + Math.Sqrt(1.0 - q * q) * v31[i];
                }

                //
                //  Select Math.Cos theta, which will sample along the edge from V2 to V4.
                //
                xsi2 = UniformRNG.r8_uniform_01(ref seed);
                z = 1.0 - xsi2 * (1.0 - typeMethods.r8vec_dot_product(3, v4, v2));
                //
                //  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
                //
                w = typeMethods.r8vec_dot_product(3, v4, v2);
                for (i = 0; i < 3; i++)
                {
                    v42[i] = v4[i] - w * v2[i];
                }

                norm = typeMethods.r8vec_norm(3, v42);
                for (i = 0; i < 3; i++)
                {
                    v42[i] = v42[i] / norm;
                }

                //
                //  Construct the point.
                //
                for (i = 0; i < 3; i++)
                {
                    x[i + j * 3] = z * v2[i] + Math.Sqrt(1.0 - z * z) * v42[i];
                }
            }

            return x;
        }

        public static void sphere01_triangle_sides_to_angles(double as_, double bs, double cs,
                ref double a, ref double b, ref double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_SIDES_TO_ANGLES: angles of spherical triangle on unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double AS, BS, CS, the (geodesic) length of the sides of the
            //    triangle.
            //
            //    Output, double *A, *B, *C, the spherical angles of the triangle.
            //    Angle A is opposite the side of length AS, and so on.
            //
        {
            double asu;
            double bsu;
            double csu;
            double ssu;
            double tan_a2;
            double tan_b2;
            double tan_c2;

            asu = as_;
            bsu = bs;
            csu = cs;
            ssu = (asu + bsu + csu) / 2.0;

            tan_a2 = Math.Sqrt((Math.Sin(ssu - bsu) * Math.Sin(ssu - csu)) /
                               (Math.Sin(ssu) * Math.Sin(ssu - asu)));

            a = 2.0 * Math.Atan(tan_a2);

            tan_b2 = Math.Sqrt((Math.Sin(ssu - asu) * Math.Sin(ssu - csu)) /
                               (Math.Sin(ssu) * Math.Sin(ssu - bsu)));

            b = 2.0 * Math.Atan(tan_b2);

            tan_c2 = Math.Sqrt((Math.Sin(ssu - asu) * Math.Sin(ssu - bsu)) /
                               (Math.Sin(ssu) * Math.Sin(ssu - csu)));

            c = 2.0 * Math.Atan(tan_c2);

            return;
        }

        public static void sphere01_triangle_vertices_to_angles(double[] v1, double[] v2,
                double[] v3, ref double a, ref double b, ref double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_VERTICES_TO_ANGLES: angles of spherical triangle on unit sphere.
            //
            //  Discussion:
            //
            //    A unit sphere centered at 0 in 3D satisfies the equation:
            //
            //      X*X + Y*Y + Z*Z = 1
            //
            //    A spherical triangle is specified by three points on the surface
            //    of the sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
            //
            //    Output, double *A, *B, *C, the angles of the spherical triangle.
        {
            double as_ = 0;
            double bs = 0;
            double cs = 0;
            //
            //  Compute the lengths of the sides of the spherical triangle.
            //
            sphere01_triangle_vertices_to_sides(v1, v2, v3, ref as_, ref bs, ref cs);
            //
            //  Get the spherical angles.
            //
            sphere01_triangle_sides_to_angles(as_, bs, cs, ref a, ref b, ref c);

            return;
        }

        public static double sphere01_triangle_vertices_to_area(double[] v1, double[] v2, double[] v3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_VERTICES_TO_AREA: area of a spherical triangle on unit sphere.
            //
            //  Discussion:
            //
            //    A unit sphere centered at 0 in 3D satisfies the equation:
            //
            //      X*X + Y*Y + Z*Z = 1
            //
            //    A spherical triangle is specified by three points on the surface
            //    of the sphere.
            //
            //    The area formula is known as Girard's formula.
            //
            //    The area of a spherical triangle on the unit sphere is:
            //
            //      AREA = A + B + C - PI
            //
            //    where A, B and C are the (surface) angles of the triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
            //
            //    Output, double STRI_VERTICES_TO_AREA, the area of the 
            //    spherical triangle.
        {
            double area = 0;
            double a = 0;
            double as_ = 0;
            double b = 0;
            double bs = 0;
            double c = 0;
            double cs = 0;
            //
            //  Compute the lengths of the sides of the spherical triangle.
            //
            sphere01_triangle_vertices_to_sides(v1, v2, v3, ref as_, ref bs, ref cs);
            //
            //  Get the spherical angles.
            //
            sphere01_triangle_sides_to_angles(as_, bs, cs, ref a, ref b, ref c);
            //
            //  Get the area
            //
            area = sphere01_triangle_angles_to_area(a, b, c);

            return area;
        }

        public static double[] sphere01_triangle_vertices_to_centroid(double[] v1, double[] v2, double[] v3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_VERTICES_TO_CENTROID: centroid of spherical triangle on unit sphere.
            //
            //  Discussion:
            //
            //    A sphere centered at 0 in 3D satisfies the equation:
            //
            //      X*X + Y*Y + Z*Z = 1
            //
            //    A spherical triangle is specified by three points on the sphere.
            //
            //    The (true) centroid of a spherical triangle is the point
            //
            //      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
            //
            //    Note that the true centroid does NOT, in general, lie on the sphere.  
            //
            //    The "flat" centroid VF is the centroid of the planar triangle defined by
            //    the vertices of the spherical triangle.
            //
            //    The "spherical" centroid VS of a spherical triangle is computed by
            //    the intersection of the geodesic bisectors of the triangle angles.
            //    The spherical centroid lies on the sphere.
            //
            //    VF, VT and VS lie on a line through the center of the sphere.  We can
            //    easily calculate VF by averaging the vertices, and from this determine
            //    VS by normalizing.
            //
            //    (Of course, we still will not have actually computed VT, which lies
            //    somewhere between VF and VS!)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
            //
            //    Output, double SPHERE01_TRIANGLE_VERTICES_TO_CENTROID[3], the coordinates of the 
            //    "spherical centroid" of the spherical triangle.
            //
        {
            int DIM_NUM = 3;

            int i;
            double norm;
            double[] vs;

            vs = new double[3];

            for (i = 0; i < DIM_NUM; i++)
            {
                vs[i] = (v1[i] + v2[i] + v3[i]) / 3.0;
            }

            norm = typeMethods.r8vec_norm(DIM_NUM, vs);

            for (i = 0; i < DIM_NUM; i++)
            {
                vs[i] = vs[i] / norm;
            }

            return vs;
        }

        public static void sphere01_triangle_vertices_to_midpoints(double[] v1, double[] v2, double[] v3,
                ref double[] m1, ref double[] m2, ref double[] m3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS gets the midsides of a spherical triangle.
            //
            //  Discussion:
            //
            //    The points are assumed to lie on the unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
            //
            //    Output, double M1[3], M2[3], M3[3], the coordinates of 
            //    the midpoints of the sides of the spherical triangle.
            //
        {
            int i;
            double norm;

            for (i = 0; i < 3; i++)
            {
                m1[i] = (v1[i] + v2[i]) / 2.0;
            }

            norm = typeMethods.r8vec_norm(3, m1);
            for (i = 0; i < 3; i++)
            {
                m1[i] = m1[i] / norm;
            }

            for (i = 0; i < 3; i++)
            {
                m2[i] = (v2[i] + v3[i]) / 2.0;
            }

            norm = typeMethods.r8vec_norm(3, m2);
            for (i = 0; i < 3; i++)
            {
                m2[i] = m2[i] / norm;
            }

            for (i = 0; i < 3; i++)
            {
                m3[i] = (v3[i] + v1[i]) / 2.0;
            }

            norm = typeMethods.r8vec_norm(3, m3);
            for (i = 0; i < 3; i++)
            {
                m3[i] = m3[i] / norm;
            }

            return;
        }

        public static void sphere01_triangle_vertices_to_sides(double[] v1, double[] v2,
                double[] v3, ref double as_, ref double bs, ref double cs)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_TRIANGLE_VERTICES_TO_SIDES_3D: sides of spherical triangle on unit sphere.
            //
            //  Discussion:
            //
            //    We can use the ACOS system call here, but the ARC_COSINE routine
            //    will automatically take care of cases where the input argument is
            //    (usually slightly) out of bounds.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V1[3], V2[3], V3[3], the vertices of the spherical
            //    triangle.
            //
            //    Output, double *AS, *BS, *CS, the (geodesic) length of the sides of the
            //    triangle.
            //
        {
            as_ = Helpers.arc_cosine(typeMethods.r8vec_dot_product(3, v2, v3));
            bs = Helpers.arc_cosine(typeMethods.r8vec_dot_product(3, v3, v1));
            cs = Helpers.arc_cosine(typeMethods.r8vec_dot_product(3, v1, v2));
        }
    }
}