using System;
using Burkardt.Types;

namespace Burkardt.SphereNS;

public static class Quad
{
    public static double sphere01_quad_icos1c(int factor,
            Func<int, double[], double[], double[]> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_ICOS1C: centroid rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over the surface of the unit sphere.
        //
        //    This function sets up an icosahedral grid, and subdivides each
        //    edge of the icosahedron into FACTOR subedges.  These edges define a grid
        //    within each triangular icosahedral face.  The centroids of these
        //    triangles can be determined.  All of these calculations are done,
        //    essentially, on the FLAT faces of the icosahedron.  Only then are
        //    the triangle vertices and centroids projected to the sphere.  
        //
        //    The resulting grid of spherical triangles and projected centroids
        //    is used to apply a centroid quadrature rule over the surface of
        //    the unit sphere.
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
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Output, int *NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_QUAD_ICOS1C, the estimated integral.
        //
    {
        double[] a_xyz = new double[3];
        double[] b_xyz = new double[3];
        double[] c_xyz = new double[3];
        int edge_num = 0;
        int face;
        int face_num = 0;
        int face_order_max = 0;
        int point_num = 0;
        double[] v = new double[1];
        //
        //  Size the icosahedron.
        //
        Icosahedron.Geometry.icos_size(ref point_num, ref edge_num, ref face_num, ref face_order_max);
        //
        //  Set the icosahedron.
        //
        double[] point_coord = new double[3 * point_num];
        int[] edge_point = new int[2 * edge_num];
        int[] face_order = new int[face_num];
        int[] face_point = new int[face_order_max * face_num];

        Icosahedron.Geometry.icos_shape(point_num, edge_num, face_num, face_order_max,
            ref point_coord, ref edge_point, ref face_order, ref face_point);
        //
        //  Initialize the integral data.
        //
        double result = 0.0;
        node_num = 0;
        //
        //  Pick a face of the icosahedron, and identify its vertices as A, B, C.
        //
        for (face = 0; face < face_num; face++)
        {
            int a = face_point[0 + face * 3];
            int b = face_point[1 + face * 3];
            int c = face_point[2 + face * 3];

            int i;
            for (i = 0; i < 3; i++)
            {
                a_xyz[i] = point_coord[i + a * 3];
                b_xyz[i] = point_coord[i + b * 3];
                c_xyz[i] = point_coord[i + c * 3];
            }

            //
            //  Some subtriangles will have the same direction as the face.
            //  Generate each in turn, by determining the barycentric coordinates
            //  of the centroid (F1,F2,F3), from which we can also work out the barycentric
            //  coordinates of the vertices of the subtriangle.
            //
            double[] a2_xyz;
            double area;
            double[] node_xyz;
            int f1;
            int f2;
            int f3;
            double[] c2_xyz;
            double[] b2_xyz;
            for (f3 = 1; f3 <= 3 * factor - 2; f3 += 3)
            {
                for (f2 = 1; f2 <= 3 * factor - f3 - 1; f2 += 3)
                {
                    f1 = 3 * factor - f3 - f2;

                    node_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1);
                    b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1);
                    c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    fun(1, node_xyz, v);

                    node_num += 1;
                    result += area * v[0];
                }
            }

            //
            //  The other subtriangles have the opposite direction from the face.
            //  Generate each in turn, by determining the barycentric coordinates
            //  of the centroid (F1,F2,F3), from which we can also work out the barycentric
            //  coordinates of the vertices of the subtriangle.
            //
            for (f3 = 2; f3 <= 3 * factor - 4; f3 += 3)
            {
                for (f2 = 2; f2 <= 3 * factor - f3 - 2; f2 += 3)
                {
                    f1 = 3 * factor - f3 - f2;

                    node_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1);
                    b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1);
                    c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    fun(1, node_xyz, v);

                    node_num += 1;
                    result += area * v[0];
                }
            }
        }

        return result;
    }

    public static double sphere01_quad_icos1m(int factor,
            Func<int, double[], double[], double[]> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_ICOS1M: midside rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over the surface of the unit sphere.
        //
        //    This function sets up an icosahedral grid, and subdivides each
        //    edge of the icosahedron into FACTOR subedges.  These edges define a grid
        //    within each triangular icosahedral face.  The midsides of these
        //    triangles can be determined.  All of these calculations are done,
        //    essentially, on the FLAT faces of the icosahedron.  Only then are
        //    the triangle vertices and midsides projected to the sphere.  
        //
        //    The resulting grid of spherical triangles and projected midsides
        //    is used to apply a midside quadrature rule over the surface of
        //    the unit sphere.
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
        //  Parameters:
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Output, int *NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_QUAD_ICOS1M, the estimated integral.
        //
    {
        double[] a_xyz = new double[3];
        double[] b_xyz = new double[3];
        double[] c_xyz = new double[3];
        int edge_num = 0;
        int face;
        int face_num = 0;
        int face_order_max = 0;
        int point_num = 0;
        double[] va = new double[1];
        double[] vb = new double[1];
        double[] vc = new double[1];
        //
        //  Size the icosahedron.
        //
        Icosahedron.Geometry.icos_size(ref point_num, ref edge_num, ref face_num, ref face_order_max);
        //
        //  Set the icosahedron.
        //
        double[] point_coord = new double[3 * point_num];
        int[] edge_point = new int[2 * edge_num];
        int[] face_order = new int[face_num];
        int[] face_point = new int[face_order_max * face_num];

        Icosahedron.Geometry.icos_shape(point_num, edge_num, face_num, face_order_max,
            ref point_coord, ref edge_point, ref face_order, ref face_point);
        //
        //  Initialize the integral data.
        //
        double result = 0.0;
        node_num = 0;
        //
        //  Pick a face of the icosahedron, and identify its vertices as A, B, C.
        //
        for (face = 0; face < face_num; face++)
        {
            int a = face_point[0 + face * 3];
            int b = face_point[1 + face * 3];
            int c = face_point[2 + face * 3];

            int i;
            for (i = 0; i < 3; i++)
            {
                a_xyz[i] = point_coord[i + a * 3];
                b_xyz[i] = point_coord[i + b * 3];
                c_xyz[i] = point_coord[i + c * 3];
            }

            //
            //  Deal with subtriangles that have same orientation as face.
            //
            double[] a3_xyz;
            int f1;
            int f2;
            int f3;
            double[] c3_xyz;
            double[] c2_xyz;
            double[] b3_xyz;
            double[] b2_xyz;
            double[] a2_xyz;
            double area;
            for (f1 = 0; f1 <= factor - 1; f1++)
            {
                for (f2 = 0; f2 <= factor - f1 - 1; f2++)
                {
                    f3 = factor - f1 - f2;

                    a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1, f2, f3 - 1);
                    b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2 + 1, f3 - 1);
                    c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    a3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2 + 1,
                        2 * f3 - 2);
                    b3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1, 2 * f2 + 1,
                        2 * f3 - 1);
                    c3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2,
                        2 * f3 - 1);

                    node_num += 3;
                    fun(1, a3_xyz, va);
                    fun(1, b3_xyz, vb);
                    fun(1, c3_xyz, vc);
                    result += area * (va[0] + vb[0] + vc[0]) / 3.0;
                }
            }

            //
            //  Deal with subtriangles that have opposite orientation as face.
            //
            for (f3 = 0; f3 <= factor - 2; f3++)
            {
                for (f2 = 1; f2 <= factor - f3 - 1; f2++)
                {
                    f1 = factor - f2 - f3;

                    a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1, f2, f3 + 1);
                    b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2 - 1, f3 + 1);
                    c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    a3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2 - 1,
                        2 * f3 + 2);
                    b3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1, 2 * f2 - 1,
                        2 * f3 + 1);
                    c3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2,
                        2 * f3 + 1);

                    node_num += 3;
                    fun(1, a3_xyz, va);
                    fun(1, b3_xyz, vb);
                    fun(1, c3_xyz, vc);
                    result += area * (va[0] + vb[0] + vc[0]) / 3.0;
                }
            }
        }

        return result;
    }

    public static double sphere01_quad_icos1v(int factor,
            Func<int, double[], double[], double[]> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_ICOS1V: vertex rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over the surface of the unit sphere.
        //
        //    This function sets up an icosahedral grid, and subdivides each
        //    edge of the icosahedron into FACTOR subedges.  These edges define a grid
        //    within each triangular icosahedral face.  The vertices of these
        //    triangles can be determined.  All of these calculations are done,
        //    essentially, on the FLAT faces of the icosahedron.  Only then are
        //    the triangle vertices projected to the sphere.  
        //
        //    The resulting grid of spherical triangles is used to apply a vertex
        //    quadrature rule over the surface of the unit sphere.
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
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Output, int *NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_QUAD_ICOS2M, the estimated integral.
        //
    {
        double[] a_xyz = new double[3];
        double[] b_xyz = new double[3];
        double[] c_xyz = new double[3];
        int edge_num = 0;
        int face;
        int face_num = 0;
        int face_order_max = 0;
        int point_num = 0;
        double[] va = new double[1];
        double[] vb = new double[1];
        double[] vc = new double[1];
        //
        //  Size the icosahedron.
        //
        Icosahedron.Geometry.icos_size(ref point_num, ref edge_num, ref face_num, ref face_order_max);
        //
        //  Set the icosahedron.
        //
        double[] point_coord = new double[3 * point_num];
        int[] edge_point = new int[2 * edge_num];
        int[] face_order = new int[face_num];
        int[] face_point = new int[face_order_max * face_num];

        Icosahedron.Geometry.icos_shape(point_num, edge_num, face_num, face_order_max,
            ref point_coord, ref edge_point, ref face_order, ref face_point);
        //
        //  Initialize the integral data.
        //
        double result = 0.0;
        node_num = 0;
        //
        //  Pick a face of the icosahedron, and identify its vertices as A, B, C.
        //
        for (face = 0; face < face_num; face++)
        {
            int a = face_point[0 + face * 3];
            int b = face_point[1 + face * 3];
            int c = face_point[2 + face * 3];

            int i;
            for (i = 0; i < 3; i++)
            {
                a_xyz[i] = point_coord[i + a * 3];
                b_xyz[i] = point_coord[i + b * 3];
                c_xyz[i] = point_coord[i + c * 3];
            }

            //
            //  Deal with subtriangles that have same orientation as face.
            //
            double[] a2_xyz;
            int f1;
            double area;
            int f3;
            int f2;
            double[] c2_xyz;
            double[] b2_xyz;
            for (f1 = 0; f1 <= factor - 1; f1++)
            {
                for (f2 = 0; f2 <= factor - f1 - 1; f2++)
                {
                    f3 = factor - f1 - f2;

                    a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1, f2, f3 - 1);
                    b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2 + 1, f3 - 1);
                    c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    node_num += 3;
                    fun(1, a2_xyz, va);
                    fun(1, b2_xyz, vb);
                    fun(1, c2_xyz, vc);
                    result += area * (va[0] + vb[0] + vc[0]) / 3.0;
                }
            }

            //
            //  Deal with subtriangles that have opposite orientation as face.
            //
            for (f3 = 0; f3 <= factor - 2; f3++)
            {
                for (f2 = 1; f2 <= factor - f3 - 1; f2++)
                {
                    f1 = factor - f2 - f3;

                    a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1, f2, f3 + 1);
                    b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2 - 1, f3 + 1);
                    c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    node_num += 3;
                    fun(1, a2_xyz, va);
                    fun(1, b2_xyz, vb);
                    fun(1, c2_xyz, vc);
                    result += area * (va[0] + vb[0] + vc[0]) / 3.0;
                }
            }
        }

        return result;
    }

    public static double sphere01_quad_icos2v(int factor,
            Func<int, double[], double[], double[]> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_ICOS2V: vertex rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over the surface of the unit sphere.
        //
        //    This function sets up an icosahedral grid, and subdivides each
        //    edge of the icosahedron into FACTOR subedges.  These edges define a grid
        //    within each triangular icosahedral face.  The vertices of these
        //    triangles can be determined.  All of these calculations are done,
        //    essentially, on the FLAT faces of the icosahedron.  Only then are
        //    the triangle vertices projected to the sphere.  
        //
        //    The resulting grid of spherical triangles is used to apply a vertex
        //    quadrature rule over the surface of the unit sphere.
        //
        //    This is a revision of SPHERE01_QUAD_ICOS2V that attempted to use a more
        //    sophisticated scheme to map points from the planar triangle to the surface
        //    of the unit sphere.  Very little improvement to the estimated integral
        //    was observed, so development of this scheme has been set aside for now.
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
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Output, int *NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_QUAD_ICOS2V, the estimated integral.
        //
    {
        double[] a_xyz = new double[3];
        double[] b_xyz = new double[3];
        double[] c_xyz = new double[3];
        int edge_num = 0;
        int face;
        int face_num = 0;
        int face_order_max = 0;
        int point_num = 0;
        double[] va = new double[1];
        double[] vb = new double[1];
        double[] vc = new double[1];
        //
        //  Size the icosahedron.
        //
        Icosahedron.Geometry.icos_size(ref point_num, ref edge_num, ref face_num, ref face_order_max);
        //
        //  Set the icosahedron.
        //
        double[] point_coord = new double[3 * point_num];
        int[] edge_point = new int[2 * edge_num];
        int[] face_order = new int[face_num];
        int[] face_point = new int[face_order_max * face_num];

        Icosahedron.Geometry.icos_shape(point_num, edge_num, face_num, face_order_max,
            ref point_coord, ref edge_point, ref face_order, ref face_point);
        //
        //  Initialize the integral data.
        //
        double result = 0.0;
        node_num = 0;
        //
        //  Pick a face of the icosahedron, and identify its vertices as A, B, C.
        //
        for (face = 0; face < face_num; face++)
        {
            int a = face_point[0 + face * 3];
            int b = face_point[1 + face * 3];
            int c = face_point[2 + face * 3];

            int i;
            for (i = 0; i < 3; i++)
            {
                a_xyz[i] = point_coord[i + a * 3];
                b_xyz[i] = point_coord[i + b * 3];
                c_xyz[i] = point_coord[i + c * 3];
            }

            //
            //  Deal with subtriangles that have same orientation as face.
            //
            double[] a2_xyz;
            double[] c2_xyz;
            double[] b2_xyz;
            int f1;
            int f3;
            double area;
            int f2;
            for (f1 = 0; f1 <= factor - 1; f1++)
            {
                for (f2 = 0; f2 <= factor - f1 - 1; f2++)
                {
                    f3 = factor - f1 - f2;

                    a2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1 + 1, f2, f3 - 1);
                    b2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1, f2 + 1, f3 - 1);
                    c2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    node_num += 3;
                    fun(1, a2_xyz, va);
                    fun(1, b2_xyz, vb);
                    fun(1, c2_xyz, vc);
                    result += area * (va[0] + vb[0] + vc[0]) / 3.0;
                }
            }

            //
            //  Deal with subtriangles that have opposite orientation as face.
            //
            for (f3 = 0; f3 <= factor - 2; f3++)
            {
                for (f2 = 1; f2 <= factor - f3 - 1; f2++)
                {
                    f1 = factor - f2 - f3;

                    a2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1 - 1, f2, f3 + 1);
                    b2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1, f2 - 1, f3 + 1);
                    c2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                    area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                    node_num += 3;
                    fun(1, a2_xyz, va);
                    fun(1, b2_xyz, vb);
                    fun(1, c2_xyz, vc);
                    result += area * (va[0] + vb[0] + vc[0]) / 3.0;
                }
            }
        }

        return result;
    }

    public static double sphere01_quad_llc(Func<int, double[], double[], double[]> f, double h,
            ref int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_LLC: Longitude/Latitude grid with centroid rule.
        //
        //  Discussion:
        //
        //    The sphere is broken up into spherical triangles, whose sides
        //    do not exceed the length H.  Then a centroid rule is used on
        //    each spherical triangle.
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
        //    Input, void F ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Input, double H, the maximum length of a side of the spherical
        //    quadrilaterals.
        //
        //    Output, int *N, the number of points used.
        //
        //    Output, double SPHERE01_QUAD_LLC, the approximate integral.
        //
    {
        double phi;

        double sphere_area;
        double theta;
        double[] v = new double[1];
        double[] x;
        //
        //  Choose PHI and THETA counts that make short sides.
        //
        int phi_num = (int)(Math.PI / h);

        switch (h * phi_num)
        {
            case < Math.PI:
                phi_num += 1;
                break;
        }

        int theta_num = (int)(2.0 * Math.PI / h);

        switch (h * theta_num)
        {
            case < Math.PI:
                theta_num += 1;
                break;
        }

        n = 0;
        double result;
        switch (theta_num)
        {
            //
            //  Only one THETA (and hence, only one PHI.)
            //
            case 1:
                sphere_area = 4.0 * Math.PI;

                theta = 0.0;
                phi = Math.PI / 2.0;
                x = typeMethods.tp_to_xyz(theta, phi);

                v = f(1, x, v);
                n += 1;
                result = sphere_area * v[0];
                break;
            //
            default:
            {
                int j;
                switch (phi_num)
                {
                    case 1:
                    {
                        sphere_area = 4.0 * Math.PI;
                        double sector_area = sphere_area / theta_num;

                        result = 0.0;

                        for (j = 1; j <= theta_num; j++)
                        {
                            theta = (j - 1) * 2 * Math.PI / theta_num;
                            phi = Math.PI / 2.0;
                            x = typeMethods.tp_to_xyz(theta, phi);
                            v = f(1, x, v);
                            n += 1;
                            result += sector_area * v[0];
                        }

                        break;
                    }
                    //
                    default:
                    {
                        result = 0.0;
                        //
                        //  Picture in top row, with V1 = north pole:
                        //
                        //        V1
                        //       .  .
                        //      .    .
                        //    V12----V22
                        //
                        double phi1 = 0.0;
                        double phi2 = Math.PI / phi_num;

                        double area;
                        double[] x22;
                        double theta1;
                        double theta2;
                        double[] x12;
                        for (j = 1; j <= theta_num; j++)
                        {
                            theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                            theta2 = j * 2.0 * Math.PI / theta_num;

                            double[] x1 = typeMethods.tp_to_xyz(theta1, phi1);
                            x12 = typeMethods.tp_to_xyz(theta1, phi2);
                            x22 = typeMethods.tp_to_xyz(theta2, phi2);

                            area = Triangle.sphere01_triangle_vertices_to_area(x1, x12, x22);
                            x = Triangle.sphere01_triangle_vertices_to_centroid(x1, x12, x22);
                            v = f(1, x, v);
                            n += 1;
                            result += area * v[0];
                        }

                        //
                        //  Picture in all intermediate rows:
                        //
                        //    V11--V21
                        //     | .  |
                        //     |  . |
                        //    V12--V22
                        //
                        int i;
                        double[] x21;
                        double[] x11;
                        for (i = 2; i <= phi_num - 1; i++)
                        {
                            phi1 = (i - 1) * Math.PI / phi_num;
                            phi2 = i * Math.PI / phi_num;

                            for (j = 1; j <= theta_num; j++)
                            {
                                theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                                theta2 = j * 2.0 * Math.PI / theta_num;

                                x11 = typeMethods.tp_to_xyz(theta1, phi1);
                                x21 = typeMethods.tp_to_xyz(theta2, phi1);
                                x12 = typeMethods.tp_to_xyz(theta1, phi2);
                                x22 = typeMethods.tp_to_xyz(theta2, phi2);

                                area = Triangle.sphere01_triangle_vertices_to_area(x11, x12, x22);
                                x = Triangle.sphere01_triangle_vertices_to_centroid(x11, x12, x22);
                                v = f(1, x, v);
                                n += 1;
                                result += area * v[0];

                                area = Triangle.sphere01_triangle_vertices_to_area(x22, x21, x11);
                                x = Triangle.sphere01_triangle_vertices_to_centroid(x22, x21, x11);
                                v = f(1, x, v);
                                n += 1;
                                result += area * v[0];
                            }
                        }

                        //
                        //  Picture in last row, with V2 = south pole:
                        //
                        //    V11----V21
                        //      .    .
                        //       .  .
                        //        V2
                        //
                        phi1 = (phi_num - 1) * Math.PI / phi_num;
                        phi2 = Math.PI;

                        for (j = 1; j <= theta_num; j++)
                        {
                            theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                            theta2 = j * 2.0 * Math.PI / theta_num;

                            x11 = typeMethods.tp_to_xyz(theta1, phi1);
                            x21 = typeMethods.tp_to_xyz(theta2, phi1);
                            double[] x2 = typeMethods.tp_to_xyz(theta2, phi2);

                            area = Triangle.sphere01_triangle_vertices_to_area(x11, x2, x21);
                            x = Triangle.sphere01_triangle_vertices_to_centroid(x11, x2, x21);
                            v = f(1, x, v);
                            n += 1;
                            result += area * v[0];
                        }

                        break;
                    }
                }

                break;
            }
        }

        return result;
    }

    public static double sphere01_quad_llm(Func<int, double[], double[], double[]> f, double h,
            ref int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_LLM: longitude/latitude grid plus midside rule.
        //
        //  Discussion:
        //
        //    The sphere is broken up into spherical triangles, whose sides
        //    do not exceed the length H.  Then the function is evaluated
        //    at the midsides, and the average is multiplied by the area.
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
        //    Input, void F ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Input, double H, the maximum length of a side of the spherical
        //    quadrilaterals.
        //
        //    Output, int *N, the number of points used.
        //
        //    Output, double SPHERE01_QUAD_LLM, the approximate integral.
        //
    {
        double[] m1 = new double[3];
        double[] m2 = new double[3];
        double[] m3 = new double[3];
        double phi;

        double result;
        double sphere_area;
        double theta;
        double[] v = new double[1];
        double[] x;
        //
        //  Choose PHI and THETA counts that make short sides.
        //
        int phi_num = (int)(Math.PI / h);

        switch (h * phi_num)
        {
            case < Math.PI:
                phi_num += 1;
                break;
        }

        int theta_num = (int)(2.0 * Math.PI / h);

        switch (h * theta_num)
        {
            case < Math.PI:
                theta_num += 1;
                break;
        }

        n = 0;
        switch (theta_num)
        {
            //
            //  Only one THETA (and hence, only one PHI.)
            //
            case 1:
                sphere_area = 4.0 * Math.PI;

                theta = 0.0;
                phi = Math.PI / 2.0;
                x = typeMethods.tp_to_xyz(theta, phi);
                v = f(1, x, v);
                n += 1;
                result = sphere_area * v[0];
                break;
            //
            default:
            {
                int j;
                switch (phi_num)
                {
                    case 1:
                    {
                        sphere_area = 4.0 * Math.PI;
                        double sector_area = sphere_area / theta_num;

                        result = 0.0;

                        for (j = 1; j <= theta_num; j++)
                        {
                            theta = (j - 1) * 2 * Math.PI / theta_num;
                            phi = Math.PI / 2.0;
                            x = typeMethods.tp_to_xyz(theta, phi);
                            v = f(1, x, v);
                            n += 1;
                            result += sector_area * v[0];
                        }

                        break;
                    }
                    //
                    default:
                    {
                        result = 0.0;
                        //
                        //  Picture:
                        //
                        //        V1
                        //       .  .
                        //      .    .
                        //    V12----V22
                        //
                        double phi1 = 0.0;
                        double phi2 = Math.PI / phi_num;

                        double area;
                        double[] x22;
                        double theta1;
                        double[] x12;
                        double theta2;
                        for (j = 1; j <= theta_num; j++)
                        {
                            theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                            theta2 = j * 2.0 * Math.PI / theta_num;

                            double[] x1 = typeMethods.tp_to_xyz(theta1, phi1);
                            x12 = typeMethods.tp_to_xyz(theta1, phi2);
                            x22 = typeMethods.tp_to_xyz(theta2, phi2);

                            area = Triangle.sphere01_triangle_vertices_to_area(x1, x12, x22);

                            Triangle.sphere01_triangle_vertices_to_midpoints(x1, x12, x22, ref m1, ref m2, ref m3);

                            v = f(1, m1, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, m2, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, m3, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                        }

                        //
                        //  Picture:
                        //
                        //    V11--V21
                        //     | .  |
                        //     |  . |
                        //    V12--V22
                        //
                        int i;
                        double[] x11;
                        double[] x21;
                        for (i = 2; i <= phi_num - 1; i++)
                        {
                            phi1 = (i - 1) * Math.PI / phi_num;
                            phi2 = i * Math.PI / phi_num;

                            for (j = 1; j <= theta_num; j++)
                            {
                                theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                                theta2 = j * 2.0 * Math.PI / theta_num;

                                x11 = typeMethods.tp_to_xyz(theta1, phi1);
                                x21 = typeMethods.tp_to_xyz(theta2, phi1);
                                x12 = typeMethods.tp_to_xyz(theta1, phi2);
                                x22 = typeMethods.tp_to_xyz(theta2, phi2);

                                area = Triangle.sphere01_triangle_vertices_to_area(x11, x12, x22);

                                Triangle.sphere01_triangle_vertices_to_midpoints(x11, x12, x22, ref m1, ref m2, ref m3);

                                v = f(1, m1, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, m2, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, m3, v);
                                n += 1;
                                result += area * v[0] / 3.0;

                                area = Triangle.sphere01_triangle_vertices_to_area(x22, x21, x11);

                                Triangle.sphere01_triangle_vertices_to_midpoints(x22, x21, x11, ref m1, ref m2, ref m3);

                                v = f(1, m1, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, m2, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, m3, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                            }
                        }

                        //
                        //  Picture:
                        //
                        //    V11----V21
                        //      .    .
                        //       .  .
                        //        V2
                        //
                        phi1 = (phi_num - 1) * Math.PI / phi_num;
                        phi2 = Math.PI;

                        for (j = 1; j <= theta_num; j++)
                        {
                            theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                            theta2 = j * 2.0 * Math.PI / theta_num;

                            x11 = typeMethods.tp_to_xyz(theta1, phi1);
                            x21 = typeMethods.tp_to_xyz(theta2, phi1);
                            double[] x2 = typeMethods.tp_to_xyz(theta2, phi2);

                            area = Triangle.sphere01_triangle_vertices_to_area(x11, x2, x21);

                            Triangle.sphere01_triangle_vertices_to_midpoints(x11, x2, x21, ref m1, ref m2, ref m3);

                            v = f(1, m1, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, m2, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, m3, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                        }

                        break;
                    }
                }

                break;
            }
        }

        return result;
    }

    public static double sphere01_quad_llv(Func<int, double[], double[], double[]> f, double h,
            ref int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_LLV: longitude/latitude grid with vertex rule.
        //
        //  Discussion:
        //
        //    The sphere is broken up into spherical triangles, whose sides
        //    do not exceed the length H.  Then the function is evaluated
        //    at the vertices, and the average is multiplied by the area.
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
        //    Input, void F ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Input, double H, the maximum length of a side of the spherical
        //    quadrilaterals.
        //
        //    Output, int *N, the number of points used.
        //
        //    Output, double SPHERE01_QUAD_LLV, the approximate integral.
        //
    {
        double phi;

        double sphere_area;
        double theta;
        double[] v = new double[1];
        double[] x;
        //
        //  Choose PHI and THETA counts that make short sides.
        //
        int phi_num = (int)(Math.PI / h);

        switch (h * phi_num)
        {
            case < Math.PI:
                phi_num += 1;
                break;
        }

        int theta_num = (int)(2.0 * Math.PI / h);

        switch (h * theta_num)
        {
            case < Math.PI:
                theta_num += 1;
                break;
        }

        n = 0;
        double result;
        switch (theta_num)
        {
            //
            //  Only one THETA (and hence, only one PHI.)
            //
            case 1:
                sphere_area = 4.0 * Math.PI;

                theta = 0.0;
                phi = Math.PI / 2.0;
                x = typeMethods.tp_to_xyz(theta, phi);
                v = f(1, x, v);
                result = sphere_area * v[0];
                break;
            //
            default:
            {
                int j;
                switch (phi_num)
                {
                    case 1:
                    {
                        sphere_area = 4.0 * Math.PI;
                        double sector_area = sphere_area / theta_num;

                        result = 0.0;

                        for (j = 1; j <= theta_num; j++)
                        {
                            theta = (j - 1) * 2 * Math.PI / theta_num;
                            phi = Math.PI / 2.0;
                            x = typeMethods.tp_to_xyz(theta, phi);
                            v = f(1, x, v);
                            n += 1;
                            result += sector_area * v[0];
                        }

                        break;
                    }
                    //
                    default:
                    {
                        result = 0.0;
                        //
                        //  Picture:
                        //
                        //        V1
                        //       .  .
                        //      .    .
                        //    V12----V22
                        //
                        double phi1 = 0.0;
                        double phi2 = Math.PI / phi_num;

                        double area;
                        double theta1;
                        double theta2;
                        double[] x22;
                        double[] x12;
                        for (j = 1; j <= theta_num; j++)
                        {
                            theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                            theta2 = j * 2.0 * Math.PI / theta_num;

                            double[] x1 = typeMethods.tp_to_xyz(theta1, phi1);
                            x12 = typeMethods.tp_to_xyz(theta1, phi2);
                            x22 = typeMethods.tp_to_xyz(theta2, phi2);

                            area = Triangle.sphere01_triangle_vertices_to_area(x1, x12, x22);

                            v = f(1, x1, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, x12, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, x22, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                        }

                        //
                        //  Picture:
                        //
                        //    V11--V21
                        //     | .  |
                        //     |  . |
                        //    V12--V22
                        //
                        int i;
                        double[] x21;
                        double[] x11;
                        for (i = 2; i <= phi_num - 1; i++)
                        {
                            phi1 = (i - 1) * Math.PI / phi_num;
                            phi2 = i * Math.PI / phi_num;

                            for (j = 1; j <= theta_num; j++)
                            {
                                theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                                theta2 = j * 2.0 * Math.PI / theta_num;

                                x11 = typeMethods.tp_to_xyz(theta1, phi1);
                                x21 = typeMethods.tp_to_xyz(theta2, phi1);
                                x12 = typeMethods.tp_to_xyz(theta1, phi2);
                                x22 = typeMethods.tp_to_xyz(theta2, phi2);

                                area = Triangle.sphere01_triangle_vertices_to_area(x11, x12, x22);

                                v = f(1, x11, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, x12, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, x22, v);
                                n += 1;
                                result += area * v[0] / 3.0;

                                area = Triangle.sphere01_triangle_vertices_to_area(x22, x21, x11);

                                v = f(1, x22, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, x21, v);
                                n += 1;
                                result += area * v[0] / 3.0;
                                v = f(1, x11, v);
                                n += 1;
                                result += area * v[0] / 3.0;

                            }
                        }

                        //
                        //  Picture:
                        //
                        //    V11----V21
                        //      \    /
                        //       \  /
                        //        V2
                        //
                        phi1 = (phi_num - 1) * Math.PI / phi_num;
                        phi2 = Math.PI;

                        for (j = 1; j <= theta_num; j++)
                        {
                            theta1 = (j - 1) * 2.0 * Math.PI / theta_num;
                            theta2 = j * 2.0 * Math.PI / theta_num;

                            x11 = typeMethods.tp_to_xyz(theta1, phi1);
                            x21 = typeMethods.tp_to_xyz(theta2, phi1);
                            double[] x2 = typeMethods.tp_to_xyz(theta2, phi2);

                            area = Triangle.sphere01_triangle_vertices_to_area(x11, x2, x21);

                            v = f(1, x11, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, x2, v);
                            n += 1;
                            result += area * v[0] / 3.0;
                            v = f(1, x21, v);
                            n += 1;
                            result += area * v[0] / 3.0;

                        }

                        break;
                    }
                }

                break;
            }
        }

        return result;
    }

    public static double sphere01_quad_mc(Func<int, double[], double[], double[]> f, double h,
            ref int seed, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_MC uses the Monte Carlo rule for sphere quadrature.
        //
        //  Discussion:
        //
        //    A number of points N are chosen at random on the sphere, with N
        //    being determined so that, if the points were laid out on a regular
        //    grid, the average spacing would be no more than H.
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
        //    Input, void F ( int n, double x[], double v[] ), evaluates the 
        //    integrand.
        //
        //    Input, double H, the maximum length of a side of the spherical
        //    quadrilaterals.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Input, int N, the number of points used.
        //
        //    Output, double SPHERE01_QUAD_MC, the approximate integral.
        //
    {
        const double sphere_area = 4.0 * Math.PI;

        double[] x = MonteCarlo.sphere01_sample(n, ref seed);

        double[] v = new double[n];

        v = f(n, x, v);

        double result = sphere_area * typeMethods.r8vec_sum(n, v) / n;

        return result;
    }

    public static int sphere01_quad_mc_size(double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_QUAD_MC_SIZE sizes a Monte Carlo rule for sphere quadrature.
        //
        //  Discussion:
        //
        //    A number of points N are chosen at random on the sphere, with N
        //    being determined so that, if the points were laid out on a regular
        //    grid, the average spacing would be no more than H.
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
        //    Input, double H, the maximum length of a side of the spherical
        //    quadrilaterals.
        //
        //    Output, int SPHERE01_QUAD_MC_SIZE, the number of points to use.
        //
    {
        //
        //  The sphere's area is 4 * PI.
        //  Choose N so that we divide this area into N subareas of PI * H * H.
        //
        const double sphere_area = 4.0 * Math.PI;

        int n = (int)(sphere_area / h / h);
        n = Math.Max(n, 1);

        return n;
    }
}