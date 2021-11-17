using System;

namespace Burkardt.SphereNS;

public static class TriangleQuad
{
    public static double sphere01_triangle_quad_00(int n, double[] v1, double[] v2,
            double[] v3, Func<double[], int, double> f, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_00: quadrature over a triangle on the unit sphere.
        //
        //  Discussion:
        //
        //    This is a Monte Carlo approach.
        //
        //    The integral is approximated by averaging the values at N random points,
        //    multiplied by the area.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, integer N, the number of sample points.
        //
        //    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
        //    the vertices of the triangle.
        //
        //    Input, double F ( double x[] ), evaluates the integrand.
        //
        //    Input/output, integer &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_00, the approximate integral.
        //
    {
        double area;
        int j;
        double quad;
        double result;
        double[] vc;

        area = Triangle.sphere01_triangle_vertices_to_area(v1, v2, v3);

        vc = Triangle.sphere01_triangle_sample(n, v1, v2, v3, ref seed);

        quad = 0.0;
        for (j = 0; j < n; j++)
        {
            quad += f(vc, +3 * j);
        }

        result = quad * area / n;


        return result;
    }

    public static double sphere01_triangle_quad_01(double[] v1, double[] v2, double[] v3,
            Func<double[], int, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_01: quadrature over a triangle on the unit sphere.
        //
        //  Discussion:
        //
        //    The integral is approximated by the value at the centroid,
        //    multiplied by the area.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
        //    the vertices of the triangle.
        //
        //    Input, double F ( double x[] ), evaluates the integrand.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_01, the approximate integral.
        //
    {
        double area;
        double quad;
        double result;
        double[] vc;

        area = Triangle.sphere01_triangle_vertices_to_area(v1, v2, v3);

        vc = Triangle.sphere01_triangle_vertices_to_centroid(v1, v2, v3);

        quad = f(vc, 0);
        result = quad * area;

        return result;
    }

    public static double sphere01_triangle_quad_02(double[] v1, double[] v2, double[] v3,
            Func<double[], int, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_02: quadrature over a triangle on the unit sphere.
        //
        //  Discussion:
        //
        //    The integral is approximated by the average of the vertex values,
        //    multiplied by the area.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
        //    the vertices of the triangle.
        //
        //    Input, double F ( double x[] ), evaluates the integrand.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_02, the approximate integral.
        //
    {
        double area;
        double quad;
        double result;

        area = Triangle.sphere01_triangle_vertices_to_area(v1, v2, v3);

        quad = (f(v1, 0) + f(v2, 0) + f(v3, 0)) / 3.0;

        result = quad * area;

        return result;
    }

    public static double sphere01_triangle_quad_03(double[] v1, double[] v2, double[] v3,
            Func<double[], int, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_03: quadrature over a triangle on the unit sphere.
        //
        //  Discussion:
        //
        //    The integral is approximated by the average of the midside values,
        //    multiplied by the area.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
        //    the vertices of the triangle.
        //
        //    Input, double F ( double x[] ), evaluates the integrand.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_03, the approximate integral.
        //
    {
        double area;
        double quad;
        double result;
        double[] v4 = new double[3];
        double[] v5 = new double[3];
        double[] v6 = new double[3];

        area = Triangle.sphere01_triangle_vertices_to_area(v1, v2, v3);

        Triangle.sphere01_triangle_vertices_to_midpoints(v1, v2, v3, ref v4, ref v5, ref v6);

        quad = (f(v4, 0) + f(v5, 0) + f(v6, 0)) / 3.0;

        result = quad * area;

        return result;
    }

    public static double sphere01_triangle_quad_icos1c(double[] a_xyz, double[] b_xyz,
            double[] c_xyz, int factor, Func<double[], int, double> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_ICOS1C: centroid rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over a spherical triangle on the
        //    unit sphere.
        //
        //    This function subdivides each edge of the triangle into FACTOR subedges.  
        //    These edges define a grid within the triangle.  The centroids of these
        //    triangles can be determined.  All of these calculations are done,
        //    essentially, on the FLAT faces of the planar triangle.  Only then are
        //    the triangle vertices and centroids projected to the sphere.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
        //    of the spherical triangle.
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, double FUN ( double x[] ), evaluates the integrand.
        //
        //    Output, int &NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_ICOS1C, the estimated integral.
        //
    {
        double[] a2_xyz;
        double area;
        double[] b2_xyz;
        double[] c2_xyz;
        int f1;
        int f2;
        int f3;
        double[] node_xyz;
        double result;
        double v;
        //
        //  Initialize the integral data.
        //
        result = 0.0;
        node_num = 0;
        //
        //  Some subtriangles will have the same direction as the face.
        //  Generate each in turn, by determining the barycentric coordinates
        //  of the centroid (F1,F2,F3), from which we can also work out the barycentric
        //  coordinates of the vertices of the subtriangle.
        //
        for (f3 = 1; f3 <= 3 * factor - 2; f3 += 3)
        {
            for (f2 = 1; f2 <= 3 * factor - f3 - 1; f2 += 3)
            {
                f1 = 3 * factor - f3 - f2;

                node_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3);

                a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 2,
                    f2 - 1, f3 - 1);
                b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1,
                    f2 + 2, f3 - 1);
                c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1,
                    f2 - 1, f3 + 2);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                v = fun(node_xyz, 0);

                node_num += 1;
                result += area * v;
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

                a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 2,
                    f2 + 1, f3 + 1);
                b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1,
                    f2 - 2, f3 + 1);
                c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1,
                    f2 + 1, f3 - 2);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                v = fun(node_xyz, 0);

                node_num += 1;
                result += area * v;
            }
        }

        return result;
    }

    public static double sphere01_triangle_quad_icos1m(double[] a_xyz, double[] b_xyz,
            double[] c_xyz, int factor, Func<double[], int, double> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_ICOS1M: midside rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over a spherical triangle on the
        //    unit sphere.
        //
        //    This function subdivides each edge of the triangle into FACTOR subedges.  
        //    These edges define a grid within the triangle.  The midsides of the
        //    edges of these triangles can be determined.  All of these calculations 
        //    are done, essentially, on the FLAT faces of the planar triangle.  Only
        //    then are the triangle vertices and midsides projected to the sphere.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
        //    of the spherical triangle.
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, double FUN ( double x[] ), evaluates the integrand.
        //
        //    Output, int &NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_ICOS1M, the estimated integral.
        //
    {
        double[] a2_xyz;
        double[] a3_xyz;
        double area;
        double[] b2_xyz;
        double[] b3_xyz;
        double[] c2_xyz;
        double[] c3_xyz;
        int f1;
        int f2;
        int f3;
        double result;
        double va;
        double vb;
        double vc;
        //
        //  Initialize the integral data.
        //
        result = 0.0;
        node_num = 0;
        //
        //  Some subtriangles will have the same direction as the face.
        //
        for (f1 = 0; f1 <= factor - 1; f1++)
        {
            for (f2 = 0; f2 <= factor - f1 - 1; f2++)
            {
                f3 = factor - f1 - f2;

                a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1,
                    f2, f3 - 1);
                b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2 + 1, f3 - 1);
                c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2, f3);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                a3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 + 1,
                    2 * f2 + 1, 2 * f3 - 2);
                b3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1,
                    2 * f2 + 1, 2 * f3 - 1);
                c3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 + 1,
                    2 * f2, 2 * f3 - 1);

                node_num += 3;
                va = fun(a3_xyz, 0);
                vb = fun(b3_xyz, 0);
                vc = fun(c3_xyz, 0);
                result += area * (va + vb + vc) / 3.0;
            }
        }

        //
        //  The other subtriangles have the opposite direction from the face.
        //
        for (f3 = 0; f3 <= factor - 2; f3++)
        {
            for (f2 = 1; f2 <= factor - f3 - 1; f2++)
            {
                f1 = factor - f2 - f3;

                a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1,
                    f2, f3 + 1);
                b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2 - 1, f3 + 1);
                c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2, f3);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                a3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 - 1,
                    2 * f2 - 1, 2 * f3 + 2);
                b3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1,
                    2 * f2 - 1, 2 * f3 + 1);
                c3_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, 2 * f1 - 1,
                    2 * f2, 2 * f3 + 1);

                node_num += 3;
                va = fun(a3_xyz, 0);
                vb = fun(b3_xyz, 0);
                vc = fun(c3_xyz, 0);
                result += area * (va + vb + vc) / 3.0;
            }
        }

        return result;
    }

    public static double sphere01_triangle_quad_icos1v(double[] a_xyz, double[] b_xyz,
            double[] c_xyz, int factor, Func<double[], int, double> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_ICOS1V: vertex rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over a spherical triangle on the
        //    unit sphere.
        //
        //    This function subdivides each edge of the triangle into FACTOR subedges.  
        //    These edges define a grid within the triangle.  The vertices of these
        //    triangles can be determined.  All of these calculations 
        //    are done, essentially, on the FLAT faces of the planar triangle.  Only
        //    then are the triangle vertices projected to the sphere.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
        //    of the spherical triangle.
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, double FUN ( double x[] ), evaluates the integrand.
        //
        //    Output, int &NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_ICOS1V, the estimated integral.
        //
    {
        double[] a2_xyz;
        double area;
        double[] b2_xyz;
        double[] c2_xyz;
        int f1;
        int f2;
        int f3;
        double result;
        double va;
        double vb;
        double vc;
        //
        //  Initialize the integral data.
        //
        result = 0.0;
        node_num = 0;
        //
        //  Some subtriangles will have the same direction as the face.
        //
        for (f1 = 0; f1 <= factor - 1; f1++)
        {
            for (f2 = 0; f2 <= factor - f1 - 1; f2++)
            {
                f3 = factor - f1 - f2;

                a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + 1,
                    f2, f3 - 1);
                b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2 + 1, f3 - 1);
                c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2, f3);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                node_num += 1;
                va = fun(a2_xyz, 0);
                vb = fun(b2_xyz, 0);
                vc = fun(c2_xyz, 0);
                result += area * (va + vb + vc) / 3.0;
            }
        }

        //
        //  The other subtriangles have the opposite direction from the face.
        //
        for (f3 = 0; f3 <= factor - 2; f3++)
        {
            for (f2 = 1; f2 <= factor - f3 - 1; f2++)
            {
                f1 = factor - f2 - f3;

                a2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 - 1,
                    f2, f3 + 1);
                b2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2 - 1, f3 + 1);
                c2_xyz = Triangle.sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1,
                    f2, f3);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                node_num += 1;
                va = fun(a2_xyz, 0);
                vb = fun(b2_xyz, 0);
                vc = fun(c2_xyz, 0);
                result += area * (va + vb + vc) / 3.0;
            }
        }

        return result;
    }

    public static double sphere01_triangle_quad_icos2v(double[] a_xyz, double[] b_xyz,
            double[] c_xyz, int factor, Func<double[], int, double> fun, ref int node_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_TRIANGLE_QUAD_ICOS2V: vertex rule, subdivide then project.
        //
        //  Discussion:
        //
        //    This function estimates an integral over a spherical triangle on the
        //    unit sphere.
        //
        //    This function subdivides each edge of the triangle into FACTOR subedges.  
        //    These edges define a grid within the triangle.  The vertices of these
        //    triangles can be determined.  All of these calculations 
        //    are done, essentially, on the FLAT faces of the planar triangle.  Only
        //    then are the triangle vertices projected to the sphere.  
        //
        //    This function uses a more sophisticated projection scheme than that
        //    used by SPHERE01_TRIANGLE_QUAD_ICOS1V.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
        //    of the spherical triangle.
        //
        //    Input, int FACTOR, the subdivision factor, which must
        //    be at least 1.
        //
        //    Input, double FUN ( double x[] ), evaluates the integrand.
        //
        //    Output, int &NODE_NUM, the number of evaluation points.
        //
        //    Output, double SPHERE01_TRIANGLE_QUAD_ICOS2V, the estimated integral.
        //
    {
        double[] a2_xyz;
        double area;
        double[] b2_xyz;
        double[] c2_xyz;
        int f1;
        int f2;
        int f3;
        double result;
        double va;
        double vb;
        double vc;
        //
        //  Initialize the integral data.
        //
        result = 0.0;
        node_num = 0;
        //
        //  Some subtriangles will have the same direction as the face.
        //
        for (f1 = 0; f1 <= factor - 1; f1++)
        {
            for (f2 = 0; f2 <= factor - f1 - 1; f2++)
            {
                f3 = factor - f1 - f2;

                a2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1 + 1,
                    f2, f3 - 1);
                b2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1,
                    f2 + 1, f3 - 1);
                c2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1,
                    f2, f3);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                node_num += 1;
                va = fun(a2_xyz, 0);
                vb = fun(b2_xyz, 0);
                vc = fun(c2_xyz, 0);
                result += area * (va + vb + vc) / 3.0;
            }
        }

        //
        //  The other subtriangles have the opposite direction from the face.
        //
        for (f3 = 0; f3 <= factor - 2; f3++)
        {
            for (f2 = 1; f2 <= factor - f3 - 1; f2++)
            {
                f1 = factor - f3 - f2;

                a2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1 - 1,
                    f2, f3 + 1);
                b2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1,
                    f2 - 1, f3 + 1);
                c2_xyz = Triangle.sphere01_triangle_project2(a_xyz, b_xyz, c_xyz, f1,
                    f2, f3);

                area = Triangle.sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz);

                node_num += 1;
                va = fun(a2_xyz, 0);
                vb = fun(b2_xyz, 0);
                vc = fun(c2_xyz, 0);
                result += area * (va + vb + vc) / 3.0;
            }
        }

        return result;
    }
}