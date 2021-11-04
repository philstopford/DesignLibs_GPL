using System;
using Burkardt.Types;

namespace Burkardt.MeasureNS
{
    public static class Measure
    {
        public static void alpha_measure(int n, double[] z, int triangle_order, int triangle_num,
                int[] triangle_node, ref double alpha_min, ref double alpha_ave,
                ref double alpha_area)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
            //
            //  Discusion:
            //
            //    The ALPHA measure evaluates the uniformity of the shapes of the triangles
            //    defined by a triangulated pointset.
            //
            //    We compute the minimum angle among all the triangles in the triangulated
            //    dataset and divide by the maximum possible value (which, in degrees,
            //    is 60).  The best possible value is 1, and the worst 0.  A good
            //    triangulation should have an ALPHA score close to 1.
            //
            //    The code has been modified to 'allow' 6-node triangulations.
            //    However, no effort is made to actually process the midside nodes.
            //    Only information from the vertices is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, real ( kind = 8 ) Z(2,N), the points.
            //
            //    Input, int TRIANGLE_ORDER, the order of the triangles.
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
            //    the triangulation.
            //
            //    Output, double *ALPHA_MIN, the minimum value of ALPHA over all
            //    triangles.
            //
            //    Output, double *ALPHA_AVE, the value of ALPHA averaged over
            //    all triangles.
            //
            //    Output, double *ALPHA_AREA, the value of ALPHA averaged over
            //    all triangles and weighted by area.
            //
        {
            double a_angle;
            int a_index;
            double a_x;
            double a_y;
            double ab_len;
            double area;
            double area_total;
            double b_angle;
            int b_index;
            double b_x;
            double b_y;
            double bc_len;
            double c_angle;
            int c_index;
            double c_x;
            double c_y;
            double ca_len;
            
            int triangle;

            alpha_min = typeMethods.r8_huge();
            alpha_ave = 0.0;
            alpha_area = 0.0;
            area_total = 0.0;

            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                a_index = triangle_node[0 + triangle * triangle_order];
                b_index = triangle_node[1 + triangle * triangle_order];
                c_index = triangle_node[2 + triangle * triangle_order];

                a_x = z[((0 + (a_index - 1) * 2) + z.Length) % z.Length];
                a_y = z[((1 + (a_index - 1) * 2) + z.Length) % z.Length];
                b_x = z[((0 + (b_index - 1) * 2) + z.Length) % z.Length];
                b_y = z[((1 + (b_index - 1) * 2) + z.Length) % z.Length];
                c_x = z[((0 + (c_index - 1) * 2) + z.Length) % z.Length];
                c_y = z[((1 + (c_index - 1) * 2) + z.Length) % z.Length];

                area = 0.5 * Math.Abs(a_x * (b_y - c_y)
                                      + b_x * (c_y - a_y)
                                      + c_x * (a_y - b_y));

                ab_len = Math.Sqrt(Math.Pow(a_x - b_x, 2) + Math.Pow(a_y - b_y, 2));
                bc_len = Math.Sqrt(Math.Pow(b_x - c_x, 2) + Math.Pow(b_y - c_y, 2));
                ca_len = Math.Sqrt(Math.Pow(c_x - a_x, 2) + Math.Pow(c_y - a_y, 2));
                //
                //  Take care of a ridiculous special case.
                //
                if (ab_len == 0.0 && bc_len == 0.0 && ca_len == 0.0)
                {
                    a_angle = 2.0 * Math.PI / 3.0;
                    b_angle = 2.0 * Math.PI / 3.0;
                    c_angle = 2.0 * Math.PI / 3.0;
                }
                else
                {
                    if (ca_len == 0.0 || ab_len == 0.0)
                    {
                        a_angle = Math.PI;
                    }
                    else
                    {
                        a_angle = Helpers.arc_cosine(
                            (ca_len * ca_len + ab_len * ab_len - bc_len * bc_len)
                            / (2.0 * ca_len * ab_len));
                    }

                    if (ab_len == 0.0 || bc_len == 0.0)
                    {
                        b_angle = Math.PI;
                    }
                    else
                    {
                        b_angle = Helpers.arc_cosine(
                            (ab_len * ab_len + bc_len * bc_len - ca_len * ca_len)
                            / (2.0 * ab_len * bc_len));
                    }

                    if (bc_len == 0.0 || ca_len == 0.0)
                    {
                        c_angle = Math.PI;
                    }
                    else
                    {
                        c_angle = Helpers.arc_cosine(
                            (bc_len * bc_len + ca_len * ca_len - ab_len * ab_len)
                            / (2.0 * bc_len * ca_len));
                    }
                }

                alpha_min = Math.Min(alpha_min, a_angle);
                alpha_min = Math.Min(alpha_min, b_angle);
                alpha_min = Math.Min(alpha_min, c_angle);

                alpha_ave = alpha_ave + alpha_min;

                alpha_area = alpha_area + area * alpha_min;

                area_total = area_total + area;
            }

            alpha_ave = alpha_ave / (double) (triangle_num);
            alpha_area = alpha_area / area_total;
            //
            //  Normalize angles from [0,pi/3] radians into qualities in [0,1].
            //
            alpha_min = alpha_min * 3.0 / Math.PI;
            alpha_ave = alpha_ave * 3.0 / Math.PI;
            alpha_area = alpha_area * 3.0 / Math.PI;
        }

        public static void area_measure(int n, double[] z, int triangle_order, int triangle_num,
                int[] triangle_node, ref double area_min, ref double area_max, ref double area_ratio,
                ref double area_ave, ref double area_std)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AREA_MEASURE determines the area ratio quality measure.
            //
            //  Discusion:
            //
            //    This measure computes the area of every triangle, and returns
            //    the ratio of the minimum to the maximum triangle.  A value of
            //    1 is "perfect", indicating that all triangles have the same area.
            //    A value of 0 is the worst possible result.
            //
            //    The code has been modified to 'allow' 6-node triangulations.
            //    However, no effort is made to actually process the midside nodes.
            //    Only information from the vertices is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double Z[2*N], the points.
            //
            //    Input, int TRIANGLE_ORDER, the order of the triangles.
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
            //    the triangulation.
            //
            //    Output, double *AREA_MIN, *AREA_MAX, the minimum and maximum
            //    areas.
            //
            //    Output, double *AREA_RATIO, the ratio of the minimum to the
            //    maximum area.
            //
            //    Output, double *AREA_AVE, the average area.
            //
            //    Output, double *AREA_STD, the standard deviation of the areas.
            //
        {
            double area;
            int triangle;
            double x1;
            double x2;
            double x3;
            double y1;
            double y2;
            double y3;

            area_max = 0.0;
            area_min = typeMethods.r8_huge();
            area_ave = 0.0;

            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                x1 = z[0 + (triangle_node[0 + triangle * triangle_order] - 1) * 2];
                y1 = z[1 + (triangle_node[0 + triangle * triangle_order] - 1) * 2];
                x2 = z[0 + (triangle_node[1 + triangle * triangle_order] - 1) * 2];
                y2 = z[1 + (triangle_node[1 + triangle * triangle_order] - 1) * 2];
                x3 = z[0 + (triangle_node[2 + triangle * triangle_order] - 1) * 2];
                y3 = z[1 + (triangle_node[2 + triangle * triangle_order] - 1) * 2];

                area = 0.5 * Math.Abs(x1 * (y2 - y3)
                                      + x2 * (y3 - y1)
                                      + x3 * (y1 - y2));

                area_min = Math.Min(area_min, area);
                area_max = Math.Max(area_max, area);
                area_ave = area_ave + area;
            }

            area_ave = area_ave / (double) (triangle_num);
            area_std = 0.0;
            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                x1 = z[0 + (triangle_node[0 + triangle * triangle_order] - 1) * 2];
                y1 = z[1 + (triangle_node[0 + triangle * triangle_order] - 1) * 2];
                x2 = z[0 + (triangle_node[1 + triangle * triangle_order] - 1) * 2];
                y2 = z[1 + (triangle_node[1 + triangle * triangle_order] - 1) * 2];
                x3 = z[0 + (triangle_node[2 + triangle * triangle_order] - 1) * 2];
                y3 = z[1 + (triangle_node[2 + triangle * triangle_order] - 1) * 2];

                area = 0.5 * Math.Abs(x1 * (y2 - y3)
                                      + x2 * (y3 - y1)
                                      + x3 * (y1 - y2));

                area_std = area_std + Math.Pow(area - area_ave, 2);
            }

            area_std = Math.Sqrt(area_std / (double) (triangle_num));

            if (0.0 < area_max)
            {
                area_ratio = area_min / area_max;
            }
            else
            {
                area_ratio = 0.0;
            }
        }

        public static void area_measure(int n, double[] z, int triangle_order, int triangle_num,
                int[] triangle_node, ref double area_min, ref double area_max, ref double area_ratio,
                ref double area_ave, ref double area_std, ref int area_negative, ref int area_zero)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AREA_MEASURE determines the area ratio quality measure.
            //
            //  Discusion:
            //
            //    This measure computes the area of every triangle, and returns
            //    the ratio of the minimum to the maximum triangle.  A value of
            //    1 is "perfect", indicating that all triangles have the same area.
            //    A value of 0 is the worst possible result.
            //
            //    The code has been modified to 'allow' 6-node triangulations.
            //    However, no effort is made to actually process the midside nodes.
            //    Only information from the vertices is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double Z[2*N], the points.
            //
            //    Input, int TRIANGLE_ORDER, the order of the triangles.
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
            //    the triangulation.
            //
            //    Output, double &AREA_MIN, &AREA_MAX, the minimum and maximum 
            //    areas.
            //
            //    Output, double &AREA_RATIO, the ratio of the minimum to the
            //    maximum area.
            //
            //    Output, double &AREA_AVE, the average area.
            //
            //    Output, double &AREA_STD, the standard deviation of the areas.
            //
            //    Output, int &AREA_NEGATIVE, the number of triangles with negative area.
            //
            //    Output, int &AREA_ZERO, the number of triangles with zero area.
        {
            double area;
            int triangle;
            double x1;
            double x2;
            double x3;
            double y1;
            double y2;
            double y3;

            area_max = 0.0;
            area_min = typeMethods.r8_huge();
            area_ave = 0.0;
            area_negative = 0;
            area_zero = 0;

            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                x1 = z[((0 + (triangle_node[0 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y1 = z[((1 + (triangle_node[0 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                x2 = z[((0 + (triangle_node[1 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y2 = z[((1 + (triangle_node[1 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                x3 = z[((0 + (triangle_node[2 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y3 = z[((1 + (triangle_node[2 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];

                area = 0.5 * Math.Abs(x1 * (y2 - y3)
                                      + x2 * (y3 - y1)
                                      + x3 * (y1 - y2));

                area_min = Math.Min(area_min, Math.Abs(area));
                area_max = Math.Max(area_max, Math.Abs(area));
                area_ave = area_ave + Math.Abs(area);

                if (area < 0.0)
                {
                    area_negative = area_negative + 1;
                }

                if (area == 0.0)
                {
                    area_zero = area_zero + 1;
                }
            }

            area_ave = area_ave / (double) (triangle_num);
            area_std = 0.0;
            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                x1 = z[((0 + (triangle_node[0 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y1 = z[((1 + (triangle_node[0 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                x2 = z[((0 + (triangle_node[1 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y2 = z[((1 + (triangle_node[1 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                x3 = z[((0 + (triangle_node[2 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y3 = z[((1 + (triangle_node[2 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];

                area = 0.5 * Math.Abs(x1 * (y2 - y3)
                                      + x2 * (y3 - y1)
                                      + x3 * (y1 - y2));

                area_std = area_std + Math.Pow(area - area_ave, 2);
            }

            area_std = Math.Sqrt(area_std / (double) (triangle_num));

            if (0.0 < area_max)
            {
                area_ratio = area_min / area_max;
            }
            else
            {
                area_ratio = 0.0;
            }
        }

        public static void q_measure(int n, double[] z, int triangle_order, int triangle_num,
        int[] triangle_node, ref double q_min, ref double q_max, ref double q_ave,
        ref double q_area )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q_MEASURE determines the triangulated pointset quality measure Q.
        //
        //  Discussion:
        //
        //    The Q measure evaluates the uniformity of the shapes of the triangles
        //    defined by a triangulated pointset.
        //
        //    For a single triangle T, the value of Q(T) is defined as follows:
        //
        //      TAU_IN = radius of the inscribed circle,
        //      TAU_OUT = radius of the circumscribed circle,
        //
        //      Q(T) = 2 * TAU_IN / TAU_OUT
        //        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
        //
        //    where A, B and C are the lengths of the sides of the triangle T.
        //
        //    The Q measure computes the value of Q(T) for every triangle T in the
        //    triangulation, and then computes the minimum of this
        //    set of values:
        //
        //      Q_MEASURE = min ( all T in triangulation ) Q(T)
        //
        //    In an ideally regular mesh, all triangles would have the same
        //    equilateral shape, for which Q = 1.  A good mesh would have
        //    0.5 < Q.
        //
        //    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
        //    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
        //    triangles.  Generally, a maximal triangulation is expected, namely,
        //    a triangulation whose image is a planar graph, but for which the
        //    addition of any new triangle would mean the graph was no longer planar.
        //    A Delaunay triangulation is a maximal triangulation which maximizes
        //    the minimum angle that occurs in any triangle.
        //
        //    The code has been modified to 'allow' 6-node triangulations.
        //    However, no effort is made to actually process the midside nodes.
        //    Only information from the vertices is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //    Per-Olof Persson and Gilbert Strang,
        //    A Simple Mesh Generator in MATLAB,
        //    SIAM Review,
        //    Volume 46, Number 2, pages 329-345, June 2004.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double Z[2*N], the points.
        //
        //    Input, int TRIANGLE_ORDER, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the triangulation.
        //
        //    Output, double *Q_MIN, *Q_MAX, the minimum and maximum values
        //    of Q over all triangles.
        //
        //    Output, double *Q_AVE, the average value of Q.
        //
        //    Output, double *Q_AREA, the average value of Q, weighted by
        //    the area of each triangle.
        //
        {
            int a_index;
            double ab_length;
            double area;
            double area_total;
            int b_index;
            double bc_length;
            int c_index;
            double ca_length;
            double q;
            int triangle;
            double x1;
            double x2;
            double x3;
            double y1;
            double y2;
            double y3;

            q_min = typeMethods.r8_huge();
            q_max = -typeMethods.r8_huge();
            q_ave = 0.0;
            q_area = 0.0;
            area_total = 0.0;

            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                a_index = triangle_node[0 + triangle * triangle_order];
                b_index = triangle_node[1 + triangle * triangle_order];
                c_index = triangle_node[2 + triangle * triangle_order];

                ab_length = Math.Sqrt(
                    Math.Pow(z[((0 + (a_index - 1) * 2) + z.Length) % z.Length] - z[((0 + (b_index - 1) * 2) + z.Length) % z.Length], 2)
                    + Math.Pow(z[((1 + (a_index - 1) * 2) + z.Length) % z.Length] - z[((1 + (b_index - 1) * 2) + z.Length) % z.Length], 2));

                bc_length = Math.Sqrt(
                    Math.Pow(z[((0 + (b_index - 1) * 2) + z.Length) % z.Length] - z[((0 + (c_index - 1) * 2) + z.Length) % z.Length], 2)
                    + Math.Pow(z[((1 + (b_index - 1) * 2) + z.Length) % z.Length] - z[((1 + (c_index - 1) * 2) + z.Length) % z.Length], 2));

                ca_length = Math.Sqrt(
                    Math.Pow(z[((0 + (c_index - 1) * 2) + z.Length) % z.Length] - z[((0 + (a_index - 1) * 2) + z.Length) % z.Length], 2)
                    + Math.Pow(z[((1 + (c_index - 1) * 2) + z.Length) % z.Length] - z[((1 + (a_index - 1) * 2) + z.Length) % z.Length], 2));

                q = (bc_length + ca_length - ab_length)
                    * (ca_length + ab_length - bc_length)
                    * (ab_length + bc_length - ca_length)
                    / (ab_length * bc_length * ca_length);

                x1 = z[((0 + (triangle_node[0 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y1 = z[((1 + (triangle_node[0 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                x2 = z[((0 + (triangle_node[1 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y2 = z[((1 + (triangle_node[1 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                x3 = z[((0 + (triangle_node[2 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];
                y3 = z[((1 + (triangle_node[2 + triangle * triangle_order] - 1) * 2) + z.Length) % z.Length];

                area = 0.5 * Math.Abs(x1 * (y2 - y3)
                                    + x2 * (y3 - y1)
                                    + x3 * (y1 - y2));

                q_min = Math.Min(q_min, q);
                q_max = Math.Max(q_max, q);
                q_ave = q_ave + q;
                q_area = q_area + q * area;

                area_total = area_total + area;
            }

            q_ave = q_ave / (double) (triangle_num);

            if (0.0 < area_total)
            {
                q_area = q_area / area_total;
            }
            else
            {
                q_area = 0.0;
            }
        }

    }
}