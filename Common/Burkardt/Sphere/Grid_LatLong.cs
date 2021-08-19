using System;

namespace Burkardt.SphereNS
{
    public static class Grid_LatLong
    {
        public static int[] sphere_ll_lines(int nlat, int nlong, int line_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_LL_LINES produces lines from a latitude/longitude grid.
            //
            //  Discussion:
            //
            //    The point numbering system is the same used in SPHERE_LL_POINTS,
            //    and that routine may be used to compute the coordinates of the points.
            //
            //    The implicit form of a sphere in 3D is:
            //
            //        pow ( P[0] - PC[0], 2 ) 
            //      + pow ( P[1] - PC[1], 2 ) 
            //      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int LINE_MAX, the maximum number of gridlines.
            //
            //    Input, int NLAT, NLONG, the number of latitude and longitude
            //    lines to draw.  The latitudes do not include the North and South
            //    poles, which will be included automatically, so NLAT = 5, for instance,
            //    will result in points along 7 lines of latitude.
            //
            //    Input, int LINE_NUM, the number of grid lines.
            //
            //    Output, int SPHERE_LL_LINE[2*LINE_NUM], contains pairs of point indices for
            //    line segments that make up the grid.
            //
        {
            int i;
            int j;
            int l;
            int[] line;
            int next;
            int newcol;
            int old;

            line = new int [2 * line_num];
            l = 0;
            //
            //  "Vertical" lines.
            //
            for (j = 0; j <= nlong - 1; j++)
            {
                old = 1;
                next = j + 2;
                line[0 + l * 2] = old;
                line[1 + l * 2] = next;
                l = l + 1;

                for (i = 1; i <= nlat - 1; i++)
                {
                    old = next;
                    next = old + nlong;
                    line[0 + l * 2] = old;
                    line[1 + l * 2] = next;
                    l = l + 1;
                }

                old = next;
                line[0 + l * 2] = old;
                line[1 + l * 2] = 1 + nlat * nlong + 1;
                l = l + 1;
            }

            //
            //  "Horizontal" lines.
            //
            for (i = 1; i <= nlat; i++)
            {
                next = 1 + (i - 1) * nlong + 1;

                for (j = 0; j <= nlong - 2; j++)
                {
                    old = next;
                    next = old + 1;
                    line[0 + l * 2] = old;
                    line[1 + l * 2] = next;
                    l = l + 1;
                }

                old = next;
                next = 1 + (i - 1) * nlong + 1;
                line[0 + l * 2] = old;
                line[1 + l * 2] = next;
                l = l + 1;
            }

            //
            //  "Diagonal" lines.
            //
            for (j = 0; j <= nlong - 1; j++)
            {
                old = 1;
                next = j + 2;
                newcol = j;

                for (i = 1; i <= nlat - 1; i++)
                {
                    old = next;
                    next = old + nlong + 1;

                    newcol = newcol + 1;
                    if (nlong - 1 < newcol)
                    {
                        newcol = 0;
                        next = next - nlong;
                    }

                    line[0 + l * 2] = old;
                    line[1 + l * 2] = next;
                    l = l + 1;
                }
            }

            return line;
        }

        public static int sphere_ll_line_num(int lat_num, int long_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_LL_LINE_NUM counts lines for a latitude/longitude grid.
            //
            //  Discussion:
            //
            //    The number returned is the number of pairs of points to be connected.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int LAT_NUM, LONG_NUM, the number of latitude and
            //    longitude lines to draw.  The latitudes do not include the North and South
            //    poles, which will be included automatically, so LAT_NUM = 5, for instance,
            //    will result in points along 7 lines of latitude.
            //
            //    Output, int SPHERE_LL_LINE_NUM, the number of grid lines.
            //
        {
            int line_num;

            line_num = long_num * (lat_num + 1)
                       + lat_num * long_num
                       + long_num * (lat_num - 1);

            return line_num;
        }

        public static double[] sphere_ll_points(double r, double[] pc, int lat_num, int lon_num,
                int point_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_LL_POINTS produces points on a latitude/longitude grid.
            //
            //  Discussion:
            //
            //    The implicit form of a sphere in 3D is:
            //
            //        pow ( P[0] - PC[0], 2 ) 
            //      + pow ( P[1] - PC[1], 2 ) 
            //      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the sphere.
            //
            //    Input, double PC[3], the coordinates of the center of the sphere.
            //
            //    Input, int LAT_NUM, LON_NUM, the number of latitude and longitude
            //    lines to draw.  The latitudes do not include the North and South
            //    poles, which will be included automatically, so LAT_NUM = 5, for instance,
            //    will result in points along 7 lines of latitude.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Output, double SPHERE_LL_POINTS[3*POINT_NUM], the coordinates 
            //    of the grid points.
            //
        {
            int lat;
            int lon;
            int n;
            double[] p;
            double phi;
            const double pi = 3.141592653589793;
            double theta;

            p = new double[3 * point_num];
            n = 0;
            //
            //  The north pole.
            //
            theta = 0.0;
            phi = 0.0;

            p[0 + n * 3] = pc[0] + r * Math.Sin(phi) * Math.Cos(theta);
            p[1 + n * 3] = pc[1] + r * Math.Sin(phi) * Math.Sin(theta);
            p[2 + n * 3] = pc[2] + r * Math.Cos(phi);
            n = n + 1;
            //
            //  Do each intermediate ring of latitude.
            //
            for (lat = 1; lat <= lat_num; lat++)
            {
                phi = (double)(lat) * pi / (double)(lat_num + 1);
                //
                //  Along that ring of latitude, compute points at various longitudes.
                //
                for (lon = 0; lon < lon_num; lon++)
                {
                    theta = (double)(lon) * 2.0 * pi / (double)(lon_num);

                    p[0 + n * 3] = pc[0] + r * Math.Sin(phi) * Math.Cos(theta);
                    p[1 + n * 3] = pc[1] + r * Math.Sin(phi) * Math.Sin(theta);
                    p[2 + n * 3] = pc[2] + r * Math.Cos(phi);
                    n = n + 1;
                }
            }

            //
            //  The south pole.
            //
            theta = 0.0;
            phi = pi;
            p[0 + n * 3] = pc[0] + r * Math.Sin(phi) * Math.Cos(theta);
            p[1 + n * 3] = pc[1] + r * Math.Sin(phi) * Math.Sin(theta);
            p[2 + n * 3] = pc[2] + r * Math.Cos(phi);
            n = n + 1;

            return p;
        }

        public static int sphere_ll_point_num(int lat_num, int long_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_LL_POINT_NUM counts points for a latitude/longitude grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int LAT_NUM, LONG_NUM, the number of latitude 
            //    and longitude lines to draw.  The latitudes do not include the North and 
            //    South poles, which will be included automatically, so LAT_NUM = 5, for 
            //    instance, will result in points along 7 lines of latitude.
            //
            //    Output, int SPHERE_LL_POINT_NUM, the number of grid points.
            //
        {
            int point_num;

            point_num = 2 + lat_num * long_num;

            return point_num;
        }

        public static int[] sphere_llq_lines(int nlat, int nlong, int line_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_LLQ_LINES: latitude/longitude quadrilateral grid lines.
            //
            //  Discussion:
            //
            //    The point numbering system is the same used in SPHERE_LL_POINTS,
            //    and that routine may be used to compute the coordinates of the points.
            //
            //    The implicit form of a sphere in 3D is:
            //
            //        pow ( P[0] - PC[0], 2 ) 
            //      + pow ( P[1] - PC[1], 2 ) 
            //      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int LINE_MAX, the maximum number of gridlines.
            //
            //    Input, int NLAT, NLONG, the number of latitude and longitude
            //    lines to draw.  The latitudes do not include the North and South
            //    poles, which will be included automatically, so NLAT = 5, for instance,
            //    will result in points along 7 lines of latitude.
            //
            //    Input, int LINE_NUM, the number of grid lines.
            //
            //    Output, int SPHERE_LLQ_LINE[2*LINE_NUM], contains pairs of point indices for
            //    line segments that make up the grid.
            //
        {
            int i;
            int j;
            int l;
            int[] line;
            int next;
            int old;

            line = new int [2 * line_num];
            l = 0;
            //
            //  "Vertical" lines.
            //
            for (j = 0; j <= nlong - 1; j++)
            {
                old = 1;
                next = j + 2;
                line[0 + l * 2] = old;
                line[1 + l * 2] = next;
                l = l + 1;

                for (i = 1; i <= nlat - 1; i++)
                {
                    old = next;
                    next = old + nlong;
                    line[0 + l * 2] = old;
                    line[1 + l * 2] = next;
                    l = l + 1;
                }

                old = next;
                line[0 + l * 2] = old;
                line[1 + l * 2] = 1 + nlat * nlong + 1;
                l = l + 1;
            }

            //
            //  "Horizontal" lines.
            //
            for (i = 1; i <= nlat; i++)
            {
                next = 1 + (i - 1) * nlong + 1;

                for (j = 0; j <= nlong - 2; j++)
                {
                    old = next;
                    next = old + 1;
                    line[0 + l * 2] = old;
                    line[1 + l * 2] = next;
                    l = l + 1;
                }

                old = next;
                next = 1 + (i - 1) * nlong + 1;
                line[0 + l * 2] = old;
                line[1 + l * 2] = next;
                l = l + 1;
            }

            return line;
        }

        public static int sphere_llq_line_num(int lat_num, int long_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_LLQ_LINE_NUM counts lines for a latitude/longitude quadrilateral grid.
            //
            //  Discussion:
            //
            //    The number returned is the number of pairs of points to be connected.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int LAT_NUM, LONG_NUM, the number of latitude and
            //    longitude lines to draw.  The latitudes do not include the North and South
            //    poles, which will be included automatically, so LAT_NUM = 5, for instance,
            //    will result in points along 7 lines of latitude.
            //
            //    Output, int SPHERE_LLQ_LINE_NUM, the number of grid lines.
            //
        {
            int line_num;

            line_num = long_num * (lat_num + 1)
                       + lat_num * long_num;

            return line_num;
        }
    }
}