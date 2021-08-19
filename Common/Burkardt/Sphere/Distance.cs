using System;
using Burkardt.Types;

namespace Burkardt.SphereNS
{
    public static class Distance
    {
        public static double sphere_distance_xyz(double[] xyz1, double[] xyz2, int index1 = 0, int index2 = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_DISTANCE_XYZ computes great circle distances on a sphere.
            //
            //  Discussion:
            //
            //    XYZ coordinates are used.
            //
            //    We assume the points XYZ1 and XYZ2 lie on the same sphere.
            //
            //    This computation is a special form of the Vincenty formula.
            //    It should be less sensitive to errors associated with very small 
            //    or very large angular separations.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    "Great-circle distance",
            //    Wikipedia.
            //
            //  Parameters:
            //
            //    Input, double XYZ1[3], the coordinates of the first point.
            //
            //    Input, double XYZ2[3], the coordinates of the second point.
            //
            //    Output, double DIST, the great circle distance between the points.
            //
        {
            double bot;
            double dist;
            double lat1;
            double lat2;
            double lon1;
            double lon2;
            double r;
            double top;

            r = typeMethods.r8vec_norm(3, xyz1, index1);

            lat1 = Helpers.arc_sine(xyz1[index1 + 2]);
            lon1 = Helpers.atan4(xyz1[index1 + 1], xyz1[index1 + 0]);

            lat2 = Helpers.arc_sine(xyz2[index2 + 2]);
            lon2 = Helpers.atan4(xyz2[index2 + 1], xyz2[index2 + 0]);

            top = Math.Pow(Math.Cos(lat2) * Math.Sin(lon1 - lon2), 2)
                  + Math.Pow(Math.Cos(lat1) * Math.Sin(lat2)
                             - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2), 2);

            top = Math.Sqrt(top);

            bot = Math.Sin(lat1) * Math.Sin(lat2)
                  + Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2);

            dist = r * Math.Atan2(top, bot);

            return dist;
        }
    }
}