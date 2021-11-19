using System;
using Burkardt.Types;

namespace Burkardt.SphereNS;

public static class Distance
{
    public static double sphere01_distance_xyz ( double[] xyz1, double[] xyz2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_DISTANCE_XYZ computes great circle distances on a unit sphere.
        //
        //  Discussion:
        //
        //    XYZ coordinates are used.
        //
        //    We assume the points XYZ1 and XYZ2 lie on the unit sphere.
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
        //    23 September 2010
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
        double lat1 = Helpers.arc_sine ( xyz1[2] );
        double lon1 = Helpers.atan4 ( xyz1[1], xyz1[0] );

        double lat2 = Helpers.arc_sine ( xyz2[2] );
        double lon2 = Helpers.atan4 ( xyz2[1], xyz2[0] );

        double top = Math.Pow ( Math.Cos ( lat2 ) * Math.Sin ( lon1 - lon2 ), 2 )
                     + Math.Pow ( Math.Cos ( lat1 ) * Math.Sin ( lat2 ) 
                                  - Math.Sin ( lat1 ) * Math.Cos ( lat2 ) * Math.Cos ( lon1 - lon2 ), 2 );

        top = Math.Sqrt ( top );

        double bot = Math.Sin ( lat1 ) * Math.Sin ( lat2 ) 
                     + Math.Cos ( lat1 ) * Math.Cos ( lat2 ) * Math.Cos ( lon1 - lon2 );

        double dist = Math.Atan2 ( top, bot );

        return dist;
    }
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
        double r = typeMethods.r8vec_norm(3, xyz1, index1);

        double lat1 = Helpers.arc_sine(xyz1[index1 + 2]);
        double lon1 = Helpers.atan4(xyz1[index1 + 1], xyz1[index1 + 0]);

        double lat2 = Helpers.arc_sine(xyz2[index2 + 2]);
        double lon2 = Helpers.atan4(xyz2[index2 + 1], xyz2[index2 + 0]);

        double top = Math.Pow(Math.Cos(lat2) * Math.Sin(lon1 - lon2), 2)
                     + Math.Pow(Math.Cos(lat1) * Math.Sin(lat2)
                                - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2), 2);

        top = Math.Sqrt(top);

        double bot = Math.Sin(lat1) * Math.Sin(lat2)
                     + Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2);

        double dist = r * Math.Atan2(top, bot);

        return dist;
    }
}