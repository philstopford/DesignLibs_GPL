using System;
using Burkardt.Types;

namespace Burkardt.Geometry;

public static class Radec
{
    public static double radec_distance_3d(double ra1, double dec1, double ra2, double dec2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RADEC_DISTANCE_3D - angular distance, astronomical units, sphere in 3D.
        //
        //  Discussion:
        //
        //    Right ascension is measured in hours, between 0 and 24, and
        //    essentially measures longitude.
        //
        //    Declination measures the angle from the equator towards the north pole,
        //    and ranges from -90 (South Pole) to 90 (North Pole).
        //
        //    On the unit sphere, the angular separation between two points is
        //    equal to their geodesic or great circle distance.  On any other
        //    sphere, multiply the angular separation by the radius of the
        //    sphere to get the geodesic or great circle distance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double RA1, DEC1, RA2, DEC2, the right ascension and declination
        //    of the two points.
        //
        //    Output, double RADEC_DISTANCE_3D, the angular separation between the points,
        //    in radians.
        //
    {
        const int DIM_NUM = 3;

        int i;
        double[] v1 = new double[DIM_NUM];
        double[] v2 = new double[DIM_NUM];

        double theta1 = Helpers.degrees_to_radians(15.0 * ra1);
        double phi1 = Helpers.degrees_to_radians(dec1);

        v1[0] = Math.Cos(theta1) * Math.Cos(phi1);
        v1[1] = Math.Sin(theta1) * Math.Cos(phi1);
        v1[2] = Math.Sin(phi1);

        double norm_v1 = typeMethods.r8vec_norm(DIM_NUM, v1);

        double theta2 = Helpers.degrees_to_radians(15.0 * ra2);
        double phi2 = Helpers.degrees_to_radians(dec2);

        v2[0] = Math.Cos(theta2) * Math.Cos(phi2);
        v2[1] = Math.Sin(theta2) * Math.Cos(phi2);
        v2[2] = Math.Sin(phi2);

        double norm_v2 = typeMethods.r8vec_norm(DIM_NUM, v2);

        double cos_theta = 0.0;
        for (i = 0; i < 3; i++)
        {
            cos_theta += v2[i] * v2[i];
        }

        cos_theta = Math.Sqrt(cos_theta);

        cos_theta /= norm_v1 * norm_v2;

        double value = typeMethods.r8_acos(cos_theta);

        return value;
    }

    public static double[] radec_to_xyz(double ra, double dec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RADEC_TO_XYZ converts right ascension/declination to (X,Y,Z) coordinates.
        //
        //  Discussion:
        //
        //    Right ascension is measured in hours, between 0 and 24, and
        //    essentially measures longitude.
        //
        //    Declination measures the angle from the equator towards the north pole,
        //    and ranges from -90 (South Pole) to 90 (North Pole).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double RA, DEC, the right ascension and declination of a point.
        //
        //    Output, double RADEC_TO_XYZ[3], the corresponding coordinates of a
        //    point with radius 1.
        //
    {
        const int DIM_NUM = 3;

        double theta = Helpers.degrees_to_radians(15.0 * ra);
        double phi = Helpers.degrees_to_radians(dec);

        double[] p = new double[DIM_NUM];

        p[0] = Math.Cos(theta) * Math.Cos(phi);
        p[1] = Math.Sin(theta) * Math.Cos(phi);
        p[2] = Math.Sin(phi);

        return p;
    }

}