using System;
using Burkardt.Types;

namespace Burkardt.Geometry;

public static class XY
{
    public static void xy_to_polar(double[] xy, ref double r, ref double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_TO_POLAR converts XY coordinates to polar coordinates.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double XY[2], the Cartesian coordinates.
        //
        //    Output, double *R, *T, the radius and angle (in radians).
        //
    {
        r = Math.Sqrt(xy[0] * xy[0] + xy[1] * xy[1]);

        t = r switch
        {
            0.0 => 0.0,
            _ => Math.Atan2(xy[0], xy[1])
        };
    }

    public static void xyz_to_radec(double[] p, ref double ra, ref double dec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYZ_TO_RADEC converts (X,Y,Z) to right ascension/declination coordinates.
        //
        //  Discussion:
        //
        //    Given an XYZ point, compute its distance R from the origin, and
        //    regard it as lying on a sphere of radius R, whose axis is the Z
        //    axis.
        //
        //    The right ascension of the point is the "longitude", measured in hours,
        //    between 0 and 24, with the X axis having right ascension 0, and the
        //    Y axis having right ascension 6.
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
        //    28 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P[3], the coordinates of a point in 3D.
        //
        //    Output, double *RA, *DEC, the corresponding right ascension
        //    and declination.
        //
    {
        const int DIM_NUM = 3;

        double norm_v = typeMethods.r8vec_norm(DIM_NUM, p);

        double phi = Math.Asin(p[2] / norm_v);

        double theta = Math.Cos(phi) switch
        {
            0.0 => 0.0,
            _ => typeMethods.r8_atan(p[1], p[0])
        };

        dec = Helpers.radians_to_degrees(phi);
        ra = Helpers.radians_to_degrees(theta) / 15.0;

    }

    public static void xyz_to_rtp(double[] xyz, ref double r, ref double theta, ref double phi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYZ_TO_RTP converts (X,Y,Z) to (R,Theta,Phi) coordinates.
        //
        //  Discussion:
        //
        //    Given an XYZ point, compute its distance R from the origin, and
        //    regard it as lying on a sphere of radius R, whose axis is the Z
        //    axis.
        //
        //    Theta measures the "longitude" of the point, between 0 and 2 PI.
        //
        //    PHI measures the angle from the "north pole", between 0 and PI.
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
        //  Parameters:
        //
        //    Input, double XYZ[3], the coordinates of a point in 3D.
        //
        //    Output, double *R, *THETA, *PHI, the radius, longitude and
        //    declination of the point.
        //
    {
        r = Math.Sqrt(Math.Pow(xyz[0], 2)
                      + Math.Pow(xyz[1], 2)
                      + Math.Pow(xyz[2], 2));

        switch (r)
        {
            case 0.0:
                theta = 0.0;
                phi = 0.0;
                return;
            default:
                phi = typeMethods.r8_acos(xyz[2] / r);

                theta = typeMethods.r8_atan(xyz[1], xyz[0]);
                break;
        }
    }

    public static void xyz_to_tp(double[] xyz, ref double theta, ref double phi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates.
        //
        //  Discussion:
        //
        //    Given an XYZ point, regard it as lying on a sphere of radius R,
        //    centered at the origin, whose axis is the Z axis.
        //
        //    We assume that the actual value of R is of no interest, and do
        //    not report it.  This is especially appropriate if the point is
        //    expected to lie on the unit sphere, for instance.
        //
        //    THETA measures the "longitude" of the point, between 0 and 2 PI.
        //
        //    PHI measures the angle from the "north pole", between 0 and PI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double XYZ[3], the coordinates of a point in 3D.
        //
        //    Output, double *THETA, *PHI, the longitude and declination of the point.
        //
    {
        double r = Math.Sqrt(Math.Pow(xyz[0], 2)
                             + Math.Pow(xyz[1], 2)
                             + Math.Pow(xyz[2], 2));

        switch (r)
        {
            case 0.0:
                theta = 0.0;
                phi = 0.0;
                return;
            default:
                phi = typeMethods.r8_acos(xyz[2] / r);

                theta = typeMethods.r8_atan(xyz[1], xyz[0]);
                break;
        }
    }


    public static void rtp_to_xyz(double r, double theta, double phi, ref double[] xyz)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RTP_TO_XYZ converts (R,Theta,Phi) to (X,Y,Z) coordinates.
        //
        //  Discussion:
        //
        //    R measures the distance of the point to the origin.
        //
        //    Theta measures the "longitude" of the point, between 0 and 2 PI.
        //
        //    PHI measures the angle from the "north pole", between 0 and PI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, THETA, PHI, the radius, longitude, and
        //    declination of a point.
        //
        //    Output, double XYZ[3], the corresponding Cartesian coordinates.
        //
    {
        xyz[0] = r * Math.Cos(theta) * Math.Sin(phi);
        xyz[1] = r * Math.Sin(theta) * Math.Sin(phi);
        xyz[2] = r * Math.Cos(phi);

    }

    public static double[] tp_to_xyz(double theta, double phi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TP_TO_XYZ converts unit spherical TP coordinates to XYZ coordinates.
        //
        //  Discussion:
        //
        //    The point is assume to lie on the unit sphere centered at the origin.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double THETA, PHI, the angular coordinates of a point
        //    on the unit sphere.
        //
        //    Output, double TP_TO_XYZ[3], the XYZ coordinates.
        //
    {
        double[] v = new double[3];

        v[0] = Math.Cos(theta) * Math.Sin(phi);
        v[1] = Math.Sin(theta) * Math.Sin(phi);
        v[2] = Math.Cos(phi);

        return v;
    }


}