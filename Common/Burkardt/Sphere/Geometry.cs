using System;
using Burkardt.Probability;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SphereNS;

public static class Geometry
{
    public static double sphere_cap_area_2d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CAP_AREA_2D computes the surface area of a spherical cap in 2D.
        //
        //  Discussion:
        //
        //    Draw any radius of the sphere and note the point P where the radius
        //    intersects the sphere.  Consider the point on the radius line which is
        //    H units from P.  Draw the circle that lies in the plane perpendicular to
        //    the radius, and which intersects the sphere.  The circle divides the sphere
        //    into two pieces, and the corresponding disk divides the solid sphere into
        //    two pieces.  The spherical cap is the part of the solid sphere that
        //    includes the point P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H, the "height" of the spherical cap.
        //    H must be between 0 and 2 * R.
        //
        //    Output, double SPHERE_CAP_AREA_2D, the area of the spherical cap.
        //
    {
        double area;
        double theta;

        switch (h)
        {
            case <= 0.0:
                area = 0.0;
                break;
            default:
            {
                if (2.0 * r <= h)
                {
                    area = 2.0 * Math.PI * r;
                }
                else
                {
                    theta = 2.0 * typeMethods.r8_asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);
                    area = r * theta;
                    if (r <= h)
                    {
                        area = 2.0 * Math.PI * r - area;
                    }
                }

                break;
            }
        }

        return area;
    }

    public static double sphere_cap_area_3d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CAP_AREA_3D computes the surface area of a spherical cap in 3D.
        //
        //  Discussion:
        //
        //    Draw any radius of the sphere and note the point P where the radius
        //    intersects the sphere.  Consider the point on the radius line which is
        //    H units from P.  Draw the circle that lies in the plane perpendicular to
        //    the radius, and which intersects the sphere.  The circle divides the sphere
        //    into two pieces, and the corresponding disk divides the solid sphere into
        //    two pieces.  The spherical cap is the part of the solid sphere that
        //    includes the point P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H, the "height" of the spherical cap.
        //    H must be between 0 and 2 * R.
        //
        //    Output, double SPHERE_CAP_AREA_3D, the area of the spherical cap.
        //
    {
        double area;

        switch (h)
        {
            case <= 0.0:
                area = 0.0;
                break;
            default:
            {
                if (2.0 * r <= h)
                {
                    area = 4.0 * Math.PI * r * r;
                }
                else
                {
                    area = 2.0 * Math.PI * r * h;
                }

                break;
            }
        }

        return area;
    }

    public static double sphere_cap_area_nd(int dim_num, double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CAP_AREA_ND computes the area of a spherical cap in ND.
        //
        //  Discussion:
        //
        //    The spherical cap is a portion of the surface of the sphere:
        //
        //      sum ( X(1:N)^2 ) = R^2
        //
        //    which is no more than H units from the uppermost point on the sphere.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Thomas Ericson, Victor Zinoviev,
        //    Codes on Euclidean Spheres,
        //    Elsevier, 2001, pages 439-441.
        //    QA166.7 E75
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H, the "thickness" of the spherical cap,
        //    which is normally between 0 and 2 * R.
        //
        //    Output, double SPHERE_CAP_AREA_ND, the area of the spherical cap.
        //
    {
        double area;
        double area2;
        double haver_sine;
        int i;
        double theta;
        double ti;
        double tj;
        double tk;

        switch (h)
        {
            case <= 0.0:
                area = 0.0;
                return area;
        }

        if (2.0 * r <= h)
        {
            area = sphere_imp_area_nd(dim_num, r);
            return area;
        }

        //
        //  For cases where R < H < 2 * R, work with the complementary region.
        //
        haver_sine = Math.Sqrt((2.0 * r - h) * h);

        theta = typeMethods.r8_asin(haver_sine / r);

        switch (dim_num)
        {
            case < 1:
                area = -1.0;
                break;
            case 1:
                area = 0.0;
                break;
            case 2:
                area = 2.0 * theta * r;
                break;
            default:
            {
                ti = theta;

                tj = ti;
                ti = 1.0 - Math.Cos(theta);

                for (i = 2; i <= dim_num - 2; i++)
                {
                    tk = tj;
                    tj = ti;
                    ti = ((i - 1) * tk
                          - Math.Cos(theta) * Math.Pow(Math.Sin(theta), i - 1))
                         / i;
                }

                area = sphere_k(dim_num - 1) * ti * Math.Pow(r, dim_num - 1);
                break;
            }
        }

        //
        //  Adjust for cases where R < H < 2R.
        //
        if (r < h)
        {
            area2 = sphere_imp_area_nd(dim_num, r);
            area = area2 - area;
        }

        return area;
    }

    public static double sphere_cap_volume_2d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CAP_VOLUME_2D computes the volume of a spherical cap in 2D.
        //
        //  Discussion:
        //
        //    Draw any radius R of the circle and denote as P the point where the
        //    radius intersects the circle.  Now consider the point Q which lies
        //    on the radius and which is H units from P.  The line which is
        //    perpendicular to the radius R and passes through Q divides the
        //    circle into two pieces.  The piece including the point P is the
        //    spherical (circular) cap of height (or thickness) H.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H, the "height" of the spherical cap.  H must
        //    be between 0 and 2 * R.
        //
        //    Output, double SPHERE_CAP_VOLUME_2D, the volume (area) of the spherical cap.
        //
    {
        double theta;
        double volume;

        switch (h)
        {
            case <= 0.0:
                volume = 0.0;
                break;
            default:
            {
                if (2.0 * r <= h)
                {
                    volume = Math.PI * r * r;
                }
                else
                {
                    theta = 2.0 * typeMethods.r8_asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);
                    volume = r * r * (theta - Math.Sin(theta)) / 2.0;

                    if (r < h)
                    {
                        volume = Math.PI * r * r - volume;
                    }
                }

                break;
            }
        }

        return volume;
    }

    public static double sphere_cap_volume_3d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CAP_VOLUME_3D computes the volume of a spherical cap in 3D.
        //
        //  Discussion:
        //
        //    Draw any radius of the sphere and note the point P where the radius
        //    intersects the sphere.  Consider the point on the radius line which is
        //    H units from P.  Draw the circle that lies in the plane perpendicular to
        //    the radius, and which intersects the sphere.  The circle divides the sphere
        //    into two pieces, and the corresponding disk divides the solid sphere into
        //    two pieces.  The spherical cap is the part of the solid sphere that
        //    includes the point P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H, the "height" of the spherical cap.  H must be between
        //    0 and 2 * R.
        //
        //    Output, double SPHERE_CAP_VOLUME_3D, the volume of the spherical cap.
        //
    {

        switch (h)
        {
            case <= 0.0:
                return 0.0;
        }

        if (2.0 * r <= h)
        {
            return 4.0 / 3.0 * Math.PI * r * r * r;
        }
        return 1.0 / 3.0 * Math.PI * h * h * (3.0 * r - h);
    }

    public static double sphere_cap_volume_nd(int dim_num, double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CAP_VOLUME_ND computes the volume of a spherical cap in ND.
        //
        //  Discussion:
        //
        //    The spherical cap is a portion of the surface and interior of the sphere:
        //
        //      sum ( X(1:N)^2 ) <= R^2
        //
        //    which is no more than H units from some point P on the sphere.
        //
        //
        //    The algorithm proceeds from the observation that the N-dimensional
        //    sphere can be parameterized by a quantity RC that runs along the
        //    radius from the center to the point P.  The value of RC at the
        //    base of the spherical cap is (R-H) and at P it is R.  We intend to
        //    use RC as our integration parameeter.
        //
        //    The volume of the spherical cap is then the integral, as RC goes
        //    from (R-H) to R, of the N-1 dimensional volume of the sphere
        //    of radius RS, where RC^2 + RS^2 = R^2.
        //
        //    The volume of the N-1 dimensional sphere of radius RS is simply
        //    some constants times RS^(N-1).
        //
        //    After factoring out the constant terms, and writing RC = R * Math.Cos ( T ),
        //    and RS = R * Math.Sin ( T ), and letting
        //      T_MAX = typeMethods.r8_asin ( Math.Sqrt ( ( 2.0 * r - h ) * h / r ) ),
        //    the "interesting part" of our integral becomes
        //
        //      constants * R^N * Integral ( T = 0 to T_MAX ) sin**N ( T ) dT
        //
        //    The integral of sin^N ( T ) dT can be handled by recursion.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H, the "thickness" of the spherical cap,
        //    which is normally between 0 and 2 * R.
        //
        //    Output, double SPHERE_CAP_VOLUME_ND, the volume of the spherical cap.
        //
    {
        double angle;
        double arg;
        double factor1;
        double factor2;
        double volume;
        double volume2;

        switch (h)
        {
            case <= 0.0:
                volume = 0.0;
                return volume;
        }

        if (2.0 * r <= h)
        {
            volume = sphere_imp_volume_nd(dim_num, r);
            return volume;
        }

        switch (dim_num)
        {
            case < 1:
                volume = -1.0;
                break;
            case 1:
                volume = h;
                break;
            default:
            {
                factor1 = sphere_unit_volume_nd(dim_num - 1);

                angle = typeMethods.r8_asin(Math.Sqrt((2.0 * r - h) * h / r));

                arg = 0.0;
                factor2 = Misc.sin_power_int(arg, angle, dim_num);

                volume = factor1 * factor2 * Math.Pow(r, dim_num);

                if (r < h)
                {
                    volume2 = sphere_imp_volume_nd(dim_num, r);
                    volume = volume2 - volume;
                }

                break;
            }
        }

        return volume;
    }

    public static void sphere_dia2imp_3d(double[] p1, double[] p2, ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_DIA2IMP_3D converts a diameter to an implicit sphere in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], are the coordinates
        //    of two points which form a diameter of the sphere.
        //
        //    Output, double *R, the computed radius of the sphere.
        //
        //    Output, double PC[3], the computed center of the sphere.
        //
    {
        r = 0.5 * Math.Sqrt(Math.Pow(p1[0] - p2[0], 2)
                            + Math.Pow(p1[1] - p2[1], 2)
                            + Math.Pow(p1[2] - p2[2], 2));

        pc[0] = 0.5 * (p1[0] + p2[0]);
        pc[1] = 0.5 * (p1[1] + p2[1]);
        pc[2] = 0.5 * (p1[2] + p2[2]);

    }

    public static double sphere_distance_xyz(double[] xyz1, double[] xyz2)

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

        r = typeMethods.r8vec_norm(3, xyz1);

        lat1 = typeMethods.r8_asin(xyz1[2]);
        lon1 = typeMethods.r8_atan(xyz1[1], xyz1[0]);

        lat2 = typeMethods.r8_asin(xyz2[2]);
        lon2 = typeMethods.r8_atan(xyz2[1], xyz2[0]);

        top = Math.Pow(Math.Cos(lat2) * Math.Sin(lon1 - lon2), 2)
              + Math.Pow(Math.Cos(lat1) * Math.Sin(lat2)
                         - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2), 2);

        top = Math.Sqrt(top);

        bot = Math.Sin(lat1) * Math.Sin(lat2)
              + Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2);

        dist = r * Math.Atan2(top, bot);

        return dist;
    }

    public static double sphere_distance1(double lat1, double lon1, double lat2, double lon2,
            double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_DISTANCE1 computes great circle distances on a sphere.
        //
        //  Discussion:
        //
        //    This computation is based on the law of cosines for spheres.
        //    This formula can suffer from rounding errors when the angular
        //    distances are small.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 February 2009
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
        //    Input, double LAT1, LON1, the latitude and longitude of
        //    the first point.
        //
        //    Input, double LAT2, LON2, the latitude and longitude of
        //    the second point.
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Output, double DIST, the great circle distance between
        //    the points, measured in the same units as R.
        //
    {
        double c;
        double dist;

        c = Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2)
            + Math.Sin(lat1) * Math.Sin(lat2);

        dist = r * Math.Acos(c);

        return dist;
    }

    public static double sphere_distance2(double lat1, double lon1, double lat2, double lon2,
            double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_DISTANCE2 computes great circle distances on a sphere.
        //
        //  Discussion:
        //
        //    This computation is written in terms of haversines, and can be more
        //    accurate when measuring small angular distances.  It can be somewhat
        //    inaccurate when the two points are antipodal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 February 2009
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
        //    Input, double LAT1, LON1, the latitude and longitude of
        //    the first point.
        //
        //    Input, double LAT2, LON2, the latitude and longitude of
        //    the second point.
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Output, double DIST, the great circle distance between
        //    the points, measured in the same units as R.
        //
    {
        double dist;
        double s;

        s = Math.Pow(Math.Sin((lat1 - lat2) / 2.0), 2)
            + Math.Cos(lat1) * Math.Cos(lat2) * Math.Pow(Math.Sin((lon1 - lon2) / 2.0), 2);
        s = Math.Sqrt(s);

        dist = 2.0 * r * Math.Asin(s);

        return dist;
    }

    public static double sphere_distance3(double lat1, double lon1, double lat2, double lon2,
            double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_DISTANCE3 computes great circle distances on a sphere.
        //
        //  Discussion:
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
        //    18 February 2009
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
        //    Input, double LAT1, LON1, the latitude and longitude of
        //    the first point.
        //
        //    Input, double LAT2, LON2, the latitude and longitude of
        //    the second point.
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Output, double DIST, the great circle distance between
        //    the points, measured in the same units as R.
        //
    {
        double bot;
        double dist;
        double top;

        top = Math.Pow(Math.Cos(lat2) * Math.Sin(lon1 - lon2), 2)
              + Math.Pow(Math.Cos(lat1) * Math.Sin(lat2)
                         - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2), 2);

        top = Math.Sqrt(top);

        bot = Math.Sin(lat1) * Math.Sin(lat2)
              + Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2);

        dist = r * Math.Atan2(top, bot);

        return dist;
    }

    public static bool sphere_exp_contains_point_3d(double[] p1, double[] p2, double[] p3,
            double[] p4, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_EXP_CONTAINS_POINT_3D determines if an explicit sphere contains a point in 3D.
        //
        //  Discussion:
        //
        //    An explicit sphere in 3D is determined by four points,
        //    which should be distinct, and not coplanar.
        //
        //    The computation checks the determinant of:
        //
        //      x1  y1  z1  x1^2+y1^2+z1^2  1
        //      x2  y2  z2  x2^2+y2^2+z2^2  1
        //      x3  y3  z3  x3^2+y3^2+z3^2  1
        //      x4  y4  z4  x4^2+y4^2+z4^2  1
        //      x   y   z   x^2 +y^2 +z^2   1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], P3[3], P4[3], the coordinates of four points
        //    that lie on a circle.
        //
        //    Input, double P[3], the coordinates of a point, whose
        //    position relative to the sphere is desired.
        //
        //    Output, bool SPHERE_EXP_CONTAINS_POINT_3D, is TRUE if the point
        //    is in the sphere, FALSE otherwise.
        //
    {
        double[] a = new double[5 * 5];
        //
        //  Set up the matrix.
        //
        a[0 + 0 * 5] = p1[0];
        a[1 + 0 * 5] = p2[0];
        a[2 + 0 * 5] = p3[0];
        a[3 + 0 * 5] = p4[0];
        a[4 + 0 * 5] = p[0];

        a[0 + 1 * 5] = p1[1];
        a[1 + 1 * 5] = p2[1];
        a[2 + 1 * 5] = p3[1];
        a[3 + 1 * 5] = p4[1];
        a[4 + 1 * 5] = p[1];

        a[0 + 2 * 5] = p1[2];
        a[1 + 2 * 5] = p2[2];
        a[2 + 2 * 5] = p3[2];
        a[3 + 2 * 5] = p4[2];
        a[4 + 2 * 5] = p[2];

        a[0 + 3 * 5] = p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2];
        a[1 + 3 * 5] = p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2];
        a[2 + 3 * 5] = p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2];
        a[3 + 3 * 5] = p4[0] * p4[0] + p4[1] * p4[1] + p4[2] * p4[2];
        a[4 + 3 * 5] = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];

        a[0 + 4 * 5] = 1.0;
        a[1 + 4 * 5] = 1.0;
        a[2 + 4 * 5] = 1.0;
        a[3 + 4 * 5] = 1.0;
        a[4 + 4 * 5] = 1.0;

        if (typeMethods.r8mat_det_5d(a) < 0.0)
        {
            return false;
        }

        return true;
    }

    public static void sphere_exp_point_near_3d(double[] p1, double[] p2, double[] p3,
            double[] p4, double[] p, ref double[] pn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_EXP_POINT_NEAR_3D finds the nearest point on an explicit sphere to a point in 3D.
        //
        //  Discussion:
        //
        //    An explicit sphere in 3D is determined by four points,
        //    which should be distinct, and not coplanar.
        //
        //    If the center of the sphere is PC, and the point is P, then
        //    the desired point lies at a positive distance R along the vector
        //    P-PC unless P = PC, in which case any
        //    point on the sphere is "nearest".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], P3[3], P4[3], the coordinates of four points
        //    that lie on a sphere.
        //
        //    Input, double P[3], the coordinates of a point whose nearest point on the
        //    sphere is desired.
        //
        //    Output, double PN[3], the nearest point on the sphere.
        //
    {
        double norm;
        double r = 0;
        double[] pc = new double[3];
        //
        //  Find the center.
        //
        sphere_exp2imp_3d(p1, p2, p3, p4, ref r, ref pc);
        //
        //  If P = PC, bail out now.
        //
        norm = Math.Sqrt(Math.Pow(p[0] - pc[0], 2)
                         + Math.Pow(p[1] - pc[1], 2)
                         + Math.Pow(p[2] - pc[2], 2));

        switch (norm)
        {
            case 0.0:
                pn[0] = pc[0] + r;
                pn[1] = pc[1];
                pn[2] = pc[2];
                return;
        }

        //
        //  Compute the nearest point.
        //
        pn[0] = pc[0] + r * (p[0] - pc[0]) / norm;
        pn[1] = pc[1] + r * (p[1] - pc[1]) / norm;
        pn[2] = pc[2] + r * (p[2] - pc[2]) / norm;
    }

    public static void sphere_exp2imp_3d(double[] p1, double[] p2, double[] p3, double[] p4,
            ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_EXP2IMP_3D converts a sphere from explicit to implicit form in 3D.
        //
        //  Discussion:
        //
        //    An explicit sphere in 3D is determined by four points,
        //    which should be distinct, and not coplanar.
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], P3[3], P4[3], the coordinates of four
        //    distinct noncoplanar points on the sphere.
        //
        //    Output, double *R, PC[3], the radius and coordinates of the
        //    center of the sphere.  If the linear system is
        //    singular, then R = -1, PC[] = 0.
        //
    {
        int DIM_NUM = 3;

        double[] tet = new double[DIM_NUM * 4];

        typeMethods.r8vec_copy(DIM_NUM, p1, ref tet, a2index: +0 * 3);
        typeMethods.r8vec_copy(DIM_NUM, p2, ref tet, a2index: +1 * 3);
        typeMethods.r8vec_copy(DIM_NUM, p3, ref tet, a2index: +2 * 3);
        typeMethods.r8vec_copy(DIM_NUM, p4, ref tet, a2index: +3 * 3);

        TetrahedronNS.Geometry.tetrahedron_circumsphere_3d(tet, ref r, ref pc);
    }

    public static void sphere_exp2imp_nd(int n, double[] p, ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_EXP2IMP_ND finds an N-dimensional sphere through N+1 points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double P[N*(N+1)], the points.
        //
        //    Output, ref double R, the radius of the sphere.
        //
        //    Output, double PC[N], the center of the sphere.
        //
    {
        double[] a;
        int i;
        int info;
        int j;
        double t;
        //
        //  Set up the linear system.
        //
        a = new double[n * (n + 1)];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a[i + j * n] = p[j + (i + 1) * n];
            }
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a[i + j * n] -= p[j + 0 * n];
            }
        }

        for (i = 0; i < n; i++)
        {
            t = 0.0;
            for (j = 0; j < n; j++)
            {
                t += a[i + j * n] * a[i + j * n];
            }

            a[i + n * n] = t;
        }

        //
        //  Solve the linear system.
        //
        info = typeMethods.r8mat_solve(n, 1, ref a);
        //
        //  If the system was singular, return a consolation prize.
        //
        if (info != 0)
        {
            r = -1.0;
            for (i = 0; i < n; i++)
            {
                pc[i] = 0.0;
            }

            return;
        }

        //
        //  Compute the radius and center.
        //
        r = 0.0;
        for (i = 0; i < n; i++)
        {
            r += a[i + n * n] * a[i + n * n];
        }

        r = 0.5 * Math.Sqrt(r);

        for (i = 0; i < n; i++)
        {
            pc[i] = p[i + 0 * n] + 0.5 * a[i + n * n];
        }

    }

    public static double sphere_imp_area_3d(double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_AREA_3D computes the surface area of an implicit sphere in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Output, double SPHERE_IMP_AREA_3D, the area of the sphere.
        //
    {
        double area;

        area = 4.0 * Math.PI * r * r;

        return area;
    }

    public static double sphere_imp_area_nd(int dim_num, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_AREA_ND computes the surface area of an implicit sphere in ND.
        //
        //  Discussion:
        //
        //    DIM_NUM   Area
        //
        //    2      2       * PI   * R
        //    3      4       * PI   * R^2
        //    4      2       * PI^2 * R^3
        //    5      (8/3)   * PI^2 * R^4
        //    6                PI^3 * R^5
        //    7      (16/15) * PI^3 * R^6
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Output, double SPHERE_IMP_AREA_ND, the area of the sphere.
        //
    {
        double area;

        area = Math.Pow(r, dim_num - 1) * sphere_unit_area_nd(dim_num);

        return area;
    }

    public static bool sphere_imp_contains_point_3d(double r, double[] pc, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_CONTAINS_POINT_3D determines if an implicit sphere contains a point in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
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
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double PC[3], the coordinates of the center of the sphere.
        //
        //    Input, double P[3], the point to be checked.
        //
        //    Output, bool SPHERE_IMP_CONTAINS_POINT_3D, is TRUE if the point is inside or
        //    on the sphere, FALSE otherwise.
        //
    {
        if (Math.Pow(p[0] - pc[0], 2)
            + Math.Pow(p[1] - pc[1], 2)
            + Math.Pow(p[2] - pc[2], 2) <= r * r)
        {
            return true;
        }

        return false;
    }

    public static void sphere_imp_grid_icos_size(int factor, ref int node_num, ref int edge_num,
            ref int triangle_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_GRID_ICOS_SIZE sizes an icosahedral grid on a sphere.
        //
        //  Discussion:
        //
        //    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
        //
        //    With FACTOR = 2, each triangle of the icosahedron is subdivided into
        //    2x2 subtriangles, resulting in 80 faces, 120 edges, and
        //    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
        //
        //    With FACTOR = 3, each triangle of the icosahedron is subdivided into
        //    3x3 subtriangles, resulting in 180 faces, 270 edges, and
        //    72 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
        //
        //    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
        //    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR*FACTOR edges, and
        //      12
        //    + 20 * 3          * (FACTOR-1) / 2
        //    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2007
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
        //    Output, int *NODE_NUM, the number of nodes.
        //
        //    Output, int *EDGE_NUM, the number of edges.
        //
        //    Output, int *TRIANGLE_NUM, the number of triangles.
        //
    {
        node_num = 12
                   + 10 * 3 * (factor - 1)
                   + 10 * (factor - 2) * (factor - 1);

        edge_num = 30 * factor * factor;

        triangle_num = 20 * factor * factor;

    }

    public static void sphere_imp_gridfaces_3d(int maxtri, int nlat, int nlong, ref int ntri,
            ref int[] tri)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_GRIDFACES_3D produces a grid of triangles on an implicit sphere in 3D.
        //
        //  Discussion:
        //
        //    The point numbering system is the same used in SPHERE_IMP_GRIDPOINTS_3D,
        //    and that routine may be used to compute the coordinates of the points.
        //
        //    The two dimensional array TRI[3,MAXTRI] is stored by columns.
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MAXTRI, the maximum number of triangles.
        //
        //    Input, int NLAT, NLONG, the number of latitude and longitude
        //    lines to draw.  The latitudes do not include the North and South
        //    poles, which will be included automatically, so NLAT = 5, for instance,
        //    will result in points along 7 lines of latitude.
        //
        //    Output, int *NTRI, the number of triangles.
        //
        //    Output, int TRI[3*MAXTRI], contains NTRI triples of point indices for
        //    the triangles that make up the grid.
        //
    {
        int i;
        int j;
        int n;
        int n_max;
        int n_min;
        int ne;
        int nw;
        int s;
        int s_max;
        int s_min;
        int se;
        int sw;

        ntri = 0;
        //
        //  The first row.
        //
        n = 1;

        sw = 2;
        se = sw + 1;

        s_min = 2;
        s_max = nlong + 1;

        for (j = 0; j <= nlong - 1; j++)
        {
            if (ntri < maxtri)
            {
                tri[0 + ntri * 3] = sw;
                tri[1 + ntri * 3] = se;
                tri[2 + ntri * 3] = n;
                ntri += 1;
            }

            sw = se;

            if (se == s_max)
            {
                se = s_min;
            }
            else
            {
                se += 1;
            }

        }

        //
        //  The intermediate rows.
        //
        for (i = 1; i <= nlat; i++)
        {
            n_max = s_max;
            n_min = s_min;

            s_max += nlong;
            s_min += nlong;

            nw = n_min;
            ne = nw + 1;
            sw = s_min;
            se = sw + 1;

            for (j = 0; j <= nlong - 1; j++)
            {
                if (ntri < maxtri)
                {
                    tri[0 + ntri * 3] = sw;
                    tri[1 + ntri * 3] = se;
                    tri[2 + ntri * 3] = nw;
                    ntri += 1;
                }

                if (ntri < maxtri)
                {
                    tri[0 + ntri * 3] = ne;
                    tri[1 + ntri * 3] = nw;
                    tri[2 + ntri * 3] = se;
                    ntri += 1;
                }

                sw = se;
                nw = ne;

                if (se == s_max)
                {
                    se = s_min;
                }
                else
                {
                    se += 1;
                }

                if (ne == n_max)
                {
                    ne = n_min;
                }
                else
                {
                    ne += 1;
                }

            }

        }

        //
        //  The last row.
        //
        n_max = s_max;
        n_min = s_min;

        s = n_max + 1;

        nw = n_min;
        ne = nw + 1;

        for (j = 0; j <= nlong - 1; j++)
        {
            if (ntri < maxtri)
            {
                tri[0 + ntri * 3] = ne;
                tri[1 + ntri * 3] = nw;
                tri[2 + ntri * 3] = s;
                ntri += 1;
            }

            nw = ne;

            if (ne == n_max)
            {
                ne = n_min;
            }
            else
            {
                ne += 1;
            }

        }

    }

    public static int sphere_imp_line_project_3d(double r, double[] pc, int n, double[] p,
            int maxpnt2, ref double[] pp, double thetamin, double thetamax)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_LINE_PROJECT_3D projects a line onto an implicit sphere in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //    The line to be projected is specified as a sequence of points.
        //    If two successive points subtend a small angle, then the second
        //    point is essentially dropped.  If two successive points subtend
        //    a large angle, then intermediate points are inserted, so that
        //    the projected line stays closer to the sphere.
        //
        //    Note that if any P coincides with the center of the sphere, then
        //    its projection is mathematically undefined.  P will
        //    be returned as PC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.  If R is
        //    zero, PP will be returned as PC, and if R is
        //    negative, points will end up diametrically opposite from where
        //    you would expect them for a positive R.
        //
        //    Input, double PC[3], the coordinates of the center of the sphere.
        //
        //    Input, int N, the number of points on the line that is
        //    to be projected.
        //
        //    Input, double P[3*N], the coordinates of the points
        //    on the line that is to be projected.
        //
        //    Input, int MAXPNT2, the maximum number of points on the projected
        //    line.  Even if the routine thinks that more points are needed,
        //    no more than MAXPNT2 will be generated.
        //
        //    Output, double PP[3*N2], the coordinates of the
        //    points representing the projected line.  The value N2 is returned
        //    as the function value of this routine.  These points lie on the
        //    sphere.  Successive points are separated by at least THETAMIN
        //    radians, and by no more than THETAMAX radians.
        //
        //    Input, double THETAMIN, THETAMAX, the minimum and maximum angular
        //    projections allowed between successive projected points.
        //    If two successive points on the original line have projections
        //    separated by more than THETAMAX radians, then intermediate points
        //    will be inserted, in an attempt to keep the line closer to the
        //    sphere.  If two successive points are separated by less than
        //    THETAMIN radians, then the second point is dropped, and the
        //    line from the first point to the next point is considered.
        //
        //    Output, int SPHERE_IMP_LINE_PROJECT_3D, the number of points on
        //    the projected line.  This value can be zero, if the line has an
        //    angular projection of less than THETAMIN radians.
        //
    {
        int DIM_NUM = 3;

        double alpha;
        double ang3d;
        double dot;
        int i;
        int j;
        int nfill;
        int n2;
        double tnorm;
        double[] p1 = new double[DIM_NUM];
        double[] p2 = new double[DIM_NUM];
        double[] pi = new double[DIM_NUM];
        switch (r)
        {
            //
            //  Check the input.
            //
            case 0.0:
                n2 = 0;
                return n2;
        }

        typeMethods.r8vec_copy(DIM_NUM, pc, ref p1);
        typeMethods.r8vec_copy(DIM_NUM, pc, ref p2);

        n2 = 0;

        for (i = 0; i < n; i++)
        {
            if (typeMethods.r8vec_eq(DIM_NUM, p, pc))
            {
            }
            else
            {
                typeMethods.r8vec_copy(DIM_NUM, p2, ref p1);

                alpha = Math.Sqrt(Math.Pow(p[0 + i * 3] - pc[0], 2)
                                  + Math.Pow(p[1 + i * 3] - pc[1], 2)
                                  + Math.Pow(p[2 + i * 3] - pc[2], 2));

                p2[0] = pc[0] + r * (p[0 + i * 3] - pc[0]) / alpha;
                p2[1] = pc[1] + r * (p[1 + i * 3] - pc[1]) / alpha;
                p2[2] = pc[2] + r * (p[2 + i * 3] - pc[2]) / alpha;
                switch (n2)
                {
                    //
                    //  If we haven't gotten any points yet, take this point as our start.
                    //
                    case 0:
                        pp[0 + n2 * 3] = p2[0];
                        pp[1 + n2 * 3] = p2[1];
                        pp[2 + n2 * 3] = p2[2];
                        n2 += 1;
                        break;
                    //
                    //  Compute the angular projection of P1 to P2.
                    //
                    case >= 1:
                    {
                        dot = (p1[0] - pc[0]) * (p2[0] - pc[0])
                              + (p1[1] - pc[1]) * (p2[1] - pc[1])
                              + (p1[2] - pc[2]) * (p2[2] - pc[2]);
                        ang3d = typeMethods.r8_acos(dot / (r * r));
                        //
                        //  If the angle is at least THETAMIN, (or it's the last point),
                        //  then we will draw a line segment.
                        //
                        if (thetamin < Math.Abs(ang3d) || i == n)
                        {
                            //
                            //  Now we check to see if the line segment is too long.
                            //
                            if (thetamax < Math.Abs(ang3d))
                            {
                                nfill = (int) (Math.Abs(ang3d) / thetamax);

                                for (j = 1; j < nfill; j++)
                                {
                                    pi[0] = (nfill - j) * (p1[0] - pc[0])
                                            + j * (p2[0] - pc[0]);
                                    pi[1] = (nfill - j) * (p1[1] - pc[1])
                                            + j * (p2[1] - pc[1]);
                                    pi[2] = (nfill - j) * (p1[2] - pc[2])
                                            + j * (p2[2] - pc[2]);

                                    tnorm = typeMethods.r8vec_norm(DIM_NUM, pi);

                                    if (tnorm != 0.0)
                                    {
                                        pi[0] = pc[0] + r * pi[0] / tnorm;
                                        pi[1] = pc[1] + r * pi[1] / tnorm;
                                        pi[2] = pc[2] + r * pi[2] / tnorm;
                                        pp[0 + n2 * 3] = pi[0];
                                        pp[1 + n2 * 3] = pi[1];
                                        pp[2 + n2 * 3] = pi[2];
                                        n2 += 1;
                                    }
                                }
                            }

                            //
                            //  Now tack on the projection of point 2.
                            //
                            pp[0 + n2 * 3] = p2[0];
                            pp[1 + n2 * 3] = p2[1];
                            pp[2 + n2 * 3] = p2[2];
                            n2 += 1;
                        }

                        break;
                    }
                }
            }
        }

        return n2;
    }

    public static void sphere_imp_local2xyz_3d(double r, double[] pc, double theta, double phi,
            ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_LOCAL2XYZ_3D converts local to XYZ coordinates on an implicit sphere in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //    The "local" spherical coordinates of a point are two angles, THETA and PHI.
        //    PHI measures the angle that the vector from the origin to the point
        //    makes with the positive Z axis.  THETA measures the angle that the
        //    projection of the vector onto the XY plane makes with the positive X axis.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 June 2005
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
        //    Input, double THETA, PHI, the local (THETA,PHI) spherical coordinates
        //    of a point on the sphere.  THETA and PHI are angles, measure in
        //    radians.  Usually, 0 <= THETA < 2 * PI, and 0 <= PHI <= PI.
        //
        //    Output, double P[3], the XYZ coordinates of the point.
        //
    {
        p[0] = pc[0] + r * Math.Sin(phi) * Math.Cos(theta);
        p[1] = pc[1] + r * Math.Sin(phi) * Math.Sin(theta);
        p[2] = pc[2] + r * Math.Cos(phi);
    }

    public static void sphere_imp_point_near_3d(double r, double[] pc, double[] p,
            ref double[] pn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_POINT_NEAR_3D finds the nearest point on an implicit sphere to a point in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //    If the center of the sphere is PC, and the point is P, then
        //    the desired point lies at a positive distance R along the vector
        //    P-PC unless P = PC, in which case any point
        //    on the sphere is "nearest".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2005
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
        //    Input, double P[3], the coordinates of a point whose
        //    nearest point on the sphere is desired.
        //
        //    Output, double PN[3], the nearest point on the sphere.
        //
    {
        double norm;
        //
        //  If P = PC, bail out now.
        //
        norm = Math.Sqrt(Math.Pow(p[0] - pc[0], 2)
                         + Math.Pow(p[1] - pc[1], 2)
                         + Math.Pow(p[2] - pc[2], 2));

        switch (norm)
        {
            case 0.0:
                pn[0] = pc[0] + r;
                pn[1] = pc[1];
                pn[2] = pc[2];
                break;
            //
            default:
                pn[0] = pc[0] + r * (p[0] - pc[0]) / norm;
                pn[1] = pc[1] + r * (p[1] - pc[1]) / norm;
                pn[2] = pc[2] + r * (p[2] - pc[2]) / norm;
                break;
        }

    }

    public static void sphere_imp_point_project_3d(double r, double[] pc, double[] p,
            ref double[] pp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_POINT_PROJECT_3D projects a point onto an implicit sphere, in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2005
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
        //    Input, double P[3], the coordinates of a point.
        //
        //    Output, double PP[3], the coordinates of the point as projected
        //    onto the sphere from the center.
        //
    {
        int DIM_NUM = 3;

        double norm;

        switch (r)
        {
            case 0.0:
                typeMethods.r8vec_copy(DIM_NUM, pc, ref pp);
                break;
            default:
            {
                if (typeMethods.r8vec_eq(DIM_NUM, p, pc))
                {
                    pp[0] = pc[0];
                    pp[1] = pc[1];
                    pp[2] = pc[2] + r;
                }
                else
                {
                    norm = Math.Sqrt(Math.Pow(p[0] - pc[0], 2)
                                     + Math.Pow(p[1] - pc[1], 2)
                                     + Math.Pow(p[2] - pc[2], 2));

                    pp[0] = pc[0] + r * (p[0] - pc[0]) / norm;
                    pp[1] = pc[1] + r * (p[1] - pc[1]) / norm;
                    pp[2] = pc[2] + r * (p[2] - pc[2]) / norm;
                }

                break;
            }
        }
    }

    public static double sphere_imp_volume_3d(double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_VOLUME_3D computes the volume of an implicit sphere in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Output, double SPHERE_IMP_VOLUME_3D, the volume of the sphere.
        //
    {
        double volume;

        volume = 4.0 * Math.PI * r * r * r / 3.0;

        return volume;
    }

    public static double sphere_imp_volume_nd(int dim_num, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_VOLUME_ND computes the volume of an implicit sphere in ND.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //    DIM_NUM  Volume
        //
        //    2             PI   * R^2
        //    3  (4/3)    * PI   * R^3
        //    4  (1/2)    * PI^2 * R^4
        //    5  (8/15)   * PI^2 * R^5
        //    6  (1/6)    * PI^3 * R^6
        //    7  (16/105) * PI^3 * R^7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Output, double SPHERE_IMP_VOLUME_ND, the volume of the sphere.
        //
    {
        double value = 0;

        value = Math.Pow(r, dim_num) * sphere_unit_volume_nd(dim_num);

        return value;
    }

    public static double sphere_imp_zone_area_3d(double r, double h1, double h2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_ZONE_AREA_3D computes the surface area of a spherical zone in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //    Draw any radius of the sphere and note the point P where the radius
        //    intersects the sphere.  Now choose two points on the radius line, a
        //    distance H1 and H2 from the point P.  Consider all the points on or within
        //    the sphere whose projection onto the radius lies between these two points.
        //    These points constitute the spherical zone, which can also be considered
        //    the difference of two spherical caps.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H1, H2, the distances that define the thickness of the zone.
        //    H1 and H2 must be between 0 and 2 * R.
        //
        //    Output, double SPHERE_IMP_ZONE_AREA_3D, the area of the spherical zone.
        //
    {
        double h;

        h = Math.Abs(h1 - h2);

        switch (h)
        {
            case <= 0.0:
                return 0.0;
        }

        if (2.0 * r <= h)
        {
            return 4.0 * Math.PI * r * r;
        }
        return 2.0 * Math.PI * r * h;
    }

    public static double sphere_imp_zone_volume_3d(double r, double h1, double h2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP_ZONE_VOLUME_3D computes the volume of a spherical zone in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a sphere in 3D is:
        //
        //        Math.Pow ( P[0] - PC[0], 2 )
        //      + Math.Pow ( P[1] - PC[1], 2 )
        //      + Math.Pow ( P[2] - PC[2], 2 ) = Math.Pow ( R, 2 )
        //
        //    Draw any radius of the sphere and note the point P where the radius
        //    intersects the sphere.  Now choose two points on the radius line, a
        //    distance H1 and H2 from the point P.  Consider all the points on or within
        //    the sphere whose projection onto the radius lies between these two points.
        //    These points constitute the spherical zone, which can also be considered
        //    the difference of two spherical caps.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double H1, H2, the distances that define the thickness of the zone.
        //    H1 and H2 must be between 0 and 2 * R.
        //
        //    Output, double SPHERE_IMP_ZONE_VOLUME_3D, the volume of the spherical zone
        //
    {
        double h11;
        double h22;

        h11 = Math.Min(h1, h2);
        h11 = Math.Max(h11, 0.0);

        if (2.0 * r <= h11)
        {
            return 0.0;
        }

        h22 = Math.Max(h1, h2);
        h22 = Math.Min(h22, 2.0 * r);

        return h22 switch
        {
            <= 0.0 => 0.0,
            _ => 1.0 / 3.0 * Math.PI * (h22 * h22 * (3.0 * r - h22) - h11 * h11 * (3.0 * r - h11))
        };
    }

    public static void sphere_imp2exp_3d(double r, double[] pc, ref double[] p1, ref double[] p2,
            ref double[] p3, ref double[] p4)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_IMP2EXP_3D converts a sphere from implicit to explicit form in 3D.
        //
        //  Discussion:
        //
        //    An implicit sphere in 3D satisfies the equation:
        //
        //      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
        //
        //    An explicit sphere in 3D is determined by four points,
        //    which should be distinct, and not coplanar.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2005
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
        //    Input, double R, PC[3], the radius and center of the sphere.
        //
        //    Output, double P1[3], P2[3], P3[3], P4[3],
        //    four distinct noncoplanar points on the sphere.
        //
    {
        double phi;
        double theta;

        theta = 0.0;
        phi = 0.0;

        p1[0] = pc[0] + r * Math.Cos(theta) * Math.Sin(phi);
        p1[1] = pc[1] + r * Math.Sin(theta) * Math.Sin(phi);
        p1[2] = pc[2] + r * Math.Cos(phi);

        theta = 0.0;
        phi = 2.0 * Math.PI / 3.0;

        p2[0] = pc[0] + r * Math.Cos(theta) * Math.Sin(phi);
        p2[1] = pc[1] + r * Math.Sin(theta) * Math.Sin(phi);
        p2[2] = pc[2] + r * Math.Cos(phi);

        theta = 2.0 * Math.PI / 3.0;
        phi = 2.0 * Math.PI / 3.0;

        p3[0] = pc[0] + r * Math.Cos(theta) * Math.Sin(phi);
        p3[1] = pc[1] + r * Math.Sin(theta) * Math.Sin(phi);
        p3[2] = pc[2] + r * Math.Cos(phi);

        theta = 4.0 * Math.PI / 3.0;
        phi = 2.0 * Math.PI / 3.0;

        p4[0] = pc[0] + r * Math.Cos(theta) * Math.Sin(phi);
        p4[1] = pc[1] + r * Math.Sin(theta) * Math.Sin(phi);
        p4[2] = pc[2] + r * Math.Cos(phi);

    }

    public static double sphere_k(int dim_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_K computes a factor useful for spherical computations.
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
        //  Reference:
        //
        //    Thomas Ericson, Victor Zinoviev,
        //    Codes on Euclidean Spheres,
        //    Elsevier, 2001, pages 439-441.
        //    QA166.7 E75
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Output, double SPHERE_K, the factor.
        //
    {
        double value = (dim_num % 2) switch
        {
            0 => Math.Pow(2.0 * Math.PI, (double)dim_num / 2),
            _ => 2.0 * Math.Pow(2.0 * Math.PI, (double)(dim_num - 1) / 2)
        };

        value /= typeMethods.i4_factorial2(dim_num - 2);

        return value;
    }

    public static double sphere_triangle_angles_to_area(double r, double a, double b, double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_TRIANGLE_ANGLES_TO_AREA computes the area of a spherical triangle.
        //
        //  Discussion:
        //
        //    A sphere centered at 0 in 3D satisfies the equation:
        //
        //      X^2 + Y^2 + Z^2 = R^2
        //
        //    A spherical triangle is specified by three points on the surface
        //    of the sphere.
        //
        //    The area formula is known as Girard's formula.
        //
        //    The area of a spherical triangle is:
        //
        //      AREA = ( A + B + C - PI ) * R^2
        //
        //    where A, B and C are the (surface) angles of the triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double A, B, C, the angles of the triangle.
        //
        //    Output, double SPHERE_TRIANGLE_ANGLES_TO_AREA, the area of the spherical triangle.
        //
    {
        double area;

        area = r * r * (a + b + c - Math.PI);

        return area;
    }

    public static void sphere_triangle_sides_to_angles(double r, double as_, double bs,
            double cs, ref double a, ref double b, ref double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_TRIANGLE_SIDES_TO_ANGLES computes spherical triangle angles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double AS, BS, CS, the (geodesic) length of the sides of the
        //    triangle.
        //
        //    Output, ref double A, &B, &C, the spherical angles of the triangle.
        //    Angle A is opposite the side of length AS, and so on.
        //
    {
        double asu = 0;
        double bsu = 0;
        double csu = 0;
        double ssu = 0;
        double tan_a2 = 0;
        double tan_b2 = 0;
        double tan_c2 = 0;

        asu = as_ / r;
        bsu = bs / r;
        csu = cs / r;
        ssu = (asu + bsu + csu) / 2.0;

        tan_a2 = Math.Sqrt(Math.Sin(ssu - bsu) * Math.Sin(ssu - csu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - asu)));

        a = 2.0 * Math.Atan(tan_a2);

        tan_b2 = Math.Sqrt(Math.Sin(ssu - asu) * Math.Sin(ssu - csu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - bsu)));

        b = 2.0 * Math.Atan(tan_b2);

        tan_c2 = Math.Sqrt(Math.Sin(ssu - asu) * Math.Sin(ssu - bsu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - csu)));

        c = 2.0 * Math.Atan(tan_c2);

    }

    public static void sphere_triangle_vertices_to_angles(double r, double[] v1, double[] v2,
            double[] v3, ref double a, ref double b, ref double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_TRIANGLE_VERTICES_TO_ANGLES computes the angles of a spherical triangle.
        //
        //  Discussion:
        //
        //    A sphere centered at 0 in 3D satisfies the equation:
        //
        //      X*X + Y*Y + Z*Z = R*R
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
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
        //
        //    Output, ref double A, &B, &C, the angles of the spherical triangle.
    {
        double as_ = 0;
        double bs = 0;
        double cs = 0;
        //
        //  Compute the lengths of the sides of the spherical triangle.
        //
        sphere_triangle_vertices_to_sides(r, v1, v2, v3, ref as_, ref bs, ref cs);
        //
        //  Get the spherical angles.
        //
        sphere_triangle_sides_to_angles(r, as_, bs, cs, ref a, ref b, ref c);

    }

    public static double sphere_triangle_vertices_to_area(double r, double[] v1, double[] v2,
            double[] v3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_TRIANGLE_VERTICES_TO_AREA computes the area of a spherical triangle.
        //
        //  Discussion:
        //
        //    A sphere centered at 0 in 3D satisfies the equation:
        //
        //      X*X + Y*Y + Z*Z = R*R
        //
        //    A spherical triangle is specified by three points on the surface
        //    of the sphere.
        //
        //    The area formula is known as Girard's formula.
        //
        //    The area of a spherical triangle is:
        //
        //      AREA = ( A + B + C - PI ) * R*R
        //
        //    where A, B and C are the (surface) angles of the triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
        //
        //    Output, double SPHERE_TRIANGLE_VERTICES_TO_AREA, the area of the
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
        sphere_triangle_vertices_to_sides(r, v1, v2, v3, ref as_, ref bs, ref cs);
        //
        //  Get the spherical angles.
        //
        sphere_triangle_sides_to_angles(r, as_, bs, cs, ref a, ref b, ref c);
        //
        //  Get the area
        //
        area = sphere_triangle_angles_to_area(r, a, b, c);

        return area;
    }

    public static void sphere_triangle_vertices_to_centroid(double r, double[] v1, double[] v2,
            double[] v3, ref double[] vs)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_TRIANGLE_VERTICES_TO_CENTROID gets a spherical triangle centroid.
        //
        //  Discussion:
        //
        //    A sphere centered at 0 in 3D satisfies the equation:
        //
        //      X*X + Y*Y + Z*Z = R*R
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
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
        //
        //    Output, double VS[3], the coordinates of the "spherical centroid"
        //    of the spherical triangle.
        //
    {
        int DIM_NUM = 3;

        int i;
        double norm;

        for (i = 0; i < DIM_NUM; i++)
        {
            vs[i] = (v1[i] + v2[i] + v3[i]) / 3.0;
        }

        norm = typeMethods.r8vec_norm(DIM_NUM, vs);

        for (i = 0; i < DIM_NUM; i++)
        {
            vs[i] = r * vs[i] / norm;
        }
    }

    public static int sphere_triangle_vertices_to_orientation(double[] a, double[] b, double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_TRIANGLE_VERTICES_TO_ORIENTATION seeks the orientation of a spherical triangle.
        //
        //  Discussion:
        //
        //    Three points on a sphere actually compute two triangles; typically
        //    we are interested in the smaller of the two.
        //
        //    As long as our triangle is "small", we can define an orientation
        //    by comparing the direction of the centroid against the normal
        //    vector (C-B) x (A-B).  If the dot product of these vectors
        //    is positive, we say the triangle has positive orientation.
        //
        //    By using information from the triangle orientation, we can correctly
        //    determine the area of a Voronoi polygon by summing up the pieces
        //    of Delaunay triangles, even in the case when the Voronoi vertex
        //    lies outside the Delaunay triangle.  In that case, the areas of
        //    some of the Delaunay triangle pieces must be formally negative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3], B[3], C[3], three points on a sphere.
        //
        //    Output, int SPHERE_TRIANGLE_VERTICES_TO_ORIENTATION, is +1 if the spherical triangle
        //    is judged to have positive orientation, and -1 otherwise.
        //
    {
        double[] cd = new double[3];
        double[] cp = new double[3];
        int i;
        int o;
        double[] v1 = new double[3];
        double[] v2 = new double[3];
        //
        //  Centroid.
        //
        for (i = 0; i < 3; i++)
        {
            cd[i] = (a[i] + b[i] + c[i]) / 3.0;
        }

        //
        //  Cross product ( C - B ) x ( A - B );
        //
        for (i = 0; i < 3; i++)
        {
            v1[i] = c[i] - b[i];
            v2[i] = a[i] - b[i];
        }

        cp[0] = v1[1] * v2[2] - v1[2] * v2[1];
        cp[1] = v1[2] * v2[0] - v1[0] * v2[2];
        cp[2] = v1[0] * v2[1] - v1[1] * v2[0];
        //
        //  Compare the directions.
        //
        if (typeMethods.r8vec_dot_product(3, cp, cd) < 0.0)
        {
            o = -1;
        }
        else
        {
            o = +1;
        }

        return o;
    }

    public static void sphere_triangle_vertices_to_sides(double r, double[] v1, double[] v2,
            double[] v3, ref double as_, ref double bs, ref double cs)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_TRIANGLE_VERTICES_TO_SIDES_3D computes spherical triangle sides.
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
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double V1[3], V2[3], V3[3], the vertices of the spherical
        //    triangle.
        //
        //    Output, ref double AS, &BS, &CS, the (geodesic) length of the sides of the
        //    triangle.
        //
    {
        as_ = r * typeMethods.r8_acos(typeMethods.r8vec_dot_product(3, v2, v3) / (r * r));
        bs = r * typeMethods.r8_acos(typeMethods.r8vec_dot_product(3, v3, v1) / (r * r));
        cs = r * typeMethods.r8_acos(typeMethods.r8vec_dot_product(3, v1, v2) / (r * r));

    }

    public static double sphere_unit_area_nd(int dim_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
        //
        //  Discussion:
        //
        //    The unit sphere in ND satisfies the equation:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
        //
        //    DIM_NUM   Area
        //
        //     2    2        * PI
        //     3    4        * PI
        //     4  ( 2 /   1) * PI^2
        //     5  ( 8 /   3) * PI^2
        //     6  ( 1 /   1) * PI^3
        //     7  (16 /  15) * PI^3
        //     8  ( 1 /   3) * PI^4
        //     9  (32 / 105) * PI^4
        //    10  ( 1 /  12) * PI^5
        //
        //    For the unit sphere, Area(DIM_NUM) = DIM_NUM * Volume(DIM_NUM)
        //
        //    Sphere_Unit_Area ( DIM_NUM ) = 2 * PI^(DIM_NUM/2) / Gamma ( DIM_NUM / 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Output, double SPHERE_UNIT_AREA_ND, the area of the sphere.
        //
    {
        double area;
        int i;
        int m;

        switch (dim_num % 2)
        {
            case 0:
            {
                m = dim_num / 2;
                area = 2.0 * Math.Pow(Math.PI, m);
                for (i = 1; i <= m - 1; i++)
                {
                    area /= i;
                }

                break;
            }
            default:
            {
                m = (dim_num - 1) / 2;
                area = Math.Pow(2.0, dim_num) * Math.Pow(Math.PI, m);
                for (i = m + 1; i <= 2 * m; i++)
                {
                    area /= i;
                }

                break;
            }
        }

        return area;
    }

    public static void sphere_unit_area_values(ref int n_data, ref int n, ref double area)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
        //
        //  Discussion:
        //
        //    The formula for the surface area of the unit sphere in N dimensions is:
        //
        //      Sphere_Unit_Area ( N ) = 2 * PI^(N/2) / Gamma ( N / 2 )
        //
        //    Some values of the function include:
        //
        //       N   Area
        //
        //       2    2        * PI
        //       3  ( 4 /    ) * PI
        //       4  ( 2 /   1) * PI^2
        //       5  ( 8 /   3) * PI^2
        //       6  ( 1 /   1) * PI^3
        //       7  (16 /  15) * PI^3
        //       8  ( 1 /   3) * PI^4
        //       9  (32 / 105) * PI^4
        //      10  ( 1 /  12) * PI^5
        //
        //    For the unit sphere, Area(N) = N * Volume(N)
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      2 * Pi^(n/2) / Gamma[n/2]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.
        //    On input, if N_DATA is 0, the first test data is returned, and
        //    N_DATA is set to the index of the test data.  On each subsequent
        //    call, N_DATA is incremented and that test data is returned.  When
        //    there is no more test data, N_DATA is set to 0.
        //
        //    Output, int &N, the spatial dimension.
        //
        //    Output, ref double AREA, the area of the unit sphere
        //    in that dimension.
        //
    {
        const int N_MAX = 20;

        double[] area_vec =
        {
            0.2000000000000000E+01,
            0.6283185307179586E+01,
            0.1256637061435917E+02,
            0.1973920880217872E+02,
            0.2631894506957162E+02,
            0.3100627668029982E+02,
            0.3307336179231981E+02,
            0.3246969701133415E+02,
            0.2968658012464836E+02,
            0.2550164039877345E+02,
            0.2072514267328890E+02,
            0.1602315322625507E+02,
            0.1183817381218268E+02,
            0.8389703410491089E+01,
            0.5721649212349567E+01,
            0.3765290085742291E+01,
            0.2396678817591364E+01,
            0.1478625959000308E+01,
            0.8858104195716824E+00,
            0.5161378278002812E+00
        };

        int[] n_vec =
        {
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            area = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            area = area_vec[n_data - 1];
        }
    }

    public static double[] sphere_unit_sample_2d(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_SAMPLE_2D picks a random point on the unit sphere (circle) in 2D.
        //
        //  Discussion:
        //
        //    The unit sphere in 2D satisfies the equation:
        //
        //      X * X + Y * Y = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SPHERE_UNIT_SAMPLE_2D[2], the random point on the unit circle.
        //
    {
        double u;
        double[] x;

        u = UniformRNG.r8_uniform_01(ref seed);

        x = new double[2];

        x[0] = Math.Cos(2.0 * Math.PI * u);
        x[1] = Math.Sin(2.0 * Math.PI * u);

        return x;
    }

    public static double[] sphere_unit_sample_3d(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_SAMPLE_3D picks a random point on the unit sphere in 3D.
        //
        //  Discussion:
        //
        //    The unit sphere in 3D satisfies the equation:
        //
        //      X * X + Y * Y + Z * Z = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SPHERE_UNIT_SAMPLE_3D[3], the sample point.
        //
    {
        double phi;
        double theta;
        double vdot;
        double[] x;
        //
        //  Pick a uniformly random VDOT, which must be between -1 and 1.
        //  This represents the dot product of the random vector with the Z unit vector.
        //
        //   this works because the surface area of the sphere between
        //  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
        //  a patch of area uniformly.
        //
        vdot = 2.0 * UniformRNG.r8_uniform_01(ref seed) - 1.0;

        phi = typeMethods.r8_acos(vdot);
        //
        //  Pick a uniformly random rotation between 0 and 2 Pi around the
        //  axis of the Z vector.
        //
        theta = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);

        x = new double[3];

        x[0] = Math.Cos(theta) * Math.Sin(phi);
        x[1] = Math.Sin(theta) * Math.Sin(phi);
        x[2] = Math.Cos(phi);

        return x;
    }

    public static double[] sphere_unit_sample_3d_2(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_SAMPLE_3D_2 is a BAD method for sampling the unit sphere in 3D.
        //
        //  Discussion:
        //
        //    The unit sphere in 3D satisfies the equation:
        //
        //      X * X + Y * Y + Z * Z = 1
        //
        //    Points on the unit sphere have coordinates ( PHI, THETA ) where
        //    PHI varies from 0 to PI, and THETA from 0 to 2 PI, so that:
        //
        //    x = Math.Cos ( theta ) * Math.Sin ( phi )
        //    y = Math.Sin ( theta ) * Math.Sin ( phi )
        //    z =                 Math.Cos ( phi )
        //
        //    This routine implements a sampling of the sphere that simply
        //    picks PHI and THETA uniformly at random from their ranges.
        //    This is a uniform sampling on the cylinder, but it is NOT
        //    a uniform sampling on the sphere.  I implement it here just
        //    so I can run some tests against the code in SPHERE_UNIT_SAMPLE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SPHERE_UNIT_SAMPLE_3D_2[3], the sample point.
        //
    {
        double phi;
        double theta;
        double[] x;

        phi = Math.PI * UniformRNG.r8_uniform_01(ref seed);
        theta = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);

        x = new double[3];

        x[0] = Math.Cos(theta) * Math.Sin(phi);
        x[1] = Math.Sin(theta) * Math.Sin(phi);
        x[2] = Math.Cos(phi);

        return x;
    }

    public static double[] sphere_unit_sample_nd(int dim_num, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_SAMPLE_ND picks a random point on the unit sphere in ND.
        //
        //  Discussion:
        //
        //    The unit sphere in ND satisfies the equation:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
        //
        //    DIM_NUM-1 random Givens rotations are applied to the point ( 1, 0, 0, ..., 0 ).
        //
        //    The I-th Givens rotation is in the plane of coordinate axes I and I+1,
        //    and has the form:
        //
        //     [ Math.Cos ( theta )  - Math.Sin ( theta ) ] * x(i)      = x'(i)
        //     [ Math.Sin ( theta )    Math.Cos ( theta ) ]   x(i+1)      x'(i+1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SPHERE_UNIT_SAMPLE_ND[DIM_NUM], the random point.
        //
    {
        int i;
        double random_cosine;
        double random_sign;
        double random_sine;
        double[] p;
        double pi;

        p = new double[dim_num];

        p[0] = 1.0;
        for (i = 1; i < dim_num; i++)
        {
            p[i] = 0.0;
        }

        for (i = 0; i < dim_num - 1; i++)
        {
            random_cosine = 2.0 * UniformRNG.r8_uniform_01(ref seed) - 1.0;

            random_sign = 2 * (int) (2.0 *
                                     UniformRNG.r8_uniform_01(ref seed)) - 1;

            random_sine = random_sign
                          * Math.Sqrt(1.0 - random_cosine * random_cosine);

            pi = p[i];
            p[i] = random_cosine * pi;
            p[i + 1] = random_sine * pi;
        }

        return p;
    }

    public static double[] sphere_unit_sample_nd_2(int dim_num, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_SAMPLE_ND_2 picks a random point on the unit sphere in ND.
        //
        //  Discussion:
        //
        //    The unit sphere in ND satisfies the equation:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
        //
        //    DIM_NUM independent normally distributed random numbers are generated,
        //    and then scaled to have unit norm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SPHERE_UNIT_SAMPLE_ND_2[DIM_NUM], the random point.
        //
    {
        int i;
        double norm;
        double[] p;
        typeMethods.r8vecNormalData data = new();

        p = typeMethods.r8vec_normal_01_new(dim_num, ref data, ref seed);

        norm = typeMethods.r8vec_norm(dim_num, p);

        for (i = 0; i < dim_num; i++)
        {
            p[i] /= norm;
        }

        return p;
    }

    public static double[] sphere_unit_sample_nd_3(int dim_num, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_SAMPLE_ND_3 picks a random point on the unit sphere in ND.
        //
        //  Discussion:
        //
        //    The unit sphere in ND satisfies the equation:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
        //
        //    Points in the [-1,1] cube are generated.  Points lying outside
        //    the sphere are rejected.  Points inside the unit sphere are normalized
        //    to lie on the sphere.
        //
        //    Because the volume of the unit sphere
        //    relative to the unit cube decreases drastically in higher dimensions,
        //    this routine becomes increasingly inefficient at higher DIM_NUM.
        //    Above DIM_NUM = 5, this problem will become significant.
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
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SPHERE_UNIT_SAMPLE_ND_3[DIM_NUM], the random point.
        //
    {
        int i;
        double norm;
        double[] p;

        for (;;)
        {
            p = UniformRNG.r8vec_uniform_01_new(dim_num, ref seed);

            for (i = 0; i < dim_num; i++)
            {
                p[i] = 2.0 * p[i] - 1.0;
            }

            norm = typeMethods.r8vec_norm(dim_num, p);

            if (norm <= 1.0)
            {
                for (i = 0; i < dim_num; i++)
                {
                    p[i] /= norm;
                }

                break;
            }
        }

        return p;
    }

    public static double sphere_unit_volume_nd(int dim_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
        //
        //  Discussion:
        //
        //    The unit sphere in ND satisfies the equation:
        //
        //      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
        //
        //     DIM_NUM  Volume
        //
        //     1    2
        //     2    1        * PI
        //     3  ( 4 /   3) * PI
        //     4  ( 1 /   2) * PI^2
        //     5  ( 8 /  15) * PI^2
        //     6  ( 1 /   6) * PI^3
        //     7  (16 / 105) * PI^3
        //     8  ( 1 /  24) * PI^4
        //     9  (32 / 945) * PI^4
        //    10  ( 1 / 120) * PI^5
        //
        //    For the unit sphere, Volume(N) = 2 * PI * Volume(N-2)/ N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Output, double SPHERE_UNIT_VOLUME_ND, the volume of the sphere.
        //
    {
        int i;
        int m;
        double volume;

        switch (dim_num % 2)
        {
            case 0:
            {
                m = dim_num / 2;
                volume = 1.0;
                for (i = 1; i <= m; i++)
                {
                    volume = volume * Math.PI / i;
                }

                break;
            }
            default:
            {
                m = (dim_num - 1) / 2;
                volume = Math.Pow(Math.PI, m) * Math.Pow(2.0, dim_num);
                for (i = m + 1; i <= 2 * m + 1; i++)
                {
                    volume /= i;
                }

                break;
            }
        }

        return volume;
    }

    public static void sphere_unit_volume_values(ref int n_data, ref int n, ref double volume)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
        //
        //  Discussion:
        //
        //    The formula for the volume of the unit sphere in N dimensions is
        //
        //      Volume(N) = 2 * PI^(N/2) / ( N * Gamma ( N / 2 ) )
        //
        //    This function satisfies the relationships:
        //
        //      Volume(N) = 2 * PI * Volume(N-2) / N
        //      Volume(N) = Area(N) / N
        //
        //    Some values of the function include:
        //
        //       N  Volume
        //
        //       1    1
        //       2    1        * PI
        //       3  ( 4 /   3) * PI
        //       4  ( 1 /   2) * PI^2
        //       5  ( 8 /  15) * PI^2
        //       6  ( 1 /   6) * PI^3
        //       7  (16 / 105) * PI^3
        //       8  ( 1 /  24) * PI^4
        //       9  (32 / 945) * PI^4
        //      10  ( 1 / 120) * PI^5
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      2 * Pi^(n/2) / ( n * Gamma[n/2] )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.
        //    On input, if N_DATA is 0, the first test data is returned, and
        //    N_DATA is set to the index of the test data.  On each subsequent
        //    call, N_DATA is incremented and that test data is returned.  When
        //    there is no more test data, N_DATA is set to 0.
        //
        //    Output, int &N, the spatial dimension.
        //
        //    Output, ref double VOLUME, the volume of the unit
        //    sphere in that dimension.
        //
    {
        const int N_MAX = 20;

        int[] n_vec =
        {
            1, 2,
            3, 4,
            5, 6,
            7, 8,
            9, 10,
            11, 12,
            13, 14,
            15, 16,
            17, 18,
            19, 20
        };

        double[] volume_vec =
        {
            0.2000000000000000E+01,
            0.3141592653589793E+01,
            0.4188790204786391E+01,
            0.4934802200544679E+01,
            0.5263789013914325E+01,
            0.5167712780049970E+01,
            0.4724765970331401E+01,
            0.4058712126416768E+01,
            0.3298508902738707E+01,
            0.2550164039877345E+01,
            0.1884103879389900E+01,
            0.1335262768854589E+01,
            0.9106287547832831E+00,
            0.5992645293207921E+00,
            0.3814432808233045E+00,
            0.2353306303588932E+00,
            0.1409811069171390E+00,
            0.8214588661112823E-01,
            0.4662160103008855E-01,
            0.2580689139001406E-01
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            volume = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            volume = volume_vec[n_data - 1];
        }
    }

    public static double sphere01_distance_xyz(double[] xyz1, double[] xyz2)

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
        double bot;
        double dist;
        double lat1;
        double lat2;
        double lon1;
        double lon2;
        double top;

        lat1 = typeMethods.r8_asin(xyz1[2]);
        lon1 = typeMethods.r8_atan(xyz1[1], xyz1[0]);

        lat2 = typeMethods.r8_asin(xyz2[2]);
        lon2 = typeMethods.r8_atan(xyz2[1], xyz2[0]);

        top = Math.Pow(Math.Cos(lat2) * Math.Sin(lon1 - lon2), 2)
              + Math.Pow(Math.Cos(lat1) * Math.Sin(lat2)
                         - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2), 2);

        top = Math.Sqrt(top);

        bot = Math.Sin(lat1) * Math.Sin(lat2)
              + Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos(lon1 - lon2);

        dist = Math.Atan2(top, bot);

        return dist;
    }

    public static double sphere01_polygon_area(int n, double[] lat, double[] lon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_POLYGON_AREA returns the area of a spherical polygon.
        //
        //  Discussion:
        //
        //    On a unit sphere, the area of a spherical polygon with N sides
        //    is equal to the spherical excess:
        //
        //      E = sum ( interior angles ) - ( N - 2 ) * pi.
        //
        //    On a sphere with radius R, the area is the spherical excess multiplied
        //    by R * R.
        //
        //    The code was revised in accordance with suggestions in Carvalho
        //    and Cavalcanti.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2005
        //
        //  Author:
        //
        //    Robert Miller
        //
        //  Reference:
        //
        //    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
        //    Point in Polyhedron Testing Using Spherical Polygons,
        //    Graphics Gems, Volume V,
        //    Edited by Alan Paeth,
        //    Academic Press, 1995, T385.G6975.
        //
        //    Robert Miller,
        //    Computing the Area of a Spherical Polygon,
        //    Graphics Gems, Volume IV, pages 132-138,
        //    Edited by Paul Heckbert,
        //    Academic Press, 1994, T385.G6974.
        //
        //    Eric Weisstein,
        //    "Spherical Polygon",
        //    CRC Concise Encyclopedia of Mathematics,
        //    CRC Press, 1999.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, double LAT[N], LON[N], the latitudes and longitudes of the vertices
        //    of the spherical polygon.
        //
        //    Output, double SPHERE01_POLYGON_AREA, the area of the spherical polygon
        //    in spherical radians.
        //
    {
        double a = 0.0;
        double area = 0.0;
        double b = 0.0;
        double beta1 = 0.0;
        double beta2 = 0.0;
        double c = 0.0;
        double cos_b1 = 0.0;
        double cos_b2 = 0.0;
        double excess = 0.0;
        double hav_a = 0.0;
        int j;
        int k;
        double lam = 0.0;
        double lam1 = 0.0;
        double lam2 = 0.0;
        double s;
        double t;

        area = 0.0;

        for (j = 0; j <= n; j++)
        {
            switch (j)
            {
                case 0:
                    lam1 = lon[j];
                    beta1 = lat[j];
                    lam2 = lon[j + 1];
                    beta2 = lat[j + 1];
                    cos_b1 = Math.Cos(beta1);
                    cos_b2 = Math.Cos(beta2);
                    break;
                default:
                    k = (j + 1) % (n + 1);
                    lam1 = lam2;
                    beta1 = beta2;
                    lam2 = lon[k % lon.Length];
                    beta2 = lat[k % lat.Length];
                    cos_b1 = cos_b2;
                    cos_b2 = Math.Cos(beta2);
                    break;
            }

            if (Math.Abs(lam1 - lam2) > double.Epsilon)
            {
                hav_a = Helpers.haversine(beta2 - beta1)
                        + cos_b1 * cos_b2 * Helpers.haversine(lam2 - lam1);
                a = 2.0 * Math.Asin(Math.Sqrt(hav_a));

                b = 0.5 * Math.PI - beta2;
                c = 0.5 * Math.PI - beta1;
                s = 0.5 * (a + b + c);
                //
                //  Given the three sides of a spherical triangle, we can use a formula
                //  to find the spherical excess.
                //
                t = Math.Tan(s / 2.0) * Math.Tan((s - a) / 2.0)
                                      * Math.Tan((s - b) / 2.0) * Math.Tan((s - c) / 2.0);

                excess = Math.Abs(4.0 * Math.Atan(Math.Sqrt(Math.Abs(t))));

                if (lam1 < lam2)
                {
                    lam = lam2 - lam1;
                }
                else
                {
                    lam = lam2 - lam1 + 2.0 * Math.PI;
                }

                excess = lam switch
                {
                    > Math.PI => -excess,
                    _ => excess
                };

                area += excess;
            }
        }

        area = area switch
        {
            < 0.0 => 4.0 * Math.PI + area,
            _ => area
        };
        //area = Math.Abs ( area );

        return area;
    }

    public static double sphere01_polygon_area_karney(int n, double[] lat, double[] lon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE01_POLYGON_AREA_KARNEY returns the area of a spherical polygon.
        //
        //  Discussion:
        //
        //    On a unit sphere, the area of a spherical polygon with N sides
        //    is equal to the spherical excess:
        //
        //      E = sum ( interior angles ) - ( N - 2 ) * pi.
        //
        //    On a sphere with radius R, the area is the spherical excess multiplied
        //    by R * R.
        //
        //    The text of this function was kindly supplied by the author.
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
        //    Charles Karney
        //
        //  References:
        //
        //    Charles Karney,
        //    GeographicLib,
        //    September 2010
        //    http://geographiclib.sf.net/html/geodesic.html#geodarea
        //
        //    Isaac Todhunter,
        //    Section 103, equation 2,
        //    Spherical Trigonometry,
        //    Macmillan, 1871,
        //    http://books.google.com/books?id=3uBHAAAAIAAJ&pg=PA71
        //
        //    This routine gives the area for spherical polygons which do not include a
        //    pole.  A solution of this problem for ellipsoidal polygons which include
        //    a pole is provided by the Planimeter utility of GeographicLib; see
        //    http://geographiclib.sf.net/html/utilities.html#planimeter
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, double LAT[N], LON[N], the latitudes and longitudes of the
        //    vertices of the spherical polygon, measured in radians.
        //
        //    Output, double SPHERE01_POLYGON_AREA_KARNEY, the signed area of the
        //    spherical polygon in spherical radians.  Clockwise traversal yields a
        //    positive result.
        //
    {
        double area;
        int j;
        double lam1;
        double lam2;
        double tbeta1;
        double tbeta2;

        area = 0.0;
        lam2 = lon[n - 1];
        tbeta2 = Math.Tan(lat[n - 1] / 2.0);

        for (j = 0; j < n; j++)
        {
            lam1 = lam2;
            lam2 = lon[j];
            tbeta1 = tbeta2;
            tbeta2 = Math.Tan(lat[j] / 2.0);
            area += Math.Atan2(Math.Tan((lam1 - lam2) / 2.0) * (tbeta1 + tbeta2),
                1.0 + tbeta1 * tbeta2);
        }

        area = 2.0 * area;

        return area;
    }

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
        //    Output, double SPHERE_TRIANGLE_ANGLES_TO_AREA, the area of the
        //    spherical triangle.
        //
    {
        double area;

        area = a + b + c - Math.PI;

        return area;
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
        //    Output, ref double A, &B, &C, the spherical angles of the triangle.
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

        tan_a2 = Math.Sqrt(Math.Sin(ssu - bsu) * Math.Sin(ssu - csu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - asu)));

        a = 2.0 * Math.Atan(tan_a2);

        tan_b2 = Math.Sqrt(Math.Sin(ssu - asu) * Math.Sin(ssu - csu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - bsu)));

        b = 2.0 * Math.Atan(tan_b2);

        tan_c2 = Math.Sqrt(Math.Sin(ssu - asu) * Math.Sin(ssu - bsu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - csu)));

        c = 2.0 * Math.Atan(tan_c2);
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
        //    Output, ref double A, &B, &C, the angles of the spherical triangle.
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
    }

    public static double sphere01_triangle_vertices_to_area(double[] v1, double[] v2,
            double[] v3)

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
        //    Output, double SPHERE_TRIANGLE_VERTICES_TO_AREA, the area of the
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

    public static double[] sphere01_triangle_vertices_to_centroid(double[] v1, double[] v2,
            double[] v3)

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
        //    Output, double SPHERE01_TRIANGLE_VERTICES_TO_CENTROID[3], the coordinates
        //    of the "spherical centroid" of the spherical triangle.
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
            vs[i] /= norm;
        }

        return vs;
    }

    public static void sphere01_triangle_vertices_to_midpoints(double[] v1, double[] v2,
            double[] v3, ref double[] m1, ref double[] m2, ref double[] m3)

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
            m1[i] /= norm;
        }

        for (i = 0; i < 3; i++)
        {
            m2[i] = (v2[i] + v3[i]) / 2.0;
        }

        norm = typeMethods.r8vec_norm(3, m2);
        for (i = 0; i < 3; i++)
        {
            m2[i] /= norm;
        }

        for (i = 0; i < 3; i++)
        {
            m3[i] = (v3[i] + v1[i]) / 2.0;
        }

        norm = typeMethods.r8vec_norm(3, m3);
        for (i = 0; i < 3; i++)
        {
            m3[i] /= norm;
        }
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
        //    Output, ref double AS, &BS, &CS, the (geodesic) length of the sides of the
        //    triangle.
        //
    {
        as_ = typeMethods.r8_acos(typeMethods.r8vec_dot_product(3, v2, v3));
        bs = typeMethods.r8_acos(typeMethods.r8vec_dot_product(3, v3, v1));
        cs = typeMethods.r8_acos(typeMethods.r8vec_dot_product(3, v1, v2));

    }

}