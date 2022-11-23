using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.CircleNS;

public static class Segment
{
    public static double circle_segment_angle_from_chord(double r, double[] c, double[] p1,
            double[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_ANGLE_FROM_CHORD computes the angle of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double P1[2], P2[2], the ends of the chord.
        //
        //    Output, double CIRCLE_SEGMENT_ANGLE_FROM_CHORD, the value of THETA,
        //    the angle of the circle segment.  0 <= THETA < 2 * PI.
        //
    {
        double[] v1 = new double[2];
        double[] v2 = new double[2];
        //
        //  Compute the radial vectors V1 and V2.
        //
        v1[0] = p1[0] - c[0];
        v1[1] = p1[1] - c[1];
        v2[0] = p2[0] - c[0];
        v2[1] = p2[1] - c[1];
        //
        //  The arc cosine will only give us an answer between 0 and PI.
        //
        double theta = typeMethods.r8_atan(v2[1], v2[0]) - typeMethods.r8_atan(v1[1], v1[0]);
        //
        //  Force 0 <= THETA < 2 * PI.
        //
        while (theta < 0.0)
        {
            theta += 2.0 * Math.PI;
        }

        while (2.0 * Math.PI <= theta)
        {
            theta -= 2.0 * Math.PI;
        }

        return theta;
    }

    public static double circle_segment_angle_from_chord_angles(double omega1, double omega2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_ANGLE_FROM_CHORD_ANGLES computes angle of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double OMEGA1, OMEGA2, the angles of the points P1 
        //    and P2.  OMEGA1 <= OMEGA2.
        //
        //    Output, double CIRCLE_SEGMENT_ANGLE_FROM_CHORD_ANGLES, the angle THETA
        //    of the circle segment.  Essentially, THETA = OMEGA2 - OMEGA1.
        //
    {
        while (omega2 < omega1)
        {
            omega2 += 2.0 * Math.PI;
        }

        double theta = omega2 - omega1;

        return theta;
    }

    public static double circle_segment_angle_from_height(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT computes the angle of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double H, the "height" of the circle segment.
        //    0 <= H <= 2 * R.
        //
        //    Output, double CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT, the angle THETA
        //    of the circle segment.
        //
    {
        double theta;

        switch (h)
        {
            case <= 0.0:
                theta = 0.0;
                break;
            default:
            {
                if (h <= r)
                {
                    theta = 2.0 * Math.Acos((r - h) / r);
                    theta = 2.0 * Math.Asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);
                }
                else if (h <= 2.0 * r)
                {
                    theta = 2.0 * Math.Acos((r - h) / r);
                    theta = 2.0 * Math.Asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);
                    theta = 2.0 * Math.PI - theta;
                }
                else
                {
                    theta = 2.0 * Math.PI;
                }

                break;
            }
        }

        return theta;
    }

    public static double circle_segment_area_from_angle(double r, double theta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_AREA_FROM_ANGLE computes the area of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double THETA, the angle of the circle segment.
        //
        //    Output, double CIRCLE_SEGMENT_AREA_FROM_ANGLE, the area of the 
        //    circle segment.
        //
    {
        double area = r * r * (theta - Math.Sin(theta)) / 2.0;

        return area;
    }

    public static double circle_segment_area_from_chord(double r, double[] c, double[] p1,
            double[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_AREA_FROM_CHORD computes the area of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double P1[2], P2[2], the ends of the chord.
        //
        //    Output, double CIRCLE_SEGMENT_AREA_FROM_CHORD, the area of the 
        //    circle segment.
        //
    {
        double theta = circle_segment_angle_from_chord(r, c, p1, p2);

        double area = r * r * (theta - Math.Sin(theta)) / 2.0;

        return area;
    }

    public static double circle_segment_area_from_height(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double H, the height of the circle segment.
        //    0 <= H <= 2 * R.
        //
        //    Output, double CIRCLE_SEGMENT_AREA_FROM_HEIGHT, the area of the 
        //    circle segment.
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
                double theta;
                if (h <= r)
                {
                    theta = 2.0 * Math.Asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);
                    area = r * r * (theta - Math.Sin(theta)) / 2.0;
                }
                else if (h <= 2.0 * r)
                {
                    theta = 2.0 * Math.Asin(Math.Sqrt(r * r - (r - h) * (r - h)) / r);
                    theta = 2.0 * Math.PI - theta;
                    area = r * r * (theta - Math.Sin(theta)) / 2.0;
                }
                else
                {
                    area = Math.PI * r * r;
                }

                break;
            }
        }

        return area;
    }

    public static double circle_segment_area_from_sample(double r, double[] c, double[] p1,
            double[] p2, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_AREA_FROM_SAMPLE computes the area of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double P1[2], P2[2], the ends of the chord.
        //
        //    Input, int N, the number of sample points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CIRCLE_SEGMENT_AREA_FROM_SAMPLE, the area of the 
        //    circle segment.
        //
    {
        int i;
        double[] p = new double[2];
        //
        //  Determine the angles of the chord endpoints.
        //
        double omega1 = typeMethods.r8_atan(p1[1] - c[1], p1[0] - c[0]);
        while (omega1 < 0.0)
        {
            omega1 += 2.0 * Math.PI;
        }

        double omega2 = typeMethods.r8_atan(p2[1] - c[1], p2[0] - c[0]);
        while (omega2 < omega1)
        {
            omega2 += 2.0 * Math.PI;
        }

        //
        //  Get N random points in the circle.
        //  To simplify angle measurement, take OMEGA1 as your smallest angle.
        //  That way, the check OMEGA1 <= ANGLE <= OMEGA2 will be legitimate.
        //
        double[] angle = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        for (i = 0; i < n; i++)
        {
            angle[i] = omega1 + 2.0 * Math.PI * angle[i];
        }

        double[] r2 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        for (i = 0; i < n; i++)
        {
            r2[i] = Math.Sqrt(r2[i]);
        }

        double[] x = new double[n];
        double[] y = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = c[0] + r2[i] * Math.Cos(angle[i]);
            y[i] = c[0] + r2[i] * Math.Sin(angle[i]);
        }

        //
        //  Determine the vector that touches the circle segment base.
        //
        p[0] = 0.5 * (p1[0] + p2[0]) - c[0];
        p[1] = 0.5 * (p1[1] + p2[1]) - c[1];

        double rmh = Math.Sqrt(p[0] * p[0] + p[1] * p[1]);

        p[0] /= rmh;
        p[1] /= rmh;

        switch (omega2 - omega1)
        {
            case > Math.PI:
                p[0] = -p[0];
                p[1] = -p[1];
                rmh = -rmh;
                break;
        }

        //
        //  Compute the projection of each point onto P.
        //
        double[] vdotp = new double[n];
        for (i = 0; i < n; i++)
        {
            vdotp[i] = (x[i] - c[0]) * p[0] + (y[i] - c[1]) * p[1];
        }

        //
        //  Points in the segment lie in the sector, and project at least RMH onto P.
        //
        int m = 0;
        for (i = 0; i < n; i++)
        {
            if (omega1 < angle[i] && angle[i] < omega2 && rmh < vdotp[i])
            {
                m += 1;
            }
        }

        //
        //  The area of the segment is its relative share of the circle area.
        //
        double area = Math.PI * r * r * m / n;

        return area;
    }

    public static double circle_segment_cdf(double r, double h, double h2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_CDF computes a CDF related to a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //    Now, suppose we want to assign a cumulative density function or CDF
        //    based on a variable H2 which measures the height of the circle segment
        //    formed by an arbitrary point in the circle segment.  CDF(H2) will
        //    measure the probability that a point drawn uniformly at random from
        //    the circle segment defines a (smaller) circle segment of height H2.
        //
        //    If we can define this CDF, then we will be able to sample uniformly
        //    from the circle segment, since our "Y" value can be determined from H2,
        //    and our X value is chosen uniformly over the horizontal chord 
        //    through (0,Y).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double H, the "height" of the circle segment.
        //    0 <= H <= 2 * R.
        //
        //    Input, double H2, the "height" of the new circle segment 
        //    defined by a given point in the circle segment.  0 <= H2 <= H.
        //
        //    Output, double CDF, the cumulative density function for H2, 
        //    the probability that a point chosen at random in the circle segment 
        //    would define a smaller circle segment of height H2 or less.
        //
    {
        double cdf;

        switch (h2)
        {
            case <= 0.0:
                cdf = 0.0;
                break;
            default:
            {
                if (h <= h2)
                {
                    cdf = 1.0;
                }
                else
                {
                    double a = circle_segment_area_from_height(r, h);
                    double a2 = circle_segment_area_from_height(r, h2);
                    cdf = a2 / a;
                }

                break;
            }
        }

        return cdf;
    }

    public static double[] circle_segment_centroid_from_chord(double r, double[] c,
            double[] p1, double[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_CENTROID_FROM_CHORD computes the centroid of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //    For this function, we assume that the center of the circle is at (0,0),
        //    that the chord is horizontal, and that the circle segment is at the top.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double P1[2], P2[2], the coordinates of the endpoints 
        //    of the chord.
        //
        //    Output, double CIRCLE_SEGMENT_CENTROID_FROM_CHORD[2], the coordinates 
        //    of the centroid.
        //
    {
        double[] v1 = new double[2];
        //
        //  Get the angle subtended by P1:P2.
        //
        double theta = circle_segment_angle_from_chord(r, c, p1, p2);
        //
        //  Construct V1, the vector from C to P1.
        //
        v1[0] = p1[0] - c[0];
        v1[1] = p1[1] - c[1];
        //
        //  Rotate V1 through THETA / 2.
        //
        double thetah = theta / 2.0;

        double[] d = new double[2];
        d[0] = Math.Cos(thetah) * v1[0] - Math.Sin(thetah) * v1[1];
        d[1] = Math.Sin(thetah) * v1[0] + Math.Cos(thetah) * v1[1];
        //
        //  Scale this vector so it represents the distance to the centroid
        //  relative to R.
        //
        double s = 4.0 * Math.Pow(Math.Sin(theta / 2.0), 3)
                   / 3.0 / (theta - Math.Sin(theta));

        d[0] = s * d[0];
        d[1] = s * d[1];
        //
        //  Add the center.
        //
        d[0] += c[0];
        d[1] += c[1];

        return d;
    }

    public static double[] circle_segment_centroid_from_height(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_CENTROID_FROM_HEIGHT computes centroid of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //    For this function, we assume that the center of the circle is at (0,0).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double H, the "height" of the circle segment.
        //    0 <= H <= 2 * R.
        //
        //    Output, double CIRCLE_SEGMENT_CENTROID_FROM_HEIGHT[2], the coordinates 
        //    of the centroid.
        //
    {
        double theta = circle_segment_angle_from_height(r, h);

        double[] d = new double[2];

        d[0] = 0.0;
        d[1] = 4.0 * r * Math.Pow(Math.Sin(theta / 2.0), 3) / 3.0
                                                            / (theta - Math.Sin(theta));

        return d;
    }

    public static double[] circle_segment_centroid_from_sample(double r, double[] c,
            double[] p1, double[] p2, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE estimates a circle segment centroid.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double P1[2], P2[2], the ends of the chord.
        //
        //    Input, int N, the number of sample points.
        //
        //    Input/output, int *&EED, a seed for the random 
        //    number generator.
        //
        //    Output, double CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE[2], the estimated 
        //    centroid of the circle segment.
        //
    {
        double[] x = new double[n];
        double[] y = new double[n];

        circle_segment_sample_from_chord(r, c, p1, p2, n, ref seed, ref x, ref y);

        double[] d = new double[2];

        d[0] = typeMethods.r8vec_sum(n, x) / n;
        d[1] = typeMethods.r8vec_sum(n, y) / n;

        return d;
    }

    public static int circle_segment_contains_point(double r, double[] c, double omega1,
            double omega2, double[] xy, int xyIndex )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_CONTAINS_POINT reports whether a point is in a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //    In this function, we allow the circle to have an arbitrary center C,
        //    arbitrary radius R, and we describe the points P1 and P2 by specifying
        //    their angles OMEGA1 and OMEGA2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double OMEGA1, OMEGA2, the angles of the two points on 
        //    the circumference of the circle that define the circle segment.
        //    OMEGA1 < OMEGA2 <= OMEGA1 + 2 * PI
        //
        //    Input, double XY[2], a point.
        //
        //    Output, int CIRCLE_SEGMENT_CONTAINS_POINT, is TRUE if the point is inside 
        //    the circle segment.
        //
    {
        double[] v = new double[2];
        int value;

        switch (r)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_SEGMENT_CONTAINS_POINT - Fatal error!");
                Console.WriteLine("  R <= 0.0.");
                return 1;
        }

        while (omega2 < omega1)
        {
            omega2 += 2.0 * Math.PI;
        }

        //
        //  Compute the vector V = XY - C:
        //
        v[0] = xy[0 + xyIndex] - c[0];
        v[1] = xy[1 + xyIndex] - c[1];
        //
        //  a: Point must be inside the circle, so ||V|| <= R.
        //
        double v_r = Math.Sqrt(v[0] * v[0] + v[1] * v[1]);

        if (r < v_r)
        {
            value = 0;
            return value;
        }

        //
        //  b: Angle made by the vector V must be between OMEGA1 and OMEGA2.
        //
        double v_omega = Math.Atan2(v[1], v[0]);

        while (omega1 <= v_omega + 2.0 * Math.PI)
        {
            v_omega -= 2.0 * Math.PI;
        }

        while (v_omega + 2.0 * Math.PI <= omega1)
        {
            v_omega += 2.0 * Math.PI;
        }

        if (omega2 < v_omega)
        {
            value = 0;
            return value;
        }

        //
        //  c: Projection of V onto unit centerline must be at least R-H.
        //
        double omegah = 0.5 * (omega1 + omega2);
        double v_project = v[0] * Math.Cos(omegah) + v[1] * Math.Sin(omegah);

        double theta = omega2 - omega1;
        double h = circle_segment_height_from_angle(r, theta);

        if (v_project < r - h)
        {
            value = 0;
            return value;
        }

        value = 1;

        return value;
    }

    public static double circle_segment_height_from_angle(double r, double angle)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE: height of a circle segment from angle.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //    This function is given the radius R and angle of the segment, and
        //    determines the corresponding height.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double ANGLE, the angle of the circle segment.
        //    0 <= ANGLE <= 2.0 * PI.
        //
        //    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE, the height of the 
        //    circle segment.
        //
    {
        double h;

        switch (angle)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!");
                Console.WriteLine("  ANGLE < 0.0.");
                return 1;
            case 0.0:
                h = 0.0;
                return h;
            case 2.0 * Math.PI:
                h = 2.0 * r;
                return h;
            case > 2.0 * Math.PI:
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!");
                Console.WriteLine("  2.0 * Math.PI < ANGLE.");
                return 1;
        }

        switch (r)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!");
                Console.WriteLine("  R <= 0.0.");
                return 1;
        }

        h = angle switch
        {
            <= Math.PI => r * (1.0 - Math.Cos(angle / 2.0)),
            _ => r * (1.0 + Math.Cos((2.0 * Math.PI - angle) / 2.0))
        };

        return h;
    }

    public static double circle_segment_height_from_area(double r, double area)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_HEIGHT_FROM_AREA: height of a circle segment from area.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //    This function is given the radius R and area of the segment, and
        //    determines the corresponding height.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double AREA, the area of the circle segment.
        //    0 <= AREA <= 2.0 * PI * R^2.
        //
        //    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_AREA, the height of the
        //    circle segment.
        //
    {
        double h = 0;

        switch (area)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!");
                Console.WriteLine("  AREA < 0.0.");
                return 1;
        }

        double area_circle = 2.0 * Math.PI * r * r;

        switch (area)
        {
            case 0.0:
                h = 0.0;
                return h;
        }

        if (Math.Abs(area - area_circle) <= typeMethods.r8_epsilon())
        {
            h = 2.0 * r;
            return h;
        }

        if (area_circle < area)
        {
            Console.WriteLine("");
            Console.WriteLine("CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!");
            Console.WriteLine("  2.0 * Math.PI * r^2 < AREA.");
            return 1;
        }

        switch (r)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!");
                Console.WriteLine("  R <= 0.0.");
                return 1;
        }

        double h1 = 0.0;
        //circle_segment_area_from_height ( r, h1 );
        double h2 = 2.0 * r;
        //circle_segment_area_from_height ( r, h2 );

        int it = 0;
        double eps = typeMethods.r8_epsilon();

        while (it < 30)
        {
            h = 0.5 * (h1 + h2);
            double a = circle_segment_area_from_height(r, h);
            it += 1;

            if (Math.Abs(a - area) < Math.Sqrt(eps) * area)
            {
                break;
            }

            if (a < area)
            {
                h1 = h;
                //    a1 = a;
            }
            else
            {
                h2 = h;
                //    a2 = a;
            }
        }

        return h;
    }

    public static double circle_segment_height_from_chord(double r, double[] c, double[] p1,
            double[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_HEIGHT_FROM_CHORD: height of a circle segment from chord.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the coordinates of the circle center.
        //
        //    Input, double P1[2], P2[2], the coordinates of the 
        //    chord endpoints.
        //
        //    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_CHORD, the height of the circle segment.
        //
    {
        double theta = circle_segment_angle_from_chord(r, c, p1, p2);

        double h = circle_segment_height_from_angle(r, theta);

        return h;
    }

    public static double circle_segment_rotation_from_chord(double r, double[] c, double[] p1,
            double[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_ROTATION_FROM_CHORD computes the rotation of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double P1[2], P2[2], the ends of the chord.
        //    Warning! If P1 = P2, we can't tell whether the segment is the whole
        //    circle or none of it!
        //
        //    Output, double CIRCLE_SEGMENT_ROTATION_FROM_CHORD, the rotation of the 
        //    circle segment.  0 <= ALPHA < 2 * PI.
        //
    {
        double[] v1 = new double[2];
        double[] v2 = new double[2];
        //
        //  Compute the radial vectors V1 and V2.
        //
        v1[0] = p1[0] - c[0];
        v1[1] = p1[1] - c[1];
        v2[0] = p2[0] - c[0];
        v2[1] = p2[1] - c[1];
        //
        //  Use R8_ATAN to guarantee that 0 <= RHO1, RHO2 <= 2 * PI.
        //
        double rho1 = typeMethods.r8_atan(v1[1], v1[0]);
        double rho2 = typeMethods.r8_atan(v2[1], v2[0]);
        //
        //  Force RHO2 to be bigger than RHO1.
        //
        while (rho2 <= rho1)
        {
            rho2 += 2.0 * Math.PI;
        }

        //
        //  Compute THETA.
        //
        double theta = rho2 - rho1;
        //
        //  ALPHA is RHO1, plus half of the angular distance between P1 and P2.
        //
        double alpha = rho1 + 0.5 * theta;

        while (2.0 * Math.PI <= alpha)
        {
            alpha -= 2.0 * Math.PI;
        }

        return alpha;
    }

    public static void circle_segment_sample_from_chord(double r, double[] c, double[] p1,
            double[] p2, int n, ref int seed, ref double[] x, ref double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_SAMPLE_FROM_CHORD samples points from a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double C[2], the center of the circle.
        //
        //    Input, double P1[2], P2[2], the endpoints of the chord.
        //
        //    Input, int N, the number of sample points.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double X[N], Y[N], the sample points.
        //
    {
        double[] c2 = new double[2];
        int i;
        double[] vc = new double[2];
        double[] vr = new double[2];
        //
        //  Determine unit vectors VR and VC.
        //  VR points to the center of the chord from the radius.
        //  VC points along the chord, from P1 to P2.
        //
        vr[0] = 0.5 * (p1[0] + p2[0]) - c[0];
        vr[1] = 0.5 * (p1[1] + p2[1]) - c[1];

        double t = Math.Sqrt(vr[0] * vr[0] + vr[1] * vr[1]);
        vr[0] /= t;
        vr[1] /= t;

        vc[0] = p2[0] - p1[0];
        vc[1] = p2[1] - p1[1];

        t = Math.Sqrt(vc[0] * vc[0] + vc[1] * vc[1]);
        vc[0] /= t;
        vc[1] /= t;
        //
        //  Get the height of the circle segment.
        //
        c2[0] = 0.0;
        c2[1] = 0.0;
        double h = circle_segment_height_from_chord(r, c2, p1, p2);
        //
        //  Sample (xi,eta) in the reference coordinates, where the chord
        //  is horizontal.
        //
        double[] xi = new double[n];
        double[] eta = new double[n];
        circle_segment_sample_from_height(r, h, n, ref seed, ref xi, ref eta);
        //
        //  XI is the left/right coordinate along VC.
        //  ETA is the distance along VR.
        //
        for (i = 0; i < n; i++)
        {
            x[i] = c[0] + eta[i] * vr[0] + xi[i] * vc[0];
            y[i] = c[1] + eta[i] * vr[1] + xi[i] * vc[1];
        }
    }

    public static void circle_segment_sample_from_height(double r, double h, int n, ref int seed,
            ref double[] x, ref double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples points from a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double H, the height of the circle segment.
        //    0 <= H <= 2 * R.
        //
        //    Input, int N, the number of sample points.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double X[N], Y[N], the sample points.
        //
    {
        int i;

        double area = circle_segment_area_from_height(r, h);
        //
        //  Pick CDF's randomly.
        //
        double[] u = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Choose points randomly by choosing ordered areas randomly.
        //
        double[] area2 = new double[n];
        for (i = 0; i < n; i++)
        {
            area2[i] = u[i] * area;
        }

        //
        //  Each area corresponds to a height H2.  Find it.
        //
        double[] h2 = new double[n];
        for (i = 0; i < n; i++)
        {
            h2[i] = circle_segment_height_from_area(r, area2[i]);
        }

        //
        //  Determine the half-width WH of the segment for each H2.
        //
        double[] wh = new double[n];
        for (i = 0; i < n; i++)
        {
            wh[i] = Math.Sqrt(h2[i] * (2.0 * r - h2[i]));
        }

        //
        //  Choose an X position randomly in [-WH,+WH].
        //
        u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        for (i = 0; i < n; i++)
        {
            x[i] = (2.0 * u[i] - 1.0) * wh[i];
        }

        //
        //  Our circle center is at (0,0).  Our height of H2 is subtracted
        //  from the height R at the peak of the circle.  Determine the Y
        //  coordinate using this fact.
        //
        for (i = 0; i < n; i++)
        {
            y[i] = r - h2[i];
        }
    }

    public static double circle_segment_width_from_height(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SEGMENT_WIDTH_FROM_HEIGHT computes the width of a circle segment.
        //
        //  Discussion:
        //
        //    Begin with a circle of radius R.  Choose two points P1 and P2 on the
        //    circle, and draw the chord P1:P2.  This chord divides the circle
        //    into two pieces, each of which is called a circle segment.
        //    Consider one of the pieces.  The "angle" of this segment is the angle 
        //    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
        //    the chord P1:P2 which is closest to C.  The "height" of the segment
        //    is the distance from Q to the perimeter of the circle.  The "width"
        //    of the circle segment is the length of P1:P2.
        //
        //    This function is given the radius R and height H of the segment, and
        //    determines the corresponding width W.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //    0 < R.
        //
        //    Input, double H, the height of the circle segment.
        //    0 <= H <= 2 * R.
        //
        //    Output, double CIRCLE_SEGMENT_WIDTH_FROM_HEIGHT, the width of the 
        //    circle segment.
        //
    {
        double w = 2.0 * Math.Sqrt(h * (2.0 * r - h));

        return w;
    }
}