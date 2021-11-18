using System;
using Burkardt.SortNS;
using Burkardt.Types;

namespace Burkardt.Geometry;

public static class Segments
{
    public static void segment_contains_point_1d(double p1, double p2, double p3, ref double u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_CONTAINS_POINT_1D reports if a line segment contains a point in 1D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1, P2, two points defining a line segment.
        //    The line segment has origin at P1, and unit at P2.
        //
        //    Input, double P3, a point to be tested.
        //
        //    Output, double *U, the coordinate of P3 in units of (P2-P1).
        //    The point P3 is contained in the line segment if 0 <= U <= 1.
        //
    {
        double unit;

        unit = p2 - p1;

        switch (unit)
        {
            case 0.0 when Math.Abs(p3 - p1) <= double.Epsilon:
                u = 0.5;
                break;
            case 0.0 when p3 < p1:
                u = -typeMethods.r8_huge();
                break;
            case 0.0:
            {
                if (p1 < p3)
                {
                    u = typeMethods.r8_huge();
                }

                break;
            }
            default:
                u = (p3 - p1) / unit;
                break;
        }

    }

    public static void segment_contains_point_2d(double[] p1, double[] p2, double[] p3,
            ref double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    In exact arithmetic, point P3 is on the line segment between
        //    P1 and P2 if and only if 0 <= U(1) <= 1 and U(2) = 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the endpoints of a line segment.
        //
        //    Input, double P3[2], a point to be tested.
        //
        //    Output, double U[2], U[0] is the coordinate of P3 along the axis from
        //    with origin at P1 and unit at P2, and U[1] is the magnitude of the off-axis
        //    portion of the  vector P3-P1, measured in units of (P2-P1).
        //
    {
        int DIM_NUM = 2;

        double t1;
        double t2;
        double unit;

        unit = Math.Sqrt((p2[0] - p1[0]) * (p2[0] - p1[0])
                         + (p2[1] - p1[1]) * (p2[1] - p1[1]));

        switch (unit)
        {
            case 0.0 when typeMethods.r8vec_eq(DIM_NUM, p1, p3):
                u[0] = 0.5;
                u[1] = 0.0;
                break;
            case 0.0:
                u[0] = 0.5;
                u[1] = typeMethods.r8_huge();
                break;
            default:
                u[0] = ((p3[0] - p1[0]) * (p2[0] - p1[0])
                        + (p3[1] - p1[1]) * (p2[1] - p1[1]))
                       / (unit * unit);

                t1 = (u[0] - 1.0) * p1[0] - u[0] * p2[0] + p3[0];
                t2 = (u[0] - 1.0) * p1[1] - u[0] * p2[1] + p3[1];

                u[1] = Math.Sqrt(t1 * t1 + t2 * t2) / unit;
                break;
        }
    }

    public static void segment_point_coords_2d(double[] p1, double[] p2, double[] p,
            ref double s, ref double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_COORDS_2D: coordinates of a point on a line segment in 2D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points P1 and P2.
        //
        //    By the coordinates of a point P with respect to a line segment [P1,P2]
        //    we mean numbers S and T such that S gives us the distance from the
        //    point P to the nearest point PN on the line (not the line segment!),
        //    and T gives us the position of PN relative to P1 and P2.
        //
        //    If S is zero, then P lies on the line.
        //
        //    If 0 <= T <= 1, then PN lies on the line segment.
        //
        //    If both conditions hold, then P lies on the line segment.
        //
        //    If E is the length of the line segment, then the distance of the
        //    point to the line segment is:
        //
        //      Math.Sqrt ( S^2 +  T^2    * E^2 )     if      T <= 0;
        //             S                         if 0 <= T <= 1
        //      Math.Sqrt ( S^2 + (T-1)^2 * E62 )     if 1 <= T
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the endpoints of the line segment.
        //
        //    Input, double P[2], the point to be considered.
        //
        //    Output, double *S, the distance of P to the nearest point PN
        //    on the line through P1 and P2.  (S will always be nonnegative.)
        //
        //    Output, double *T, the relative position of the point PN
        //    to the points P1 and P2.
        //
    {
        int DIM_NUM = 2;

        double bot;
        int i;
        double[] pn = new double[DIM_NUM];
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2))
        {
            t = 0.0;
        }
        else
        {
            bot = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                bot += Math.Pow(p2[i] - p1[i], 2);
            }

            t = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                t += (p[i] - p1[i]) * (p2[i] - p1[i]);
            }

            t /= bot;
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[i] + t * (p2[i] - p1[i]);
        }

        s = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            s += Math.Pow(p[i] - pn[i], 2);
        }

        s = Math.Sqrt(s);
    }

    public static void segment_point_coords_3d(double[] p1, double[] p2, double[] p,
            ref double s, ref double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_COORDS_3D: coordinates of a point on a line segment in 3D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points P1 and P2.
        //
        //    By the coordinates of a point P with respect to a line segment [P1,P2]
        //    we mean numbers S and T such that S gives us the distance from the
        //    point P to the nearest point PN on the line (not the line segment!),
        //    and T gives us the position of PN relative to P1 and P2.
        //
        //    If S is zero, then P lies on the line.
        //
        //    If 0 <= T <= 1, then PN lies on the line segment.
        //
        //    If both conditions hold, then P lies on the line segment.
        //
        //    If E is the length of the line segment, then the distance of the
        //    point to the line segment is:
        //
        //      Math.Sqrt ( S^2 +  T^2    * E^2 )     if      T <= 0;
        //             S                         if 0 <= T <= 1
        //      Math.Sqrt ( S^2 + (T-1)^2 * E^2 )     if 1 <= T
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the endpoints of the line segment.
        //
        //    Input, double P[3], the point to be considered.
        //
        //    Output, double *S, the distance of P to the nearest point PN
        //    on the line through P1 and P2.  (S will always be nonnegative.)
        //
        //    Output, double *T, the relative position of the point PN
        //    to the points P1 and P2.
        //
    {
        int DIM_NUM = 3;

        double bot;
        int i;
        double[] pn = new double[DIM_NUM];
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2))
        {
            t = 0.0;
        }
        else
        {
            bot = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                bot += Math.Pow(p2[i] - p1[i], 2);
            }

            t = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                t += (p[i] - p1[i]) * (p2[i] - p1[i]);
            }

            t /= bot;
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[i] + t * (p2[i] - p1[i]);
        }

        s = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            s += Math.Pow(p[i] - pn[i], 2);
        }

        s = Math.Sqrt(s);
    }

    public static double segment_point_dist_2d(double[] p1, double[] p2, double[] p, int p1Index = 0, int p2Index = 0, int pIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_DIST_2D: distance ( line segment, point ) in 2D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    The nearest point will satisfy the condition
        //
        //      PN = (1-T) * P1 + T * P2.
        //
        //    T will always be between 0 and 1.
        //
        //    Thanks to Kirill Speransky for pointing out that a previous version
        //    of this routine was incorrect, 02 May 2006.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the endpoints of the line segment.
        //
        //    Input, double P[2], the point whose nearest neighbor on the line
        //    segment is to be determined.
        //
        //    Output, double SEGMENT_POINT_DIST_2D, the distance from the point
        //    to the line segment.
        //
    {
        int DIM_NUM = 2;

        double bot;
        double dist;
        int i;
        double t;
        double[] pn = new double[DIM_NUM];
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2, p1Index, p2Index))
        {
            t = 0.0;
        }
        else
        {
            bot = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                bot += Math.Pow(p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length], 2);
            }

            t = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                t += (p[(i + pIndex) % p.Length] - p1[(i + p1Index) % p1.Length]) * (p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length]);
            }

            t /= bot;
            t = Math.Max(t, 0.0);
            t = Math.Min(t, 1.0);
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[(i + p1Index) % p1.Length] + t * (p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length]);
        }

        dist = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dist += Math.Pow(p[(i + pIndex) % p.Length] - pn[i], 2);
        }

        dist = Math.Sqrt(dist);

        return dist;
    }

    public static double segment_point_dist_3d(double[] p1, double[] p2, double[] p, int p1Index = 0, int p2Index = 0, int pIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_DIST_3D: distance ( line segment, point ) in 3D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    Thanks to Kirill Speransky for pointing out that a previous version
        //    of this routine was incorrect, 02 May 2006.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the endpoints of the line segment.
        //
        //    Input, double P[3], the point whose nearest neighbor on the line
        //    segment is to be determined.
        //
        //    Output, double SEGMENT_POINT_DIST_3D, the distance from the point
        //    to the line segment.
        //
    {
        int DIM_NUM = 3;

        double bot;
        double dist;
        int i;
        double t;
        double[] pn = new double[DIM_NUM];
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2, p1Index, p2Index))
        {
            t = 0.0;
        }
        else
        {
            bot = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                bot += Math.Pow(p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length], 2);
            }

            t = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                t += (p[(i + pIndex) % p.Length] - p1[(i + p1Index) % p1.Length]) * (p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length]);
            }

            t /= bot;
            t = Math.Max(t, 0.0);
            t = Math.Min(t, 1.0);
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[(i + p1Index) % p1.Length] + t * (p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length]);
        }

        dist = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dist += Math.Pow(p[(i + pIndex) % p.Length] - pn[i], 2);
        }

        dist = Math.Sqrt(dist);

        return dist;
    }

    public static void segment_point_near_2d(double[] p1, double[] p2, double[] p,
            ref double[] pn, ref double dist, ref double t, int p1Index = 0, int p2Index = 0, int pIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_NEAR_2D finds the point on a line segment nearest a point in 2D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    The nearest point will satisfy the condition:
        //
        //      PN = (1-T) * P1 + T * P2.
        //
        //    and T will always be between 0 and 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the two endpoints of the line segment.
        //
        //    Input, double P[2], the point whose nearest neighbor
        //    on the line segment is to be determined.
        //
        //    Output, double PN[2], the point on the line segment which is nearest P.
        //
        //    Output, double *DIST, the distance from the point to the nearest point
        //    on the line segment.
        //
        //    Output, double *T, the relative position of the point Pn to the
        //    points P1 and P2.
        //
    {
        int DIM_NUM = 2;

        double bot;
        int i;
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2, p1Index, p2Index))
        {
            t = 0.0;
        }
        else
        {
            bot = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                bot += Math.Pow(p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length], 2);
            }

            t = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                t += (p[(i + pIndex) % p.Length] - p1[(i + p1Index) % p1.Length]) * (p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length]);
            }

            t /= bot;
            t = Math.Max(t, 0.0);
            t = Math.Min(t, 1.0);
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[(i + p1Index) % p1.Length] + t * (p2[(i + p2Index) % p2.Length] - p1[(i + p1Index) % p1.Length]);
        }

        dist = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dist += Math.Pow(p[(i + pIndex) % p.Length] - pn[i], 2);
        }

        dist = Math.Sqrt(dist);

    }

    public static void segment_point_near_3d(double[] p1, double[] p2, double[] p,
            ref double[] pn, ref double dist, ref double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_NEAR_3D finds the point on a line segment nearest a point in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the two endpoints of the line segment.
        //
        //    Input, double P[3], the point whose nearest neighbor
        //    on the line segment is to be determined.
        //
        //    Output, double PN[3], the point on the line segment which is nearest to P.
        //
        //    Output, double *DIST, the distance from the point to the nearest point
        //    on the line segment.
        //
        //    Output, double *T, the relative position of the nearest point
        //    PN to the defining points P1 and P2.
        //
        //      PN = (1-T)*P1 + T*P2.
        //
        //    T will always be between 0 and 1.
        //
        //
    {
        int DIM_NUM = 3;

        double bot;
        int i;
        //
        //  If the line segment is actually a point, then the answer is easy.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2))
        {
            t = 0.0;
        }
        else
        {
            bot = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                bot += Math.Pow(p2[i] - p1[i], 2);
            }

            t = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                t += (p[i] - p1[i]) * (p2[i] - p1[i]);
            }

            t /= bot;
            t = Math.Max(t, 0.0);
            t = Math.Min(t, 1.0);
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[i] + t * (p2[i] - p1[i]);
        }

        dist = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dist += Math.Pow(p[i] - pn[i], 2);
        }

        dist = Math.Sqrt(dist);

    }

    public static double segments_curvature_2d(double[] p1, double[] p2, double[] p3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENTS_CURVATURE_2D computes the curvature of two line segments in 2D.
        //
        //  Discussion:
        //
        //    We assume that the segments [P1,P2] and [P2,P3] are given.
        //
        //    We compute the circle that passes through P1, P2 and P3.
        //
        //    The inverse of the radius of this circle is the local "curvature".
        //
        //    If curvature is 0, the two line segments have the same slope,
        //    and the three points are collinear.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], P3[2], the points.
        //
        //    Output, double SEGMENTS_CURVATURE_2D, the local curvature.
        //
    {
        int DIM_NUM = 2;

        double curvature;
        double[] pc = new double[DIM_NUM];
        double r = 0;

        CircleNS.Geometry.circle_exp2imp_2d(p1, p2, p3, ref r, ref pc);

        curvature = r switch
        {
            > 0.0 => 1.0 / r,
            _ => 0.0
        };

        return curvature;
    }

    public static double segments_dist_2d(double[] p1, double[] p2, double[] q1,
            double[] q2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENTS_DIST_2D computes the distance between two line segments in 2D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    If the lines through [P1,P2] and [Q1,Q2] intersect, and both
        //    line segments include the point of intersection, then the distance
        //    is zero and we are done.
        //
        //    Therefore, we compute the intersection of the two lines, and
        //    find the coordinates of that intersection point on each line.
        //    This will tell us if the zero distance case has occurred.
        //
        //    Otherwise, let PN and QN be points in [P1,P2] and [Q1,Q2] for which
        //    the distance is minimal.  If the lines do not intersect, then it
        //    cannot be the case that both PN and QN are strictly interior to their
        //    line segments, aside from the exceptional singular case when
        //    the line segments overlap or are parallel.  Even then, one of PN
        //    and QN may be taken to be a segment endpoint.
        //
        //    Therefore, our second computation finds the minimum of:
        //
        //      Distance ( P1, [Q1,Q2] );
        //      Distance ( P2, [Q1,Q2] );
        //      Distance ( Q1, [P1,P2] );
        //      Distance ( Q2, [P1,P2] );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the endpoints of the first segment.
        //
        //    Input, double Q1[2], Q2[2], the endpoints of the second segment.
        //
        //    Output, double SEGMENTS_DIST_2D, the distance between the line segments.
        //
    {
        int DIM_NUM = 2;

        double dist;
        double dist2;
        int ival = 0;
        double[] r = new double[DIM_NUM];
        double rps = 0;
        double rpt = 0;
        double rqs = 0;
        double rqt = 0;
        //
        //  Determine whether and where the underlying lines intersect.
        //
        LineNS.Geometry.lines_exp_int_2d(p1, p2, q1, q2, ref ival, ref r);
        switch (ival)
        {
            //
            //  If there is exactly one intersection point part of both lines,
            //  check that it is part of both line segments.
            //
            case 1:
            {
                segment_point_coords_2d(p1, p2, r, ref rps, ref rpt);
                segment_point_coords_2d(q1, q2, r, ref rqs, ref rqt);

                switch (rpt)
                {
                    case >= 0.0 and <= 1.0 when 0.0 <= rqt && rqt <= 1.0:
                        dist = 0.0;
                        return dist;
                }

                break;
            }
        }

        //
        //  If there is no intersection, or the intersection point is
        //  not part of both line segments, then an endpoint of one
        //  line segment achieves the minimum distance.
        //
        dist2 = segment_point_dist_2d(q1, q2, p1);
        dist = dist2;
        dist2 = segment_point_dist_2d(q1, q2, p2);
        dist = Math.Min(dist, dist2);
        dist2 = segment_point_dist_2d(p1, p2, q1);
        dist = Math.Min(dist, dist2);
        dist2 = segment_point_dist_2d(p1, p2, q2);
        dist = Math.Min(dist, dist2);

        return dist;
    }

    public static double segments_dist_3d(double[] p1, double[] p2, double[] q1,
            double[] q2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENTS_DIST_3D computes the distance between two line segments in 3D.
        //
        //  Discussion:
        //
        //
        //    NOTE: The special cases for identical and parallel lines have not been
        //    worked out yet; those cases are exceptional, and so this code
        //    is made available in a slightly unfinished form!
        //
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points P1 and P2.
        //
        //    Given two line segments, consider the underlying lines on which
        //    they lie.
        //
        //    A) If the lines are identical, then the distance between the line segments
        //    is 0, if the segments overlap, or otherwise is attained by the
        //    minimum of the distances between each endpoint and the opposing
        //    line segment.
        //
        //    B) If the lines are parallel, then the distance is either the distance
        //    between the lines, if the projection of one line segment onto
        //    the other overlaps, or otherwise is attained by the
        //    minimum of the distances between each endpoint and the opposing
        //    line segment.
        //
        //    C) If the lines are not identical, and not parallel, then there are
        //    unique points PN and QN which are the closest pair of points on the lines.
        //    If PN is interior to [P1,P2] and QN is interior to [Q1,Q2],
        //    then the distance between the two line segments is the distance
        //    between PN and QN.  Otherwise, the nearest distance can be computed
        //    by taking the minimum of the distance from each endpoing to the
        //    opposing line segment.
        //
        //    Therefore, our computation first checks whether the lines are
        //    identical, parallel, or other, and checks for the special case
        //    where the minimum occurs in the interior.
        //
        //    If that case is ruled out, it computes and returns the minimum of:
        //
        //      Distance ( P1, [Q1,Q2] );
        //      Distance ( P2, [Q1,Q2] );
        //      Distance ( Q1, [P1,P2] );
        //      Distance ( Q2, [P1,P2] );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the endpoints of the first
        //    segment.
        //
        //    Input, double Q1[3], Q2[3], the endpoints of the second
        //    segment.
        //
        //    Output, double SEGMENTS_DIST_3D, the distance between the line segments.
        //
    {
        int DIM_NUM = 3;

        double a;
        double b;
        double c;
        double d;
        double det;
        double dist;
        double dist2;
        double e;
        int i;
        double[] pn = new double[DIM_NUM];
        double[] qn = new double[DIM_NUM];
        double sn;
        double tn;
        double[] u = new double[DIM_NUM];
        double[] v = new double[DIM_NUM];
        double[] w0 = new double[DIM_NUM];
        //
        //  The lines are identical.
        //  THIS CASE NOT SET UP YET
        //
        // if ( lines_exp_equal_3d ( p1, p2, q1, q2 ) ) then
        // end if
        //
        //  The lines are not identical, but parallel
        //  THIS CASE NOT SET UP YET.
        //
        // if ( lines_exp_parallel_3d ( p1, p2, q1, q2 ) ) then
        // end if
        //
        //  C: The lines are not identical, not parallel.
        //

        //
        //  Let U = (P2-P1) and V = (Q2-Q1) be the direction vectors on
        //  the two lines.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            u[i] = p2[i] - p1[i];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v[i] = q2[i] - q1[i];
        }

        //
        //  Let SN be the unknown coordinate of the nearest point PN on line 1,
        //  so that PN = P(SN) = P1 + SN * (P2-P1).
        //
        //  Let TN be the unknown coordinate of the nearest point QN on line 2,
        //  so that QN = Q(TN) = Q1 + TN * (Q2-Q1).
        //
        //  Let W0 = (P1-Q1).
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            w0[i] = p1[i] - q1[i];
        }

        //
        //  The vector direction WC = P(SN) - Q(TC) is unique (among directions)
        //  perpendicular to both U and V, so
        //
        //    U dot WC = 0
        //    V dot WC = 0
        //
        //  or, equivalently:
        //
        //    U dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
        //    V dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
        //
        //  or, equivalently:
        //
        //    (u dot u ) * sn - (u dot v ) tc = -u * w0
        //    (v dot u ) * sn - (v dot v ) tc = -v * w0
        //
        //  or, equivalently:
        //
        //   ( a  -b ) * ( sn ) = ( -d )
        //   ( b  -c )   ( tc )   ( -e )
        //
        a = typeMethods.r8vec_dot_product(DIM_NUM, u, u);
        b = typeMethods.r8vec_dot_product(DIM_NUM, u, v);
        c = typeMethods.r8vec_dot_product(DIM_NUM, v, v);
        d = typeMethods.r8vec_dot_product(DIM_NUM, u, w0);
        e = typeMethods.r8vec_dot_product(DIM_NUM, v, w0);
        //
        //  Check the determinant.
        //
        det = -a * c + b * b;

        switch (det)
        {
            case 0.0:
            {
                sn = 0.0;
                if (Math.Abs(b) < Math.Abs(c))
                {
                    tn = e / c;
                }
                else
                {
                    tn = d / b;
                }

                break;
            }
            default:
                sn = (c * d - b * e) / det;
                tn = (b * d - a * e) / det;
                break;
        }

        switch (sn)
        {
            //
            //  Now if both nearest points on the lines
            //  also happen to lie inside their line segments,
            //  then we have found the nearest points on the line segments.
            //
            case >= 0.0 and <= 1.0 when 0.0 <= tn && tn <= 1.0:
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    pn[i] = p1[i] + sn * (p2[i] - p1[i]);
                }

                for (i = 0; i < DIM_NUM; i++)
                {
                    qn[i] = q1[i] + tn * (q2[i] - q1[i]);
                }

                dist = 0.0;
                for (i = 0; i < DIM_NUM; i++)
                {
                    dist += Math.Pow(pn[i] - qn[i], 2);
                }

                dist = Math.Sqrt(dist);

                return dist;
            }
        }

        //
        //  The nearest point did not occur in the interior.
        //  Therefore it must be achieved at an endpoint.
        //
        dist2 = segment_point_dist_3d(q1, q2, p1);
        dist = dist2;
        dist2 = segment_point_dist_3d(q1, q2, p2);
        dist = Math.Min(dist, dist2);
        dist2 = segment_point_dist_3d(p1, p2, q1);
        dist = Math.Min(dist, dist2);
        dist2 = segment_point_dist_3d(p1, p2, q2);
        dist = Math.Min(dist, dist2);

        return dist;
    }

    public static double segments_dist_3d_old(double[] p1, double[] p2, double[] p3,
            double[] p4)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENTS_DIST_3D_OLD computes the distance between two line segments in 3D.
        //
        //  Discussion:
        //
        //    A line segment is the portion of an infinite line that lies between
        //    two given points.  The behavior of the distance function is a bit
        //    complicated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 November 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the endpoints of the first segment.
        //
        //    Input, double P3[3], P4[3], the endpoints of the second segment.
        //
        //    Output, double SEGMENTS_DIST_3D, the distance between the line segments.
        //
    {
        int DIM_NUM = 3;

        double d1 = 0;
        double d2 = 0;
        double dist = 0;
        double dl;
        double dm;
        double dr;
        double[] pn1 = new double[DIM_NUM];
        double[] pn2 = new double[DIM_NUM];
        double[] pt = new double[DIM_NUM];
        bool result;
        double t1 = 0;
        double t2 = 0;
        double tl;
        double tm;
        double tmin = 0;
        double tr;
        //
        //  Find the nearest points on line 2 to the endpoints of line 1.
        //
        segment_point_near_3d(p3, p4, p1, ref pn1, ref d1, ref t1);

        segment_point_near_3d(p3, p4, p2, ref pn2, ref d2, ref t2);

        if (Math.Abs(t1 - t2) <= double.Epsilon)
        {
            dist = segment_point_dist_3d(p1, p2, pn1);
            return dist;
        }

        //
        //  On line 2, over the interval between the points nearest to line 1,
        //  the square of the distance of any point to line 1 is a quadratic function.
        //  Evaluate it at three points, and seek its local minimum.
        //
        dl = segment_point_dist_3d(p1, p2, pn1);

        pt[0] = 0.5 * (pn1[0] + pn2[0]);
        pt[1] = 0.5 * (pn1[1] + pn2[1]);
        pt[2] = 0.5 * (pn1[2] + pn2[2]);

        dm = segment_point_dist_3d(p1, p2, pt);

        dr = segment_point_dist_3d(p1, p2, pn2);

        tl = 0.0;
        tm = 0.5;
        tr = 1.0;

        dl *= dl;
        dm *= dm;
        dr *= dr;

        result = LocalMinimum.minquad(tl, dl, tm, dm, tr, dr, ref tmin, ref dist);

        switch (result)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("SEGMENTS_DIST_3D - Fatal error!");
                Console.WriteLine("  MINQUAD returned error condition.");
                return 1;
            default:
                dist = Math.Sqrt(dist);

                return dist;
        }
    }

    public static double segments_int_1d(double p1, double p2, double q1, double q2,
            ref double r1, ref double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENTS_INT_1D computes the intersection of two line segments in 1D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    In 1D, two line segments "intersect" if they overlap.
        //
        //    Using a real number DIST to report overlap is preferable to
        //    returning a TRUE/FALSE flag, since DIST is better able to
        //    handle cases where the segments "almost" interlap.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1, P2, the endpoints of the first segment.
        //
        //    Input, double Q1, Q2, the endpoints of the second segment.
        //
        //    Output, double *R1, *R2, the endpoints of the intersection
        //    segment.
        //    If DIST < 0, then the interval [R1,R2] is the common intersection
        //    of the two segments.
        //    If DIST = 0, then R1 = R2 is the single common point of the two segments.
        //    If DIST > 0, then (R1,R2) is an open interval separating the two
        //    segments, which do not overlap at all.
        //
        //    Output, double SEGMENTS_INT_1D, the "distance" DIST between the segments.
        //    < 0, the segments overlap, and the overlap is DIST units long;
        //    = 0, the segments overlap at a single point;
        //    > 0, the segments do not overlap.  The distance between the nearest
        //    points is DIST units.
        //
    {
        double dist;

        r1 = Math.Max(Math.Min(p1, p2),
            Math.Min(q1, q2));

        r2 = Math.Min(Math.Max(p1, p2),
            Math.Max(q1, q2));

        dist = r1 - r2;

        return dist;
    }

    public static void segments_int_2d(double[] p1, double[] p2, double[] p3,
            double[] p4, ref int flag, ref double[] p5)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENTS_INT_2D computes the intersection of two line segments in 2D.
        //
        //  Discussion:
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //    In 2D, two line segments might not intersect, even though the
        //    lines, of which they are portions, intersect.
        //
        //    Thanks to Siavosh Bahrami for pointing out an error involving incorrect
        //    indexing of the U array, 17 August 2005.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the endpoints of the first segment.
        //
        //    Input, double P3[2], P4[2], the endpoints of the second segment.
        //
        //    Output, int *FLAG, records the results.
        //    0, the line segments do not intersect.
        //    1, the line segments intersect.
        //
        //    Output, double *P5[2].
        //    If FLAG = 0, P5 = 0.
        //    If FLAG = 1, then P5 is a point of intersection.
        //
    {
        int DIM_NUM = 2;

        int ival = 0;
        double tol = 0.001;
        double[] u = new double[DIM_NUM];
        //
        //  Find the intersection of the two lines.
        //
        LineNS.Geometry.lines_exp_int_2d(p1, p2, p3, p4, ref ival, ref p5);

        switch (ival)
        {
            case 0:
                flag = 0;
                p5[0] = 0.0;
                p5[1] = 0.0;
                return;
        }

        //
        //  Is the intersection point on the first line segment?
        //
        segment_contains_point_2d(p1, p2, p5, ref u);

        if (u[0] < 0.0 || 1.0 < u[0] || tol < u[1])
        {
            flag = 0;
            p5[0] = 0.0;
            p5[1] = 0.0;
            return;
        }

        //
        //  Is the intersection point on the second line segment?
        //
        segment_contains_point_2d(p3, p4, p5, ref u);

        if (u[0] < 0.0 || 1.0 < u[0] || tol < u[1])
        {
            flag = 0;
            p5[0] = 0.0;
            p5[1] = 0.0;
            return;
        }

        flag = 1;

    }

    public static void string_2d(int vec_num, double[] p1, double[] p2, ref int string_num,
            ref int[] order, ref int[] string_)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRING_2D groups line segments into connected lines in 2D.
        //
        //  Discussion:
        //
        //    The routine receives an unordered set of line segments, described by
        //    pairs of coordinates P1 and P2, and tries to group them
        //    into ordered lists that constitute connected jagged lines.
        //
        //    This routine will not match two endpoints unless they are exactly equal.
        //
        //    On input, line segment I has endpoints PI(I) and P2(I).
        //
        //    On output, the order of the components may have been
        //    switched.  That is, for some I, P1(I) and P2(I) may have been swapped.
        //
        //    More importantly, all the entries P1(I) and P2(I)
        //    may have been swapped with another index J.
        //
        //    The resulting coordinates will have been sorted in order
        //    of the string to which they belong, and then by the order
        //    of their traversal within that string.
        //
        //    The array STRING(I) identifies the string to which segment I belongs.
        //
        //    If two segments I and J have the same value of STRING, then
        //    ORDER(I) and ORDER(J) give the relative order of the two segments
        //    in the string.  Thus if ORDER(I) = -3 and ORDER(J) = 2, then when
        //    the string is traversed, segment I is traversed first, then four other
        //    segments are traversed, and then segment J is traversed.
        //
        //    For each string, the segment with ORDER(I) = 0 is the initial segment
        //    from which the entire string was "grown" (with growth possible to both the
        //    left and the right).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int VEC_NUM, the number of line segments to be analyzed.
        //
        //    Input/output, double P1[2*VEC_NUM], P2[2*VEC_NUM], the line segments.
        //
        //    Output, int *STRING_NUM, the number of strings created.
        //
        //    Output, int ORDER[VEC_NUM], the order vector.
        //
        //    Output, int STRING[VEC_NUM], the string to which each segment I belongs.
        //
    {
        int i;
        int indx;
        int isgn;
        int itemp;
        int j;
        int jval;
        int kval;
        int match;
        int seed;
        double temp;
        double x1val;
        double x2val;
        double y1val;
        double y2val;
        //
        //  Mark STRING so that each segment is alone.
        //
        for (i = 0; i < vec_num; i++)
        {
            order[i] = 0;
            string_[i] = vec_num + i + 1;
        }

        //
        //  Starting with the lowest numbered group of line segments,
        //  see if any higher numbered groups belong.
        //
        seed = 0;
        string_num = 1;
        string_[seed] = string_num;

        for (;;)
        {
            x1val = p1[0 + seed * 2];
            x2val = p2[0 + seed * 2];
            y1val = p1[1 + seed * 2];
            y2val = p2[1 + seed * 2];
            jval = order[seed];
            kval = order[seed];

            for (;;)
            {
                match = 0;

                for (j = 0; j < vec_num; j++)
                {
                    if (string_num < string_[j])
                    {
                        if (Math.Abs(x1val - p1[0 + j * 2]) <= double.Epsilon && Math.Abs(y1val - p1[1 + j * 2]) <= double.Epsilon)
                        {
                            jval -= 1;
                            order[j] = jval;
                            string_[j] = string_num;
                            x1val = p2[0 + j * 2];
                            y1val = p2[1 + j * 2];
                            match += 1;

                            temp = p1[0 + j * 2];
                            p1[0 + j * 2] = p2[0 + j * 2];
                            p2[0 + j * 2] = temp;

                            temp = p1[1 + j * 2];
                            p1[1 + j * 2] = p2[1 + j * 2];
                            p2[1 + j * 2] = temp;
                        }
                        else if (Math.Abs(x1val - p2[0 + j * 2]) <= double.Epsilon && Math.Abs(y1val - p2[1 + j * 2]) <= double.Epsilon)
                        {
                            jval -= 1;
                            order[j] = jval;
                            string_[j] = string_num;
                            x1val = p1[0 + j * 2];
                            y1val = p1[1 + j * 2];
                            match += 1;
                        }
                        else if (Math.Abs(x2val - p1[0 + j * 2]) <= double.Epsilon && Math.Abs(y2val - p1[1 + j * 2]) <= double.Epsilon)
                        {
                            kval += 1;
                            order[j] = kval;
                            string_[j] = string_num;
                            x2val = p2[0 + j * 2];
                            y2val = p2[1 + j * 2];
                            match += 1;
                        }
                        else if (Math.Abs(x2val - p2[0 + j * 2]) <= double.Epsilon && Math.Abs(y2val - p2[1 + j * 2]) <= double.Epsilon)
                        {
                            kval += 1;
                            order[j] = kval;
                            string_[j] = string_num;
                            x2val = p1[0 + j * 2];
                            y2val = p1[1 + j * 2];
                            match += 1;

                            temp = p1[0 + j * 2];
                            p1[0 + j * 2] = p2[0 + j * 2];
                            p2[0 + j * 2] = temp;

                            temp = p1[1 + j * 2];
                            p1[1 + j * 2] = p2[1 + j * 2];
                            p2[1 + j * 2] = temp;
                        }
                    }
                }

                //
                //  If the string has closed on itself, then we don't want to
                //  look for any more matches for this string.
                //
                if (Math.Abs(x1val - x2val) <= double.Epsilon && Math.Abs(y1val - y2val) <= double.Epsilon)
                {
                    break;
                }

                //
                //  If we made no matches this pass, we're done.
                //
                if (match <= 0)
                {
                    break;
                }
            }

            //
            //  This string is "exhausted".  Are there any line segments we
            //  haven't looked at yet?
            //
            seed = 0;

            for (i = 0; i < vec_num; i++)
            {
                if (string_num < string_[i])
                {
                    seed = i;
                    string_num += 1;
                    string_[i] = string_num;
                    break;
                }
            }

            if (seed == 0)
            {
                break;
            }
        }

        //
        //  There are no more line segments to look at.  Renumber the
        //  isolated segments.
        //
        //  Question: Can this ever happen?
        //
        for (i = 0; i < vec_num; i++)
        {
            if (vec_num < string_[i])
            {
                string_num += 1;
                string_[i] = string_num;
            }
        }

        //
        //  Now sort the line segments by string and by order of traversal.
        //
        i = 0;
        isgn = 0;
        j = 0;

        indx = 0;

        SortHeapExternalData data = new();

        for (;;)
        {
            Sort.sort_heap_external(ref data, vec_num, ref indx, ref i, ref j, isgn);

            if (0 < indx)
            {
                itemp = order[i - 1];
                order[i - 1] = order[j - 1];
                order[j - 1] = itemp;

                itemp = string_[i - 1];
                string_[i - 1] = string_[j - 1];
                string_[j - 1] = itemp;

                temp = p1[0 + (i - 1) * 2];
                p1[0 + (i - 1) * 2] = p1[0 + (j - 1) * 2];
                p1[0 + (j - 1) * 2] = temp;

                temp = p1[1 + (i - 1) * 2];
                p1[1 + (i - 1) * 2] = p1[1 + (j - 1) * 2];
                p1[1 + (j - 1) * 2] = temp;

                temp = p2[0 + (i - 1) * 2];
                p2[0 + (i - 1) * 2] = p2[0 + (j - 1) * 2];
                p2[0 + (j - 1) * 2] = temp;

                temp = p2[1 + (i - 1) * 2];
                p2[1 + (i - 1) * 2] = p2[1 + (j - 1) * 2];
                p2[1 + (j - 1) * 2] = temp;
            }
            else if (indx < 0)
            {
                if (string_[i - 1] < string_[j - 1] ||
                    string_[i - 1] == string_[j - 1] && order[i - 1] < order[j - 1])
                {
                    isgn = -1;
                }
                else
                {
                    isgn = +1;
                }
            }
            else
            {
                break;
            }
        }
    }


}