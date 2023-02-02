using System;
using Burkardt.Types;

namespace Burkardt.LineNS;

public static class Geometry
{
    public static bool line_exp_is_degenerate_nd(int dim_num, double[] p1, double[] p2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
        //
        //  Discussion:
        //
        //    The explicit form of a line in ND is:
        //
        //      the line through the points P1 and P2.
        //
        //    An explicit line is degenerate if the two defining points are equal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double P1[DIM_NUM], P2[DIM_NUM], two points on the line.
        //
        //    Output, bool LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
        //    is degenerate.
        //
    {
        bool value = typeMethods.r8vec_eq(dim_num, p1, p2);

        return value;
    }

    public static double[] line_exp_normal_2d(double[] p1, double[] p2, int p1Index = 0, int p2Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_NORMAL_2D computes the unit normal vector to a line in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through the points P1 and P2.
        //
        //    The sign of the normal vector N is chosen so that the normal vector
        //    points "to the left" of the direction of the line.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], two distinct points on the line.
        //
        //    Output, double LINE_EXP_NORMAL_2D[2], a unit normal vector to the line.
        //
    {
        const int DIM_NUM = 2;

        double[] normal = new double[DIM_NUM];

        double norm = Math.Sqrt((p2[(0 + p2Index) % p2.Length] - p1[(0 + p1Index) % p1.Length]) * (p2[(0 + p2Index) % p2.Length] - p1[(0 + p1Index) % p1.Length])
                                + (p2[(1 + p2Index) % p2.Length] - p1[(1 + p1Index) % p1.Length]) * (p2[(1 + p2Index) % p2.Length] - p1[(1 + p1Index) % p1.Length]));

        switch (norm)
        {
            case 0.0:
                normal[0] = Math.Sqrt(2.0);
                normal[1] = Math.Sqrt(2.0);
                break;
            default:
                normal[0] = -(p2[(1 + p2Index) % p2.Length] - p1[(1 + p1Index) % p1.Length]) / norm;
                normal[1] = (p2[(0 + p2Index) % p2.Length] - p1[(0 + p1Index) % p1.Length]) / norm;
                break;
        }

        return normal;
    }

    public static double[] line_exp_perp_2d(double[] p1, double[] p2, double[] p3,
            ref bool flag, int p1Index = 0, int p2Index = 0, int p3Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_PERP_2D computes a line perpendicular to a line and through a point.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], two points on the given line.
        //
        //    Input, double P3[2], a point not on the given line, through which the
        //    perpendicular must pass.
        //
        //    Output, bool *FLAG, is TRUE if the value could not be computed.
        //
        //    Output, double LINE_EXP_PERP_2D[2], a point on the given line, such that
        //    the line through P3 and P4 is perpendicular to the given line.
        //
    {
        const int DIM_NUM = 2;

        double[] p4 = new double[DIM_NUM];
        flag = false;

        double bot = Math.Pow(p2[(0 + p2Index) % p2.Length] - p1[(0 + p1Index) % p1.Length], 2) + Math.Pow(p2[(1 + p2Index) % p2.Length] - p1[(1 + p1Index) % p1.Length], 2);

        switch (bot)
        {
            case 0.0:
                flag = true;
                p4[0] = typeMethods.r8_huge();
                p4[1] = typeMethods.r8_huge();
                return p4;
        }

        //
        //  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
        //
        //  (P3-P1) dot (P2-P1) / Norm(P3-P1)^2 = normalized coordinate T
        //  of the projection of (P3-P1) onto (P2-P1).
        //
        double t = ((p1[(0 + p1Index) % p1.Length] - p3[(0 + p3Index) % p3.Length]) * (p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length])
                    + (p1[(1 + p1Index) % p1.Length] - p3[(1 + p3Index) % p3.Length]) * (p1[(1 + p1Index) % p1.Length] - p2[(1 + p2Index) % p2.Length])) / bot;

        p4[0] = p1[(0 + p1Index) % p1.Length] + t * (p2[(0 + p2Index) % p2.Length] - p1[(0 + p1Index) % p1.Length]);
        p4[1] = p1[(1 + p1Index) % p1.Length] + t * (p2[(1 + p2Index) % p2.Length] - p1[(1 + p1Index) % p1.Length]);

        return p4;
    }

    public static double line_exp_point_dist_2d(double[] p1, double[] p2, double[] p, int pIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], two points on the line.
        //
        //    Input, double P[2], the point whose distance from the line is
        //    to be measured.
        //
        //    Output, double LINE_EXP_DIST_2D, the distance from the point to the line.
        //
    {
        double[] pn = new double[2];

        double bot = Math.Pow(p2[0] - p1[0], 2)
                     + Math.Pow(p2[1] - p1[1], 2);

        switch (bot)
        {
            case 0.0:
                pn[0] = p1[0];
                pn[1] = p1[1];
                break;
            //
            default:
                double dot = (p[(0 + pIndex) % p.Length] - p1[0]) * (p2[0] - p1[0])
                             + (p[(1 + pIndex) % p.Length] - p1[1]) * (p2[1] - p1[1]);

                double t = dot / bot;

                pn[0] = p1[0] + t * (p2[0] - p1[0]);
                pn[1] = p1[1] + t * (p2[1] - p1[1]);
                break;
        }

        double dist = Math.Sqrt(Math.Pow(p[(0 + pIndex) % p.Length] - pn[0], 2)
                                + Math.Pow(p[(1 + pIndex) % p.Length] - pn[1], 2));

        return dist;
    }

    public static double line_exp_point_dist_3d(double[] p1, double[] p2, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_POINT_DIST_3D: distance ( explicit line, point ) in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two points on a line.
        //
        //    Input, double P[3], the point whose distance from the line is
        //    to be measured.
        //
        //    Output, double LINE_EXP_POINT_DIST_3D, the distance from the point
        //    to the line.
        //
    {
        const int DIM_NUM = 3;

        double[] pn = new double[DIM_NUM];

        double bot = Math.Pow(p2[0] - p1[0], 2)
                     + Math.Pow(p2[1] - p1[1], 2)
                     + Math.Pow(p2[2] - p1[2], 2);

        switch (bot)
        {
            case 0.0:
                typeMethods.r8vec_copy(DIM_NUM, p1, ref pn);
                break;
            //
            default:
                double t = (
                    (p[0] - p1[0]) * (p2[0] - p1[0]) +
                    (p[1] - p1[1]) * (p2[1] - p1[1]) +
                    (p[2] - p1[2]) * (p2[2] - p1[2])) / bot;

                pn[0] = p1[0] + t * (p2[0] - p1[0]);
                pn[1] = p1[1] + t * (p2[1] - p1[1]);
                pn[2] = p1[2] + t * (p2[2] - p1[2]);
                break;
        }

        //
        //  Now compute the distance between the projection point and P.
        //
        double dist = Math.Sqrt(Math.Pow(p[0] - pn[0], 2)
                                + Math.Pow(p[1] - pn[1], 2)
                                + Math.Pow(p[2] - pn[2], 2));

        return dist;
    }

    public static double line_exp_point_dist_signed_2d(double[] p1, double[] p2, double[] p, int p1Index = 0, int p2Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( explicit line, point ) in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //    The signed distance has two interesting properties:
        //
        //    *  The absolute value of the signed distance is the
        //       usual (Euclidean) distance.
        //
        //    *  Points with signed distance 0 lie on the line,
        //       points with a negative signed distance lie on one side
        //         of the line,
        //       points with a positive signed distance lie on the
        //         other side of the line.
        //
        //    Assuming that C is nonnegative, then if a point is a positive
        //    distance away from the line, it is on the same side of the
        //    line as the point (0,0), and if it is a negative distance
        //    from the line, it is on the opposite side from (0,0).
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
        //    Input, double P1[2], P2[2], two points that determine the line.
        //
        //    Input, double P[2], the point whose signed distance is desired.
        //
        //    Output, double LINE_EXP_DIST_SIGNED_2D, the signed distance from the
        //    point to the line.
        //
    {
        double a = 0;
        double b = 0;
        double c = 0;
        //
        //  Convert the line to A*x+B*y+C form.
        //
        line_exp2imp_2d(p1, p2, ref a, ref b, ref c, p1Index, p2Index);
        //
        //  Compute the signed distance from the point to the line.
        //
        double dist_signed = (a * p[0] + b * p[1] + c) / Math.Sqrt(a * a + b * b);

        return dist_signed;
    }

    public static void line_exp_point_near_2d(double[] p1, double[] p2, double[] p,
            ref double[] pn, ref double dist, ref double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_POINT_NEAR_2D computes the point on an explicit line nearest a point in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //    The nearest point PN will have the form:
        //
        //      PN = (1-T) * P1 + T * P2.
        //
        //    If T is less than 0, then PN is furthest away from P2.
        //    If T is between 0 and 1, PN is between P1 and P2.
        //    If T is greater than 1, PN is furthest away from P1.
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
        //    Input, double P1[2], P2[2], two points that define a line.
        //
        //    Input, double P[2], the point whose nearest neighbor on the line is
        //    to be determined.
        //
        //    Output, double PN[2], the nearest point on the line to P.
        //
        //    Output, double DIST, the distance from the point to the line.
        //
        //    Output, double[] T, the relative position of the point PN to the points P1 and P2.
        //
    {
        double bot = Math.Pow(p2[0] - p1[0], 2) + Math.Pow(p2[1] - p1[1], 2);

        switch (bot)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("LINE_EXP_POINT_NEAR_2D - Fatal error!");
                Console.WriteLine("  The points P1 and P2 are identical.");
                return;
        }

        //
        //  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
        //
        //  (P-P1) dot (P2-P1) / Norm(P-P1)^2 = normalized coordinate T
        //  of the projection of (P-P1) onto (P2-P1).
        //
        t = ((p1[0] - p[0]) * (p1[0] - p2[0])
             + (p1[1] - p[1]) * (p1[1] - p2[1])) / bot;

        pn[0] = p1[0] + t * (p2[0] - p1[0]);
        pn[1] = p1[1] + t * (p2[1] - p1[1]);

        dist = Math.Sqrt(Math.Pow(p[0] - pn[0], 2) + Math.Pow(p[1] - pn[1], 2));

    }

    public static void line_exp_point_near_3d(double[] p1, double[] p2, double[] p,
            ref double[] pn, ref double dist, ref double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP_POINT_NEAR_3D: nearest point on explicit line to point in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //    The nearest point PN will have the form:
        //
        //      PN = (1-T) * P1 + T * P2.
        //
        //    If T is less than 0, then PN is furthest away from P2.
        //    If T is between 0 and 1, PN is between P1 and P2.
        //    If T is greater than 1, PN is furthest away from P1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two points that define a line.
        //
        //    Input, double P[3], the point whose nearest neighbor on the line is
        //    to be determined.
        //
        //    Output, double PN[3], the nearest point on the line to P.
        //
        //    Output, double DIST, the distance from the point to the line.
        //
        //    Output, double[] T, the relative position of the point PN to the points P1 and P2.
        //
    {
        double bot = Math.Pow(p2[0] - p1[0], 2)
                     + Math.Pow(p2[1] - p1[1], 2)
                     + Math.Pow(p2[2] - p1[2], 2);

        switch (bot)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("LINE_EXP_POINT_NEAR_3D - Fatal error!");
                Console.WriteLine("  The points P1 and P2 are identical.");
                return;
        }

        //
        //  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
        //
        //  (P-P1) dot (P2-P1) / Norm(P-P1)^2 = normalized coordinate T
        //  of the projection of (P-P1) onto (P2-P1).
        //
        t = ((p1[0] - p[0]) * (p1[0] - p2[0])
             + (p1[1] - p[1]) * (p1[1] - p2[1])
             + (p1[2] - p[2]) * (p1[2] - p2[2])) / bot;
        //
        //  Now compute the location of the projection point.
        //
        pn[0] = p1[0] + t * (p2[0] - p1[0]);
        pn[1] = p1[1] + t * (p2[1] - p1[1]);
        pn[2] = p1[2] + t * (p2[2] - p1[2]);
        //
        //  Now compute the distance between the projection point and P.
        //
        dist = Math.Sqrt(Math.Pow(p[0] - pn[0], 2)
                         + Math.Pow(p[1] - pn[1], 2)
                         + Math.Pow(p[2] - pn[2], 2));
    }

    public static void line_exp2imp_2d(double[] p1, double[] p2, ref double a, ref double b,
            ref double c, int p1Index = 0, int p2Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
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
        //    Input, double P1[2], P2[2], two distinct points on the line.
        //
        //    Output, double[] A, *B, *C, three coefficients which describe
        //    the line that passes through P1 and P2.
        //
    {
        //
        //  Take care of degenerate cases.
        //
        if (typeMethods.r8vec_eq(2, p1, p2, p1Index, p2Index))
        {
            Console.WriteLine("");
            Console.WriteLine("LINE_EXP2IMP_2D - Fatal error!");
            Console.WriteLine("  P1 = P2");
            Console.WriteLine("  P1 = " + p1[(0 + p1Index) % p1.Length] + " " + p1[(1 + p1Index) % p1.Length] + "");
            Console.WriteLine("  P2 = " + p2[(0 + p2Index) % p2.Length] + " " + p2[(1 + p2Index) % p2.Length] + "");
            return;
        }

        a = p2[(1 + p2Index) % p2.Length] - p1[(1 + p1Index) % p1.Length];
        b = p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length];
        c = p2[(0 + p2Index) % p2.Length] * p1[(1 + p1Index) % p1.Length] - p1[(0 + p1Index) % p1.Length] * p2[(1 + p2Index) % p2.Length];

    }

    public static void line_exp2par_2d(double[] p1, double[] p2, ref double f, ref double g,
            ref double x0, ref double y0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP2PAR_2D converts a line from explicit to parametric form in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], two points on the line.
        //
        //    Output, double[] F, *G, *X0, *Y0, the parametric parameters of the line.
        //
    {
        x0 = p1[0];
        y0 = p1[1];

        double norm = Math.Sqrt((p2[0] - p1[0]) * (p2[0] - p1[0])
                                + (p2[1] - p1[1]) * (p2[1] - p1[1]));

        switch (norm)
        {
            case 0.0:
                f = 0.0;
                g = 0.0;
                break;
            default:
                f = (p2[0] - p1[0]) / norm;
                g = (p2[1] - p1[1]) / norm;
                break;
        }

        switch (f)
        {
            case < 0.0:
                f = -f;
                g = -g;
                break;
        }
    }

    public static void line_exp2par_3d(double[] p1, double[] p2, ref double f, ref double g,
            ref double h, ref double x0, ref double y0, ref double z0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_EXP2PAR_3D converts an explicit line into parametric form in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through P1 and P2.
        //
        //    The parametric form of a line in 3D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //      Z = Z0 + H * T
        //
        //    For normalization, we choose F*F+G*G+H*H = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two points on the line.
        //
        //    Output, double[] F, *G, *H, the components of the direction vector.
        //
        //    Output, double[] X0, *Y0, *Z0, the base vector.
        //
    {
        f = p2[0] - p1[0];
        g = p2[1] - p1[1];
        h = p2[2] - p1[2];

        double norm = Math.Sqrt(Math.Pow(f, 2) + Math.Pow(g, 2) + Math.Pow(h, 2));

        switch (norm)
        {
            case 0.0:
                f = 0.0;
                g = 0.0;
                h = 0.0;
                break;
            default:
                f /= norm;
                g /= norm;
                h /= norm;
                break;
        }

        switch (f)
        {
            case < 0.0:
                f = -f;
                g = -g;
                h = -h;
                break;
        }

        x0 = p1[0];
        y0 = p1[1];
        z0 = p1[2];

    }

    public static bool line_imp_is_degenerate_2d(double a, double b, double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
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
        //  Parameters:
        //
        //    Input, double A, B, C, the implicit line parameters.
        //
        //    Output, bool LINE_IMP_IS_DEGENERATE_2D, is true if the
        //    line is degenerate.
        //
    {
        bool value = a * a + b * b == 0.0;

        return value;
    }

    public static double line_imp_point_dist_2d(double a, double b, double c, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_IMP_POINT_DIST_2D: distance ( implicit line, point ) in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the implicit line parameters.
        //
        //    Input, double P[2], the point whose distance from the line is
        //    to be measured.
        //
        //    Output, double LINE_IMP_POINT_DIST_2D, the distance from the
        //    point to the line.
        //
    {
        switch (a * a + b * b)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("LINE_IMP_POINT_DIST_2D - Fatal error!");
                Console.WriteLine("  A * A + B * B = 0.");
                return 1;
            default:
                return Math.Abs(a * p[0] + b * p[1] + c) / Math.Sqrt(a * a + b * b);
        }
    }

    public static double line_imp_point_dist_signed_2d(double a, double b, double c,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_IMP_POINT_DIST_SIGNED_2D: signed distance ( implicit line, point ) in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C * Z + D = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the equation of the line is A*X + B*Y + C = 0.
        //
        //    Input, double P[2], the coordinates of the point.
        //
        //    Output, double LINE_IMP_POINT_DIST_SIGNED_2D, the signed distance
        //    from the point to the line.
        //
    {
        switch (a * a + b * b)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("LINE_IMP_POINT_DIST_SIGNED_2D - Fatal error!");
                Console.WriteLine("  A * A + B * B = 0.");
                return 1;
            default:
                double dist = -typeMethods.r8_sign(c) * (a * p[0] + b * p[1] + c) / Math.Sqrt(a * a + b * b);

                return dist;
        }
    }

    public static void line_imp2exp_2d(double a, double b, double c, ref double[] p1,
            ref double[] p2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_IMP2EXP_2D converts an implicit line to explicit form in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of line in 2D is:
        //
        //      A * X + B * Y + C = 0
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through the points P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2005
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
        //    Input, double A, B, C, the implicit line parameters.
        //
        //    Output, double P1[2], P2[2], two points on the line.
        //
    {
        if (line_imp_is_degenerate_2d(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("LINE_IMP2EXP_2D - Fatal error!");
            Console.WriteLine("  The line is degenerate.");
            return;
        }

        double normsq = a * a + b * b;

        p1[0] = -a * c / normsq;
        p1[1] = -b * c / normsq;

        if (Math.Abs(b) < Math.Abs(a))
        {
            p2[0] = -(a - b / a) * c / normsq;
            p2[1] = -(b + 1.0) * c / normsq;
        }
        else
        {
            p2[0] = -(a + 1.0) * c / normsq;
            p2[1] = -(b - a / b) * c / normsq;
        }

    }

    public static void line_imp2par_2d(double a, double b, double c, ref double f, ref double g,
            ref double x0, ref double y0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_IMP2PAR_2D converts an implicit line to parametric form in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of line in 2D is:
        //
        //      A * X + B * Y + C = 0
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2003
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
        //    Input, double A, B, C, the implicit parameters of the line.
        //
        //    Output, double[] F, *G, *X0, *Y0, the parametric parameters of the line.
        //
    {
        double test = a * a + b * b;

        switch (test)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("LINE_IMP2PAR_2D - Fatal error!");
                Console.WriteLine("  A * A + B * B = 0.");
                return;
        }

        x0 = -a * c / test;
        y0 = -b * c / test;

        f = b / Math.Sqrt(test);
        g = -a / Math.Sqrt(test);

        switch (f)
        {
            case < 0.0:
                f = -f;
                g = -g;
                break;
        }

    }

    public static double line_par_point_dist_2d(double f, double g, double x0, double y0,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_PAR_POINT_DIST_2D: distance ( parametric line, point ) in 2D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we may choose F*F+G*G = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
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
        //    Input, double F, G, X0, Y0, the parametric line parameters.
        //
        //    Input, double P[2], the point whose distance from the line is
        //    to be measured.
        //
        //    Output, double LINE_PAR_POINT_DIST_2D, the distance from the
        //    point to the line.
        //
    {
        double dx = g * g * (p[0] - x0) - f * g * (p[1] - y0);
        double dy = -f * g * (p[0] - x0) + f * f * (p[1] - y0);

        double value = Math.Sqrt(dx * dx + dy * dy) / (f * f + g * g);

        return value;
    }

    public static double line_par_point_dist_3d(double f, double g, double h, double x0,
            double y0, double z0, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_PAR_POINT_DIST_3D: distance ( parametric line, point ) in 3D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 3D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //      Z = Z0 + H * T
        //
        //    For normalization, we may choose F*F+G*G+H*H = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
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
        //    Input, double F, G, H, X0, Y0, Z0, the parametric line parameters.
        //
        //    Input, double P[3], the point whose distance from the line is
        //    to be measured.
        //
        //    Output, double LINE_PAR_POINT_DIST_3D, the distance from the point
        //    to the line.
        //
    {
        double dx = g * (f * (p[1] - y0) - g * (p[0] - x0))
                    + h * (f * (p[2] - z0) - h * (p[0] - x0));

        double dy = h * (g * (p[2] - z0) - h * (p[1] - y0))
                    - f * (f * (p[1] - y0) - g * (p[0] - x0));

        double dz = -f * (f * (p[2] - z0) - h * (p[0] - x0))
                    - g * (g * (p[2] - z0) - h * (p[1] - y0));

        double value = Math.Sqrt(dx * dx + dy * dy + dz * dz)
                       / (f * f + g * g + h * h);

        return value;
    }

    public static double[] line_par_point_near_2d(double f, double g, double x0, double y0,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_PAR_POINT_NEAR_2D: nearest point on parametric line to point in 2D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we may choose F*F+G*G = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2013
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
        //    Input, double F, G, X0, Y0, the parametric line parameters.
        //
        //    Input, double P[2], the point whose distance from the line is
        //    to be measured.
        //
        //    Output, double LINE_PAR_POINT_DIST_2D[2], the nearest point.
        //
    {
        double t = (f * (p[0] - x0) + g * (p[1] - y0)) / (f * f + g * g);

        double[] pn = new double[2];

        pn[0] = x0 + t * f;
        pn[1] = y0 + t * g;

        return pn;
    }

    public static double[] line_par_point_near_3d(double f, double g, double h, double x0,
            double y0, double z0, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_PAR_POINT_DIST_3D: distance ( parametric line, point ) in 3D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 3D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //      Z = Z0 + H * T
        //
        //    For normalization, we may choose F*F+G*G+H*H = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2013
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
        //    Input, double F, G, H, X0, Y0, Z0, the parametric line parameters.
        //
        //    Input, double P[3], the point whose distance from the line is
        //    to be measured.
        //
        //    Output, double LINE_PAR_POINT_NEAR_3D[3], the nearest point.
        //
    {
        double t = (f * (p[0] - x0) + g * (p[1] - y0) + h * (p[2] - z0))
                   / (f * f + g * g + h * h);

        double[] pn = new double[3];

        pn[0] = x0 + t * f;
        pn[1] = y0 + t * g;
        pn[2] = z0 + t * h;

        return pn;
    }

    public static void line_par2exp_2d(double f, double g, double x0, double y0,
            ref double[] p1, ref double[] p2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_PAR2EXP_2D converts a parametric line to explicit form in 2D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through the points P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2005
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
        //    Input, double F, G, X0, Y0, the parametric line parameters.
        //
        //    Output, double P1[2], P2[2], two points on the line.
        //
    {
        p1[0] = x0;
        p1[1] = y0;

        p2[0] = p1[0] + f;
        p2[1] = p1[1] + g;

    }

    public static void line_par2exp_3d(double f, double g, double h, double x0, double y0,
            double z0, ref double[] p1, ref double[] p2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_PAR2EXP_2D converts a parametric line to explicit form in 3D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 3D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //      Z = Z0 + H * T
        //
        //    For normalization, we choose F*F+G*G+H*H = 1 and 0 <= F.
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through the points P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2013
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
        //    Input, double F, G, H, X0, Y0, Z0, the parametric line parameters.
        //
        //    Output, double P1[3], P2[3], two points on the line.
        //
    {
        p1[0] = x0;
        p1[1] = y0;
        p1[2] = z0;

        p2[0] = p1[0] + f;
        p2[1] = p1[1] + g;
        p2[2] = p1[2] + h;

    }

    public static void line_par2imp_2d(double f, double g, double x0, double y0, ref double a,
            ref double b, ref double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_PAR2IMP_2D converts a parametric line to implicit form in 2D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double F, G, X0, Y0, the parametric parameters of the line.
        //
        //    Output, double[] A, *B, *C, the implicit parameters of the line.
        //
    {
        a = -g;
        b = f;
        c = g * x0 - f * y0;

    }

    public static double lines_exp_angle_3d(double[] p1, double[] p2, double[] p3,
            double[] p4)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two distince points on the first line.
        //
        //    Input, double P3[3], P4[3], two distinct points on the second line.
        //
        //    Output, double LINES_EXP_ANGLE_3D, the angle in radians between the
        //    two lines.  The angle is computed using the ACOS function, and so
        //    lies between 0 and PI.  But if one of the lines is degenerate,
        //    the angle is returned as -1.0.
        //
    {
        double angle;

        double pnorm = Math.Sqrt(Math.Pow(p2[0] - p1[0], 2)
                                 + Math.Pow(p2[1] - p1[1], 2)
                                 + Math.Pow(p2[2] - p1[2], 2));

        double qnorm = Math.Sqrt(Math.Pow(p4[0] - p3[0], 2)
                                 + Math.Pow(p4[1] - p3[1], 2)
                                 + Math.Pow(p4[2] - p3[2], 2));

        double pdotq = (p2[0] - p1[0]) * (p4[0] - p3[0])
                       + (p2[1] - p1[1]) * (p4[1] - p3[1])
                       + (p2[2] - p1[2]) * (p4[2] - p3[2]);

        if (pnorm <= 0.0 || qnorm <= 0.0)
        {
            Console.WriteLine("");
            Console.WriteLine("LINES_EXP_ANGLE_3D - Warning!");
            Console.WriteLine("  One of the lines is degenerate!");
            angle = typeMethods.r8_huge();
        }
        else
        {
            double ctheta = pdotq / (pnorm * qnorm);
            angle = typeMethods.r8_acos(ctheta);
        }

        return angle;
    }

    public static double lines_exp_angle_nd(double[] p1, double[] p2, double[] q1, double[] q2,
            int dim_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_ANGLE_ND returns the angle between two explicit lines in ND.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[DIM_NUM], P2[DIM_NUM], two points on the first line.
        //
        //    Input, double Q1[DIM_NUM], Q2[DIM_NUM], two points on the second line.
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Output, double LINES_EXP_ANGLE_ND, the angle in radians between the two lines.
        //    The angle is computed using the ACOS function, and so lies between 0 and PI.
        //    But if one of the lines is degenerate, the angle is returned as -1.0.
        //
    {
        int i;

        double pnorm = 0.0;
        for (i = 0; i < dim_num; i++)
        {
            pnorm += Math.Pow(p2[i] - p1[i], 2);
        }

        pnorm = Math.Sqrt(pnorm);

        double qnorm = 0.0;
        for (i = 0; i < dim_num; i++)
        {
            qnorm += Math.Pow(q2[i] - q1[i], 2);
        }

        qnorm = Math.Sqrt(qnorm);

        double pdotq = 0.0;
        for (i = 0; i < dim_num; i++)
        {
            pdotq += (p2[i] - p1[i]) * (q2[i] - q1[i]);
        }

        if (pnorm == 0.0 || qnorm == 0.0)
        {
            Console.WriteLine("");
            Console.WriteLine("LINES_EXP_ANGLE_ND - Fatal error!");
            Console.WriteLine("  One of the lines is degenerate!");
            return 1;
        }

        double ctheta = pdotq / (pnorm * qnorm);
        double angle = typeMethods.r8_acos(ctheta);

        return angle;
    }

    public static double lines_exp_dist_3d(double[] p1, double[] p2, double[] q1,
            double[] q2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_DIST_3D computes the distance between two explicit lines in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two distinct points on the first line.
        //
        //    Input, double Q1[3], Q2[3], two distinct points on the second line.
        //
        //    Output, double LINES_EXP_DIST_3D, the distance between the lines.
        //
    {
        const int DIM_NUM = 3;

        double[] a1 = new double[DIM_NUM];
        double[] a2 = new double[DIM_NUM];
        double[] a3 = new double[DIM_NUM];
        double dist;
        //
        //  The distance is found by computing the volume of a parallelipiped,
        //  and dividing by the area of its base.
        //
        //  But if the lines are parallel, we compute the distance by
        //  finding the distance between the first line and any point
        //  on the second line.
        //
        a1[0] = q1[0] - p1[0];
        a1[1] = q1[1] - p1[1];
        a1[2] = q1[2] - p1[2];

        a2[0] = p2[0] - p1[0];
        a2[1] = p2[1] - p1[1];
        a2[2] = p2[2] - p1[2];

        a3[0] = q2[0] - q1[0];
        a3[1] = q2[1] - q1[1];
        a3[2] = q2[2] - q1[2];

        double[] cr = typeMethods.r8vec_cross_product_3d(a2, a3);

        double bot = typeMethods.r8vec_norm(3, cr);

        switch (bot)
        {
            case 0.0:
                dist = line_exp_point_dist_3d(p1, p2, q1);
                break;
            default:
                double top = Math.Abs(a1[0] * (a2[1] * a3[2] - a2[2] * a3[1])
                                      - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0])
                                      + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]));

                dist = top / bot;
                break;
        }

        return dist;
    }

    public static double lines_exp_dist_3d_2(double[] p1, double[] p2, double[] q1,
            double[] q2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_DIST_3D_2 computes the distance between two explicit lines in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through the points P1 and P2.
        //
        //    This routine uses a method that is essentially independent of dimension.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two points on the first line.
        //
        //    Input, double Q1[3], Q2[3], two points on the second line.
        //
        //    Output, double LINES_EXP_DIST_3D_2, the distance between the lines.
        //
    {
        const int DIM_NUM = 3;

        int i;
        double[] pn = new double[DIM_NUM];
        double[] qn = new double[DIM_NUM];
        double sn;
        double tn;
        double[] u = new double[DIM_NUM];
        double[] v = new double[DIM_NUM];
        double[] w0 = new double[DIM_NUM];
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
        double a = typeMethods.r8vec_dot_product(DIM_NUM, u, u);
        double b = typeMethods.r8vec_dot_product(DIM_NUM, u, v);
        double c = typeMethods.r8vec_dot_product(DIM_NUM, v, v);
        double d = typeMethods.r8vec_dot_product(DIM_NUM, u, w0);
        double e = typeMethods.r8vec_dot_product(DIM_NUM, v, w0);
        //
        //  Check the determinant.
        //
        double det = -a * c + b * b;

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

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[i] + sn * (p2[i] - p1[i]);
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            qn[i] = q1[i] + tn * (q2[i] - q1[i]);
        }

        double dist = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dist += Math.Pow(pn[i] - qn[i], 2);
        }

        dist = Math.Sqrt(dist);

        return dist;
    }

    public static bool lines_exp_equal_2d(double[] p1, double[] p2, double[] q1,
            double[] q2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_EQUAL_2D determines if two explicit lines are equal in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through the points P1 and P2.
        //
        //    It is essentially impossible to accurately determine whether two
        //    explicit lines are equal in 2D.  However, for form's sake, and
        //    because occasionally the correct result can be determined, we
        //    provide this routine.  Since divisions are avoided, if the
        //    input data is exactly representable, the result should be
        //    correct.
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
        //    Input, double P1[2], P2[2], two points on the first line.
        //
        //    Input, double Q1[2], Q2[2], two points on the second line.
        //
        //    Output, bool LINES_EXP_EQUAL_2D, is TRUE if the two lines are
        //    determined to be identical.
        //
    {
        //
        //  Slope (P1,P2) = Slope (P2,Q1).
        //
        double test1 = (p2[1] - p1[1]) * (q1[0] - p2[0])
                       - (p2[0] - p1[0]) * (q1[1] - p2[1]);

        if (test1 != 0.0)
        {
            return false;
        }

        //
        //  Slope (Q1,Q2) = Slope (P2,Q1).
        //
        double test2 = (q2[1] - q1[1]) * (q1[0] - p2[0])
                       - (q2[0] - q1[0]) * (q1[1] - p2[1]);

        if (test2 != 0.0)
        {
            return false;
        }

        //
        //  Slope (P1,P2) = Slope (P1,Q2).
        //
        double test3 = (p2[1] - p1[1]) * (q2[0] - p1[0])
                       - (p2[0] - p1[0]) * (q2[1] - p1[1]);

        if (test3 != 0.0)
        {
            return false;
        }

        //
        //  Slope (Q1,Q2) = Slope (P1,Q2).
        //
        double test4 = (q2[1] - q1[1]) * (q2[0] - p1[0])
                       - (q2[0] - q1[0]) * (q2[1] - p1[1]);

        return test4 == 0.0;
    }

    public static void lines_exp_int_2d(double[] p1, double[] p2, double[] p3, double[] p4,
            ref int ival, ref double[] p, int p1Index = 0, int p2Index = 0, int p3Index = 0, int p4Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
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
        //    Input, double P1[2], P2[2], define the first line.
        //
        //    Input, double P3[2], P4[2], define the second line.
        //
        //    Output, int *IVAL, reports on the intersection:
        //    0, no intersection, the lines may be parallel or degenerate.
        //    1, one intersection point, returned in P.
        //    2, infinitely many intersections, the lines are identical.
        //
        //    Output, double P[2], if IVAl = 1, then P contains
        //    the intersection point.  Otherwise, P = 0.
        //
    {
        const int DIM_NUM = 2;

        double a1 = 0.0;
        double a2 = 0.0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c1 = 0.0;
        double c2 = 0.0;

        ival = 0;
        p[0] = 0.0;
        p[1] = 0.0;
        //
        //  Check whether either line is a point.
        //
        bool point_1 = typeMethods.r8vec_eq(DIM_NUM, p1, p2, p1Index, p2Index);

        bool point_2 = typeMethods.r8vec_eq(DIM_NUM, p3, p4, p3Index, p4Index);

        switch (point_1)
        {
            //
            //  Convert the lines to ABC format.
            //
            case false:
                line_exp2imp_2d(p1, p2, ref a1, ref b1, ref c1, p1Index, p2Index);
                break;
        }

        switch (point_2)
        {
            case false:
                line_exp2imp_2d(p3, p4, ref a2, ref b2, ref c2, p3Index, p4Index);
                break;
        }

        switch (point_1)
        {
            //
            //  Search for intersection of the lines.
            //
            case true when point_2:
            {
                if (typeMethods.r8vec_eq(DIM_NUM, p1, p3, p1Index, p3Index))
                {
                    ival = 1;
                    typeMethods.r8vec_copy(DIM_NUM, p1, ref p, p1Index);
                }

                break;
            }
            case true:
            {
                if (Math.Abs(a2 * p1[(0 + p1Index) % p1.Length] + b2 * p1[(1 + p1Index) % p1.Length] - c2) <= typeMethods.r8_epsilon())
                {
                    ival = 1;
                    typeMethods.r8vec_copy(DIM_NUM, p1, ref p, p1Index);
                }

                break;
            }
            default:
            {
                switch (point_2)
                {
                    case true:
                    {
                        if (Math.Abs(a1 * p3[(0 + p3Index) % p3.Length] + b1 * p3[(1 + p3Index) % p3.Length] - c1) <= typeMethods.r8_epsilon())
                        {
                            ival = 1;
                            typeMethods.r8vec_copy(DIM_NUM, p3, ref p, p3Index);
                        }

                        break;
                    }
                    default:
                        lines_imp_int_2d(a1, b1, c1, a2, b2, c2, ref ival, ref p);
                        break;
                }

                break;
            }
        }

    }

    public static void lines_exp_near_3d(double[] p1, double[] p2, double[] q1,
            double[] q2, ref double[] pn, ref double[] qn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_NEAR_3D computes nearest points on two explicit lines in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through the points P1 and P2.
        //
        //    This routine uses a method that is essentially independent of dimension.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two points on the first line.
        //
        //    Input, double Q1[3], Q2[3], two points on the second line.
        //
        //    Output, double PN[3], QN[3], the nearest points on the lines.
        //
    {
        const int DIM_NUM = 3;

        int i;
        double sn;
        double tn;
        double[] u = new double[DIM_NUM];
        double[] v = new double[DIM_NUM];
        double[] w0 = new double[DIM_NUM];
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
        double a = typeMethods.r8vec_dot_product(DIM_NUM, u, u);
        double b = typeMethods.r8vec_dot_product(DIM_NUM, u, v);
        double c = typeMethods.r8vec_dot_product(DIM_NUM, v, v);
        double d = typeMethods.r8vec_dot_product(DIM_NUM, u, w0);
        double e = typeMethods.r8vec_dot_product(DIM_NUM, v, w0);
        //
        //  Check the determinant.
        //
        double det = -a * c + b * b;

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

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = p1[i] + sn * (p2[i] - p1[i]);
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            qn[i] = q1[i] + tn * (q2[i] - q1[i]);
        }
    }

    public static bool lines_exp_parallel_2d(double[] p1, double[] p2, double[] q1,
            double[] q2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_PARALLEL_2D determines if two lines are parallel in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 2D is:
        //
        //      the line through P1 and P2.
        //
        //    The test is essentially a comparison of slopes, but should be
        //    more accurate than an explicit slope comparison, and unfazed
        //    by degenerate cases.
        //
        //    On the other hand, there is NO tolerance for error.  If the
        //    slopes differ by a single digit in the last place, then the
        //    lines are judged to be nonparallel.  A more robust test would
        //    be to compute the angle between the lines, because then it makes
        //    sense to say the lines are "almost" parallel: the angle is small.
        //
        //    If the lines are determined to be parallel, then you can
        //    determine whether they are identical or distinct by evaluating:
        //
        //      lines_exp_parallel_2d ( p1, q2, q1, p2 )
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
        //    Input, double P1[2], P2[2], define the first line.
        //
        //    Input, double Q1[2], Q2[2], define the second line.
        //
        //    Output, bool LINES_EXP_PARALLEL_2D is TRUE if the lines are parallel.
        //
    {
        bool value = Math.Abs((p2[1] - p1[1]) * (q2[0] - q1[0]) - (q2[1] - q1[1]) * (p2[0] - p1[0])) > typeMethods.r8_epsilon();

        return !value;
    }

    public static bool lines_exp_parallel_3d(double[] p1, double[] p2, double[] q1,
            double[] q2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_EXP_PARALLEL_3D determines if two lines are parallel in 3D.
        //
        //  Discussion:
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through P1 and P2.
        //
        //    The points P1, P2 define a direction (P2-P1).  Similarly, the
        //    points (Q1,Q2) define a direction (Q2-Q1).  The quantity
        //
        //      (P2-P1) dot (Q2-Q1) = norm(P2-P1) * norm(Q2-Q1) * cos ( angle )
        //
        //    Therefore, the following value is between 0 and 1;
        //
        //      abs ( (P2-P1) dot (Q2-Q1) / ( norm(P2-P1) * norm(Q2-Q1) ) )
        //
        //    and the lines are parallel if
        //
        //      abs ( (P2-P1) dot (Q2-Q1) / ( norm(P2-P1) * norm(Q2-Q1) ) ) = 1
        //
        //    We can rephrase this as requiring:
        //
        //      ( (P2-P1)dot(Q2-Q1) )^2 = (P2-P1)dot(P2-P1) * (Q2-Q1)dot(Q2-Q1)
        //
        //    which avoids division and square roots.
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
        //    Input, double P1[3], P2[3], define the first line.
        //
        //    Input, double Q1[3], Q2[3, define the second line.
        //
        //    Output, bool LINES_EXP_PARALLEL_3D is TRUE if the lines are parallel.
        //
    {
        const int DIM_NUM = 3;

        int i;

        double[] p = new double[DIM_NUM];
        double[] q = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            p[i] = p2[i] - p1[i];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            q[i] = q2[i] - q1[i];
        }

        double pdotq = typeMethods.r8vec_dot_product(DIM_NUM, p, q);
        double pdotp = typeMethods.r8vec_dot_product(DIM_NUM, p, p);
        double qdotq = typeMethods.r8vec_dot_product(DIM_NUM, q, q);

        bool value = Math.Abs(pdotq * pdotq - pdotp * qdotq) > typeMethods.r8_epsilon();

        return !value;
    }

    public static double lines_imp_angle_2d(double a1, double b1, double c1,
            double a2, double b2, double c2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_IMP_ANGLE_2D finds the angle between two implicit lines in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 1999
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
        //    Input, double A1, B1, C1, the implicit parameters of the first line.
        //
        //    Input, double A2, B2, C2, the implicit parameters of the second line.
        //
        //    Output, double LINES_IMP_ANGLE_2D, the angle between the two lines.
        //
    {
        double pdotq = a1 * a2 + b1 * b2;
        double pnorm = Math.Sqrt(a1 * a1 + b1 * b1);
        double qnorm = Math.Sqrt(a2 * a2 + b2 * b2);

        double ctheta = pdotq / (pnorm * qnorm);

        double theta = Math.Acos(ctheta);

        return theta;
    }

    public static double lines_imp_dist_2d(double a1, double b1, double c1, double a2,
            double b2, double c2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_IMP_DIST_2D determines the distance between two implicit lines in 2D.
        //
        //  Discussion:
        //
        //    If the lines are not parallel, then they must intersect, so their
        //    distance is zero.
        //
        //    If the two lines are parallel, then they may have a nonzero distance.
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
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
        //    Input, double A1, B1, C1, define the first line.
        //    At least one of A1 and B1 must be nonzero.
        //
        //    Input, double A2, B2, C2, define the second line.
        //    At least one of A2 and B2 must be nonzero.
        //
        //    Output, double LINES_IMP_DIST_2D, the distance between the two lines.
        //
    {
        double value = 0;
        switch (a1)
        {
            //
            //  Refuse to handle degenerate lines.
            //
            case 0.0 when b1 == 0.0:
                Console.WriteLine("");
                Console.WriteLine("LINES_IMP_DIST_2D - Fatal error!");
                Console.WriteLine("  Line 1 is degenerate.");
                return 1;
        }

        switch (a2)
        {
            case 0.0 when b2 == 0.0:
                Console.WriteLine("");
                Console.WriteLine("LINES_IMP_DIST_2D - Fatal error!");
                Console.WriteLine("  Line 2 is degenerate.");
                return 1;
        }

        //
        //  If the lines are not parallel, they intersect, and have distance 0.
        //
        if (Math.Abs(a1 * b2 - a2 * b1) > typeMethods.r8_epsilon())
        {
            value = 0.0;
            return value;
        }

        //
        //  Determine the distance between the parallel lines.
        //
        value = Math.Abs(c2 / Math.Sqrt(a2 * a2 + b2 * b2)
                         - c1 / Math.Sqrt(a1 * a1 + b1 * b1));

        return value;
    }

    public static void lines_imp_int_2d(double a1, double b1, double c1, double a2, double b2,
            double c2, ref int ival, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
        //
        //  Discussion:
        //
        //    The implicit form of a line in 2D is:
        //
        //      A * X + B * Y + C = 0
        //
        //    22 May 2004: Thanks to John Asmuth for pointing out that the
        //    B array was not being deallocated on exit.
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
        //    Input, double A1, B1, C1, define the first line.
        //    At least one of A1 and B1 must be nonzero.
        //
        //    Input, double A2, B2, C2, define the second line.
        //    At least one of A2 and B2 must be nonzero.
        //
        //    Output, int *IVAL, reports on the intersection.
        //    -1, both A1 and B1 were zero.
        //    -2, both A2 and B2 were zero.
        //     0, no intersection, the lines are parallel.
        //     1, one intersection point, returned in P.
        //     2, infinitely many intersections, the lines are identical.
        //
        //    Output, double P[2], if IVAL = 1, then P contains
        //    the intersection point.  Otherwise, P = 0.
        //
    {
        const int DIM_NUM = 2;

        double[] a = new double[DIM_NUM * 2];

        p[0] = 0.0;
        p[1] = 0.0;
        switch (a1)
        {
            //
            //  Refuse to handle degenerate lines.
            //
            case 0.0 when b1 == 0.0:
                ival = -1;
                return;
        }

        switch (a2)
        {
            case 0.0 when b2 == 0.0:
                ival = -2;
                return;
        }

        //
        //  Set up a linear system, and compute its inverse.
        //
        a[0 + 0 * 2] = a1;
        a[0 + 1 * 2] = b1;
        a[1 + 0 * 2] = a2;
        a[1 + 1 * 2] = b2;

        double[] b = typeMethods.r8mat_inverse_2d(a);
        //
        //  If the inverse exists, then the lines intersect.
        //  Multiply the inverse times -C to get the intersection point.
        //
        if (b != null)
        {

            ival = 1;
            p[0] = -b[0 + 0 * 2] * c1 - b[0 + 1 * 2] * c2;
            p[1] = -b[1 + 0 * 2] * c1 - b[1 + 1 * 2] * c2;
        }
        //
        //  If the inverse does not exist, then the lines are parallel
        //  or coincident.  Check for parallelism by seeing if the
        //  C entries are in the same ratio as the A or B entries.
        //
        else
        {
            ival = 0;

            switch (a1)
            {
                case 0.0:
                {
                    if (Math.Abs(b2 * c1 - c2 * b1) <= typeMethods.r8_epsilon())
                    {
                        ival = 2;
                    }

                    break;
                }
                default:
                {
                    if (Math.Abs(a2 * c1 - c2 * a1) <= typeMethods.r8_epsilon())
                    {
                        ival = 2;
                    }

                    break;
                }
            }
        }
    }

    public static double lines_par_angle_2d(double f1, double g1, double x01, double y01,
            double f2, double g2, double x02, double y02)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_PAR_ANGLE_2D finds the angle between two parametric lines in 2D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double F1, G1, X01, Y01, the parametric parameters of the
        //    first line.
        //
        //    Input, double F2, G2, X02, Y02, the parametric parameters of the
        //    second line.
        //
        //    Output, double LINES_PAR_ANGLE_2D, the angle between the two lines.
        //
    {
        double pdotq = f1 * f2 + g1 * g2;
        double pnorm = Math.Sqrt(f1 * f1 + g1 * g1);
        double qnorm = Math.Sqrt(f2 * f2 + g2 * g2);

        double value = typeMethods.r8_acos(pdotq / (pnorm * qnorm));

        return value;
    }

    public static double lines_par_angle_3d(double f1, double g1, double h1, double x01,
            double y01, double z01, double f2, double g2, double h2, double x02,
            double y02, double z02)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_PAR_ANGLE_3D finds the angle between two parametric lines in 3D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 3D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //      Z = Z0 + H * T
        //
        //    For normalization, we choose F*F+G*G+H*H = 1 and 0 <= F.
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double F1, G1, H1, X01, Y01, Z01, the parametric parameters
        //    of the first line.
        //
        //    Input, double F2, G2, H2, X02, Y02, Z02, the parametric parameters
        //    of the second line.
        //
        //    Output, double LINES_PAR_ANGLE_3D, the angle between the two lines.
        //
    {
        double pdotq = f1 * f2 + g1 * g2 + h1 * h2;
        double pnorm = Math.Sqrt(f1 * f1 + g1 * g1 + h1 * h1);
        double qnorm = Math.Sqrt(f2 * f2 + g2 * g2 + h2 * h2);

        double value = typeMethods.r8_acos(pdotq / (pnorm * qnorm));

        return value;
    }

    public static double lines_par_dist_3d(double f1, double g1, double h1, double x01,
            double y01, double z01, double f2, double g2, double h2, double x02,
            double y02, double z02)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_PAR_DIST_3D finds the distance between two parametric lines in 3D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 3D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //      Z = Z0 + H * T
        //
        //    For normalization, we choose F*F+G*G+H*H = 1 and 0 <= F.
        //
        //    This code does not work for parallel or near parallel lines.
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double F1, G1, H1, X01, Y01, Z01, the parametric parameters
        //    of the first line.
        //
        //    Input, double F2, G2, H2, X02, Y02, Z02, the parametric parameters
        //    of the second line.
        //
        //    Output, double LINES_PAR_DIST_3D, the distance between the two lines.
        //
    {
        double value = 0;

        value = Math.Abs((x02 - x01) * (g1 * h2 - g2 * h1)
                         + (y02 - y01) * (h1 * f2 - h2 * f1)
                         + (z02 - z01) * (f1 * g2 - f2 * g1)) /
                ((f1 * g2 - f2 * g1) * (f1 * g2 - f2 * g1)
                 + (g1 * h2 - g2 * h1) * (g1 * h2 - g2 * h1)
                 + (h1 * f2 - h2 * f1) * (h1 * f2 - h2 * f1));

        return value;
    }

    public static void lines_par_int_2d(double f1, double g1, double x1, double y1, double f2,
            double g2, double x2, double y2, ref double t1, ref double t2, ref double[] pint)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINES_PAR_INT_2D determines where two parametric lines intersect in 2D.
        //
        //  Discussion:
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
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
        //    Input, double F1, G1, X1, Y1, define the first parametric line.
        //
        //    Input, double F2, G2, X2, Y2, define the second parametric line.
        //
        //    Output, double[] T1, *T2, the T parameters on the first and second
        //    lines of the intersection point.
        //
        //    Output, double PINT[2], the intersection point.
        //
    {
        double det = f2 * g1 - f1 * g2;

        switch (det)
        {
            case 0.0:
                t1 = 0.0;
                t2 = 0.0;
                typeMethods.r8vec_zero(2, ref pint);
                break;
            default:
                t1 = (f2 * (y2 - y1) - g2 * (x2 - x1)) / det;
                t2 = (f1 * (y2 - y1) - g1 * (x2 - x1)) / det;
                pint[0] = x1 + f1 * t1;
                pint[1] = y1 + g1 * t1;
                break;
        }
    }
}