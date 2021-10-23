using System;
using Burkardt.Geometry;
using Burkardt.Types;

namespace Burkardt.Plane
{
    public static class Geometry
    {
        public static void plane_exp_grid_3d(double[] p1, double[] p2, double[] p3, ref int ncor3,
                ref int line_num, ref double[] cor3, ref int[] lines, int maxcor3, int line_max,
                ref int ierror)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_GRID_3D computes points and lines making up a planar grid in 3D.
            //
            //  Discussion:
            //
            //    The data format used is that of SGI Inventor.
            //
            //    On input, if NCOR3 is zero (or negative), then the data computed by
            //    this routine will be stored normally in COR3.  But if NCOR3 is
            //    positive, it is assumed that COR3 already contains NCOR3 items
            //    of useful data.  The new data is appended to COR3.  On output, NCOR3
            //    is increased by the number of points computed by this routine.
            //
            //    On input, if LINE_NUM is zero (or negative), then the data computed by
            //    this routine will be stored normally in LINES.  But if LINE_NUM is positive,
            //    it is assumed that LINES already contains some useful data.  The
            //    new data is appended to LINES.  On output, LINE_NUM is increased by the
            //    number of line data items computed by this routine.
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
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Input/output, int *NCOR3, the number of points stored in COR3.
            //
            //    Input/output, int *LINE_NUM, the number of line data items.
            //
            //    Output, double COR3[3*MAXCOR3], the coordinates of points
            //    used in the grid.
            //
            //    Output, int LINES[LINE_MAX], the indices of points used in
            //    the lines of the grid.  Successive entries of LINES are joined
            //    by a line, unless an entry equals -1.  Note that indices begin
            //    with 0.
            //
            //    Input, int MAXCOR3, the maximum number of points.
            //
            //    Input, int LINE_MAX, the maximum number of lines.
            //
            //    Output, int *IERROR, error indicator.
            //    0, no error.
            //    1, more space for point coordinates is needed.
            //    2, more space for line data is needed.
            //
        {
            int DIM_NUM = 3;
            int NX = 5;
            int NY = 5;

            double a = 0.0;
            double amax = 0.0;
            double amin = 0.0;
            double b = 0.0;
            double bmax = 0.0;
            double bmin = 0.0;
            double dot = 0.0;
            int i;
            int j;
            int k;
            int nbase;
            double[] v1 = new double[DIM_NUM];
            double[] v2 = new double[DIM_NUM];

            ierror = 0;

            if (ncor3 <= 0)
            {
                ncor3 = 0;
            }

            if (line_num <= 0)
            {
                line_num = 0;
            }

            nbase = ncor3;
            //
            //  Compute the two basis vectors for the affine plane.
            //
            v1[0] = p2[0] - p1[0];
            v1[1] = p2[1] - p1[1];
            v1[2] = p2[2] - p1[2];

            Burkardt.Vector.Geometry.vector_unit_nd(DIM_NUM, ref v1);

            v2[0] = p3[0] - p1[0];
            v2[1] = p3[1] - p1[1];
            v2[2] = p3[2] - p1[2];

            dot = typeMethods.r8vec_dot_product(3, v1, v2);
            //
            //  Remove the component of V1 from V2, and give the
            //  resulting vector unit norm.  V1 and V2 are now orthogonal
            //  and of unit length, and represent the two direction vectors
            //  of our plane.
            //
            for (i = 0; i < DIM_NUM; i++)
            {
                v2[i] = v2[i] - dot * v1[i];
            }

            Burkardt.Vector.Geometry.vector_unit_nd(DIM_NUM, ref v2);
            //
            //  Compute the (V1,V2) coordinate range of the input data, if any.
            //
            if (ncor3 == 0)
            {
                amin = 0.0;
                amax = 1.0;
                bmin = 0.0;
                bmax = 1.0;
            }
            else
            {
                for (i = 0; i < ncor3; i++)
                {
                    a = 0.0;
                    b = 0.0;
                    for (j = 0; j < 3; j++)
                    {
                        a = a + v1[j] * cor3[j + i * 3];
                        b = b + v2[j] * cor3[j + i * 3];
                    }

                    if (i == 0)
                    {
                        amin = a;
                        amax = a;
                        bmin = b;
                        bmax = b;
                    }
                    else
                    {
                        amin = Math.Min(amin, a);
                        amax = Math.Max(amax, a);
                        bmin = Math.Min(bmin, b);
                        bmax = Math.Max(bmax, b);
                    }
                }
            }

            //
            //  Generate the points we will use.
            //
            if (maxcor3 < ncor3 + NX * NY)
            {
                ierror = 1;
                return;
            }

            for (j = 1; j <= NY; j++)
            {
                b = ((double) (NY - j) * bmin
                     + (double) (j - 1) * bmax)
                    / (double) (NY - 1);

                for (i = 1; i <= NX; i++)
                {
                    a = ((double) (NX - i) * amin
                         + (double) (i - 1) * amax)
                        / (double) (NX - 1);

                    for (k = 0; k < 3; k++)
                    {
                        cor3[k + (ncor3) * 3] = a * v1[k] + b * v2[k];
                    }

                    ncor3 = ncor3 + 1;
                }
            }

            //
            //  Do the "horizontals".
            //
            for (i = 1; i <= NX; i++)
            {
                for (j = 1; j <= NY; j++)
                {
                    if (line_max <= line_num)
                    {
                        ierror = 2;
                        return;
                    }

                    lines[line_num] = nbase + (j - 1) * NX + i;
                    line_num = line_num + 1;
                }

                if (line_max <= line_num)
                {
                    ierror = 2;
                    return;
                }

                lines[line_num] = -1;
                line_num = line_num + 1;
            }

            //
            //  Do the "verticals".
            //
            for (j = 1; j <= NY; j++)
            {
                for (i = 1; i <= NX; i++)
                {
                    if (line_max <= line_num)
                    {
                        ierror = 2;
                        return;
                    }

                    lines[line_num] = nbase + (j - 1) * NX + i;
                    line_num = line_num + 1;
                }

                if (line_max <= line_num)
                {
                    ierror = 2;
                    return;
                }

                lines[line_num] = -1;
                line_num = line_num + 1;
            }

        }

        public static double plane_exp_point_dist_3d(double[] p1, double[] p2, double[] p3,
                double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_POINT_DIST_3D: distance ( explicit plane, point ) in 3D.
            //
            //  Discussion:
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
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
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Input, double P[3], the coordinates of the point.
            //
            //    Output, double PLANE_EXP_POINT_DIST_3D, the distance from the
            //    point to the plane.
            //
        {
            double a = 0;
            double b = 0;
            double c = 0;
            double d = 0;
            double dist;

            plane_exp2imp_3d(p1, p2, p3, ref a, ref b, ref c, ref d);

            dist = plane_imp_point_dist_3d(a, b, c, d, p);

            return dist;
        }

        public static void plane_exp_normal_3d(double[] p1, double[] p2, double[] p3,
                ref double[] pn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
            //
            //  Discussion:
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
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
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Output, double PN[3], the unit normal vector to the plane.
            //
        {
            double norm;
            //
            //  The cross product (P2-P1) x (P3-P1) is a vector normal to
            //  (P2-P1) and (P3-P1).
            //
            pn[0] = (p2[1] - p1[1]) * (p3[2] - p1[2])
                    - (p2[2] - p1[2]) * (p3[1] - p1[1]);

            pn[1] = (p2[2] - p1[2]) * (p3[0] - p1[0])
                    - (p2[0] - p1[0]) * (p3[2] - p1[2]);

            pn[2] = (p2[0] - p1[0]) * (p3[1] - p1[1])
                    - (p2[1] - p1[1]) * (p3[0] - p1[0]);

            norm = Math.Sqrt(Math.Pow(pn[0], 2) + Math.Pow(pn[1], 2) + Math.Pow(pn[2], 2));

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_EXP_NORMAL_3D - Fatal error!");
                Console.WriteLine("  The plane is poorly defined.");
            }
            else
            {
                pn[0] = pn[0] / norm;
                pn[1] = pn[1] / norm;
                pn[2] = pn[2] / norm;
            }

        }

        public static void plane_exp_pro2(double[] p1, double[] p2, double[] p3, int n,
                double[] pp, ref double[] alpha, ref double[] beta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_PRO2 produces 2D coordinates of points that lie in a plane, in 3D.
            //
            //  Discussion:
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
            //
            //    The first thing to do is to compute two orthonormal vectors V1 and
            //    V2, so that any point P that lies in the plane may be written as
            //
            //      P = P1 + alpha * V1 + beta * V2
            //
            //    The vector V1 lies in the direction P2-P1, and V2 lies in
            //    the plane, is orthonormal to V1, and has a positive component
            //    in the direction of P3-P1.
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
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Input, int N, the number of points to project.
            //
            //    Input, double PP[3*N], the Cartesian coordinates of points which lie on
            //    the plane spanned by the three points.  These points are not checked to
            //    ensure that they lie on the plane.
            //
            //    Output, double ALPHA[N], BETA[N], the "in-plane" coordinates of
            //    the points.
            //
        {
            double dot;
            int i;
            double[] v1 = new double[3];
            double[] v2 = new double[3];
            //
            //  Compute the two basis vectors for the affine plane.
            //
            v1[0] = p2[0] - p1[0];
            v1[1] = p2[1] - p1[1];
            v1[2] = p2[2] - p1[2];

            Burkardt.Vector.Geometry.vector_unit_nd(3, ref v1);

            v2[0] = p3[0] - p1[0];
            v2[1] = p3[1] - p1[1];
            v2[2] = p3[2] - p1[2];

            dot = typeMethods.r8vec_dot_product(3, v1, v2);

            for (i = 0; i < 3; i++)
            {
                v2[i] = v2[i] - dot * v1[i];
            }

            Burkardt.Vector.Geometry.vector_unit_nd(3, ref v2);
            //
            //  Now decompose each point.
            //
            for (i = 0; i < n; i++)
            {
                alpha[i] = (pp[0 + i * 3] - p1[0]) * v1[0]
                           + (pp[1 + i * 3] - p1[1]) * v1[1]
                           + (pp[2 + i * 3] - p1[2]) * v1[2];

                beta[i] = (pp[0 + i * 3] - p1[0]) * v2[0]
                          + (pp[1 + i * 3] - p1[1]) * v2[1]
                          + (pp[2 + i * 3] - p1[2]) * v2[2];
            }

        }

        public static void plane_exp_pro3(double[] p1, double[] p2, double[] p3, int n,
                double[] po, ref double[] pp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_PRO3 projects points orthographically onto a plane, in 3D.
            //
            //  Discussion:
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
            //
            //    PP may share the same memory as PO, in
            //    which case the projections will overwrite the original data.
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
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Input, int N, the number of points to project.
            //
            //    Input, double PO[3*N], the object points.
            //
            //    Output, double PP[3*N], the projections of the object points.
            //
        {
            double a = 0;
            double b = 0;
            double c = 0;
            double d = 0;
            int i;
            //
            //  Put the plane into ABCD form.
            //
            plane_exp2imp_3d(p1, p2, p3, ref a, ref b, ref c, ref d);
            //
            //  For each point, its image in the plane is the nearest point
            //  in the plane.
            //
            for (i = 0; i < n; i++)
            {
                plane_imp_point_near_3d(a, b, c, d, po, ref pp, pIndex: +i * 3, pnIndex: +i * 3);
            }
        }

        public static void plane_exp_project_3d(double[] p1, double[] p2, double[] p3,
                double[] pf, int n, double[] po, ref double[] pp, ref int[] ivis)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_PROJECT_3D projects points through a point onto a plane in 3D.
            //
            //  Discussion:
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
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
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Input, double PF[3], the focus point.
            //
            //    Input, int N, the number of points to project.
            //
            //    Input, double PO[3*N], the object points.
            //
            //    Output, double PP[3*N], the projections of the object points through the
            //    focus point onto the plane.  PP may share the same memory as PO,
            //    in which case the projections will overwrite the original data.
            //
            //    Output, int IVIS[N], visibility indicator:
            //    3, the object was behind the plane;
            //    2, the object was already on the plane;
            //    1, the object was between the focus and the plane;
            //    0, the line from the object to the focus is parallel to the plane,
            //    so the object is "invisible".
            //    -1, the focus is between the object and the plane.  The object
            //    might be considered invisible.
            //
        {
            int DIM_NUM = 3;

            double a = 0;
            double alpha;
            double b = 0;
            double beta;
            double c = 0;
            double d = 0;
            double disfo;
            double disfn;
            int i;
            double[] pn = new double[DIM_NUM];
            //
            //  Put the plane into ABCD form.
            //
            plane_exp2imp_3d(p1, p2, p3, ref a, ref b, ref c, ref d);
            //
            //  Get the nearest point on the plane to the focus.
            //
            plane_imp_point_near_3d(a, b, c, d, pf, ref pn);
            //
            //  Get the distance from the focus to the plane.
            //
            disfn = Burkardt.PointsNS.Geometry.points_dist_3d(pf, pn);
            //
            //  If the focus lies in the plane, this is bad.  We could still
            //  project points that actually lie in the plane, but we'll
            //  just bail out.
            //
            if (disfn == 0.0)
            {
                for (i = 0; i < n; i++)
                {
                    ivis[i] = 0;
                    pp[0 + i * 3] = pf[0];
                    pp[1 + i * 3] = pf[1];
                    pp[2 + i * 3] = pf[2];
                }

                return;
            }

            //
            //  Process the points.
            //
            for (i = 0; i < n; i++)
            {
                //
                //  Get the distance from the focus to the object.
                //
                disfo = Burkardt.PointsNS.Geometry.points_dist_3d(pf, po, p2Index: +i * 3);

                if (disfo == 0.0)
                {
                    ivis[i] = 0;
                    pp[0 + i * 3] = pn[0];
                    pp[1 + i * 3] = pn[1];
                    pp[2 + i * 3] = pn[2];
                }
                else
                {
                    //
                    //  Compute ALPHA, the angle between (OBJECT-FOCUS) and (NEAREST-FOCUS).
                    //
                    alpha = Angle.angle_rad_3d(po, pf, pn, p1Index: +i * 3);

                    if (Math.Cos(alpha) == 0.0)
                    {
                        ivis[i] = 0;
                        pp[0 + i * 3] = pn[0];
                        pp[1 + i * 3] = pn[1];
                        pp[2 + i * 3] = pn[2];
                    }
                    else
                    {
                        //
                        //  BETA is Dist(NEAREST-FOCUS) / ( Cos(ALPHA)*Dist(OBJECT-FOCUS) )
                        //
                        beta = disfn / (Math.Cos(alpha) * disfo);

                        if (1.0 < beta)
                        {
                            ivis[i] = 1;
                        }
                        else if (beta == 1.0)
                        {
                            ivis[i] = 2;
                        }
                        else if (0.0 < beta)
                        {
                            ivis[i] = 3;
                        }
                        else
                        {
                            ivis[i] = -1;
                        }

                        //
                        //  Set the projected point.
                        //
                        pp[0 + i * 3] = pf[0] + beta * (po[0 + i * 3] - pf[0]);
                        pp[1 + i * 3] = pf[1] + beta * (po[1 + i * 3] - pf[1]);
                        pp[2 + i * 3] = pf[2] + beta * (po[2 + i * 3] - pf[2]);
                    }
                }
            }
        }

        public static void plane_exp2imp_3d(double[] p1, double[] p2, double[] p3, ref double a,
                ref double b, ref double c, ref double d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
            //
            //  Discussion:
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
            //
            //    The implicit form of a plane in 3D is
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
            //  Reference:
            //
            //    Adrian Bowyer, John Woodwark,
            //    A Programmer's Geometry,
            //    Butterworths, 1983.
            //
            //  Parameters:
            //
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Output, double *A, *B, *C, *D, coefficients which describe the plane.
            //
        {
            a = (p2[1] - p1[1]) * (p3[2] - p1[2])
                - (p2[2] - p1[2]) * (p3[1] - p1[1]);

            b = (p2[2] - p1[2]) * (p3[0] - p1[0])
                - (p2[0] - p1[0]) * (p3[2] - p1[2]);

            c = (p2[0] - p1[0]) * (p3[1] - p1[1])
                - (p2[1] - p1[1]) * (p3[0] - p1[0]);

            d = -p2[0] * (a) - p2[1] * (b) - p2[2] * (c);

        }

        public static void plane_exp2normal_3d(double[] p1, double[] p2, double[] p3,
                ref double[] pp, ref double[] pn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP2NORMAL_3D converts an explicit plane to normal form in 3D.
            //
            //  Discussion;
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
            //
            //    The normal form of a plane in 3D is
            //
            //      PP, a point on the plane, and
            //      PN, the unit normal to the plane.
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
            //    Input, double P1[3], P2[3], P3[3], three points on the plane.
            //
            //    Output, double PP[3], a point on the plane.
            //
            //    Output, double PN[3], the unit normal vector to the plane.
            //
        {
            int DIM_NUM = 3;

            double norm;

            typeMethods.r8vec_copy(DIM_NUM, p1, ref pp);

            pn[0] = (p2[1] - p1[1]) * (p3[2] - p1[2])
                    - (p2[2] - p1[2]) * (p3[1] - p1[1]);
            pn[1] = (p2[2] - p1[2]) * (p3[0] - p1[0])
                    - (p2[0] - p1[0]) * (p3[2] - p1[2]);
            pn[2] = (p2[0] - p1[0]) * (p3[1] - p1[1])
                    - (p2[1] - p1[1]) * (p3[0] - p1[0]);

            norm = Math.Sqrt(pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2]);

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_EXP2NORMAL_3D - Fatal error!");
                Console.WriteLine("  The normal vector is null.");
                Console.WriteLine("  Two points coincide, or nearly so.");
                return;
            }

            pn[0] = pn[0] / norm;
            pn[1] = pn[1] / norm;
            pn[2] = pn[2] / norm;

        }

        public static bool plane_imp_is_degenerate_3d(double a, double b, double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_IS_DEGENERATE_3D is TRUE if an implicit plane is degenerate.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is:
            //
            //      A * X + B * Y + C * Z + D = 0
            //
            //    The implicit plane is degenerate if A = B = C = 0.
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
            //    Input, double A, B, C, the implicit plane coefficients.
            //
            //    Output, bool PLANE_IMP_IS_DEGENERATE_3D, is TRUE if the plane
            //    is degenerate.
            //
        {
            if (a == 0.0 && b == 0.0 && c == 0.0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public static bool plane_imp_line_par_int_3d(double a, double b, double c, double d,
                double x0, double y0, double z0, double f, double g, double h, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_LINE_PAR_INT_3D: intersection ( implicit plane, parametric line ) in 3D.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is:
            //
            //      A * X + B * Y + C * Z + D = 0
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
            //    09 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Adrian Bowyer, John Woodwark,
            //    A Programmer's Geometry,
            //    Butterworths, 1983, page 111.
            //
            //  Parameters:
            //
            //    Input, double A, B, C, D, parameters that define the implicit plane.
            //
            //    Input, double X0, Y0, Z0, F, G, H, parameters that define the
            //    parametric line.
            //
            //    Output, double P[3], is a point of intersection of the line
            //    and the plane, if the line and plane intersect.
            //
            //    Output, bool PLANE_IMP_LINE_PAR_INT_3D, is TRUE if the line and
            //    the plane intersect, and false otherwise.
            //
        {
            double denom;
            double norm1;
            double norm2;
            double t;
            double TOL = 0.00001;
            //
            //  Check.
            //
            norm1 = Math.Sqrt(a * a + b * b + c * c);

            if (norm1 == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_IMP_LINE_PAR_INT_3D - Fatal error!");
                Console.WriteLine("  The plane normal vector is null.");
                return (true);
            }

            norm2 = Math.Sqrt(f * f + g * g + h * h);

            if (norm2 == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_IMP_LINE_PAR_INT_3D - Fatal error!");
                Console.WriteLine("  The line direction vector is null.");
                return (true);
            }

            denom = a * f + b * g + c * h;
            //
            //  The line and the plane may be parallel.
            //
            if (Math.Abs(denom) < TOL * norm1 * norm2)
            {
                if (a * x0 + b * y0 + c * z0 + d == 0.0)
                {
                    p[0] = x0;
                    p[1] = y0;
                    p[2] = z0;
                    return true;
                }
                else
                {
                    typeMethods.r8vec_zero(3, ref p);
                    return false;
                }
            }
            //
            //  If they are not parallel, they must intersect.
            //
            else
            {
                t = -(a * x0 + b * y0 + c * z0 + d) / denom;
                p[0] = x0 + t * f;
                p[1] = y0 + t * g;
                p[2] = z0 + t * h;
                return true;
            }
        }

        public static double plane_imp_point_dist_3d(double a, double b, double c, double d,
                double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_POINT_DIST_3D: distance ( point, implicit plane ) in 3D.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is:
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
            //  Reference:
            //
            //    Adrian Bowyer, John Woodwark,
            //    A Programmer's Geometry,
            //    Butterworths, 1983.
            //
            //  Parameters:
            //
            //    Input, double A, B, C, D, coefficients that define the plane as
            //    the set of points for which A*X+B*Y+C*Z+D = 0.
            //
            //    Input, double P[3], the coordinates of the point.
            //
            //    Output, double PLANE_IMP_POINT_DIST_3D, the distance from the point to
            //    the plane.
            //
        {
            double dist;

            dist =
                Math.Abs(a * p[0] + b * p[1] + c * p[2] + d) /
                Math.Sqrt(a * a + b * b + c * c);

            return dist;
        }

        public static double plane_imp_point_dist_signed_3d(double a, double b, double c, double d,
                double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_POINT_DIST_SIGNED_3D: signed distance ( implicit plane, point) in 3
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is:
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
            //    Input, double A, B, C, D, determine the equation of the
            //    plane, which is:
            //
            //      A*X + B*Y + C*Z + D = 0.
            //
            //    Input, double P[3], the coordinates of the point.
            //
            //    Output, double PLANE_IMP_POINT_DIST_SIGNED_3D, the signed distance from
            //    the point to the plane.
            //
        {
            double dist;

            dist = -(a * p[0] + b * p[1] + c * p[2] + d)
                   / Math.Sqrt(a * a + b * b + c * c);

            if (d < 0.0)
            {
                dist = -dist;
            }

            return dist;
        }

        public static void plane_imp_point_near_3d(double a, double b, double c, double d,
                double[] p, ref double[] pn, int pIndex = 0, int pnIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_POINT_NEAR_3D: nearest point on a implicit plane to a point in 3D.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is:
            //
            //      A * X + B * Y + C * Z + D = 0
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
            //    Input, double A, B, C, D, coefficients that define the plane as
            //    the set of points for which A*X+B*Y+C*Z+D = 0.
            //
            //    Input, double P[3], the coordinates of the point.
            //
            //    Output, double PN[3], the coordinates of the nearest point on
            //    the plane.
            //
        {
            double t;

            if (plane_imp_is_degenerate_3d(a, b, c))
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_IMP_POINT_NEAR_3D - Fatal error!");
                Console.WriteLine("  A = B = C = 0.");
                return;
            }

            //
            //  The normal N to the plane is (A,B,C).
            //
            //  The line defined by (XN-X)/A = (YN-Y)/B = (ZN-Z)/C = T
            //  goes through (X,Y,Z) and is parallel to N.
            //
            //  Solving for the point PN we get
            //
            //    XN = A*T+X
            //    YN = B*T+Y
            //    ZN = C*T+Z
            //
            //  Now place these values in the equation for the plane:
            //
            //    A*(A*T+X) + B*(B*T+Y) + C*(C*T+Z) + D = 0
            //
            //  and solve for T:
            //
            //    T = (-A*X-B*Y-C*Z-D) / (A * A + B * B + C * C )
            //
            t = -(a * p[(0 + pIndex) % p.Length] + b * p[(1 + pIndex) % p.Length] + c * p[(2 + pIndex) % p.Length] +
                  d) / (a * a + b * b + c * c);

            pn[(0 + pnIndex) % pn.Length] = p[(0 + pIndex) % p.Length] + a * t;
            pn[(1 + pnIndex) % pn.Length] = p[(1 + pIndex) % p.Length] + b * t;
            pn[(2 + pnIndex) % pn.Length] = p[(2 + pIndex) % p.Length] + c * t;

        }

        public static void plane_imp_segment_near_3d(double[] p1, double[] p2, double a, double b,
                double c, double d, ref double dist, ref double[] pnp, ref double[] pnl, int p1Index = 0,
                int p2Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_SEGMENT_NEAR_3D: nearest ( implicit plane, line segment ) in 3D
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
            //    Input, double P1[3], P2[3], the endpoints of the line
            //    segment.
            //
            //    Input, double A, B, C, D, the parameters that define the implicit
            //    plane.
            //
            //    Output, double *DIST, the distance between the line segment and
            //    the plane.
            //
            //    Output, double PNP[3], the nearest point on the plane.
            //
            //    Output, double PNL[3], the nearest point on the line segment
            //    to the plane.  If DIST is zero, PNL is a point of
            //    intersection of the plane and the line segment.
            //
        {
            int DIM_NUM = 3;

            double alpha;
            double an;
            double bn;
            double cn;
            double dn;
            double idiocy;
            double norm;
            double t1;
            double t2;

            typeMethods.r8vec_zero(DIM_NUM, ref pnl);
            typeMethods.r8vec_zero(DIM_NUM, ref pnp);

            norm = Math.Sqrt(a * a + b * b + c * c);

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_IMP_SEGMENT_NEAR_3D - Fatal error!");
                Console.WriteLine("  Plane normal vector is null.");
                return;
            }

            //
            //  The normalized coefficients allow us to compute the (signed) distance.
            //
            an = a / norm;
            bn = b / norm;
            cn = c / norm;
            dn = d / norm;
            //
            //  If the line segment is actually a point, then the answer is easy.
            //
            if (typeMethods.r8vec_eq(DIM_NUM, p1, p2, startIndexA1: p1Index))
            {
                t1 = an * p1[(0 + p1Index) % p1.Length] + bn * p1[(1 + p1Index) % p1.Length] +
                     cn * p1[(2 + p1Index) % p1.Length] + dn;
                dist = Math.Abs(t1);
                typeMethods.r8vec_copy(DIM_NUM, p1, ref pnl, a1index: p1Index);

                pnp[0] = p1[(0 + p1Index) % p1.Length] - an * t1;
                pnp[1] = p1[(1 + p1Index) % p1.Length] - bn * t1;
                pnp[2] = p1[(2 + p1Index) % p1.Length] - cn * t1;

                return;
            }

            //
            //  Compute the projections of the two points onto the normal vector.
            //
            t1 = an * p1[(0 + p1Index) % p1.Length] + bn * p1[(1 + p1Index) % p1.Length] +
                 cn * p1[(2 + p1Index) % p1.Length] + dn;
            t2 = an * p2[(0 + p2Index) % p2.Length] + bn * p2[(1 + p2Index) % p2.Length] +
                 cn * p2[(2 + p2Index) % p2.Length] + dn;
            //
            //  If these have the same sign, then the line segment does not
            //  cross the plane, and one endpoint is the nearest point.
            //
            idiocy = t1 * t2;
            if (0.0 < idiocy)
            {
                t1 = Math.Abs(t1);
                t2 = Math.Abs(t2);

                if (t1 < t2)
                {
                    typeMethods.r8vec_copy(DIM_NUM, p1, ref pnl, a1index: p1Index);
                    pnp[0] = p1[(0 + p1Index) % p1.Length] - an * t1;
                    pnp[1] = p1[(1 + p1Index) % p1.Length] - bn * t1;
                    pnp[2] = p1[(2 + p1Index) % p1.Length] - cn * t1;
                    dist = t1;
                }
                else
                {
                    typeMethods.r8vec_copy(DIM_NUM, p2, ref pnl, a1index: p2Index);
                    dist = t2;
                    pnp[0] = p2[(0 + p2Index) % p2.Length] - an * t2;
                    pnp[1] = p2[(1 + p2Index) % p2.Length] - bn * t2;
                    pnp[2] = p2[(2 + p2Index) % p2.Length] - cn * t2;
                }
                //
                //  If the projections differ in sign, the line segment crosses the plane.
                //
            }
            else
            {

                if (t1 == 0.0)
                {
                    alpha = 0.0;
                }
                else if (t2 == 0.0)
                {
                    alpha = 1.0;
                }
                else
                {
                    alpha = t2 / (t2 - t1);
                }

                pnl[0] = alpha * p1[(0 + p1Index) % p1.Length] + (1.0 - alpha) * p2[(0 + p2Index) % p2.Length];
                pnl[1] = alpha * p1[(1 + p1Index) % p1.Length] + (1.0 - alpha) * p2[(1 + p2Index) % p2.Length];
                pnl[2] = alpha * p1[(2 + p1Index) % p1.Length] + (1.0 - alpha) * p2[(2 + p2Index) % p2.Length];
                typeMethods.r8vec_copy(DIM_NUM, pnl, ref pnp);

                dist = 0.0;
            }

        }

        public static void plane_imp_triangle_int_3d(double a, double b, double c, double d,
                double[] t, ref int int_num, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_TRIANGLE_INT_3D: intersection ( implicit plane, triangle ) in 3D.
            //
            //  Discussion:
            //
            //    An implicit plane in 3D is the set of points satisfying
            //      A * X + B * Y + C * Z + D = 0,
            //    for a given set of parameters A, B, C, D.  At least one of
            //    A, B and C must be nonzero.
            //
            //    There may be 0, 1, 2 or 3 points of intersection return;ed.
            //
            //    If two intersection points are return;ed, then the entire line
            //    between them comprises points of intersection.
            //
            //    If three intersection points are return;ed, then all points of
            //    the triangle intersect the plane.
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
            //    Input, double A, B, C, D, the parameters that define the implicit plane.
            //
            //    Input, double T[3*3], the vertices of the triangle.
            //
            //    Output, double DIST, the distance between the triangle and the plane.
            //
            //    Output, int *INT_NUM, the number of intersection points return;ed.
            //
            //    Output, double P[3*3], the coordinates of the intersection points.
            //
        {
            double dist1;
            double dist2;
            double dist3;
            int n;

            n = 0;
            //
            //  Compute the signed distances between the vertices and the plane.
            //
            dist1 = a * t[0 + 0 * 3] + b * t[1 + 0 * 3] + c * t[2 + 0 * 3] + d;
            dist2 = a * t[0 + 1 * 3] + b * t[1 + 1 * 3] + c * t[2 + 1 * 3] + d;
            dist3 = a * t[0 + 2 * 3] + b * t[1 + 2 * 3] + c * t[2 + 2 * 3] + d;
            //
            //  Consider any zero distances.
            //
            if (dist1 == 0.0)
            {
                p[0 + n * 3] = t[0 + 0 * 3];
                p[1 + n * 3] = t[1 + 0 * 3];
                p[2 + n * 3] = t[2 + 0 * 3];
                n = n + 1;
            }

            if (dist2 == 0.0)
            {
                p[0 + n * 3] = t[0 + 1 * 3];
                p[1 + n * 3] = t[1 + 1 * 3];
                p[2 + n * 3] = t[2 + 1 * 3];
                n = n + 1;
            }

            if (dist3 == 0.0)
            {
                p[0 + n * 3] = t[0 + 2 * 3];
                p[1 + n * 3] = t[1 + 2 * 3];
                p[2 + n * 3] = t[2 + 2 * 3];
                n = n + 1;
            }

            //
            //  If 2 or 3 of the nodes intersect, we're already done.
            //
            if (2 <= n)
            {
                int_num = n;
                return;
            }

            //
            //  If one node intersects, then we're done unless the other two
            //  are of opposite signs.
            //
            if (n == 1)
            {
                if (dist1 == 0.0)
                {
                    plane_imp_triangle_int_add_3d(t, t, dist2, dist3, ref n, p, p1Index: +1 * 3, p2Index: +2 * 3);
                }
                else if (dist2 == 0.0)
                {
                    plane_imp_triangle_int_add_3d(t, t, dist1, dist3, ref n, p, p1Index: +0 * 3, p2Index: +2 * 3);
                }
                else if (dist3 == 0.0)
                {
                    plane_imp_triangle_int_add_3d(t, t, dist1, dist2, ref n, p, p1Index: +0 * 3, p2Index: +1 * 3);
                }

                return;
            }

            //
            //  All nodal distances are nonzero, and there is at least one
            //  positive and one negative.
            //
            if (dist1 * dist2 < 0.0 && dist1 * dist3 < 0.0)
            {
                plane_imp_triangle_int_add_3d(t, t, dist1, dist2, ref n, p, p1Index: +0 * 3, p2Index: +1 * 3);
                plane_imp_triangle_int_add_3d(t, t, dist1, dist3, ref n, p, p1Index: +0 * 3, p2Index: +2 * 3);
            }
            else if (dist2 * dist1 < 0.0 && dist2 * dist3 < 0.0)
            {
                plane_imp_triangle_int_add_3d(t, t, dist2, dist1, ref n, p, p1Index: +1 * 3, p2Index: +0 * 3);
                plane_imp_triangle_int_add_3d(t, t, dist2, dist3, ref n, p, p1Index: +1 * 3, p2Index: +2 * 3);
            }
            else if (dist3 * dist1 < 0.0 && dist3 * dist2 < 0.0)
            {
                plane_imp_triangle_int_add_3d(t, t, dist3, dist1, ref n, p, p1Index: +2 * 3, p2Index: +0 * 3);
                plane_imp_triangle_int_add_3d(t, t, dist3, dist2, ref n, p, p1Index: +2 * 3, p2Index: +1 * 3);
            }

            int_num = n;

        }

        public static void plane_imp_triangle_int_add_3d(double[] p1, double[] p2, double dist1,
                double dist2, ref int int_num, double[] p, int p1Index = 0, int p2Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_TRIANGLE_INT_ADD_3D is a utility for PLANE_IMP_TRIANGLE_INT_3D.
            //
            //  Discussion:
            //
            //    This routine is called to consider the value of the signed distance
            //    from a plane of two nodes of a triangle.  If the two values
            //    have opposite signs, then there is a point of intersection between
            //    them.  The routine computes this point and adds it to the list.
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
            //    Input, double P1[3], P2[3], two vertices of a triangle.
            //
            //    Input, double DIST1, DIST2, the signed distances of the two vertices
            //    from a plane.
            //
            //    Input/output, int *INT_NUM, the number of intersection points.
            //
            //    Input/output, double P[3*(*INT_NUM)], the coordinates
            //    of the intersection points.
            //
        {
            double alpha;
            int n;

            n = int_num;

            if (dist1 == 0.0)
            {
                p[0 + n * 3] = p1[(0 + p1Index) % p1.Length];
                p[1 + n * 3] = p1[(1 + p1Index) % p1.Length];
                p[2 + n * 3] = p1[(2 + p1Index) % p1.Length];
                n = n + 1;
            }
            else if (dist2 == 0.0)
            {
                p[0 + n * 3] = p2[(0 + p2Index) % p2.Length];
                p[1 + n * 3] = p2[(1 + p2Index) % p2.Length];
                p[2 + n * 3] = p2[(2 + p2Index) % p2.Length];
                n = n + 1;
            }
            else if (dist1 * dist2 < 0.0)
            {
                alpha = dist2 / (dist2 - dist1);
                p[0 + n * 3] = alpha * p1[(0 + p1Index) % p1.Length] + (1.0 - alpha) * p2[(0 + p2Index) % p2.Length];
                p[1 + n * 3] = alpha * p1[(1 + p1Index) % p1.Length] + (1.0 - alpha) * p2[(1 + p2Index) % p2.Length];
                p[2 + n * 3] = alpha * p1[(2 + p1Index) % p1.Length] + (1.0 - alpha) * p2[(2 + p2Index) % p2.Length];
                n = n + 1;
            }

            int_num = n;
        }

        public static int plane_imp_triangle_near_3d(double[] t, double a, double b, double c,
                double d, ref double dist, ref double[] pn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_TRIANGLE_NEAR_3D: nearest ( implicit plane, triangle ) in 3D.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is:
            //
            //      A * X + B * Y + C * Z + D = 0
            //
            //    If DIST = 0, then each point is a point of intersection, and there
            //    will be at most 3 such points returned.
            //
            //    If 0 < DIST, then the points are listed in pairs, with the first
            //    being on the triangle, and the second on the plane.  Two points will
            //    be listed in the most common case, but possibly 4 or 6.
            //
            //    Please see to it that the underlying distance routine always returns
            //    one of the endpoints if the entire line segment is at zero distance.
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
            //    Input, double T[3*3], the vertices of the triangle.
            //
            //    Input, double A, B, C, D, the parameters that define the implicit plane.
            //
            //    Output, double *DIST, the distance between the triangle and the plane.
            //
            //    Output, double PN[3*6], a collection of nearest points.
            //
            //    Output, int PLANE_IMP_TRIANGLE_NEAR_3D, the number of nearest points
            //    returned.
            //
        {
            int DIM_NUM = 3;

            double dist12 = 0;
            double dist23 = 0;
            double dist31 = 0;
            int near_num;
            double[] pp = new double[DIM_NUM];
            double[] pt = new double[DIM_NUM];

            near_num = 0;
            //
            //  Consider the line segment P1 - P2.
            //
            plane_imp_segment_near_3d(t, t, a, b, c, d, ref dist12, ref pp, ref pt, p1Index: +0 * 3, p2Index: +1 * 3);

            dist = dist12;
            typeMethods.r8vec_copy(DIM_NUM, pt, ref pn, a2index: +near_num * 3);
            near_num = near_num + 1;

            if (0.0 < dist12)
            {
                typeMethods.r8vec_copy(DIM_NUM, pp, ref pn, a2index: +near_num * 3);
                near_num = near_num + 1;
            }

            //
            //  Consider the line segment P2 - P3.
            //
            plane_imp_segment_near_3d(t, t, a, b, c, d, ref dist23, ref pp, ref pt, p1Index: +1 * 3, p2Index: +2 * 3);

            if (dist23 < dist)
            {
                near_num = 0;
                dist = dist23;

                typeMethods.r8vec_copy(DIM_NUM, pt, ref pn, a2index: +near_num * 3);
                near_num = near_num + 1;

                if (0.0 < dist23)
                {
                    typeMethods.r8vec_copy(DIM_NUM, pp, ref pn, a2index: +near_num * 3);
                    near_num = near_num + 1;
                }
            }
            else if (dist23 == dist)
            {
                typeMethods.r8vec_copy(DIM_NUM, pt, ref pn, a2index: +near_num * 3);
                near_num = near_num + 1;

                if (0.0 < dist23)
                {
                    typeMethods.r8vec_copy(DIM_NUM, pp, ref pn, a2index: +near_num * 3);
                    near_num = near_num + 1;
                }
            }

            //
            //  Consider the line segment P3 - P1.
            //
            plane_imp_segment_near_3d(t, t, a, b, c, d, ref dist31, ref pp, ref pt, p1Index: +2 * 3, p2Index: +0 * 3);

            if (dist31 < dist)
            {
                near_num = 0;
                dist = dist31;

                typeMethods.r8vec_copy(DIM_NUM, pt, ref pn, a2index: +near_num * 3);
                near_num = near_num + 1;

                if (0.0 < dist31)
                {
                    typeMethods.r8vec_copy(DIM_NUM, pp, ref pn, a2index: +near_num * 3);
                    near_num = near_num + 1;
                }
            }
            else if (dist31 == dist)
            {
                typeMethods.r8vec_copy(DIM_NUM, pt, ref pn, a2index: +near_num * 3);
                near_num = near_num + 1;

                if (0.0 < dist31)
                {
                    typeMethods.r8vec_copy(DIM_NUM, pp, ref pn, a2index: +near_num * 3);
                    near_num = near_num + 1;
                }
            }

            return near_num;
        }

        public static void plane_imp2exp_3d(double a, double b, double c, double d, ref double[] p1,
                ref double[] p2, ref double[] p3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP2EXP_3D converts an implicit plane to explicit form in 3D.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is
            //
            //      A * X + B * Y + C * Z + D = 0.
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
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
            //    Input, double A, B, C, D, parameters that define the implicit plane.
            //
            //    Output, double P1[3], P2[3], P3[3], three points on the plane.
            //
        {
            int DIM_NUM = 3;

            double[] pn = new double[DIM_NUM];
            double[] pp = new double[DIM_NUM];

            plane_imp2normal_3d(a, b, c, d, ref pp, ref pn);

            plane_normal2exp_3d(pp, pn, ref p1, ref p2, ref p3);

        }

        public static void plane_imp2normal_3d(double a, double b, double c, double d,
                ref double[] pp, ref double[] pn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP2NORMAL_3D converts an implicit plane to normal form in 3D.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is
            //
            //      A * X + B * Y + C * Z + D = 0.
            //
            //    The normal form of a plane in 3D is
            //
            //      PP, a point on the plane, and
            //      PN, the unit normal to the plane.
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
            //    Input, double A, B, C, D, parameters that define the implicit plane.
            //
            //    Output, double PP[3] point on the plane.
            //
            //    Output, double PN[3], the unit normal vector to the plane.
            //
        {
            double norm;

            norm = Math.Sqrt(a * a + b * b + c * c);

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_IMP2NORMAL_3D - Fatal error!");
                Console.WriteLine("  The normal vector is null.");
                Console.WriteLine("  Two points coincide, or nearly so.");
                return;
            }

            pn[0] = a / norm;
            pn[1] = b / norm;
            pn[2] = c / norm;

            if (a != 0.0)
            {
                pp[0] = -d / a;
                pp[1] = 0.0;
                pp[2] = 0.0;
            }
            else if (b != 0.0)
            {
                pp[0] = 0.0;
                pp[1] = -d / b;
                pp[2] = 0.0;
            }
            else if (c != 0.0)
            {
                pp[0] = 0.0;
                pp[1] = 0.0;
                pp[2] = -d / c;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_IMP2NORMAL_3D - Fatal error!");
                Console.WriteLine("  The (A,B,C) vector is null.");
            }

        }

        public static void plane_normal_basis_3d(double[] pp, double[] pn, ref double[] pq,
                ref double[] pr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is:
            //
            //      PP is a point on the plane,
            //      N is a normal vector to the plane.
            //
            //    The two vectors to be computed, PQ and PR, can be regarded as
            //    the basis of a Cartesian coordinate system for points in the plane.
            //    Any point in the plane can be described in terms of the "origin"
            //    point PP plus a weighted sum of the two vectors PQ and PR:
            //
            //      P = PP + a * PQ + b * PR.
            //
            //    The vectors PQ and PR have unit length, and are perpendicular to N
            //    and to each other.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double PN[3], a normal vector to the plane.  The
            //    vector must not have zero length, but it is not necessary for PN
            //    to have unit length.
            //
            //    Output, double PQ[3], a vector of unit length, perpendicular
            //    to the vector PN and the vector PR.
            //
            //    Output, double PR[3], a vector of unit length, perpendicular
            //    to the vector PN and the vector PQ.
            //
        {
            int DIM_NUM = 3;

            int i;
            double normal_norm;
            double pr_norm;
            double[] temp;
            //
            //  Compute the length of NORMAL.
            //
            normal_norm = typeMethods.r8vec_norm(DIM_NUM, pn);

            if (normal_norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_NORMAL_BASIS_3D - Fatal error!");
                Console.WriteLine("  The normal vector is 0.");
                return;
            }

            //
            //  Find a vector PQ that is normal to PN and has unit length.
            //
            temp = typeMethods.r8vec_any_normal(DIM_NUM, pn);
            typeMethods.r8vec_copy(DIM_NUM, temp, ref pq);
            //
            //  Now just take the cross product PN x PQ to get the PR vector.
            //
            temp = typeMethods.r8vec_cross_product_3d(pn, pq);

            pr_norm = typeMethods.r8vec_norm(DIM_NUM, temp);

            for (i = 0; i < DIM_NUM; i++)
            {
                pr[i] = temp[i] / pr_norm;
            }
        }

        public static int plane_normal_line_exp_int_3d(double[] pp, double[] normal,
                double[] p1, double[] p2, ref double[] pint)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_LINE_EXP_INT_3D: intersection of plane and line in 3D.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is:
            //
            //      PP is a point on the plane,
            //      N is a normal vector to the plane.
            //
            //    The explicit form of a line in 3D is:
            //
            //      P1, P2 are two points on the line.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double NORMAL[3], a normal vector to the plane.
            //
            //    Input, double P1[3], P2[3], two distinct points on the line.
            //
            //    Output, double PINT[3], the coordinates of a
            //    common point of the plane and line, when IVAL is 1 or 2.
            //
            //    Output, integer PLANE_NORMAL_LINE_EXP_INT_3D, the kind of intersection;
            //    0, the line and plane seem to be parallel and separate;
            //    1, the line and plane intersect at a single point;
            //    2, the line and plane seem to be parallel and joined.
            //
        {
            int DIM_NUM = 3;

            double[] direction = new double[DIM_NUM];
            int i;
            int ival;
            double temp;
            double temp2;
            //
            //  Make sure the line is not degenerate.
            //
            if (Burkardt.LineNS.Geometry.line_exp_is_degenerate_nd(DIM_NUM, p1, p2))
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_NORMAL_LINE_EXP_INT_3D - Fatal error!");
                Console.WriteLine("  The line is degenerate.");
                return (1);
            }

            //
            //  Make sure the plane normal vector is a unit vector.
            //
            temp = typeMethods.r8vec_norm(DIM_NUM, normal);

            if (temp == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_NORMAL_LINE_EXP_INT_3D - Fatal error!");
                Console.WriteLine("  The normal vector of the plane is degenerate.");
                return (1);
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                normal[i] = normal[i] / temp;
            }

            //
            //  Determine the unit direction vector of the line.
            //
            for (i = 0; i < DIM_NUM; i++)
            {
                direction[i] = p2[i] - p1[i];
            }

            temp = typeMethods.r8vec_norm(DIM_NUM, direction);

            for (i = 0; i < DIM_NUM; i++)
            {
                direction[i] = direction[i] / temp;
            }

            //
            //  If the normal and direction vectors are orthogonal, then
            //  we have a special case to deal with.
            //
            if (typeMethods.r8vec_dot_product(DIM_NUM, normal, direction) == 0.0)
            {
                temp = 0.0;
                for (i = 0; i < DIM_NUM; i++)
                {
                    temp = temp + normal[i] * (p1[i] - pp[i]);
                }

                if (temp == 0.0)
                {
                    ival = 2;
                    typeMethods.r8vec_copy(DIM_NUM, p1, ref pint);
                }
                else
                {
                    ival = 0;
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        pint[i] = typeMethods.r8_huge();
                    }
                }

                return ival;
            }

            //
            //  Determine the distance along the direction vector to the intersection point.
            //
            temp = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                temp = temp + normal[i] * (pp[i] - p1[i]);
            }

            temp2 = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                temp2 = temp2 + normal[i] * direction[i];
            }

            ival = 1;
            for (i = 0; i < DIM_NUM; i++)
            {
                pint[i] = p1[i] + temp * direction[i] / temp2;
            }

            return ival;
        }

        public static double[] plane_normal_qr_to_xyz(double[] pp, double[] normal, double[] pq,
                double[] pr, int n, double[] qr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_QR_TO_XYZ: QR_TO_XYZ coordinates for a normal form plane.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is:
            //
            //      PP is a point on the plane,
            //      NORMAL is a normal vector to the plane.
            //
            //    Two vectors PQ and PR can be computed with the properties that
            //    * NORMAL, PQ and PR are pairwise orthogonal;
            //    * PQ and PR have unit length;
            //    * every point P in the plane has a "QR" representation
            //      as P = PP + q * PQ + r * PR.
            //
            //    This function is given the QR coordinates of a set of points on the
            //    plane, and returns the XYZ coordinates.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double NORMAL[3], a normal vector N to the plane.  The
            //    vector must not have zero length, but it is not necessary for N
            //    to have unit length.
            //
            //    Input, double PQ[3], a vector of unit length,
            //    perpendicular to the vector N and the vector PR.
            //
            //    Input, double PR[3], a vector of unit length,
            //    perpendicular to the vector N and the vector PQ.
            //
            //    Input, integer N, the number of points on the plane.
            //
            //    Input, double QR[2*N], the QR coordinates of the points.
            //
            //    Output, double PLANE_NORMAL_QR_TO_XYZ[3*N], the XYZ coordinates of the points.
            //
        {
            int i;
            int j;
            double[] xyz;

            xyz = new double[3 * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    xyz[i + j * 3] = pp[i] + pq[i] * qr[0 + j * 2] + pr[i] * qr[1 + j * 2];
                }
            }

            return xyz;
        }

        public static void plane_normal_tetrahedron_intersect(double[] pp,
                double[] normal, double[] t, ref int int_num, ref double[] pint)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_TETRAHEDRON_INTERSECT intersects a plane and a tetrahedron.
            //
            //  Discussion:
            //
            //    The intersection of a plane and a tetrahedron is one of:
            //    0) empty
            //    1) a single point
            //    2) a single line segment
            //    3) a triangle
            //    4) a quadrilateral.
            //
            //    In each case, the region of intersection can be described by the
            //    corresponding number of points.  In particular, cases 2, 3 and 4
            //    are described by the vertices that bound the line segment, triangle,
            //    or quadrilateral.
            //
            //    The normal form of a plane is:
            //
            //      PP is a point on the plane,
            //      N is a normal vector to the plane.
            //
            //    The form of a tetrahedron is
            //
            //      T(1:3,1:4) contains the coordinates of the vertices.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double NORMAL[3], a normal vector to the plane.
            //
            //    Input, double T[3*4], the tetrahedron vertices.
            //
            //    Output, int *INT_NUM, the number of intersection
            //    points returned.  This will be 0, 1, 2, 3 or 4.
            //
            //    Output, double PINT[3*4], the coordinates of the
            //    intersection points.
            //
        {
            double area1;
            double area2;
            double[] d = new double[4];
            double dn;
            double dpp;
            int i;
            int j;
            int j1;
            int j2;
            double temp;

            int_num = 0;
            for (j = 0; j < 4; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    pint[i + j * 3] = 0.0;
                }
            }

            //
            //  DN is the length of the normal vector.
            //
            dn = Math.Sqrt(typeMethods.r8vec_dot_product(3, normal, normal));
            //
            //  DPP is the distance between the origin and the projection of the
            //  point PP onto the normal vector.
            //
            dpp = dn - typeMethods.r8vec_dot_product(3, normal, pp) / dn;
            //
            //  D(I) is positive, zero, or negative if vertex I is above,
            //  on, or below the plane.
            //
            for (j = 0; j < 4; j++)
            {
                d[j] = dn - dpp;
                for (i = 0; i < 3; i++)
                {
                    d[j] = d[j] - normal[i] * t[i + j * 3];
                }
            }

            //
            //  If all D are positive or negative, no intersection.
            //
            if (typeMethods.r8vec_negative_strict(4, d) || typeMethods.r8vec_positive_strict(4, d))
            {
                int_num = 0;
                return;
            }

            //
            //  Points with zero distance are automatically added to the list.
            //
            //  For each point with nonzero distance, seek another point
            //  with opposite sign and higher index, and compute the intersection
            //  of the line between those points and the plane.
            //
            for (j1 = 0; j1 < 4; j1++)
            {
                if (d[j1] == 0.0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        pint[i + (int_num) * 3] = t[i + j1 * 3];
                    }

                    int_num = int_num + 1;
                }
                else
                {
                    for (j2 = j1 + 1; j2 < 4; j2++)
                    {
                        if (typeMethods.r8_sign_opposite_strict(d[j1], d[j2]))
                        {
                            for (i = 0; i < 3; i++)
                            {
                                pint[i + (int_num) * 3] = (d[j1] * t[i + j2 * 3]
                                                           - d[j2] * t[i + j1 * 3])
                                                          / (d[j1] - d[j2]);
                            }

                            int_num = int_num + 1;
                        }
                    }
                }
            }

            //
            //  If four points were found, try to order them properly.
            //
            if (int_num == 4)
            {
                area1 = Burkardt.Quadrilateral.Geometry.quad_area_3d(pint);
                for (i = 0; i < 3; i++)
                {
                    temp = pint[i + 3 * 3];
                    pint[i + 3 * 3] = pint[i + 4 * 3];
                    pint[i + 4 * 3] = temp;
                }

                area2 = Burkardt.Quadrilateral.Geometry.quad_area_3d(pint);
                if (area2 < area1)
                {
                    temp = pint[i + 3 * 3];
                    pint[i + 3 * 3] = pint[i + 4 * 3];
                    pint[i + 4 * 3] = temp;
                }
            }
        }

        public static int plane_normal_triangle_int_3d(double[] pp, double[] pn, double[] t,
                double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_TRIANGLE_INT_3D: intersection ( normal plane, triangle ) in 3D.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is:
            //
            //      PP is a point on the plane,
            //      PN is a normal vector to the plane.
            //
            //    There may be 0, 1, 2 or 3 points of intersection returned.
            //
            //    If two intersection points are returned, then the entire line
            //    between them comprises points of intersection.
            //
            //    If three intersection points are returned, then all points of
            //    the triangle intersect the plane.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double PN[3], a normal vector to the plane.  The
            //    vector must not have zero length, but it is not necessary for PN
            //    to have unit length.
            //
            //    Input, double T[3*3], the vertices of the triangle.
            //
            //    Output, double P[3*3], the coordinates of the intersection points.
            //
            //    Output, int PLANE_NORMAL_TRIANGLE_INT_3D, the number of intersection
            //    points returned.
            //
        {
            int DIM_NUM = 3;

            double d;
            double dist1;
            double dist2;
            double dist3;
            int int_num;

            int_num = 0;
            //
            //  Compute the signed distances between the vertices and the plane.
            //
            d = -typeMethods.r8vec_dot_product(DIM_NUM, pn, pp);

            dist1 = typeMethods.r8vec_dot_product(DIM_NUM, pn, t, a2Index: +0 * 3) + d;
            dist2 = typeMethods.r8vec_dot_product(DIM_NUM, pn, t, a2Index: +1 * 3) + d;
            dist3 = typeMethods.r8vec_dot_product(DIM_NUM, pn, t, a2Index: +2 * 3) + d;
            //
            //  Consider any zero distances.
            //
            if (dist1 == 0.0)
            {
                typeMethods.r8vec_copy(DIM_NUM, t, ref p, a1index: +0 * 3, a2index: +int_num * 3);
                int_num = int_num + 1;
            }

            if (dist2 == 0.0)
            {
                typeMethods.r8vec_copy(DIM_NUM, t, ref p, a1index: +1 * 3, a2index: +int_num * 3);
                int_num = int_num + 1;
            }

            if (dist3 == 0.0)
            {
                typeMethods.r8vec_copy(DIM_NUM, t, ref p, a1index: +2 * 3, a2index: +int_num * 3);
                int_num = int_num + 1;
            }

            //
            //  If 2 or 3 of the nodes intersect, we're already done.
            //
            if (2 <= int_num)
            {
                return int_num;
            }

            //
            //  If one node intersects, then we're done unless the other two
            //  are of opposite signs.
            //
            if (int_num == 1)
            {
                if (dist1 == 0.0)
                {
                    plane_imp_triangle_int_add_3d(t, t, dist2, dist3, ref int_num, p, p1Index: +1 * 3, p2Index: +2 * 3);
                }
                else if (dist2 == 0.0)
                {
                    plane_imp_triangle_int_add_3d(t, t, dist1, dist3, ref int_num, p, p1Index: +0 * 3, p2Index: +2 * 3);
                }
                else if (dist3 == 0.0)
                {
                    plane_imp_triangle_int_add_3d(t, t, dist1, dist2, ref int_num, p, p1Index: +0 * 3, p2Index: +1 * 3);
                }

                return int_num;
            }

            //
            //  All nodal distances are nonzero, and there is at least one
            //  positive and one negative.
            //
            if (dist1 * dist2 < 0.0 && dist1 * dist3 < 0.0)
            {
                plane_imp_triangle_int_add_3d(t, t, dist1, dist2, ref int_num, p, p1Index: +0 * 3, p2Index: +1 * 3);

                plane_imp_triangle_int_add_3d(t, t, dist1, dist3, ref int_num, p, p1Index: +0 * 3, p2Index: +2 * 3);
            }
            else if (dist2 * dist1 < 0.0 && dist2 * dist3 < 0.0)
            {
                plane_imp_triangle_int_add_3d(t, t, dist2, dist1, ref int_num, p, p1Index: +1 * 3, p2Index: +0 * 3);

                plane_imp_triangle_int_add_3d(t, t, dist2, dist3, ref int_num, p, p1Index: +1 * 3, p2Index: +2 * 3);
            }
            else if (dist3 * dist1 < 0.0 && dist3 * dist2 < 0.0)
            {
                plane_imp_triangle_int_add_3d(t, t, dist3, dist1, ref int_num, p, p1Index: +2 * 3, p2Index: +0 * 3);

                plane_imp_triangle_int_add_3d(t, t, dist3, dist2, ref int_num, p, p1Index: +2 * 3, p2Index: +1 * 3);
            }

            return int_num;
        }

        public static void plane_normal_uniform_3d(ref int seed, double[] pp, double[] normal)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_UNIFORM_3D generates a random normal plane in 3D.
            //
            //  Discussion:
            //
            //    The normal form of a plane is:
            //
            //      PP is a point on the plane,
            //      N is a normal vector to the plane.
            //
            //    The point PP will be chosen at random inside the unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 November 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double PP[3], a point on the plane.
            //
            //    Output, double NORMAL[3], the unit normal vector.
            //
        {
            int DIM_NUM = 3;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            int i;
            double norm;
            double[] v;
            //
            //  Pick PP as a random point inside the unit sphere in ND.
            //
            v = Burkardt.Ball.Geometry.ball01_sample_3d(ref seed);
            typeMethods.r8vec_copy(DIM_NUM, v, ref pp);
            //
            //  Get values from a standard normal distribution.
            //
            v = typeMethods.r8vec_normal_01_new(DIM_NUM, ref data, ref seed);
            typeMethods.r8vec_copy(DIM_NUM, v, ref normal);
            //
            //  Compute the length of the vector.
            //
            norm = typeMethods.r8vec_norm(DIM_NUM, normal);
            //
            //  Normalize the vector.
            //
            for (i = 0; i < DIM_NUM; i++)
            {
                normal[i] = normal[i] / norm;
            }
        }

        public static void plane_normal_uniform_nd(int dim_num, ref int seed, double[] pp,
                double[] normal)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_UNIFORM_ND generates a random normal plane in ND.
            //
            //  Discussion:
            //
            //    The normal form of a plane is:
            //
            //      PP is a point on the plane,
            //      N is a normal vector to the plane.
            //
            //    The point PP will be chosen at random inside the unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 November 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double PP[DIM_NUM], a point on the plane.
            //
            //    Output, double NORMAL[DIM_NUM], the unit normal vector.
            //
        {
            int i;
            double norm;
            double[] v;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();
            //
            //  Pick PP as a random point inside the unit sphere in ND.
            //
            v = Burkardt.Ball.Geometry.ball01_sample_nd(dim_num, ref seed);
            typeMethods.r8vec_copy(dim_num, v, ref pp);
            //
            //  Get values from a standard normal distribution.
            //
            v = typeMethods.r8vec_normal_01_new(dim_num, ref data, ref seed);
            typeMethods.r8vec_copy(dim_num, v, ref normal);
            //
            //  Compute the length of the vector.
            //
            norm = typeMethods.r8vec_norm(dim_num, normal);
            //
            //  Normalize the vector.
            //
            for (i = 0; i < dim_num; i++)
            {
                normal[i] = normal[i] / norm;
            }

        }

        public static double[] plane_normal_xyz_to_qr(double[] pp, double[] normal, double[] pq,
                double[] pr, int n, double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_XYZ_TO_QR: XYZ to QR coordinates for a normal form plane.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is:
            //
            //      PP is a point on the plane,
            //      NORMAL is a normal vector to the plane.
            //
            //    Two vectors PQ and PR can be computed with the properties that
            //    * NORMAL, PQ and PR are pairwise orthogonal;
            //    * PQ and PR have unit length;
            //    * every point P in the plane has a "QR" representation
            //      as P = PP + q * PQ + r * PR.
            //
            //    This function is given the XYZ coordinates of a set of points on the
            //    plane, and returns the QR coordinates.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double NORMAL[3], a normal vector N to the plane.  The
            //    vector must not have zero length, but it is not necessary for N
            //    to have unit length.
            //
            //    Input, double PQ[3], a vector of unit length,
            //    perpendicular to the vector N and the vector PR.
            //
            //    Input, double PR[3], a vector of unit length,
            //    perpendicular to the vector N and the vector PQ.
            //
            //    Input, int N, the number of points on the plane.
            //
            //    Input, double XYZ[3*N], the XYZ coordinates of the points.
            //
            //    Output, double PLANE_NORMAL_XYZ_TO_QR[2*N], the QR coordinates
            //    of the points.
            //
        {
            int j;
            double[] qr;

            qr = new double[2 * n];

            for (j = 0; j < n; j++)
            {
                qr[0 + j * 2] = pq[0] * (xyz[0 + j * 3] - pp[0])
                                + pq[1] * (xyz[1 + j * 3] - pp[1])
                                + pq[2] * (xyz[2 + j * 3] - pp[2]);
                qr[1 + j * 2] = pr[0] * (xyz[0 + j * 2] - pp[0])
                                + pr[1] * (xyz[1 + j * 3] - pp[1])
                                + pr[2] * (xyz[2 + j * 3] - pp[2]);
            }

            return qr;
        }

        public static void plane_normal2exp_3d(double[] pp, double[] pn, ref double[] p1,
                ref double[] p2, ref double[] p3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL2EXP_3D converts a normal plane to explicit form in 3D.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is
            //
            //      PP, a point on the plane, and
            //      PN, the unit normal to the plane.
            //
            //    The explicit form of a plane in 3D is
            //
            //      the plane through P1, P2 and P3.
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
            //    Input, double PP(3), a point on the plane.
            //
            //    Input, double PN[3], a normal vector N to the plane.  The
            //    vector must not have zero length, but it is not necessary for N
            //    to have unit length.
            //
            //    Output, double P1[3], P2[3], P3[3], three points that lie on the plane.
            //
        {
            int DIM_NUM = 3;

            double[] pq = new double[DIM_NUM];
            double[] pr = new double[DIM_NUM];

            plane_normal_basis_3d(pp, pn, ref pq, ref pr);

            p1[0] = pp[0];
            p1[1] = pp[1];
            p1[2] = pp[2];

            p2[0] = pp[0] + pq[0];
            p2[1] = pp[1] + pq[1];
            p2[2] = pp[2] + pq[2];

            p3[0] = pp[0] + pr[0];
            p3[1] = pp[1] + pr[1];
            p3[2] = pp[2] + pr[2];

        }

        public static void plane_normal2imp_3d(double[] pp, double[] pn, ref double a, ref double b,
                ref double c, ref double d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL2IMP_3D converts a normal form plane to implicit form in 3D.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is
            //
            //      PP, a point on the plane, and
            //      PN, the unit normal to the plane.
            //
            //    The implicit form of a plane in 3D is
            //
            //      A * X + B * Y + C * Z + D = 0.
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
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double PN[3], the unit normal vector to the plane.
            //
            //    Output, double *A, *B, *C, *D, parameters that define the implicit plane.
            //
        {
            int DIM_NUM = 3;

            a = pn[0];
            b = pn[1];
            c = pn[2];

            d = -typeMethods.r8vec_dot_product(DIM_NUM, pn, pp);

        }

        public static double planes_imp_angle_3d(double a1, double b1, double c1, double d1,
                double a2, double b2, double c2, double d2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANES_IMP_ANGLE_3D: dihedral angle between implicit planes in 3D.
            //
            //  Discussion:
            //
            //    The implicit form of a plane in 3D is:
            //
            //      A * X + B * Y + C * Z + D = 0
            //
            //    If two planes P1 and P2 intersect in a nondegenerate way, then there is a
            //    line of intersection L0.  Consider any plane perpendicular to L0.  The
            //    dihedral angle of P1 and P2 is the angle between the lines L1 and L2, where
            //    L1 is the intersection of P1 and P0, and L2 is the intersection of P2
            //    and P0.
            //
            //    The dihedral angle may also be calculated as the angle between the normal
            //    vectors of the two planes.  Note that if the planes are parallel or
            //    coincide, the normal vectors are identical, and the dihedral angle is 0.
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
            //  Reference:
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Math Tables and Formulae, 30th edition,
            //    Section 4.13, "Planes",
            //    CRC Press, 1996, pages 305-306.
            //
            //  Parameters:
            //
            //    Input, double A1, B1, C1, D1, coefficients that define the first plane.
            //
            //    Input, double A2, B2, C2, D2, coefficients that define the second plane.
            //
            //    Output, double PLANES_IMP_ANGLE_3D, the dihedral angle, in radians,
            //    defined by the two planes.  If either plane is degenerate, or they do
            //    not intersect, or they coincide, then the angle is set to R8_HUGE().
            //    Otherwise, the angle is between 0 and PI.
            //
        {
            double cosine;
            double norm1;
            double norm2;
            double value;

            norm1 = Math.Sqrt(a1 * a1 + b1 * b1 + c1 * c1);

            if (norm1 == 0.0)
            {
                value = typeMethods.r8_huge();
                return value;
            }

            norm2 = Math.Sqrt(a2 * a2 + b2 * b2 + c2 * c2);

            if (norm2 == 0.0)
            {
                value = typeMethods.r8_huge();
                return value;
            }

            cosine = (a1 * a2 + b1 * b2 + c1 * c2) / (norm1 * norm2);

            value = Math.Acos(cosine);
            return value;
        }

    }
}