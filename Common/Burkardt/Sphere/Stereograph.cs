using System;

namespace Burkardt.SphereNS
{
    public static class Stereograph
    {
        public static double[] sphere_stereograph(int m, int n, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_STEREOGRAPH computes the stereographic image of points on a sphere.
            //
            //  Discussion:
            //
            //    We start with a sphere of radius 1 and center (0,0,0).
            //
            //    The north pole N = (0,0,1) is the point of tangency to the sphere
            //    of a plane, and the south pole S = (0,0,-1) is the focus for the
            //    stereographic projection.
            //
            //    For any point P on the sphere, the stereographic projection Q of the
            //    point is defined by drawing the line from S through P, and computing
            //    Q as the intersection of this line with the plane.
            //
            //    Actually, we allow the spatial dimension M to be arbitrary.  Values
            //    of M make sense starting with 2.  The north and south poles are
            //    selected as the points (0,0,...,+1) and (0,0,...,-1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    C F Marcus,
            //    The stereographic projection in vector notation,
            //    Mathematics Magazine,
            //    Volume 39, Number 2, March 1966, pages 100-102.
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points.
            //
            //    Input, double P[M*N], a set of points on the unit sphere.
            //
            //    Output, double SPHERE_STEREOGRAPH[M*N], the coordinates of the
            //    image points.
            //
        {
            int i;
            int j;
            double[] q;

            q = new double[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m - 1; i++)
                {
                    q[i + j * m] = 2.0 * p[i + j * m] / (1.0 + p[m - 1 + j * m]);
                }

                q[m - 1 + j * m] = 1.0;
            }

            return q;
        }

        public static double[] sphere_stereograph_inverse(int m, int n, double[] q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_STEREOGRAPH_INVERSE computes stereographic preimages of points.
            //
            //  Discussion:
            //
            //    We start with a sphere of radius 1 and center (0,0,0).
            //
            //    The north pole N = (0,0,1) is the point of tangency to the sphere
            //    of a plane, and the south pole S = (0,0,-1) is the focus for the
            //    stereographic projection.
            //
            //    For any point Q on the plane, the stereographic inverse projection
            //    P of the point is defined by drawing the line from S through Q, and
            //    computing P as the intersection of this line with the sphere.
            //
            //    Actually, we allow the spatial dimension M to be arbitrary.  Values
            //    of M make sense starting with 2.  The north and south poles are
            //    selected as the points (0,0,...,+1) and (0,0,...,-1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    C F Marcus,
            //    The stereographic projection in vector notation,
            //    Mathematics Magazine,
            //    Volume 39, Number 2, March 1966, pages 100-102.
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points.
            //
            //    Input, double Q[M*N], the points, which are presumed to lie
            //    on the plane Z = 1.
            //
            //    Output, double SPHERE_STEREOGRAPH_INVERSE[M*N], the stereographic
            //    inverse projections of the points.
            //
        {
            int i;
            int j;
            double[] p;
            double qn;

            p = new double[m * n];

            for (j = 0; j < n; j++)
            {
                qn = 0.0;
                for (i = 0; i < m - 1; i++)
                {
                    qn = qn + Math.Pow(q[i + j * m], 2);
                }

                for (i = 0; i < m - 1; i++)
                {
                    p[i + j * m] = 4.0 * q[i + j * m] / (4.0 + qn);
                }

                p[m - 1 + j * m] = (4.0 - qn) / (4.0 + qn);
            }

            return p;
        }

        public static double[] sphere_stereograph2(int m, int n, double[] p, double[] focus,
                double[] center)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_STEREOGRAPH2 computes the stereographic image of points on a sphere.
            //
            //  Discussion:
            //
            //    We start with a sphere of center C.
            //
            //    F is a point on the sphere which is the focus of the mapping,
            //    and the antipodal point 2*C-F is the point of tangency
            //    to the sphere of a plane.
            //
            //    For any point P on the sphere, the stereographic projection Q of the
            //    point is defined by drawing the line from F through P, and computing
            //    Q as the intersection of this line with the plane.
            //
            //    The spatial dimension M is arbitrary, but should be at least 2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    C F Marcus,
            //    The stereographic projection in vector notation,
            //    Mathematics Magazine,
            //    Volume 39, Number 2, March 1966, pages 100-102.
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points.
            //
            //    Input, double P[M*N], a set of points on the unit sphere.
            //
            //    Input, double FOCUS[M], the coordinates of the focus point.
            //
            //    Input, double CENTER[M], the coordinates of the center of the sphere.
            //
            //    Output, double SPHERE_STEREOGRAPH2[M*N], the coordinates of the
            //    image points,
            //
        {
            double cf_dot_pf;
            double cf_normsq;
            int i;
            int j;
            double[] q;
            double s;

            q = new double[m * n];

            for (j = 0; j < n; j++)
            {
                cf_normsq = 0.0;
                cf_dot_pf = 0.0;
                for (i = 0; i < m; i++)
                {
                    cf_normsq = cf_normsq + Math.Pow(center[i] - focus[i], 2);
                    cf_dot_pf = cf_dot_pf + (center[i] - focus[i]) * (p[i + j * m] - focus[i]);
                }

                s = 2.0 * cf_normsq / cf_dot_pf;
                for (i = 0; i < m; i++)
                {
                    q[i + j * m] = s * p[i + j * m] + (1.0 - s) * focus[i];
                }
            }

            return q;
        }

        public static double[] sphere_stereograph2_inverse(int m, int n, double[] q, double[] focus,
                double[] center)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_STEREOGRAPH2_INVERSE computes stereographic preimages of points.
            //
            //  Discussion:
            //
            //    We start with a sphere of center C.
            //
            //    F is a point on the sphere which is the focus of the mapping,
            //    and the antipodal point 2*C-F is the point of tangency
            //    to the sphere of a plane.
            //
            //    For any point Q on the plane, the stereographic inverse projection
            //    P of the point is defined by drawing the line from F through Q, and
            //    computing P as the intersection of this line with the sphere.
            //
            //    The spatial dimension M is arbitrary, but should be at least 2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    C F Marcus,
            //    The stereographic projection in vector notation,
            //    Mathematics Magazine,
            //    Volume 39, Number 2, March 1966, pages 100-102.
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points.
            //
            //    Input, double Q[M*N], the points, which are presumed to lie
            //    on the plane.
            //
            //    Input, double FOCUS[M], the coordinates of the focus point.
            //
            //    Input, double CENTER[M], the coordinates of the center of the sphere.
            //
            //    Output, double SPHERE_STEREOGRAPH2_INVERSE[M*N], the stereographic
            //    inverse projections of the points.
            //
        {
            double cf_dot_qf;
            int i;
            int j;
            double[] p;
            double qf_normsq;
            double s;

            p = new double[m * n];

            for (j = 0; j < n; j++)
            {
                cf_dot_qf = 0.0;
                qf_normsq = 0.0;
                for (i = 0; i < m; i++)
                {
                    cf_dot_qf = cf_dot_qf + (center[i] - focus[i]) * (q[i + j * m] - focus[i]);
                    qf_normsq = qf_normsq + Math.Pow(q[i + j * m] - focus[i], 2);
                }

                s = 2.0 * cf_dot_qf / qf_normsq;
                for (i = 0; i < m; i++)
                {
                    p[i + j * m] = s * q[i + j * m] + (1.0 - s) * focus[i];
                }
            }

            return p;
        }
    }
}