using System;

namespace Burkardt.NavierStokesNS
{
    public static class Ethier
    {
        public static void uvwp_ethier(double a, double d, int n, double[] x, double[] y,
                double[] z, double[] t, ref double[] u, ref double[] v, ref double[] w, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UVWP_ETHIER evaluates the Ethier solution.
            //
            //  Discussion:
            //
            //    The reference asserts that the given velocity and pressure fields
            //    are exact solutions for the 3D incompressible time-dependent
            //    Navier Stokes equations over any region.
            //
            //    To define a typical problem, one chooses a bounded spatial region 
            //    and a starting time, and then imposes boundary and initial conditions
            //    by referencing the exact solution appropriately.
            //
            //    In the reference paper, a calculation is made for the cube centered
            //    at (0,0,0) with a "radius" of 1 unit, and over the time interval
            //    from t = 0 to t = 0.1, with parameters a = PI/4 and d = PI/2,
            //    and with Dirichlet boundary conditions on all faces of the cube.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    C Ross Ethier, David Steinman,
            //    Exact fully 3D Navier-Stokes solutions for benchmarking,
            //    International Journal for Numerical Methods in Fluids,
            //    Volume 19, Number 5, March 1994, pages 369-375.
            //
            //  Parameters:
            //
            //    Input, double A, D, the parameters.  Sample values are A = PI/4 and
            //    D = PI/2.
            //
            //    Input, int N, the number of points at which the solution is to
            //    be evaluated.
            //
            //    Input, double X[N], Y[N], Z[N], the coordinates of the points.
            //
            //    Input, double T[N], the time coordinates.
            //
            //    Output, double U[N], V[N], W[N], P[N], the velocity components and
            //    pressure at each of the points.
            //
        {
            double cxy;
            double cyz;
            double czx;
            double e2t;
            double ex;
            double exy;
            double ey;
            double eyz;
            double ez;
            double ezx;
            int i;
            double sxy;
            double syz;
            double szx;

            for (i = 0; i < n; i++)
            {
                ex = Math.Exp(a * x[i]);
                ey = Math.Exp(a * y[i]);
                ez = Math.Exp(a * z[i]);

                e2t = Math.Exp(-d * d * t[i]);

                exy = Math.Exp(a * (x[i] + y[i]));
                eyz = Math.Exp(a * (y[i] + z[i]));
                ezx = Math.Exp(a * (z[i] + x[i]));

                sxy = Math.Sin(a * x[i] + d * y[i]);
                syz = Math.Sin(a * y[i] + d * z[i]);
                szx = Math.Sin(a * z[i] + d * x[i]);

                cxy = Math.Cos(a * x[i] + d * y[i]);
                cyz = Math.Cos(a * y[i] + d * z[i]);
                czx = Math.Cos(a * z[i] + d * x[i]);

                u[i] = -a * (ex * syz + ez * cxy) * e2t;
                v[i] = -a * (ey * szx + ex * cyz) * e2t;
                w[i] = -a * (ez * sxy + ey * czx) * e2t;
                p[i] = 0.5 * a * a * e2t * e2t * (
                    +ex * ex + 2.0 * sxy * czx * eyz
                             + ey * ey + 2.0 * syz * cxy * ezx
                             + ez * ez + 2.0 * szx * cyz * exy);
            }
        }

        public static void resid_ethier(double a, double d, int n, double[] x, double[] y,
                double[] z, double[] t, ref double[] ur, ref double[] vr, ref double[] wr, ref double[] pr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RESID_ETHIER evaluates the Ethier residual.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    C Ross Ethier, David Steinman,
            //    Exact fully 3D Navier-Stokes solutions for benchmarking,
            //    International Journal for Numerical Methods in Fluids,
            //    Volume 19, Number 5, March 1994, pages 369-375.
            //
            //  Parameters:
            //
            //    Input, double A, D, the parameters.  Sample values are A = PI/4 
            //    and D = PI/2.
            //
            //    Input, int N, the number of points at which the solution 
            //    is to be evaluated.
            //
            //    Input, double X[N], Y[N], Z[N], the coordinates of the points.
            //
            //    Input, double T[N], the time coordinates.
            //
            //    Output, double UR[N], VR[N], WR[N], PR[N], the residuals.
            //
        {
            double cxy;
            double cyz;
            double czx;
            double e2t;
            double e2x;
            double e2y;
            double e2z;
            double e4t;
            double ex;
            double exy;
            double ey;
            double eyz;
            double ez;
            double ezx;
            int i;
            double px;
            double py;
            double pz;
            double sxy;
            double syz;
            double szx;
            double u;
            double ut;
            double ux;
            double uxx;
            double uy;
            double uyy;
            double uz;
            double uzz;
            double v;
            double vt;
            double vx;
            double vxx;
            double vy;
            double vyy;
            double vz;
            double vzz;
            double w;
            double wt;
            double wx;
            double wxx;
            double wy;
            double wyy;
            double wz;
            double wzz;
            //
            //  Make some temporaries.
            //
            for (i = 0; i < n; i++)
            {
                ex = Math.Exp(a * x[i]);
                ey = Math.Exp(a * y[i]);
                ez = Math.Exp(a * z[i]);

                e2x = Math.Exp(2.0 * a * x[i]);
                e2y = Math.Exp(2.0 * a * y[i]);
                e2z = Math.Exp(2.0 * a * z[i]);

                e2t = Math.Exp(-d * d * t[i]);
                e4t = Math.Exp(-2.0 * d * d * t[i]);

                exy = Math.Exp(a * (x[i] + y[i]));
                eyz = Math.Exp(a * (y[i] + z[i]));
                ezx = Math.Exp(a * (z[i] + x[i]));

                sxy = Math.Sin(a * x[i] + d * y[i]);
                syz = Math.Sin(a * y[i] + d * z[i]);
                szx = Math.Sin(a * z[i] + d * x[i]);

                cxy = Math.Cos(a * x[i] + d * y[i]);
                cyz = Math.Cos(a * y[i] + d * z[i]);
                czx = Math.Cos(a * z[i] + d * x[i]);
                //
                //  Form the functions and derivatives.
                //
                u = -a * (ex * syz
                          + ez * cxy) * e2t;
                ux = -a * (a * ex * syz
                           - a * ez * sxy) * e2t;
                uxx = -a * (a * a * ex * syz
                            - a * a * ez * cxy) * e2t;
                uy = -a * (a * ex * cyz
                           - d * ez * sxy) * e2t;
                uyy = -a * (-a * a * ex * syz
                            - d * d * ez * cxy) * e2t;
                uz = -a * (d * ex * cyz
                           + a * ez * cxy) * e2t;
                uzz = -a * (-d * d * ex * syz
                            + a * a * ez * cxy) * e2t;
                ut = +d * d * a * (ex * syz
                                   + ez * cxy) * e2t;

                v = -a * (ey * szx
                          + ex * cyz) * e2t;
                vx = -a * (d * ey * czx
                           + a * ex * cyz) * e2t;
                vxx = -a * (-d * d * ey * szx
                            + a * a * ex * cyz) * e2t;
                vy = -a * (a * ey * szx
                           - a * ex * syz) * e2t;
                vyy = -a * (a * a * ey * szx
                            - a * a * ex * cyz) * e2t;
                vz = -a * (a * ey * czx
                           - d * ex * syz) * e2t;
                vzz = -a * (-a * a * ey * szx
                            - d * d * ex * cyz) * e2t;
                vt = +d * d * a * (ey * szx
                                   + ex * cyz) * e2t;

                w = -a * (ez * sxy
                          + ey * czx) * e2t;
                wx = -a * (a * ez * cxy
                           - d * ey * szx) * e2t;
                wxx = -a * (-a * a * ez * sxy
                            - d * d * ey * czx) * e2t;
                wy = -a * (d * ez * cxy
                           + a * ey * czx) * e2t;
                wyy = -a * (-d * d * ez * sxy
                            + a * a * ey * czx) * e2t;
                wz = -a * (a * ez * sxy
                           - a * ey * szx) * e2t;
                wzz = -a * (a * a * ez * sxy
                            - a * a * ey * czx) * e2t;
                wt = +d * d * a * (ez * sxy
                                   + ey * czx) * e2t;

                //  p = - 0.5 * a * a * e4t * ( 
                //    + e2x + 2.0 * sxy * czx * eyz 
                //    + e2y + 2.0 * syz * cxy * ezx 
                //    + e2z + 2.0 * szx * cyz * exy );

                px = -0.5 * a * a * e4t * (
                    +2.0 * a * e2x
                    + 2.0 * a * cxy * czx * eyz
                    - 2.0 * d * sxy * szx * eyz
                    - 2.0 * a * syz * sxy * ezx
                    + 2.0 * a * syz * cxy * ezx
                    + 2.0 * d * czx * cyz * exy
                    + 2.0 * a * szx * cyz * exy);

                py = -0.5 * a * a * e4t * (
                    +2.0 * d * cxy * czx * eyz
                    + 2.0 * a * sxy * czx * eyz
                    + 2.0 * a * e2y
                    + 2.0 * a * cyz * cxy * ezx
                    - 2.0 * d * syz * sxy * ezx
                    - 2.0 * a * szx * syz * exy
                    + 2.0 * a * szx * cyz * exy);

                pz = -0.5 * a * a * e4t * (
                    -2.0 * a * sxy * szx * eyz
                    + 2.0 * a * sxy * czx * eyz
                    + 2.0 * d * cyz * cxy * ezx
                    + 2.0 * a * syz * cxy * ezx
                    + 2.0 * a * e2z
                    + 2.0 * a * czx * cyz * exy
                    - 2.0 * d * szx * syz * exy);
                //
                //  Evaluate the residuals.
                //
                ur[i] = ut
                        + u * ux + v * uy + w * uz + px
                        - (uxx + uyy + uzz);

                vr[i] = vt
                        + u * vx + v * vy + w * vz + py
                        - (vxx + vyy + vzz);

                wr[i] = wt
                        + u * wx + v * wy + w * wz + pz
                        - (wxx + wyy + wzz);

                pr[i] = ux + vy + wz;
            }
        }
    }
}