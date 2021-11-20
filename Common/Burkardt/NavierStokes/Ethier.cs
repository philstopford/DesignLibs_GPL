using System;

namespace Burkardt.NavierStokesNS;

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
        int i;

        for (i = 0; i < n; i++)
        {
            double ex = Math.Exp(a * x[i]);
            double ey = Math.Exp(a * y[i]);
            double ez = Math.Exp(a * z[i]);

            double e2t = Math.Exp(-d * d * t[i]);

            double exy = Math.Exp(a * (x[i] + y[i]));
            double eyz = Math.Exp(a * (y[i] + z[i]));
            double ezx = Math.Exp(a * (z[i] + x[i]));

            double sxy = Math.Sin(a * x[i] + d * y[i]);
            double syz = Math.Sin(a * y[i] + d * z[i]);
            double szx = Math.Sin(a * z[i] + d * x[i]);

            double cxy = Math.Cos(a * x[i] + d * y[i]);
            double cyz = Math.Cos(a * y[i] + d * z[i]);
            double czx = Math.Cos(a * z[i] + d * x[i]);

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
        int i;
        //
        //  Make some temporaries.
        //
        for (i = 0; i < n; i++)
        {
            double ex = Math.Exp(a * x[i]);
            double ey = Math.Exp(a * y[i]);
            double ez = Math.Exp(a * z[i]);

            double e2x = Math.Exp(2.0 * a * x[i]);
            double e2y = Math.Exp(2.0 * a * y[i]);
            double e2z = Math.Exp(2.0 * a * z[i]);

            double e2t = Math.Exp(-d * d * t[i]);
            double e4t = Math.Exp(-2.0 * d * d * t[i]);

            double exy = Math.Exp(a * (x[i] + y[i]));
            double eyz = Math.Exp(a * (y[i] + z[i]));
            double ezx = Math.Exp(a * (z[i] + x[i]));

            double sxy = Math.Sin(a * x[i] + d * y[i]);
            double syz = Math.Sin(a * y[i] + d * z[i]);
            double szx = Math.Sin(a * z[i] + d * x[i]);

            double cxy = Math.Cos(a * x[i] + d * y[i]);
            double cyz = Math.Cos(a * y[i] + d * z[i]);
            double czx = Math.Cos(a * z[i] + d * x[i]);
            //
            //  Form the functions and derivatives.
            //
            double u = -a * (ex * syz
                             + ez * cxy) * e2t;
            double ux = -a * (a * ex * syz
                              - a * ez * sxy) * e2t;
            double uxx = -a * (a * a * ex * syz
                               - a * a * ez * cxy) * e2t;
            double uy = -a * (a * ex * cyz
                              - d * ez * sxy) * e2t;
            double uyy = -a * (-a * a * ex * syz
                               - d * d * ez * cxy) * e2t;
            double uz = -a * (d * ex * cyz
                              + a * ez * cxy) * e2t;
            double uzz = -a * (-d * d * ex * syz
                               + a * a * ez * cxy) * e2t;
            double ut = +d * d * a * (ex * syz
                                      + ez * cxy) * e2t;

            double v = -a * (ey * szx
                             + ex * cyz) * e2t;
            double vx = -a * (d * ey * czx
                              + a * ex * cyz) * e2t;
            double vxx = -a * (-d * d * ey * szx
                               + a * a * ex * cyz) * e2t;
            double vy = -a * (a * ey * szx
                              - a * ex * syz) * e2t;
            double vyy = -a * (a * a * ey * szx
                               - a * a * ex * cyz) * e2t;
            double vz = -a * (a * ey * czx
                              - d * ex * syz) * e2t;
            double vzz = -a * (-a * a * ey * szx
                               - d * d * ex * cyz) * e2t;
            double vt = +d * d * a * (ey * szx
                                      + ex * cyz) * e2t;

            double w = -a * (ez * sxy
                             + ey * czx) * e2t;
            double wx = -a * (a * ez * cxy
                              - d * ey * szx) * e2t;
            double wxx = -a * (-a * a * ez * sxy
                               - d * d * ey * czx) * e2t;
            double wy = -a * (d * ez * cxy
                              + a * ey * czx) * e2t;
            double wyy = -a * (-d * d * ez * sxy
                               + a * a * ey * czx) * e2t;
            double wz = -a * (a * ez * sxy
                              - a * ey * szx) * e2t;
            double wzz = -a * (a * a * ez * sxy
                               - a * a * ey * czx) * e2t;
            double wt = +d * d * a * (ez * sxy
                                      + ey * czx) * e2t;

            //  p = - 0.5 * a * a * e4t * ( 
            //    + e2x + 2.0 * sxy * czx * eyz 
            //    + e2y + 2.0 * syz * cxy * ezx 
            //    + e2z + 2.0 * szx * cyz * exy );

            double px = -0.5 * a * a * e4t * (
                +2.0 * a * e2x
                + 2.0 * a * cxy * czx * eyz
                - 2.0 * d * sxy * szx * eyz
                - 2.0 * a * syz * sxy * ezx
                + 2.0 * a * syz * cxy * ezx
                + 2.0 * d * czx * cyz * exy
                + 2.0 * a * szx * cyz * exy);

            double py = -0.5 * a * a * e4t * (
                +2.0 * d * cxy * czx * eyz
                + 2.0 * a * sxy * czx * eyz
                + 2.0 * a * e2y
                + 2.0 * a * cyz * cxy * ezx
                - 2.0 * d * syz * sxy * ezx
                - 2.0 * a * szx * syz * exy
                + 2.0 * a * szx * cyz * exy);

            double pz = -0.5 * a * a * e4t * (
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