using System;

namespace Burkardt.NavierStokesNS;

public static class Taylor
{
    public static void resid_taylor(double nu, double rho, int n, double[] x, double[] y,
            double t, ref double[] ur, ref double[] vr, ref double[] pr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESID_TAYLOR returns Taylor residuals.
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
        //    Geoffrey Taylor,
        //    On the decay of vortices in a viscous fluid,
        //    Philosophical Magazine,
        //    Volume 46, 1923, pages 671-674.
        //
        //    Geoffrey Taylor, A E Green,
        //    Mechanism for the production of small eddies from large ones,
        //    Proceedings of the Royal Society of London, 
        //    Series A, Volume 158, 1937, pages 499-521.
        //
        //  Input:
        //
        //    double NU, the kinematic viscosity.
        //
        //    double RHO, the density.
        //
        //    int N, the number of evaluation points.
        //
        //    double X[N], Y[N], the coordinates of the points.
        //
        //    double T, the time coordinate or coordinates.
        //
        //  Output:
        //
        //    double UR[N], VR[N], PR[N], the residuals in the U, 
        //    V and P equations.
        //
    {
        double dpdx;
        double dpdy;
        double dudt;
        double dudx;
        double dudxx;
        double dudy;
        double dudyy;
        double dvdt;
        double dvdx;
        double dvdxx;
        double dvdy;
        double dvdyy;
        double[] f;
        double[] g;
        double[] h;
        int i;

        double u;
        double v;
        //
        //  Get the right hand sides.
        //
        f = new double[n];
        g = new double[n];
        h = new double[n];

        rhs_taylor(nu, rho, n, x, y, t, ref f, ref g, ref h);
        //
        //  Form the functions and derivatives for the left hand side.
        //
        for (i = 0; i < n; i++)
        {
            u = -
                Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dudx = Math.PI
                   * Math.Sin(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dudxx = Math.PI * Math.PI
                            * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dudy = -Math.PI
                   * Math.Cos(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            dudyy = Math.PI * Math.PI
                            * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dudt = +2.0 * nu * Math.PI * Math.PI
                   * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);

            v =
                Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            dvdx = Math.PI
                   * Math.Cos(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            dvdxx = -Math.PI * Math.PI
                             * Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            dvdy = -Math.PI
                   * Math.Sin(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dvdyy = -Math.PI * Math.PI
                             * Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            dvdt = -2.0 * nu * Math.PI * Math.PI
                   * Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);

            dpdx = +0.5 * rho * Math.PI * Math.Sin(2.0 * Math.PI * x[i]);
            dpdy = +0.5 * rho * Math.PI * Math.Sin(2.0 * Math.PI * y[i]);
            //
            //  Time scaling.
            //
            u *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dudx *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dudxx *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dudy *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dudyy *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dudt *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);

            v *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dvdx *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dvdxx *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dvdy *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dvdyy *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            dvdt *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);

            dpdx *= Math.Exp(-4.0 * Math.PI * Math.PI * nu * t);
            dpdy *= Math.Exp(-4.0 * Math.PI * Math.PI * nu * t);
            //
            //  Evaluate the residuals.
            //
            ur[i] = dudt + u * dudx + v * dudy
                + 1.0 / rho * dpdx - nu * (dudxx + dudyy) - f[i];

            vr[i] = dvdt + u * dvdx + v * dvdy
                + 1.0 / rho * dpdy - nu * (dvdxx + dvdyy) - g[i];

            pr[i] = dudx + dvdy - h[i];
        }
    }

    public static void rhs_taylor(double nu, double rho, int n, double[] x, double[] y,
            double t, ref double[] f, ref double[] g, ref double[] h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RHS_TAYLOR returns Taylor right hand sides.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Geoffrey Taylor,
        //    On the decay of vortices in a viscous fluid,
        //    Philosophical Magazine,
        //    Volume 46, 1923, pages 671-674.
        //
        //    Geoffrey Taylor, A E Green,
        //    Mechanism for the production of small eddies from large ones,
        //    Proceedings of the Royal Society of London, 
        //    Series A, Volume 158, 1937, pages 499-521.
        //
        //  Input:
        //
        //    double NU, the kinematic viscosity.
        //
        //    double RHO, the density.
        //
        //    int N, the number of evaluation points.
        //
        //    double X[N], Y[N], the coordinates of the points.
        //
        //    double T, the time coordinate or coordinates.
        //
        //  Output:
        //
        //    double F[N], G[N], H[N], the right hand sides in the U, 
        //    V and P equations.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            f[i] = 0.0;
            g[i] = 0.0;
            h[i] = 0.0;
        }
    }

    public static void uvp_taylor(double nu, double rho, int n, double[] x, double[] y,
            double t, ref double[] u, ref double[] v, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UVP_TAYLOR evaluates Taylor solutions.
        //
        //  Discussion:
        //
        //    This flow is known as a Taylor-Green vortex.
        //
        //    The given velocity and pressure fields are exact solutions for the 2D 
        //    incompressible time-dependent Navier Stokes equations over any region.
        //
        //    To define a typical problem, one chooses a bounded spatial region 
        //    and a starting time, and then imposes boundary and initial conditions
        //    by referencing the exact solution appropriately.
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
        //    Geoffrey Taylor,
        //    On the decay of vortices in a viscous fluid,
        //    Philosophical Magazine,
        //    Volume 46, 1923, pages 671-674.
        //
        //    Geoffrey Taylor, A E Green,
        //    Mechanism for the production of small eddies from large ones,
        //    Proceedings of the Royal Society of London, 
        //    Series A, Volume 158, 1937, pages 499-521.
        //
        //  Input:
        //
        //    double NU, the kinematic viscosity.
        //
        //    double RHO, the density.
        //
        //    int N, the number of evaluation points.
        //
        //    double X[N], Y[N], the coordinates of the points.
        //
        //    double T, the time coordinate or coordinates.
        //
        //  Output:
        //
        //    double U[N], V[N], P[N], the velocity components and
        //    pressure at each of the points.
        //
    {
        int i;
            

        for (i = 0; i < n; i++)
        {
            u[i] = -Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            v[i] = Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            p[i] = -0.25 * rho
                         * (Math.Cos(2.0 * Math.PI * x[i]) + Math.Cos(2.0 * Math.PI * y[i]));

            u[i] *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            v[i] *= Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
            p[i] *= Math.Exp(-4.0 * Math.PI * Math.PI * nu * t);
        }

    }
}