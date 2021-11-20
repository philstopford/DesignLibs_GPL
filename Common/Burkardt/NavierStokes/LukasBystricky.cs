using System;

namespace Burkardt.NavierStokesNS;

public static class LukasBystricky
{
    public static void resid_lukas(double nu, double rho, int n, double[] x, double[] y,
            double t, ref double[] ur, ref double[] vr, ref double[] pr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESID_LUKAS returns Lukas Bystricky residuals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
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
        int i;

        //
        //  Get the right hand sides.
        //
        double[] f = new double[n];
        double[] g = new double[n];
        double[] h = new double[n];

        rhs_lukas(nu, rho, n, x, y, t, ref f, ref g, ref h);
        //
        //  Form the functions and derivatives of the left hand side.
        //
        for (i = 0; i < n; i++)
        {
            double u = -Math.Cos(Math.PI * x[i]) / Math.PI;
            const double dudt = 0.0;
            double dudx = Math.Sin(Math.PI * x[i]);
            double dudxx = Math.PI * Math.Cos(Math.PI * x[i]);
            const double dudy = 0.0;
            const double dudyy = 0.0;

            double v = -y[i] * Math.Sin(Math.PI * x[i]);
            const double dvdt = 0.0;
            double dvdx = -Math.PI * y[i] * Math.Cos(Math.PI * x[i]);
            double dvdxx = +Math.PI * Math.PI * y[i] * Math.Sin(Math.PI * x[i]);
            double dvdy = -Math.Sin(Math.PI * x[i]);
            const double dvdyy = 0.0;

            const double dpdx = 0.0;
            const double dpdy = 0.0;
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

    public static void rhs_lukas(double nu, double rho, int n, double[] x, double[] y, double t,
            ref double[] f, ref double[] g, ref double[] h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RHS_LUKAS evaluates Lukas Bystricky right hand sides.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double NU, the kinematic viscosity.
        //
        //    double RHO, the fluid density.
        //
        //    int N, the number of nodes.
        //
        //    double X[N], Y[N], the coordinates of nodes.
        //
        //    double T, the current time.
        //
        //  Output:
        //
        //    double F[N], G[N], H[N], the right hand sides.
        //
    {
        int i;

        double[] dpdx = new double[n];
        double[] dpdy = new double[n];
        double[] dudt = new double[n];
        double[] dudx = new double[n];
        double[] dudxx = new double[n];
        double[] dudy = new double[n];
        double[] dudyy = new double[n];
        double[] dvdt = new double[n];
        double[] dvdx = new double[n];
        double[] dvdxx = new double[n];
        double[] dvdy = new double[n];
        double[] dvdyy = new double[n];
        double[] p = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        for (i = 0; i < n; i++)
        {
            u[i] = -Math.Cos(Math.PI * x[i]) / Math.PI;
            dudt[i] = 0.0;
            dudx[i] = Math.Sin(Math.PI * x[i]);
            dudxx[i] = Math.PI * Math.Cos(Math.PI * x[i]);
            dudy[i] = 0.0;
            dudyy[i] = 0.0;

            v[i] = -y[i] * Math.Sin(Math.PI * x[i]);
            dvdt[i] = 0.0;
            dvdx[i] = -Math.PI * y[i] * Math.Cos(Math.PI * x[i]);
            dvdxx[i] = +Math.PI * Math.PI * y[i] * Math.Sin(Math.PI * x[i]);
            dvdy[i] = -Math.Sin(Math.PI * x[i]);
            dvdyy[i] = 0.0;

            p[i] = 0.0;
            dpdx[i] = 0.0;
            dpdy[i] = 0.0;

            f[i] = dudt[i] - nu * (dudxx[i] + dudyy[i])
                   + u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho;

            g[i] = dvdt[i] - nu * (dvdxx[i] + dvdyy[i])
                   + u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho;

            h[i] = dudx[i] + dvdy[i];
        }
    }
        
    public static void uvp_lukas ( double nu, double rho, int n, double[] x, double[] y, 
            double t, ref double[] u, ref double[] v, ref double[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UVP_LUKAS evaluates Lukas Bystricky's solution.
        //
        //  Discussion:
        //
        //    There is no time dependence.
        //
        //    The pressure is 0.
        //
        //    The preferred domain is the unit square.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 March 2015
        //
        //  Author:
        //
        //    John Burkardt
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
            

        for ( i = 0; i < n; i++ )
        {
            u[i] = - Math.Cos ( Math.PI * x[i] ) / Math.PI;
            v[i] = - y[i] *  Math.Sin ( Math.PI * x[i] );
            p[i] = - 0.0;
        }
    }
}