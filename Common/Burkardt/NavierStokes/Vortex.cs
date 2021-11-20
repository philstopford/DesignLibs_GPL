using System;

namespace Burkardt.NavierStokesNS;

public static class Vortex
{
    public static void resid_vortex(double nu, double rho, int n, double[] x, double[] y,
            double t, ref double[] ur, ref double[] vr, ref double[] pr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESID_VORTEX returns Vortex residuals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
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

        rhs_vortex(nu, rho, n, x, y, t, ref f, ref g, ref h);
        //
        //  Form the functions and derivatives for the left hand side.
        //
        for (i = 0; i < n; i++)
        {
            double u = -
                Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            double dudx = Math.PI
                          * Math.Sin(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            double dudxx = Math.PI * Math.PI
                                   * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            double dudy = -Math.PI
                          * Math.Cos(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            double dudyy = Math.PI * Math.PI
                                   * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            double dudt = +2.0 * nu * Math.PI * Math.PI
                          * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);

            double v = Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            double dvdx = Math.PI
                          * Math.Cos(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            double dvdxx = -Math.PI * Math.PI
                                    * Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            double dvdy = -Math.PI
                          * Math.Sin(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            double dvdyy = -Math.PI * Math.PI
                                    * Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            double dvdt = -2.0 * nu * Math.PI * Math.PI
                          * Math.Sin(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);

            double dpdx = +0.5 * rho * Math.PI * Math.Sin(2.0 * Math.PI * x[i]);
            double dpdy = +0.5 * rho * Math.PI * Math.Sin(2.0 * Math.PI * y[i]);
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
        
    public static void rhs_vortex ( double nu, double rho, int n, double[] x, double[] y, 
            double t, ref double[] f, ref double[] g, ref double[] h )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RHS_VORTEX returns Vortex right hand sides.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
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
        //    double F[N], G[N], H[N], the right hand sides in the U, 
        //    V and P equations.
        //
    {
        int i;
            

        for ( i = 0; i < n; i++ )
        {
            f[i] = - 2.0 * nu * Math.PI * Math.PI 
                   * Math.Cos ( Math.PI * x[i] ) * Math.Sin ( Math.PI * y[i] );
            g[i] =   2.0 * nu * Math.PI * Math.PI 
                     * Math.Sin ( Math.PI * x[i] ) * Math.Cos ( Math.PI * y[i] );
            h[i] = 0.0;
        }
    }
        
    public static void uvp_vortex ( double nu, double rho, int n, double[] x, double[] y, 
            double t, ref double[] u, ref double[] v, ref double[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UVP_VORTEX evaluates Vortex solutions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2015
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
            u[i] = - Math.Cos ( Math.PI * x[i] ) * Math.Sin ( Math.PI * y[i] );
            v[i] =   Math.Sin ( Math.PI * x[i] ) * Math.Cos ( Math.PI * y[i] );
            p[i] = - 0.25 * rho 
                          * ( Math.Cos ( 2.0 * Math.PI * x[i] ) + Math.Cos ( 2.0 * Math.PI * y[i] ) );
        }
    }
}