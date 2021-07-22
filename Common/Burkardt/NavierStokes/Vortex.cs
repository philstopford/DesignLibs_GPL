using System;

namespace Burkardt.NavierStokesNS
{
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
            double p;
            const double r8_pi = 3.141592653589793;
            double u;
            double v;
            //
            //  Get the right hand sides.
            //
            f = new double[n];
            g = new double[n];
            h = new double[n];

            rhs_vortex(nu, rho, n, x, y, t, ref f, ref g, ref h);
            //
            //  Form the functions and derivatives for the left hand side.
            //
            for (i = 0; i < n; i++)
            {
                u = -
                    Math.Cos(r8_pi * x[i]) * Math.Sin(r8_pi * y[i]);
                dudx = r8_pi
                       * Math.Sin(r8_pi * x[i]) * Math.Sin(r8_pi * y[i]);
                dudxx = r8_pi * r8_pi
                              * Math.Cos(r8_pi * x[i]) * Math.Sin(r8_pi * y[i]);
                dudy = -r8_pi
                       * Math.Cos(r8_pi * x[i]) * Math.Cos(r8_pi * y[i]);
                dudyy = r8_pi * r8_pi
                              * Math.Cos(r8_pi * x[i]) * Math.Sin(r8_pi * y[i]);
                dudt = +2.0 * nu * r8_pi * r8_pi
                       * Math.Cos(r8_pi * x[i]) * Math.Sin(r8_pi * y[i]);

                v =
                    Math.Sin(r8_pi * x[i]) * Math.Cos(r8_pi * y[i]);
                dvdx = r8_pi
                       * Math.Cos(r8_pi * x[i]) * Math.Cos(r8_pi * y[i]);
                dvdxx = -r8_pi * r8_pi
                               * Math.Sin(r8_pi * x[i]) * Math.Cos(r8_pi * y[i]);
                dvdy = -r8_pi
                       * Math.Sin(r8_pi * x[i]) * Math.Sin(r8_pi * y[i]);
                dvdyy = -r8_pi * r8_pi
                               * Math.Sin(r8_pi * x[i]) * Math.Cos(r8_pi * y[i]);
                dvdt = -2.0 * nu * r8_pi * r8_pi
                       * Math.Sin(r8_pi * x[i]) * Math.Cos(r8_pi * y[i]);

                p = -0.25 * rho *
                    (Math.Cos(2.0 * r8_pi * x[i]) + Math.Cos(2.0 * r8_pi * y[i]));
                dpdx = +0.5 * rho * r8_pi * Math.Sin(2.0 * r8_pi * x[i]);
                dpdy = +0.5 * rho * r8_pi * Math.Sin(2.0 * r8_pi * y[i]);
                //
                //  Time scaling.
                //
                u = u * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dudx = dudx * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dudxx = dudxx * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dudy = dudy * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dudyy = dudyy * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dudt = dudt * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);

                v = v * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dvdx = dvdx * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dvdxx = dvdxx * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dvdy = dvdy * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dvdyy = dvdyy * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);
                dvdt = dvdt * Math.Exp(-2.0 * r8_pi * r8_pi * nu * t);

                p = p * Math.Exp(-4.0 * r8_pi * r8_pi * nu * t);
                dpdx = dpdx * Math.Exp(-4.0 * r8_pi * r8_pi * nu * t);
                dpdy = dpdy * Math.Exp(-4.0 * r8_pi * r8_pi * nu * t);
                //
                //  Evaluate the residuals.
                //
                ur[i] = dudt + u * dudx + v * dudy
                    + (1.0 / rho) * dpdx - nu * (dudxx + dudyy) - f[i];

                vr[i] = dvdt + u * dvdx + v * dvdy
                    + (1.0 / rho) * dpdy - nu * (dvdxx + dvdyy) - g[i];

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
            const double r8_pi = 3.141592653589793;

            for ( i = 0; i < n; i++ )
            {
                f[i] = - 2.0 * nu * r8_pi * r8_pi 
                       * Math.Cos ( r8_pi * x[i] ) * Math.Sin ( r8_pi * y[i] );
                g[i] =   2.0 * nu * r8_pi * r8_pi 
                         * Math.Sin ( r8_pi * x[i] ) * Math.Cos ( r8_pi * y[i] );
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
            const double r8_pi = 3.141592653589793;

            for ( i = 0; i < n; i++ )
            {
                u[i] = - Math.Cos ( r8_pi * x[i] ) * Math.Sin ( r8_pi * y[i] );
                v[i] =   Math.Sin ( r8_pi * x[i] ) * Math.Cos ( r8_pi * y[i] );
                p[i] = - 0.25 * rho 
                              * ( Math.Cos ( 2.0 * r8_pi * x[i] ) + Math.Cos ( 2.0 * r8_pi * y[i] ) );
            }
        }
    }
}