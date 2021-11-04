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

                p = -0.25 * rho *
                    (Math.Cos(2.0 * Math.PI * x[i]) + Math.Cos(2.0 * Math.PI * y[i]));
                dpdx = +0.5 * rho * Math.PI * Math.Sin(2.0 * Math.PI * x[i]);
                dpdy = +0.5 * rho * Math.PI * Math.Sin(2.0 * Math.PI * y[i]);
                //
                //  Time scaling.
                //
                u = u * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dudx = dudx * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dudxx = dudxx * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dudy = dudy * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dudyy = dudyy * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dudt = dudt * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);

                v = v * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dvdx = dvdx * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dvdxx = dvdxx * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dvdy = dvdy * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dvdyy = dvdyy * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);
                dvdt = dvdt * Math.Exp(-2.0 * Math.PI * Math.PI * nu * t);

                p = p * Math.Exp(-4.0 * Math.PI * Math.PI * nu * t);
                dpdx = dpdx * Math.Exp(-4.0 * Math.PI * Math.PI * nu * t);
                dpdy = dpdy * Math.Exp(-4.0 * Math.PI * Math.PI * nu * t);
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
}