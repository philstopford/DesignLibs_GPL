using System;

namespace Burkardt.NavierStokesNS
{
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

            rhs_lukas(nu, rho, n, x, y, t, ref f, ref g, ref h);
            //
            //  Form the functions and derivatives of the left hand side.
            //
            for (i = 0; i < n; i++)
            {
                u = -Math.Cos(Math.PI * x[i]) / Math.PI;
                dudt = 0.0;
                dudx = Math.Sin(Math.PI * x[i]);
                dudxx = Math.PI * Math.Cos(Math.PI * x[i]);
                dudy = 0.0;
                dudyy = 0.0;

                v = -y[i] * Math.Sin(Math.PI * x[i]);
                dvdt = 0.0;
                dvdx = -Math.PI * y[i] * Math.Cos(Math.PI * x[i]);
                dvdxx = +Math.PI * Math.PI * y[i] * Math.Sin(Math.PI * x[i]);
                dvdy = -Math.Sin(Math.PI * x[i]);
                dvdyy = 0.0;

                dpdx = 0.0;
                dpdy = 0.0;
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
            double[] dpdx;
            double[] dpdy;
            double[] dudt;
            double[] dudx;
            double[] dudxx;
            double[] dudy;
            double[] dudyy;
            double[] dvdt;
            double[] dvdx;
            double[] dvdxx;
            double[] dvdy;
            double[] dvdyy;
            int i;
            double[] p;
            
            double[] u;
            double[] v;

            dpdx = new double[n];
            dpdy = new double[n];
            dudt = new double[n];
            dudx = new double[n];
            dudxx = new double[n];
            dudy = new double[n];
            dudyy = new double[n];
            dvdt = new double[n];
            dvdx = new double[n];
            dvdxx = new double[n];
            dvdy = new double[n];
            dvdyy = new double[n];
            p = new double[n];
            u = new double[n];
            v = new double[n];

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
}