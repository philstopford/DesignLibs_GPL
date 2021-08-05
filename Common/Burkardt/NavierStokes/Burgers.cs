using System;
using Burkardt.FullertonFnLib;
using Burkardt.Types;

namespace Burkardt.NavierStokesNS
{
    public static class Burgers
    {
        public static void uvwp_burgers ( double nu, int n, double[] x, double[] y, 
        double[] z, double[] t, ref double[] u, ref double[] v, ref double[] w, ref double[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UVWP_BURGERS evaluates the Burgers solution.
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
        //    Martin Bazant, Henry Moffatt,
        //    Exact solutions of the Navier-Stokes equations having steady vortex 
        //    structures,
        //    Journal of Fluid Mechanics,
        //    Volume 541, pages 55-64, 2005.
        //
        //    Johannes Burgers,
        //    A mathematical model illustrating the theory of turbulence,
        //    Advances in Applied Mechanics,
        //    Volume 1, pages 171-199, 1948.
        //
        //  Parameters:
        //
        //    Input, double NU, the kinematic viscosity.
        //
        //    Input, int N, the number of points at which the solution 
        //    is to be evaluated.
        //
        //    Input, double X[N], Y[N], Z[N], the coordinates of the points.
        //
        //    Input, double T[N], the time coordinates.
        //
        //    Output, double U[N], V[N], W[N], P[N], the solution values.
        //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                //
                //  Form the functions and derivatives.
                //
                u[i] =   2.0 * x[i];
                v[i] =   - 2.0 * y[i];
                w[i] =   Helpers.Erf( y[i] / Math.Sqrt ( nu ) );
                p[i] = - 2.0 * ( x[i] * x[i] + y[i] * y[i] );
            }
        }
        
        public static void resid_burgers(double nu, int n, double[] x, double[] y,
        double[] z, double[] t, ref double[] ur, ref double[] vr, ref double[] wr, ref double[] pr )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESID_BURGERS evaluates the Burgers residual.
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
        //    Martin Bazant, Henry Moffatt,
        //    Exact solutions of the Navier-Stokes equations having steady vortex 
        //    structures,
        //    Journal of Fluid Mechanics,
        //    Volume 541, pages 55-64, 2005.
        //
        //    Johannes Burgers,
        //    A mathematical model illustrating the theory of turbulence,
        //    Advances in Applied Mechanics,
        //    Volume 1, pages 171-199, 1948.
        //
        //  Parameters:
        //
        //    Input, double NU, the kinematic viscosity.
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
            double px;
            double py;
            double pz;
            const double r8_pi = 3.141592653589793;
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

            for (i = 0; i < n; i++)
            {
                //
                //  Form the functions and derivatives.
                //
                u = 2.0 * x[i];
                ux = 2.0;
                uxx = 0.0;
                uy = 0.0;
                uyy = 0.0;
                uz = 0.0;
                uzz = 0.0;
                ut = 0.0;

                v = -2.0 * y[i];
                vx = 0.0;
                vxx = 0.0;
                vy = -2.0;
                vyy = 0.0;
                vz = 0.0;
                vzz = 0.0;
                vt = 0.0;

                w = Helpers.Erf(y[i] / Math.Sqrt(nu));
                wx = 0.0;
                wxx = 0.0;
                wy = 2.0 * Math.Sqrt(1.0 / nu / r8_pi) * Math.Exp(-y[i] * y[i] / nu);
                wyy = -4.0 * Math.Sqrt(1.0 / nu / r8_pi) * y[i] * Math.Exp(-y[i] * y[i] / nu) / nu;
                wz = 0.0;
                wzz = 0.0;
                wt = 0.0;

                //  p = - 2.0 * ( x[i] * x[i] + y[i] * y[i] );
                px = -4.0 * x[i];
                py = -4.0 * y[i];
                pz = 0.0;
                //
                //  Evaluate the residuals.
                //
                ur[i] = ut
                        + u * ux + v * uy + w * uz + px
                        - nu * (uxx + uyy + uzz);

                vr[i] = vt
                        + u * vx + v * vy + w * vz + py
                        - nu * (vxx + vyy + vzz);

                wr[i] = wt
                        + u * wx + v * wy + w * wz + pz
                        - nu * (wxx + wyy + wzz);

                pr[i] = ux + vy + wz;
            }

            return;
        }
    }
}