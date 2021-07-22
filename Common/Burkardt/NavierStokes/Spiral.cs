using System;

namespace Burkardt.NavierStokesNS
{
    public static class Spiral
    {
        public static void resid_spiral(double nu, double rho, int n, double[] x, double[] y,
                double t, ref double[] ur, ref double[] vr, ref double[] pr)

            //****************************************************************************//
            //
            //  Purpose:
            //
            //    RESID_SPIRAL evaluates Spiral residuals.
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
            //    Maxim Olshanskii, Leo Rebholz,
            //    Application of barycenter refined meshes in linear elasticity
            //    and incompressible fluid dynamics,
            //    ETNA: Electronic Transactions in Numerical Analysis, 
            //    Volume 38, pages 258-274, 2011.
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
            //    double UR[N], VR[N], PR[N], the right hand sides.
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
            double[] f;
            double[] g;
            double[] h;
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
            f = new double[n];
            g = new double[n];
            h = new double[n];
            p = new double[n];
            u = new double[n];
            v = new double[n];
            //
            //  Get the right hand side functions.
            //
            rhs_spiral(nu, rho, n, x, y, t, ref f, ref g, ref h);
            //
            //  Form the functions and derivatives for the left hand side.
            //
            for (i = 0; i < n; i++)
            {
                u[i] = (1.0 + nu * t) * 2.0
                                      * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                                      * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudt[i] = nu * 2.0
                             * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                             * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudx[i] = (1.0 + nu * t) * 2.0
                                         * (4.0 * Math.Pow(x[i], 3) - 6.0 * Math.Pow(x[i], 2) + 2.0 * x[i])
                                         * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudxx[i] = (1.0 + nu * t) * 2.0
                                          * (12.0 * Math.Pow(x[i], 2) - 12.0 * x[i] + 2.0)
                                          * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudy[i] = (1.0 + nu * t) * 2.0
                                         * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                                         * (6.0 * Math.Pow(y[i], 2) - 6.0 * y[i] + 1.0);

                dudyy[i] = (1.0 + nu * t) * 2.0
                                          * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                                          * (12.0 * y[i] - 6.0);

                v[i] = -(1.0 + nu * t) * 2.0
                                       * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                                       * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdt[i] = -nu * 2.0
                              * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                              * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdx[i] = -(1.0 + nu * t) * 2.0
                                          * (6.0 * Math.Pow(x[i], 2) - 6.0 * x[i] + 1.0)
                                          * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdxx[i] = -(1.0 + nu * t) * 2.0
                                           * (12.0 * x[i] - 6.0)
                                           * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdy[i] = -(1.0 + nu * t) * 2.0
                                          * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                                          * (4.0 * Math.Pow(y[i], 3) - 6.0 * Math.Pow(y[i], 2) + 2.0 * y[i]);

                dvdyy[i] = -(1.0 + nu * t) * 2.0
                                           * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                                           * (12.0 * Math.Pow(y[i], 2) - 12.0 * y[i] + 2.0);

                p[i] = rho * y[i];
                dpdx[i] = 0.0;
                dpdy[i] = rho;

                ur[i] = dudt[i] - nu * (dudxx[i] + dudyy[i])
                    + u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho - f[i];

                vr[i] = dvdt[i] - nu * (dvdxx[i] + dvdyy[i])
                    + u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho - g[i];

                pr[i] = dudx[i] + dvdy[i] - h[i];
            }
        }

        public static void rhs_spiral(double nu, double rho, int n, double[] x, double[] y,
                double t, ref double[] f, ref double[] g, ref double[] h)

            //****************************************************************************//
            //
            //  Purpose:
            //
            //    RHS_SPIRAL evaluates Spiral right hand sides.
            //
            //  Discussion:
            //
            //    The right hand side is artificially determined by the requirement
            //    that the specified values of U, V and P satisfy the discretized
            //    Navier Stokes and continuity equations.
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
            //    Maxim Olshanskii, Leo Rebholz,
            //    Application of barycenter refined meshes in linear elasticity
            //    and incompressible fluid dynamics,
            //    ETNA: Electronic Transactions in Numerical Analysis, 
            //    Volume 38, pages 258-274, 2011.
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
                u[i] = (1.0 + nu * t) * 2.0
                                      * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                                      * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudt[i] = nu * 2.0
                             * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                             * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudx[i] = (1.0 + nu * t) * 2.0
                                         * (4.0 * Math.Pow(x[i], 3) - 6.0 * Math.Pow(x[i], 2) + 2.0 * x[i])
                                         * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudxx[i] = (1.0 + nu * t) * 2.0
                                          * (12.0 * Math.Pow(x[i], 2) - 12.0 * x[i] + 2.0)
                                          * (2.0 * Math.Pow(y[i], 3) - 3.0 * Math.Pow(y[i], 2) + y[i]);

                dudy[i] = (1.0 + nu * t) * 2.0
                                         * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                                         * (6.0 * Math.Pow(y[i], 2) - 6.0 * y[i] + 1.0);

                dudyy[i] = (1.0 + nu * t) * 2.0
                                          * (Math.Pow(x[i], 4) - 2.0 * Math.Pow(x[i], 3) + Math.Pow(x[i], 2))
                                          * (12.0 * y[i] - 6.0);

                v[i] = -(1.0 + nu * t) * 2.0
                                       * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                                       * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdt[i] = -nu * 2.0
                              * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                              * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdx[i] = -(1.0 + nu * t) * 2.0
                                          * (6.0 * Math.Pow(x[i], 2) - 6.0 * x[i] + 1.0)
                                          * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdxx[i] = -(1.0 + nu * t) * 2.0
                                           * (12.0 * x[i] - 6.0)
                                           * (Math.Pow(y[i], 4) - 2.0 * Math.Pow(y[i], 3) + Math.Pow(y[i], 2));

                dvdy[i] = -(1.0 + nu * t) * 2.0
                                          * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                                          * (4.0 * Math.Pow(y[i], 3) - 6.0 * Math.Pow(y[i], 2) + 2.0 * y[i]);

                dvdyy[i] = -(1.0 + nu * t) * 2.0
                                           * (2.0 * Math.Pow(x[i], 3) - 3.0 * Math.Pow(x[i], 2) + x[i])
                                           * (12.0 * Math.Pow(y[i], 2) - 12.0 * y[i] + 2.0);

                p[i] = rho * y[i];
                dpdx[i] = 0.0;
                dpdy[i] = rho;

                f[i] = dudt[i] - nu * (dudxx[i] + dudyy[i])
                       + u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho;

                g[i] = dvdt[i] - nu * (dvdxx[i] + dvdyy[i])
                       + u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho;

                h[i] = dudx[i] + dvdy[i];
            }
        }
        
        public static void uvp_spiral ( double nu, double rho, int n, double[] x, double[] y, 
        double t, ref double[] u, ref double[] v, ref double[] p )

        //****************************************************************************//
        //
        //  Purpose:
        //
        //    UVP_SPIRAL evaluates Spiral solutions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Maxim Olshanskii, Leo Rebholz,
        //    Application of barycenter refined meshes in linear elasticity
        //    and incompressible fluid dynamics,
        //    ETNA: Electronic Transactions in Numerical Analysis, 
        //    Volume 38, pages 258-274, 2011.
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
        //    double U[N], V[N], the X and Y velocity.
        //
        //    double P[N], the pressure.
        //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                u[i] = ( 1.0 + nu * t ) * 2.0 
                                        * Math.Pow ( x[i], 2 ) * Math.Pow ( x[i] - 1.0, 2 )
                                        * y[i] * ( 2.0 * y[i] - 1.0 ) * ( y[i] - 1.0 );

                v[i] = - ( 1.0 + nu * t ) * 2.0 
                                          * x[i] * ( 2.0 * x[i] - 1.0 ) * ( x[i] - 1.0 )  
                                          * Math.Pow ( y[i], 2 ) * Math.Pow ( y[i] - 1.0, 2 );

                p[i] = rho * y[i];
            }
        }
    }
}