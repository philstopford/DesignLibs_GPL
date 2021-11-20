using System;
using Burkardt.Quadrature;

namespace Burkardt.NavierStokesNS;

public static class Burgers
{
    public static double[] burgers_viscous_time_exact1(double nu, int vxn, double[] vx, int vtn,
            double[] vt )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURGERS_VISCOUS_TIME_EXACT1 evaluates solution #1 to the Burgers equation.
        //
        //  Discussion:
        //
        //    The form of the Burgers equation considered here is
        //
        //      du       du        d^2 u
        //      -- + u * -- = nu * -----
        //      dt       dx        dx^2
        //
        //    for -1.0 < x < +1.0, and 0 < t.
        //
        //    Initial conditions are u(x,0) = - sin(pi*x).  Boundary conditions
        //    are u(-1,t) = u(+1,t) = 0.  The viscosity parameter nu is taken
        //    to be 0.01 / pi, although this is not essential.
        //
        //    The authors note an integral representation for the solution u(x,t),
        //    and present a better version of the formula that is amenable to
        //    approximation using Hermite quadrature.  
        //
        //    This program library does little more than evaluate the exact solution
        //    at a user-specified set of points, using the quadrature rule.
        //    Internally, the order of this quadrature rule is set to 8, but the
        //    user can easily modify this value if greater accuracy is desired.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 November 2011
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Claude Basdevant, Michel Deville, Pierre Haldenwang, J Lacroix, 
        //    J Ouazzani, Roger Peyret, Paolo Orlandi, Anthony Patera,
        //    Spectral and finite difference solutions of the Burgers equation,
        //    Computers and Fluids,
        //    Volume 14, Number 1, 1986, pages 23-41.
        //
        //  Parameters:
        //
        //    Input, double NU, the viscosity.
        //
        //    Input, int VXN, the number of spatial grid points.
        //
        //    Input, double VX[VXN], the spatial grid points.
        //
        //    Input, int VTN, the number of time grid points.
        //
        //    Input, double VT[VTN], the time grid points.
        //
        //    Output, double BURGERS_VISCOUS_TIME_EXACT1[VXN*VTN], the solution of 
        //    the Burgers equation at each space and time grid point.
        //
    {
        const int qn = 8;
        int vti;
        //
        //  Compute the rule.
        //
        double[] qx = new double[qn];
        double[] qw = new double[qn];

        HermiteQuadrature.hermite_ek_compute(qn, ref qx, ref qw);
        //
        //  Evaluate U(X,T) for later times.
        //
        double[] vu = new double[vxn * vtn];

        for (vti = 0; vti < vtn; vti++)
        {
            int vxi;
            switch (vt[vti])
            {
                case 0.0:
                {
                    for (vxi = 0; vxi < vxn; vxi++)
                    {
                        vu[vxi + vti * vxn] = -Math.Sin(Math.PI * vx[vxi]);
                    }

                    break;
                }
                default:
                {
                    for (vxi = 0; vxi < vxn; vxi++)
                    {
                        double top = 0.0;
                        double bot = 0.0;
                        int qi;
                        for (qi = 0; qi < qn; qi++)
                        {
                            double c = 2.0 * Math.Sqrt(nu * vt[vti]);

                            top -= qw[qi] * c * Math.Sin(Math.PI * (vx[vxi] - c * qx[qi]))
                                   * Math.Exp(-Math.Cos(Math.PI * (vx[vxi] - c * qx[qi]))
                                              / (2.0 * Math.PI * nu));

                            bot += qw[qi] * c
                                          * Math.Exp(-Math.Cos(Math.PI * (vx[vxi] - c * qx[qi]))
                                                     / (2.0 * Math.PI * nu));

                            vu[vxi + vti * vxn] = top / bot;
                        }
                    }

                    break;
                }
            }
        }

        return vu;
    }

    public static double[] burgers_viscous_time_exact2(double nu, int xn, double[] x, int tn,
            double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURGERS_VISCOUS_TIME_EXACT2 evaluates solution #2 to the Burgers equation.
        //
        //  Discussion:
        //
        //    The form of the Burgers equation considered here is
        //
        //      du       du        d^2 u
        //      -- + u * -- = nu * -----
        //      dt       dx        dx^2
        //
        //   for 0.0 < x < 2 Pi and 0 < t.
        //
        //    The initial condition is
        //
        //      u(x,0) = 4 - 2 * nu * dphi(x,0)/dx / phi(x,0)
        //
        //    where
        //
        //      phi(x,t) = exp ( - ( x-4*t      ) / ( 4*nu*(t+1) ) )
        //               + exp ( - ( x-4*t-2*pi ) / ( 4*nu*(t+1) ) )
        //
        //    The boundary conditions are periodic:
        //
        //      u(0,t) = u(2 Pi,t)
        //
        //    The viscosity parameter nu may be taken to be 0.01, but other values
        //    may be chosen.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2015
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Claude Basdevant, Michel Deville, Pierre Haldenwang, J Lacroix, 
        //    J Ouazzani, Roger Peyret, Paolo Orlandi, Anthony Patera,
        //    Spectral and finite difference solutions of the Burgers equation,
        //    Computers and Fluids,
        //    Volume 14, Number 1, 1986, pages 23-41.
        //
        //  Parameters:
        //
        //    Input, double NU, the viscosity.
        //
        //    Input, int XN, the number of spatial grid points.
        //
        //    Input, double X[XN], the spatial grid points.
        //
        //    Input, int TN, the number of time grid points.
        //
        //    Input, double T[TN], the time grid points.
        //
        //    Output, double BURGERS_VISCOUS_TIME_EXACT2[XN*TN], the solution of the 
        //    Burgers equation at each space and time grid point.
        //
    {
        int j;

        double[] u = new double[xn * tn];

        for (j = 0; j < tn; j++)
        {
            int i;
            for (i = 0; i < xn; i++)
            {
                double a = x[i] - 4.0 * t[j];
                double b = x[i] - 4.0 * t[j] - 2.0 * Math.PI;
                double c = 4.0 * nu * (t[j] + 1.0);
                double phi = Math.Exp(-a * a / c) + Math.Exp(-b * b / c);
                double dphi = -2.0 * a * Math.Exp(-a * a / c) / c
                              - 2.0 * b * Math.Exp(-b * b / c) / c;
                u[i + j * xn] = 4.0 - 2.0 * nu * dphi / phi;
            }
        }

        return u;
    }

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

        for (i = 0; i < n; i++)
        {
            //
            //  Form the functions and derivatives.
            //
            double u = 2.0 * x[i];
            const double ux = 2.0;
            const double uxx = 0.0;
            const double uy = 0.0;
            const double uyy = 0.0;
            const double uz = 0.0;
            const double uzz = 0.0;
            const double ut = 0.0;

            double v = -2.0 * y[i];
            const double vx = 0.0;
            const double vxx = 0.0;
            const double vy = -2.0;
            const double vyy = 0.0;
            const double vz = 0.0;
            const double vzz = 0.0;
            const double vt = 0.0;

            double w = Helpers.Erf(y[i] / Math.Sqrt(nu));
            const double wx = 0.0;
            const double wxx = 0.0;
            double wy = 2.0 * Math.Sqrt(1.0 / nu / Math.PI) * Math.Exp(-y[i] * y[i] / nu);
            double wyy = -4.0 * Math.Sqrt(1.0 / nu / Math.PI) * y[i] * Math.Exp(-y[i] * y[i] / nu) / nu;
            const double wz = 0.0;
            const double wzz = 0.0;
            const double wt = 0.0;

            //  p = - 2.0 * ( x[i] * x[i] + y[i] * y[i] );
            double px = -4.0 * x[i];
            double py = -4.0 * y[i];
            const double pz = 0.0;
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
    }
}