namespace Burkardt.NavierStokesNS;

public static class Poiseuille
{
    public static void resid_poiseuille ( double nu, double rho, int n, double[] x, double[] y, 
            double t, ref double[] ur, ref double[] vr, ref double[] pr )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESID_POISEUILLE returns Poiseuille residuals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
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

        rhs_poiseuille ( nu, rho, n, x, y, t, ref f, ref g, ref h );
        //
        //  Form the functions and derivatives of the left hand side.
        //
        for ( i = 0; i < n; i++ )
        {
            double u = 1.0 - y[i] * y[i];
            const double dudt = 0.0;
            const double dudx = 0.0;
            const double dudxx = 0.0;
            double dudy = - 2.0 * y[i];
            const double dudyy = - 2.0;

            const double v = 0.0;
            const double dvdt = 0.0;
            const double dvdx = 0.0;
            const double dvdxx = 0.0;
            const double dvdy = 0.0;
            const double dvdyy = 0.0;

            double dpdx = - 2.0 * rho * nu;
            const double dpdy = 0.0;
            //
            //  Evaluate the residuals.
            //
            ur[i] = dudt + u * dudx + v * dudy 
                + 1.0 / rho * dpdx - nu * ( dudxx + dudyy ) - f[i];

            vr[i] = dvdt + u * dvdx + v * dvdy 
                + 1.0 / rho * dpdy - nu * ( dvdxx + dvdyy ) - g[i];

            pr[i] = dudx + dvdy - h[i];
        }
    }
        
    public static void rhs_poiseuille ( double nu, double rho, int n, double[] x, double[] y, 
            double t, ref double[] f, ref double[] g, ref double[] h )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RHS_POISEUILLE evaluates Poiseuille right hand sides.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
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

        for ( i = 0; i < n; i++ )
        {
            f[i] = 0.0;
            g[i] = 0.0;
            h[i] = 0.0;
        }
    }
        
    public static void uvp_poiseuille ( double nu, double rho, int n, double[] x, double[] y, 
            double t, ref double[] u, ref double[] v, ref double[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UVP_POISEUILLE evaluate Poiseuille solutions.
        //
        //  Discussion:
        //
        //    There is no time dependence.
        //
        //    The vertical velocity is zero.
        //
        //    The preferred domain is a channel bounded by y = -1 and y = +1,
        //    along which the boundary condition u = 0 and v = 0 will be satisfied.
        //    A parabolic inflow may then be imposed along some line, such as x=0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
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
            u[i] = 1.0 - y[i] * y[i];
            v[i] = 0.0;
            p[i] = - 2.0 * rho * nu;
        }
    }
}