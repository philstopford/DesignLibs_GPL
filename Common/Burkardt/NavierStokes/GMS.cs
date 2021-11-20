using System;

namespace Burkardt.NavierStokesNS;

public static class GMS
{
    public static void all_gms(double nu, double rho, int n, double[] x, double[] y, double t,
            ref double[] u, ref double[] dudt, ref double[] dudx, ref double[] dudxx, ref double[] dudy,
            ref double[] dudyy, ref double[] us,
            ref double[] v, ref double[] dvdt, ref double[] dvdx, ref double[] dvdxx, ref double[] dvdy,
            ref double[] dvdyy, ref double[] vs,
            ref double[] p, ref double[] dpdt, ref double[] dpdx, ref double[] dpdxx, ref double[] dpdy,
            ref double[] dpdyy, ref double[] ps)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    all_gms evaluates the variables of the GMS flow.
        //
        //  Discussion:
        //
        //    The flow has been modified by a sign change that makes it slightly
        //    more plausible physically.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 August 2020
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
        //    double U[N], DUDT[N], DUDX[N], DUDXX[N], DUDY[N], DUDYY[N], US[N],
        //    the horizontal velocity values.
        //
        //    double V[N], DVDT[N], DVDX[N], DVDXX[N], DVDY[N], DVDYY[N], VS[N],
        //    the vertical velocity values.
        //
        //    double P[N], DPDT[N], DPDX[N], DPDXX[N], DPDY[N], DPDYY[N], PS[N],
        //    the pressure values.
        //
    {
        int i;
        const double pi2 = Math.PI * Math.PI;
        const double pi3 = Math.PI * Math.PI * Math.PI;

        for (i = 0; i < n; i++)
        {
            double s = Math.Pow(-1.0, Math.Ceiling(x[i]) + Math.Ceiling(y[i]));

            u[i] = Math.PI * s * Math.Sin(t) * Math.Pow(Math.Sin(Math.PI * x[i]), 2) * Math.Sin(2.0 * Math.PI * y[i]);
            dudt[i] = Math.PI * s * Math.Cos(t) * Math.Pow(Math.Sin(Math.PI * x[i]), 2) * Math.Sin(2.0 * Math.PI * y[i]);
            dudx[i] = pi2 * s * Math.Sin(t) * Math.Sin(2.0 * Math.PI * x[i]) * Math.Sin(2.0 * Math.PI * y[i]);
            dudxx[i] = 2.0 * pi3 * s * Math.Sin(t) * Math.Cos(2.0 * Math.PI * x[i]) * Math.Sin(2.0 * Math.PI * y[i]);
            dudy[i] = 2.0 * pi2 * s * Math.Sin(t) * Math.Pow(Math.Sin(Math.PI * x[i]), 2) * Math.Cos(2.0 * Math.PI * y[i]);
            dudyy[i] = -4.0 * pi3 * s * Math.Sin(t) * Math.Pow(Math.Sin(Math.PI * x[i]), 2) * Math.Sin(2.0 * Math.PI * y[i]);

            v[i] = -Math.PI * s * Math.Sin(t) * Math.Sin(2.0 * Math.PI * x[i]) * Math.Pow(Math.Sin(Math.PI * y[i]), 2);
            dvdt[i] = -Math.PI * s * Math.Cos(t) * Math.Sin(2.0 * Math.PI * x[i]) * Math.Pow(Math.Sin(Math.PI * y[i]), 2);
            dvdx[i] = -2.0 * pi2 * s * Math.Sin(t) * Math.Cos(2.0 * Math.PI * x[i]) * Math.Pow(Math.Sin(Math.PI * y[i]), 2);
            dvdxx[i] = 4.0 * pi3 * s * Math.Sin(t) * Math.Sin(2.0 * Math.PI * x[i]) * Math.Pow(Math.Sin(Math.PI * y[i]), 2);
            dvdy[i] = -pi2 * s * Math.Sin(t) * Math.Sin(2.0 * Math.PI * x[i]) * Math.Sin(2.0 * Math.PI * y[i]);
            dvdyy[i] = -2.0 * pi3 * s * Math.Sin(t) * Math.Sin(2.0 * Math.PI * x[i]) * Math.Cos(2.0 * Math.PI * y[i]);

            p[i] = rho * s * Math.Sin(t) * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dpdt[i] = rho * s * Math.Cos(t) * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dpdx[i] = -Math.PI * rho * s * Math.Sin(t) * Math.Sin(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dpdxx[i] = -pi2 * rho * s * Math.Sin(t) * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);
            dpdy[i] = Math.PI * rho * s * Math.Sin(t) * Math.Cos(Math.PI * x[i]) * Math.Cos(Math.PI * y[i]);
            dpdyy[i] = -pi2 * rho * s * Math.Sin(t) * Math.Cos(Math.PI * x[i]) * Math.Sin(Math.PI * y[i]);

            us[i] = dudt[i] + u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho - nu * (dudxx[i] + dudyy[i]);
            vs[i] = dvdt[i] + u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho - nu * (dvdxx[i] + dvdyy[i]);
            ps[i] = dudx[i] + dvdy[i];
        }
    }

    public static void resid_gms(double nu, double rho, int n, double[] x, double[] y, double t,
            ref double[] ur, ref double[] vr, ref double[] pr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    resid_gms evaluates the residual of the GMS flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2020
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
        //    double UR[N], VR[N], PR[N], the residuals.
        //
    {
        int i;

        double[] u = new double[n];
        double[] dudt = new double[n];
        double[] dudx = new double[n];
        double[] dudxx = new double[n];
        double[] dudy = new double[n];
        double[] dudyy = new double[n];
        double[] us = new double[n];
        double[] v = new double[n];
        double[] dvdt = new double[n];
        double[] dvdx = new double[n];
        double[] dvdxx = new double[n];
        double[] dvdy = new double[n];
        double[] dvdyy = new double[n];
        double[] vs = new double[n];
        double[] p = new double[n];
        double[] dpdt = new double[n];
        double[] dpdx = new double[n];
        double[] dpdxx = new double[n];
        double[] dpdy = new double[n];
        double[] dpdyy = new double[n];
        double[] ps = new double[n];

        all_gms(nu, rho, n, x, y, t,
            ref u, ref dudt, ref dudx, ref dudxx, ref dudy, ref dudyy, ref us,
            ref v, ref dvdt, ref dvdx, ref dvdxx, ref dvdy, ref dvdyy, ref vs,
            ref p, ref dpdt, ref dpdx, ref dpdxx, ref dpdy, ref dpdyy, ref ps);

        for (i = 0; i < n; i++)
        {
            ur[i] = dudt[i] + u[i] * dudx[i] + v[i] * dudy[i]
                + dpdx[i] / rho - nu * (dudxx[i] + dudyy[i]) - us[i];
            vr[i] = dvdt[i] + u[i] * dvdx[i] + v[i] * dvdy[i]
                + dpdy[i] / rho - nu * (dvdxx[i] + dvdyy[i]) - vs[i];
            pr[i] = dudx[i] + dvdy[i] - ps[i];
        }
    }
        
    public static void rhs_gms ( double nu, double rho, int n, double[] x, double[] y, double t, 
            ref double[] us, ref double[] vs, ref double[] ps )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    rhs_gms evaluates the source terms of the GMS flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2020
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
        //    double US[N], VS[N], PS[N], the source terms.
        //
    {
        double[] u = new double[n];
        double[] dudt = new double[n];
        double[] dudx = new double[n];
        double[] dudxx = new double[n];
        double[] dudy = new double[n];
        double[] dudyy = new double[n];
        double[] v = new double[n];
        double[] dvdt = new double[n];
        double[] dvdx = new double[n];
        double[] dvdxx = new double[n];
        double[] dvdy = new double[n];
        double[] dvdyy = new double[n];
        double[] p = new double[n];
        double[] dpdt = new double[n];
        double[] dpdx = new double[n];
        double[] dpdxx = new double[n];
        double[] dpdy = new double[n];
        double[] dpdyy = new double[n];

        all_gms ( nu, rho, n, x, y, t,
            ref u, ref dudt, ref dudx, ref dudxx, ref dudy, ref dudyy, ref us, 
            ref v, ref dvdt, ref dvdx, ref dvdxx, ref dvdy, ref dvdyy, ref vs, 
            ref p, ref dpdt, ref dpdx, ref dpdxx, ref dpdy, ref dpdyy, ref ps );

    }
        
    public static void uvp_gms ( double nu, double rho, int n, double[] x, double[] y, double t, 
            ref double[] u, ref double[] v, ref double[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    uvp_gms evaluates the state variables of the GMS flow.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2020
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
        //    double U[N], V[N], P[N], the state variables.
        //
    {
        double[] dudt = new double[n];
        double[] dudx = new double[n];
        double[] dudxx = new double[n];
        double[] dudy = new double[n];
        double[] dudyy = new double[n];
        double[] us = new double[n];
        double[] dvdt = new double[n];
        double[] dvdx = new double[n];
        double[] dvdxx = new double[n];
        double[] dvdy = new double[n];
        double[] dvdyy = new double[n];
        double[] vs = new double[n];
        double[] dpdt = new double[n];
        double[] dpdx = new double[n];
        double[] dpdxx = new double[n];
        double[] dpdy = new double[n];
        double[] dpdyy = new double[n];
        double[] ps = new double[n];

        all_gms ( nu, rho, n, x, y, t,
            ref u, ref dudt, ref dudx, ref dudxx, ref dudy, ref dudyy, ref us, 
            ref v, ref dvdt, ref dvdx, ref dvdxx, ref dvdy, ref dvdyy, ref vs, 
            ref p, ref dpdt, ref dpdx, ref dpdxx, ref dpdy, ref dpdyy, ref ps );
    }
}