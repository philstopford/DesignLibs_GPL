using System;
using Burkardt.Types;

namespace Burkardt.Polygon;

public static class Integrals
{
    public static double moment(int n, double[] x, double[] y, int p, int q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENT computes an unnormalized moment of a polygon.
        //
        //  Discussion:
        //
        //    Nu(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carsten Steger,
        //    On the calculation of arbitrary moments of polygons,
        //    Technical Report FGBV-96-05,
        //    Forschungsgruppe Bildverstehen, Informatik IX,
        //    Technische Universitaet Muenchen, October 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double X[N], Y[N], the vertex coordinates.
        //
        //    Input, int P, Q, the indices of the moment.
        //
        //    Output, double MOMENT, the unnormalized moment Nu(P,Q).
        //
    {
        int i;

        double nu_pq = 0.0;

        double xj = x[n - 1];
        double yj = y[n - 1];

        for (i = 0; i < n; i++)
        {
            double xi = x[i];
            double yi = y[i];

            double s_pq = 0.0;
            int k;
            for (k = 0; k <= p; k++)
            {
                int l;
                for (l = 0; l <= q; l++)
                {
                    s_pq += typeMethods.r8_choose(k + l, l) * typeMethods.r8_choose(p + q - k - l, q - l)
                                                            * Math.Pow(xi, k) * Math.Pow(xj, p - k)
                                                            * Math.Pow(yi, l) * Math.Pow(yj, q - l);
                }
            }

            nu_pq += (xj * yi - xi * yj) * s_pq;

            xj = xi;
            yj = yi;
        }

        nu_pq = nu_pq / (p + q + 2) / (p + q + 1)
                / typeMethods.r8_choose(p + q, p);

        return nu_pq;
    }

    public static double moment_central(int n, double[] x, double[] y, int p, int q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENT_CENTRAL computes central moments of a polygon.
        //
        //  Discussion:
        //
        //    The central moment Mu(P,Q) is defined by
        //
        //      Mu(P,Q) = Integral ( polygon ) (x-Alpha(1,0))^p (y-Alpha(0,1))^q dx dy
        //              / Area ( polygon )
        //
        //    where 
        //
        //      Alpha(1,0) = Integral ( polygon ) x dx dy / Area ( polygon )
        //      Alpha(0,1) = Integral ( polygon ) y dx dy / Area ( polygon )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carsten Steger,
        //    On the calculation of arbitrary moments of polygons,
        //    Technical Report FGBV-96-05,
        //    Forschungsgruppe Bildverstehen, Informatik IX,
        //    Technische Universitaet Muenchen, October 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double X[N], Y[N], the vertex coordinates.
        //
        //    Input, int P, Q, the indices of the moment.
        //
        //    Output, double MOMENT_CENTRAL, the unnormalized moment Mu(P,Q).
        //
    {
        int i;

        double alpha_10 = moment_normalized(n, x, y, 1, 0);
        double alpha_01 = moment_normalized(n, x, y, 0, 1);

        double mu_pq = 0.0;

        for (i = 0; i <= p; i++)
        {
            int j;
            for (j = 0; j <= q; j++)
            {
                double alpha_ij = moment_normalized(n, x, y, i, j);

                mu_pq += typeMethods.r8_mop(p + q - i - j)
                         * typeMethods.r8_choose(p, i) * typeMethods.r8_choose(q, j)
                         * Math.Pow(alpha_10, p - i) * Math.Pow(alpha_01, q - j) * alpha_ij;
            }
        }

        return mu_pq;
    }

    public static double moment_normalized(int n, double[] x, double[] y, int p, int q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENT_NORMALIZED computes a normalized moment of a polygon.
        //
        //  Discussion:
        //
        //    Alpha(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy / Area ( polygon )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carsten Steger,
        //    On the calculation of arbitrary moments of polygons,
        //    Technical Report FGBV-96-05,
        //    Forschungsgruppe Bildverstehen, Informatik IX,
        //    Technische Universitaet Muenchen, October 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices of the polygon.
        //
        //    Input, double X[N], Y[N], the vertex coordinates.
        //
        //    Input, int P, Q, the indices of the moment.
        //
        //    Output, double MOMENT_NORMALIZED, the normalized moment Alpha(P,Q).
        //
    {
        double nu_pq = moment(n, x, y, p, q);
        double nu_00 = moment(n, x, y, 0, 0);

        double alpha_pq = nu_pq / nu_00;

        return alpha_pq;
    }
}