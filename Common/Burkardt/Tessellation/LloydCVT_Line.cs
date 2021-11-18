using System;
using Burkardt.PlotNS;
using Burkardt.Types;

namespace Burkardt.Tessellation;

public static class LloydCVT_Line
{
    public static void line_ccvt_lloyd(int n, double a, double b, int it_num, string header,
            ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CCVT_LLOYD carries out the constrained Lloyd algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of generators.
        //
        //    Input, double A, B, the left and right endpoints.
        //
        //    Input, int IT_NUM, the number of iterations to take.
        //
        //    Input, string HEADER, an identifying string.
        //
        //    Input/output, double X[N], the point locations.
        //
    {
        int i;
        int it;

        double[] e_plot = new double[it_num + 1];

        double e = line_cvt_energy(n, a, b, x);
        e_plot[0] = e;

        double[] x_plot = new double[n * (it_num + 1)];
        for (i = 0; i < n; i++)
        {
            x_plot[i + 0 * n] = x[i];
        }

        double[] xm_plot = new double[it_num];

        double[] x_old = new double[n];

        for (it = 1; it <= it_num; it++)
        {
            typeMethods.r8vec_copy(n, x, ref x_old);

            line_ccvt_lloyd_step(n, a, b, ref x);

            for (i = 0; i < n; i++)
            {
                x_plot[i + it * n] = x[i];
            }

            e = line_cvt_energy(n, a, b, x);
            e_plot[it] = e;

            double xm = 0.0;
            for (i = 0; i < n; i++)
            {
                xm += Math.Pow(x_old[i] - x[i], 2);
            }

            xm /= n;
            xm_plot[it - 1] = xm;
        }

        Plot.energy_plot(it_num, e_plot, header);
        Plot.motion_plot(it_num, xm_plot, header);
        Plot.evolution_plot(n, it_num, x_plot, header);
    }

    public static void line_ccvt_lloyd_step(int n, double a, double b, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CCVT_LLOYD_STEP takes one step of Lloyd"s constrained CVT algorithm.
        //
        //  Discussion:
        //
        //    Each step of Lloyd"s algorithm replaces a point by the center of mass
        //    of the associated region.  For points on a line, with a uniform
        //    density, the associated region is demarcated by the midways between 
        //    successive points.
        //
        //    Here, we include the additional constraint that we want the first and last
        //    points to be fixed at the endpoints of the line, that is, X(1) = A
        //    and X(2) = B.  In that case, the calculation of the updates for the
        //    first two and last two points must be handled differently.
        //
        //    For points away from the boundary, a step of Lloyd"s method can be 
        //    regarded as replacing each point by the average of the left and right
        //    midways.  The midways, of course, are the average of two points.
        //    So for point J, we have:
        //
        //      M(J-1,J) = ( X(J-1) + X(J) ) / 2
        //      M(J,J+1) = ( X(J) + X(J+1) ) / 2
        //      X*(J) = ( M(J-1,J) + M(J,J+1) ) / 2 = ( X(J-1) + 2 X(J) + X(J+1) ) / 4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //    1 <= N.
        //
        //    Input, double A, B, the left and right endpoints.
        //
        //    Input/output, double X[N], the point locations.
        //
    {
        double[] x_old = typeMethods.r8vec_copy_new(n, x);

        switch (n)
        {
            case 1:
                x[0] = (a + b) / 2.0;
                break;
            case 2:
                x[0] = a;
                x[1] = b;
                break;
            default:
            {
                x[0] = a;

                int j;
                for (j = 1; j < n - 1; j++)
                {
                    x[j] = (0.5 * (x_old[j - 1] + x_old[j])
                            + 0.5 * (x_old[j] + x_old[j + 1])) / 2.0;
                }

                x[n - 1] = b;
                break;
            }
        }
    }

    public static double line_cvt_energy(int n, double a, double b, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_ENERGY computes the CVT energy for a given set of generators.
        //
        //  Discussion:
        //
        //    Given a set of generators G over the line [A,B], then the energy
        //    is defined as
        //      E = integral ( a <= x <= b ) ( x - g(x) )^2 dx
        //    where g(x) is the nearest generator to the point x.
        //
        //    For the 1D case, this integral can be evaluated exactly as the
        //    sum of integrals over each subinterval:
        //
        //      E(i) = integral ( xl <= x <= xr ) ( x - x(i) )^2 dx
        //           = ( ( x(i) - xl )^3 + ( xr - x(i) )^3 ) / 3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of generators.
        //
        //    Input, double A, B, the left and right endpoints.
        //
        //    Input, double X[N], the generator locations.
        //
        //    Output, double LINE_CVT_ENERGY, the energy of the generator distribution.
        //
    {
        int j;

        double e = 0.0;

        for (j = 0; j < n; j++)
        {
            double xl = j switch
            {
                0 => a,
                _ => (x[j - 1] + x[j]) / 2.0
            };

            double xr;
            if (j == n - 1)
            {
                xr = b;
            }
            else
            {
                xr = (x[j] + x[j + 1]) / 2.0;
            }

            e += (Math.Pow(x[j] - xl, 3) + Math.Pow(xr - x[j], 3)) / 3.0;
        }

        return e;
    }

    public static void line_cvt_lloyd(int n, double a, double b, int it_num, string header,
            double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_LLOYD carries out the Lloyd algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of generators.
        //
        //    Input, double A, B, the left and right endpoints.
        //
        //    Input, int IT_NUM, the number of iterations to take.
        //
        //    Input, string HEADER, an identifying string.
        //
        //    Input/output, double X[N], the point locations.
        //
    {
        int i;
        int it;

        double[] e_plot = new double[it_num + 1];

        double e = line_cvt_energy(n, a, b, x);
        e_plot[0] = e;

        double[] x_plot = new double[n * (it_num + 1)];
        for (i = 0; i < n; i++)
        {
            x_plot[i + 0 * n] = x[i];
        }

        double[] xm_plot = new double[it_num];

        double[] x_old = new double[n];

        for (it = 1; it <= it_num; it++)
        {
            typeMethods.r8vec_copy(n, x, ref x_old);

            line_cvt_lloyd_step(n, a, b, ref x);

            for (i = 0; i < n; i++)
            {
                x_plot[i + it * n] = x[i];
            }

            e = line_cvt_energy(n, a, b, x);
            e_plot[it] = e;

            double xm = 0.0;
            for (i = 0; i < n; i++)
            {
                xm += Math.Pow(x_old[i] - x[i], 2);
            }

            xm /= n;
            xm_plot[it - 1] = xm;
        }

        Plot.energy_plot(it_num, e_plot, header);
        Plot.motion_plot(it_num, xm_plot, header);
        Plot.evolution_plot(n, it_num, x_plot, header);
    }

    public static void line_cvt_lloyd_step(int n, double a, double b, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_LLOYD_STEP takes one step of Lloyd"s unconstrained CVT algorithm.
        //
        //  Discussion:
        //
        //    Each step of Lloyd"s algorithm replaces a point by the center of mass
        //    of the associated region.  For points on a line, with a uniform
        //    density, the associated region is demarcated by the midways between 
        //    successive points.
        //
        //    For points away from the boundary, a step of Lloyd"s method can be 
        //    regarded as replacing each point by the average of the left and right
        //    midways.  The midways, of course, are the average of two points.
        //    So for point J, we have:
        //
        //      M(J-1,J) = ( X(J-1) + X(J) ) / 2
        //      M(J,J+1) = ( X(J) + X(J+1) ) / 2
        //      X*(J) = ( M(J-1,J) + M(J,J+1) ) / 2 = ( X(J-1) + 2 X(J) + X(J+1) ) / 4
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //    1 <= N.
        //
        //    Input, double A, B, the left and right endpoints.
        //
        //    Input, double X[N], the point locations.
        //
    {
        double[] x_old = typeMethods.r8vec_copy_new(n, x);

        switch (n)
        {
            case 1:
                x[0] = (a + b) / 2.0;
                break;
            default:
            {
                int j = 0;
                x[j] = (a
                        + 0.5 * (x_old[j] + x_old[1])) / 2.0;

                for (j = 1; j < n - 1; j++)
                {
                    x[j] = (0.5 * (x_old[j - 1] + x_old[j])
                            + 0.5 * (x_old[j] + x_old[j + 1])) / 2.0;
                }

                j = n - 1;
                x[j] = (0.5 * (x_old[j - 1] + x_old[j])
                        + b) / 2.0;
                break;
            }
        }
    }
}