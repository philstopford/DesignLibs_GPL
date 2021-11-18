using System;

namespace Burkardt.Stroud;

public static class Octahedron
{
    public static double octahedron_unit_nd(int setting, Func<int, int, double[], double> func, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OCTAHEDRON_UNIT_ND approximates integrals in the unit octahedron in ND.
        //
        //  Integration region:
        //
        //    sum ( abs ( X(1:N) ) ) <= 1.
        //
        //  Discussion:
        //
        //    A 2*N point 3rd degree formula is used, Stroud number GN:3-1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, Func < int, double[], double> func, the name of the
        //    user supplied function to be integrated.
        //
        //    Input, int N, the dimension of the octahedron.
        //
        //    Output, double OCTAHEDRON_UNIT_ND, the approximate integral of the function.
        //
    {
        int i;

        double[] x = new double[n];

        double w = 1.0 / (2 * n);

        double r = Math.Sqrt(2 * n
                             / (double)((n + 1) * (n + 2)));

        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        double quad = 0.0;
        for (i = 0; i < n; i++)
        {
            x[i] = r;
            int j;
            for (j = 0; j < 2; j++)
            {
                quad += w * func(setting, n, x);
                x[i] = -x[i];
            }

            x[i] = 0.0;
        }

        double volume = octahedron_unit_volume_nd(n);
        double result = quad * volume;

        return result;
    }

    public static double octahedron_unit_volume_nd(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OCTAHEDRON_UNIT_VOLUME_ND returns the volume of the unit octahedron in ND.
        //
        //  Integration region:
        //
        //    sum ( abs ( X(1:N) ) ) <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the space.
        //
        //    Output, double OCTAHEDRON_UNIT_VOLUME_ND, the volume of
        //    the unit octahedron.
        //
    {
        int i;

        double volume = 1.0;
        for (i = 1; i <= n; i++)
        {
            volume = volume * 2.0 / i;
        }

        return volume;
    }

}