namespace Burkardt.TetrahedronNS;

public static class Grid
{
    public static double[] tetrahedron_grid(int n, double[] t, int ng)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_GRID computes points on a tetrahedral grid.
        //
        //  Discussion:
        //
        //    The grid is defined by specifying the coordinates of an enclosing
        //    tetrahedron T, and the number of subintervals each edge of the 
        //    tetrahedron should be divided into.
        //
        //    Choosing N = 10, for instance, breaks each side into 10 subintervals,
        //    and produces a grid of ((10+1)*(10+2)*(10+3))/6 = 286 points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, double T[3*4], the vertices of the tetrahedron.
        //
        //    Input, int NG, the number of grid points.
        //
        //    Output, double TETRAHEDRON_GIRD[3*NG], the tetrahedron grid points.
        //
    {
        int i;
        int ii;
        int j;
        int k;
        int l;
        int p;
        double[] tg;

        tg = new double[3 * ng];

        p = 0;

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n - i; j++)
            {
                for (k = 0; k <= n - i - j; k++)
                {
                    l = n - i - j - k;
                    for (ii = 0; ii < 3; ii++)
                    {
                        tg[ii + p * 3] = (i * t[ii + 0 * 3]
                                          + j * t[ii + 1 * 3]
                                          + k * t[ii + 2 * 3]
                                          + l * t[ii + 3 * 3]) / n;
                    }

                    p += 1;
                }
            }
        }

        return tg;
    }

    public static int tetrahedron_grid_count(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_GRID_COUNT counts the grid points inside a tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Output, int TETRAHEDRON_GRID_COUNT, the number of grid points.
        //
    {
        int ng;

        ng = (n + 1) * (n + 2) * (n + 3) / 6;

        return ng;
    }
}