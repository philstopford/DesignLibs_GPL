namespace Burkardt.PointsNS;

public static partial class Points
{
    public static double[] equidistant1(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIDISTANT1 returns the Type 1 Equidistant points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2018
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double EQUIDISTANT1[N], the points.
        //
    {
        int i;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = (-n + 1 + 2 * i) / (double) (n + 1);
        }

        return x;
    }

    public static double[] equidistant2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIDISTANT2 returns the Type 2 Equidistant points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2018
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double EQUIDISTANT2[N], the points.
        //
    {
        double[] x = new double[n];

        switch (n)
        {
            case 1:
                x[0] = 0.0;
                break;
            default:
            {
                int i;
                for (i = 0; i < n; i++)
                {
                    x[i] = (-n + 1 + 2 * i) / (double) (n - 1);
                }

                break;
            }
        }

        return x;
    }

    public static double[] equidistant3(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIDISTANT3 returns the Type 3 Equidistant points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2018
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double EQUIDISTANT3[N], the points.
        //
    {
        int i;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = (-n + 1 + 2 * i) / (double) n;
        }

        return x;
    }
}