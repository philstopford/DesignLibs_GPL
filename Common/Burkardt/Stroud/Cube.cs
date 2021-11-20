using System;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.Stroud;

public static class Cube
{
    public static double cube_shell_nd(int setting, Func<int, int, double[], double> func, int n, double r1,
            double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_SHELL_ND approximates an integral inside a cubic shell in N dimensions.
        //
        //  Integration region:
        //
        //    R1 <= abs ( X(1:N) ) <= R2
        //
        //  Discussion:
        //
        //    An N*2^N point third degree formula is used, Stroud number CNSHELL:3-4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 March 2008
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
        //  user supplied function.
        //
        //    Input, int N, the dimension of the space.
        //
        //    Input, double R1, R2, the inner and outer radii of the cubical
        //    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
        //    2*R1.
        //
        //    Output, double CUBE_SHELL_ND, the approximate integral of the function.
        //
    {
        int i;
        double result;

        if (Math.Abs(r1 - r2) <= double.Epsilon)
        {
            result = 0.0;
            return result;
        }

        double rmax = Math.Max(r1, r2);
        double rmin = Math.Min(r1, r2);

        double u = Math.Sqrt(n * (Math.Pow(rmax, n + 2) - Math.Pow(rmin, n + 2))
                             / ((n + 2) * (Math.Pow(rmax, n) - Math.Pow(rmin, n))));

        double v = u / Math.Sqrt(3.0);

        double[] x = new double[n];

        double quad = 0.0;
        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                x[j] = v;
            }

            x[i] = u;

            for (;;)
            {
                quad += func(setting, n, x);

                bool done = typeMethods.r8vec_mirror_next(n, ref x);

                if (done)
                {
                    break;
                }
            }
        }

        quad /= n * (int)Math.Pow(2, n);

        double volume = cube_shell_volume_nd(n, r1, r2);
        result = quad * volume;

        return result;
    }

    public static double cube_shell_volume_nd(int n, double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_SHELL_VOLUME_ND computes the volume of a cubic shell in ND.
        //
        //  Integration region:
        //
        //    R1 <= abs ( X(1:N) ) <= R2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the space.
        //
        //    Input, double R1, R2, the inner and outer radii of the cubic
        //    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
        //    2*R1.
        //
        //    Output, double CUBE_SHELL_VOLUME_ND, the volume of the cubic
        //    shell.
        //
    {
        double value = (Math.Pow(r2, n) - Math.Pow(r1, n)) * (int)Math.Pow(2, n);

        return value;
    }

    public static double cube_unit_3d(int settings, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_UNIT_3D approximates an integral inside the unit cube in 3D.
        //
        //  Integration region:
        //
        //      -1 <= X <= 1,
        //    and
        //      -1 <= Y <= 1,
        //    and
        //      -1 <= Z <= 1.
        //
        //  Discussion:
        //
        //    An 8 point third degree formula is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the
        //    user supplied routine to evaluate F(X,Y,Z).
        //
        //    Output, double CUBE_UNIT_3D, the approximate integral of the function.
        //
    {
        double s = 1.0 / Math.Sqrt(3.0);
        const double w = 1.0 / 8.0;

        double quad = w * (
            func(settings, s, s, s) + func(settings, s, s, -s)
                                    + func(settings, s, -s, s) + func(settings, s, -s, -s)
                                    + func(settings, -s, s, s) + func(settings, -s, s, -s)
                                    + func(settings, -s, -s, s) + func(settings, -s, -s, -s));

        double volume = cube_unit_volume_nd(3);
        double result = quad * volume;

        return result;
    }

    public static void cube_unit_nd(int setting, Func<int, int, double[], double> func, ref double[] qa,
            ref double[] qb, int n, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_UNIT_ND approximates an integral inside the unit cube in ND.
        //
        //  Integration region:
        //
        //    -1 <= X(1:N) <= 1
        //
        //  Discussion:
        //
        //    A K**N point product formula is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, BJJ McHugh,
        //    Integration Over Multidimensional Hypercubes, 
        //    A Progressive Procedure,
        //    The Computer Journal,
        //    Volume 6, 1963, pages 264-270.
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
        //    Output, double QA[K], QB[K], two sets of estimates for
        //    the integral.  The QB entries are obtained from the
        //    QA entries by Richardson extrapolation, and QB(K) is
        //    the best estimate for the integral.
        //
        //    Input, int N, the dimension of the cube.
        //
        //    Input, int K, the highest order of integration, and the order
        //    of Richardson extrapolation.  K can be no greater than 10.
        //
    {
        double[] g = new double[10 * 10];
        int i;
        const int kmax = 10;

        g[0 + 0 * 10] = 1.0E+00;
        g[1 + 0 * 10] = -0.3333333333333E+00;
        g[1 + 1 * 10] = 0.1333333333333E+01;
        g[2 + 0 * 10] = 0.4166666666667E-01;
        g[2 + 1 * 10] = -0.1066666666667E+01;
        g[2 + 2 * 10] = 0.2025000000000E+01;
        g[3 + 0 * 10] = -0.2777777777778E-02;
        g[3 + 1 * 10] = 0.3555555555556E+00;
        g[3 + 2 * 10] = -0.2603571428571E+01;
        g[3 + 3 * 10] = 0.3250793650794E+01;
        g[4 + 0 * 10] = 0.1157407407407E-03;
        g[4 + 1 * 10] = -0.6772486772487E-01;
        g[4 + 2 * 10] = 0.1464508928571E+01;
        g[4 + 3 * 10] = -0.5779188712522E+01;
        g[4 + 4 * 10] = 0.5382288910935E+01;
        g[5 + 0 * 10] = -0.3306878306878E-05;
        g[5 + 1 * 10] = 0.8465608465608E-02;
        g[5 + 2 * 10] = -0.4881696428571E+00;
        g[5 + 3 * 10] = 0.4623350970018E+01;
        g[5 + 4 * 10] = -0.1223247479758E+02;
        g[5 + 5 * 10] = 0.9088831168831E+01;
        g[6 + 0 * 10] = 0.6889329805996E-07;
        g[6 + 1 * 10] = -0.7524985302763E-03;
        g[6 + 2 * 10] = 0.1098381696429E+00;
        g[6 + 3 * 10] = -0.2241624712736E+01;
        g[6 + 4 * 10] = 0.1274216124748E+02;
        g[6 + 5 * 10] = -0.2516907092907E+02;
        g[6 + 6 * 10] = 0.1555944865432E+02;
        g[7 + 0 * 10] = -0.1093544413650E-08;
        g[7 + 1 * 10] = 0.5016656868509E-04;
        g[7 + 2 * 10] = -0.1797351866883E-01;
        g[7 + 3 * 10] = 0.7472082375786E+00;
        g[7 + 4 * 10] = -0.8168052081717E+01;
        g[7 + 5 * 10] = 0.3236023405166E+02;
        g[7 + 6 * 10] = -0.5082753227079E+02;
        g[7 + 7 * 10] = 0.2690606541646E+02;
        g[8 + 0 * 10] = 0.1366930517063E-10;
        g[8 + 1 * 10] = -0.2606055516108E-05;
        g[8 + 2 * 10] = 0.2246689833604E-02;
        g[8 + 3 * 10] = -0.1839281815578E+00;
        g[8 + 4 * 10] = 0.3646451822195E+01;
        g[8 + 5 * 10] = -0.2588818724133E+02;
        g[8 + 6 * 10] = 0.7782965878964E+02;
        g[8 + 7 * 10] = -0.1012934227443E+03;
        g[8 + 8 * 10] = 0.4688718347156E+02;
        g[9 + 0 * 10] = -0.1380737896023E-12;
        g[9 + 1 * 10] = 0.1085856465045E-06;
        g[9 + 2 * 10] = -0.2222000934334E-03;
        g[9 + 3 * 10] = 0.3503393934435E-01;
        g[9 + 4 * 10] = -0.1215483940732E+01;
        g[9 + 5 * 10] = 0.1456210532325E+02;
        g[9 + 6 * 10] = -0.7477751530769E+02;
        g[9 + 7 * 10] = 0.1800771959898E+03;
        g[9 + 8 * 10] = -0.1998874663788E+03;
        g[9 + 9 * 10] = 0.8220635246624E+02;

        if (kmax < k)
        {
            Console.WriteLine("");
            Console.WriteLine("CUBE_UNIT_ND - Fatal error!");
            Console.WriteLine("  K must be no greater than KMAX = " + kmax + "");
            Console.WriteLine("  but the input K is " + k + "");
            return;
        }

        for (i = 0; i < k; i++)
        {
            qa[i] = QuadratureRule.qmdpt(setting, func, n, i + 1);
        }

        qb[0] = qa[0];

        for (i = 1; i < k; i++)
        {
            qb[i] = 0.0;
            int j;
            for (j = 0; j <= i; j++)
            {
                qb[i] += g[i + j * 10] * qa[j];
            }
        }
    }

    public static double cube_unit_volume_nd(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_UNIT_VOLUME_ND returns the volume of the unit cube in ND.
        //
        //  Integration region:
        //
        //    -1 <= X(1:N) <= 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the space.
        //
        //    Output, double CUBE_UNIT_VOLUME_ND, the volume of the unit
        //    cube in ND.
        //
    {
        double value = Math.Pow(2.0, n);

        return value;
    }


}