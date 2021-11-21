using System;
using Burkardt.Quadrature;

namespace Burkardt.IntegralNS;

public static class QMult
{
    public static double qmult_1d(int setting, Func< int, double, double > func, double a, double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QMULT_1D approximates an integral over an interval in 1D.
        //
        //  Integration region:
        //
        //    A <= X <= B.
        //
        //  Discussion:
        //
        //    A 16 point 31-st degree Gauss-Legendre formula is used.
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
        //    Input, double FUNC ( double x ), the name of the user supplied
        //    function which evaluates F(X).
        //
        //    Input, double A, B, the lower and upper limits of integration.
        //
        //    Output, double QMULT_1D, the approximate integral of 
        //    the function.
        //
    {
        int i;
        const int order = 16;

        double[] xtab = new double[order];
        double[] weight = new double[order];

        LegendreQuadrature.legendre_set(order, ref xtab, ref weight);

        double quad = 0.0;
        for (i = 0; i < order; i++)
        {
            double x = 0.5 * (b - a) * xtab[i] + 0.5 * (a + b);
            quad += 0.5 * weight[i] * func(setting, x);
        }

        double volume = b - a;
        double result = quad * volume;

        return result;
    }

    public static double qmult_2d(Func<double, double, double> func, double a, double b,
            Func< double, double > fup, Func<double, double> flo )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QMULT_2D approximates an integral with varying Y dimension in 2D.
        //
        //  Integration region:
        //
        //      A <= X <= B
        //
        //    and
        //
        //      FLO(X) <= Y <= FHI(X).
        //
        //  Discussion:
        //
        //    A 256 point product of two 16 point 31-st degree Gauss-Legendre
        //    quadrature formulas is used.
        //
        //    This routine could easily be modified to use a different
        //    order product rule by changing the value of ORDER.
        //
        //    Another easy change would allow the X and Y directions to
        //    use quadrature rules of different orders.
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
        //    Input, Func < double, double, double > func, the name of the 
        //    user supplied function which evaluates F(X,Y).
        //
        //    Input, double A, B, the lower and upper limits of X integration.
        //
        //    Input, double FUP ( double x ), double FLO ( double x ), 
        //    the names of the user supplied functions which evaluate the upper 
        //    and lower limits of the Y integration.
        //
    {
        int i;
        const int order = 16;

        double[] xtab = new double[order];
        double[] weight = new double[order];

        LegendreQuadrature.legendre_set(order, ref xtab, ref weight);

        double quad = 0.0;
        for (i = 0; i < order; i++)
        {
            double w1 = 0.5 * (b - a) * weight[i];
            double x = 0.5 * (b - a) * xtab[i] + 0.5 * (b + a);
            double c = flo(x);
            double d = fup(x);

            int j;
            for (j = 0; j < order; j++)
            {
                double w2 = 0.5 * (d - c) * weight[j];
                double y = 0.5 * (d - c) * xtab[j] + 0.5 * (d + c);
                quad += w1 * w2 * func(x, y);
            }
        }

        return quad;
    }

    public static double qmult_3d(int settings, Func<int, double, double, double, double> func, double a,
            double b, Func< double, double > fup1, Func<double, double> flo1,
            Func<double, double, double> fup2, Func<double, double, double> flo2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QMULT_3D approximates an integral with varying Y and Z dimension in 3D.
        //
        //  Integration region:
        //
        //      A         <= X <= B,
        //    and
        //      FLO(X)    <= Y <= FHI(X),
        //    and
        //      FLO2(X,Y) <= Z <= FHI2(X,Y).
        //
        //  Discussion:
        //
        //    A 4096 point product of three 16 point 31-st degree Gauss-Legendre
        //    quadrature formulas is used.
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
        //    user supplied unction which evaluates F(X,Y,Z).
        //
        //    Input, double A, B, the lower and upper limits of X integration.
        //
        //    Input, double FUP1 ( double x ), double FLO1 ( double x ), the names 
        //    of the user supplied functions which evaluate the upper and lower
        //    limits of the Y integration.
        //
        //    Input, double FUP2 ( double x, double y ), 
        //    double FLO2 ( double x, double y ), the names of the user
        //    supplied functions which evaluate the upper and lower
        //    limits of the Z integration.
        //
        //    Output, double QMULT_3D, the approximate integral of 
        //    the function.
        //
    {
        int i;
        const int order = 16;

        double[] xtab = new double[order];
        double[] weight = new double[order];

        LegendreQuadrature.legendre_set(order, ref xtab, ref weight);

        double quad = 0.0;

        for (i = 0; i < order; i++)
        {
            double x = 0.5 * (b - a) * xtab[i] + 0.5 * (b + a);
            double w1 = 0.5 * weight[i];
            double c = flo1(x);
            double d = fup1(x);

            int j;
            for (j = 0; j < order; j++)
            {
                double w2 = 0.5 * (d - c) * weight[j];
                double y = 0.5 * (d - c) * xtab[j] + 0.5 * (d + c);
                double e = flo2(x, y);
                double f = fup2(x, y);

                int k;
                for (k = 0; k < order; k++)
                {
                    double w3 = 0.5 * (f - e) * weight[k];
                    double z = 0.5 * (f - e) * xtab[k] + 0.5 * (f + e);
                    quad += w1 * w2 * w3 * func(settings, x, y, z);
                }
            }
        }

        double volume = b - a;
        double result = quad * volume;

        return result;
    }

}