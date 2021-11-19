using System;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class NewtonCotesQuadrature
{
    public static void nc_rule(int norder, double a, double b, double[] xtab, ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NC_RULE computes the weights of a Newton-Cotes quadrature rule.
        //
        //  Discussion:
        //
        //    For the interval [A,B], the Newton-Cotes quadrature rule estimates
        //
        //      Integral ( A <= X <= B ) F(X) dX
        //
        //    using NORDER equally spaced abscissas XTAB(I) and a weight vector
        //    WEIGHT(I):
        //
        //      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
        //
        //    For the CLOSED rule, the abscissas include the points A and B.
        //    For the OPEN rule, the abscissas do not include A and B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NORDER, the order of the rule.
        //
        //    Input, double A, B, the left and right endpoints of the interval
        //    over which the quadrature rule is to be applied.
        //
        //    Input, double XTAB[NORDER], the abscissas of the rule.
        //
        //    Output, double WEIGHT[NORDER], the weights of the rule.
        //
    {
        int i;
        //
        //  Allocate temporary space for POLY_COF.
        //
        double[] poly_cof = new double[norder];

        for (i = 1; i <= norder; i++)
        {
            //
            //  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
            //  and zero at the other nodes.
            //
            typeMethods.r8poly_basis_1(i, norder, xtab, ref poly_cof);
            //
            //  Evaluate the antiderivative of the polynomial at the left and
            //  right endpoints.
            //
            double yvala = typeMethods.r8poly_ant_val(norder - 1, poly_cof, a);

            double yvalb = typeMethods.r8poly_ant_val(norder - 1, poly_cof, b);

            weight[i - 1] = yvalb - yvala;
        }
    }

    public static void ncc_rule(int norder, ref double[] xtab, ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NCC_RULE computes the coefficients of a Newton-Cotes closed quadrature rule.
        //
        //  Discussion:
        //
        //    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
        //
        //      Integral ( -1 <= X <= 1 ) F(X) dX
        //
        //    using NORDER equally spaced abscissas XTAB(I) and a weight vector
        //    WEIGHT(I):
        //
        //      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
        //
        //    For the CLOSED rule, the abscissas include A and B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 September 1998
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NORDER, the order of the rule.
        //
        //    Output, double XTAB[NORDER], the abscissas of the rule.
        //
        //    Output, double WEIGHT[NORDER], the weights of the rule.
        //
    {
        int i;
        //
        //  Compute a closed quadrature rule.
        //
        const double a = -1.0;
        const double b = 1.0;

        for (i = 1; i <= norder; i++)
        {
            xtab[i - 1] = ((norder - i) * a + (i - 1) * b)
                          / (norder - 1);
        }

        nc_rule(norder, a, b, xtab, ref weight);
    }

    public static double nco_abscissa ( int order, int i )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NCO_ABSCISSA returns the I-th abscissa for the Newton Cotes open rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to
        //    right.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //    1 <= ORDER.
        //
        //    Input, int I, the index of the desired abscissa.  
        //    1 <= I <= ORDER.
        //
        //    Output, double NCO_ABSCISSA, the value of the I-th 
        //    abscissa in the Newton Cotes open rule of order ORDER.
        //
    {
        double value;
        const double x_max = +1.0;
        const double x_min = -1.0;

        switch (order)
        {
            case < 1:
                value = - typeMethods.r8_huge ( );
                return value;
        }

        if ( i < 1 || order < i )
        {
            Console.WriteLine("");
            Console.WriteLine("NCO_ABSCISSA - Fatal error!");
            Console.WriteLine("  1 <= I <= ORDER is required.");
            return 1;
        }

        switch (order)
        {
            case 1:
                value = ( x_min + x_max ) / 2.0;
                return value;
            default:
                value = ( (order - i + 1) * x_min
                          + i * x_max )
                        / (order     + 1);

                return value;
        }
    }
        
    public static void nco_rule(int norder, ref double[] xtab, ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NCO_RULE computes the coefficients of a Newton-Cotes open quadrature rule.
        //
        //  Discussion:
        //
        //    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
        //
        //      Integral ( -1 <= X <= 1 ) F(X) dX
        //
        //    using NORDER equally spaced abscissas XTAB(I) and a weight vector
        //    WEIGHT(I):
        //
        //      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
        //
        //    For the OPEN rule, the abscissas do not include A and B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NORDER, the order of the rule.
        //
        //    Output, double XTAB[NORDER], the abscissas of the rule.
        //
        //    Output, double WEIGHT[NORDER], the weights of the  rule.
        //
    {
        int i;

        const double a = -1.0;
        const double b = 1.0;

        for (i = 1; i <= norder; i++)
        {
            xtab[i - 1] = ((norder + 1 - i) * a + i * b)
                          / (norder + 1);
        }

        nc_rule(norder, a, b, xtab, ref weight);
    }

    public static double[] nco_weights(int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NCO_WEIGHTS computes weights for a Newton-Cotes Open rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        int i;
        const double x_max = +1.0;
        const double x_min = -1.0;

        double[] diftab = new double[order];
        double[] w = new double[order];
        double[] x = new double[order];

        for (i = 1; i <= order; i++)
        {
            x[i - 1] = ((order + 1 - i) * x_min
                        + i * x_max)
                       / (order + 1);
        }

        for (i = 1; i <= order; i++)
        {
            //
            //  Compute the Lagrange basis polynomial which is 1 at X(I),
            //  and zero at the other nodes.
            //
            int j;
            for (j = 0; j < order; j++)
            {
                diftab[j] = 0.0;
            }

            diftab[i - 1] = 1.0;

            int k;
            for (j = 2; j <= order; j++)
            {
                for (k = j; k <= order; k++)
                {
                    diftab[order + j - k - 1] = (diftab[order + j - k - 1 - 1] - diftab[order + j - k - 1])
                                                / (x[order + 1 - k - 1] - x[order + j - k - 1]);
                }
            }

            for (j = 1; j < order; j++)
            {
                for (k = 1; k <= order - j; k++)
                {
                    diftab[order - k - 1] -= x[order - k - j] * diftab[order - k];
                }
            }

            //
            //  Evaluate the antiderivative of the polynomial at the left and
            //  right endpoints.
            //
            double yvala = diftab[order - 1] / order;
            for (j = order - 1; 1 <= j; j--)
            {
                yvala = yvala * x_min + diftab[j - 1] / j;
            }

            yvala *= x_min;

            double yvalb = diftab[order - 1] / order;
            for (j = order - 1; 1 <= j; j--)
            {
                yvalb = yvalb * x_max + diftab[j - 1] / j;
            }

            yvalb *= x_max;

            w[i - 1] = yvalb - yvala;
        }

        return w;
    }
}