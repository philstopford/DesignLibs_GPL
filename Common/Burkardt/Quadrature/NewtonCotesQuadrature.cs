using Burkardt.Types;

namespace Burkardt.Quadrature
{
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
            double[] poly_cof;
            double yvala;
            double yvalb;
            //
            //  Allocate temporary space for POLY_COF.
            //
            poly_cof = new double[norder];

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
                yvala = typeMethods.r8poly_ant_val(norder - 1, poly_cof, a);

                yvalb = typeMethods.r8poly_ant_val(norder - 1, poly_cof, b);

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
            double a;
            double b;
            int i;
            //
            //  Compute a closed quadrature rule.
            //
            a = -1.0;
            b = 1.0;

            for (i = 1; i <= norder; i++)
            {
                xtab[i - 1] = ((double) (norder - i) * a + (double) (i - 1) * b)
                              / (double) (norder - 1);
            }

            nc_rule(norder, a, b, xtab, ref weight);
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
            double a;
            double b;
            int i;

            a = -1.0;
            b = 1.0;

            for (i = 1; i <= norder; i++)
            {
                xtab[i - 1] = ((double) (norder + 1 - i) * a + (double) (i) * b)
                              / (double) (norder + 1);
            }

            nc_rule(norder, a, b, xtab, ref weight);
        }
    }
}