using System;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class TanhSinh
{
    public static double ts_abscissa(int order, int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TS_ABSCISSA returns the I-th abscissa for the tanh-sinh rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to
        //    right.
        //
        //    This rule is defined on [-1,1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, int I, the index of the desired abscissa. 
        //    1 <= I <= ORDER.
        //
        //    Output, double TS_ABSCISSA, the value of the I-th abscissa 
        //    in the rule of order ORDER.
        //
    {
        double h;
        double st;
        double t;
        double value = 0;

        switch (order)
        {
            case < 1:
                value = -typeMethods.r8_huge();
                break;
            default:
            {
                if (i < 1 || order < i)
                {
                    value = -typeMethods.r8_huge();
                }
                else
                {
                    switch (order)
                    {
                        case 1:
                            value = 0.0;
                            break;
                        default:
                        {
                            switch (2 * i - order - 1)
                            {
                                case 0:
                                    value = 0.0;
                                    break;
                                default:
                                    h = 4.0 / (order + 1);

                                    t = (2 * i - order - 1) * h / 2.0;

                                    st = Math.Sinh(t);
                                    value = Math.Tanh(0.5 * Math.PI * st);
                                    break;
                            }

                            break;
                        }
                    }
                }

                break;
            }
        }

        return value;
    }

    public static double[] ts_weights(int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TS_WEIGHTS computes weights for a tanh-sinh rule.
        //
        //  Discussion:
        //
        //    In the 1D case, a sequence of rules is used of increasing order.
        //    For low order, the weights do not sum to 2, but with increasing 
        //    order, the sum quickly converges to 2.
        //
        //    However, for sparse grid applications, the lowest order rules are
        //    involved in every grid, so it seems it might be useful to force
        //    the weights to sum to 2 immediately.  This addresses only one very
        //    obvious defect of the lower order rules.  I am not sure what to do
        //    about the fact the none of the rules have a definable precision,
        //    and the family of rules has not precision but asymptotic accuracy.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Output, double W[ORDER], the weights of the rule.
        //
    {
        double ct;
        double ct2;
        double h;
        int i;
        double st;
        double t;
        double[] w;
        double w_sum;

        switch (order)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("TS_WEIGHTS - Fatal error!");
                Console.WriteLine("  ORDER < 1.");
                return null;
        }

        w = new double[order];

        h = 4.0 / (order + 1);

        for (i = 0; i < order; i++)
        {
            t = (2 * i - order + 1) * h / 2.0;

            ct = Math.Cosh(t);
            st = Math.Sinh(t);
            ct2 = Math.Cosh(0.5 * Math.PI * st);
            ;

            w[i] = 0.5 * Math.PI * h * ct / ct2 / ct2;
        }

        //
        //  Normalize the weights so that they sum to 2.0.
        //
        w_sum = 0.0;
        for (i = 0; i < order; i++)
        {
            w_sum += w[i];
        }

        for (i = 0; i < order; i++)
        {
            w[i] = 2.0 * w[i] / w_sum;
        }

        return w;
    }
}