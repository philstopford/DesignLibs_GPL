using System;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class Fejer1
{
    public static double f1_abscissa(int order, int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1_ABSCISSA returns the I-th abscissa for the Fejer type 1 rule.
        //
        //  Discussion:
        //
        //    Our convention is that the abscissas are numbered from left to
        //    right.
        //
        //    This rule is defined on [-1,1].
        //
        //  Modified:
        //
        //    31 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the Fejer type 1 rule.
        //    1 <= ORDER.
        //
        //    Input, int I, the index of the desired abscissa.  
        //    1 <= I <= ORDER.
        //
        //    Output, double F1_ABSCISSA, the value of the I-th 
        //    abscissa in the Fejer type 1 rule of order ORDER.
        //
    {
        double value;

        switch (order)
        {
            case < 1:
                value = -typeMethods.r8_huge();
                return value;
        }

        if (i < 1 || order < i)
        {
            Console.WriteLine("");
            Console.WriteLine("F1_ABSCISSA - Fatal error!");
            Console.WriteLine("  1 <= I <= ORDER is required.");
            Console.WriteLine("  I = " + i + "");
            Console.WriteLine("  ORDER = " + order + "");
            return 1;
        }

        switch (order)
        {
            case 1:
                value = 0.0;
                break;
            default:
            {
                if (2 * (2 * order + 1 - 2 * i) == 2 * order)
                {
                    value = 0.0;
                }
                else
                {
                    value = Math.Cos((2 * order + 1 - 2 * i) * Math.PI
                                     / (2 * order));
                }

                break;
            }
        }

        return value;
    }

    public static double[] f1_weights(int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1_WEIGHTS computes weights for a Fejer type 1 rule.
        //
        //  Modified:
        //
        //    28 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //    Walter Gautschi,
        //    Numerical Quadrature in the Presence of a Singularity,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 4, Number 3, 1967, pages 357-362.
        //
        //    Joerg Waldvogel,
        //    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
        //    BIT Numerical Mathematics,
        //    Volume 43, Number 1, 2003, pages 1-18.
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Output, double F1_WEIGHTS[ORDER], the weights.
        //
    {
        int i;

        switch (order)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("F1_WEIGHTS - Fatal error!");
                Console.WriteLine("  ORDER < 1.");
                return null;
        }

        double[] w = new double[order];

        switch (order)
        {
            case 1:
                w[0] = 2.0;
                return w;
        }

        double[] theta = new double[order];

        for (i = 1; i <= order; i++)
        {
            theta[i - 1] = (2 * (order + 1 - i) - 1) * Math.PI
                           / (2 * order);
        }

        for (i = 1; i <= order; i++)
        {
            w[i - 1] = 1.0;
            int j;
            for (j = 1; j <= order / 2; j++)
            {
                w[i - 1] -= 2.0 * Math.Cos(2.0 * j * theta[i - 1])
                            / (4 * j * j - 1);
            }
        }

        for (i = 0; i < order; i++)
        {
            w[i] = 2.0 * w[i] / order;
        }

        return w;
    }
}