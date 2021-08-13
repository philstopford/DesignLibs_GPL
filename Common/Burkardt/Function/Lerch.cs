using System;

namespace Burkardt.Function
{
    public static class Lerch
    {
        public static double lerch(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LERCH estimates the Lerch transcendent function.
        //
        //  Discussion:
        //
        //    The Lerch transcendent function is defined as:
        //
        //      LERCH ( A, B, C ) = Sum ( 0 <= K < Infinity ) A^K / ( C + K )^B
        //
        //    excluding any term with ( C + K ) = 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein, editor,
        //    CRC Concise Encylopedia of Mathematics,
        //    CRC Press, 1998.
        //
        //  Thanks:
        //
        //    Oscar van Vlijmen
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the function.
        //
        //    Output, double LERCH, an approximation to the Lerch
        //    transcendent function.
        //
        {
            double sum2 = 0.0;
            int k = 0;
            double a_k = 1.0;

            for (;;)
            {
                double sum2_old = sum2;

                if (c + (double) (k) == 0.0)
                {
                    k = k + 1;
                    a_k = a_k * a;
                    continue;
                }

                sum2 = sum2 + a_k / Math.Pow(c + (double) (k), b);

                if (sum2 <= sum2_old)
                {
                    break;
                }

                k = k + 1;
                a_k = a_k * a;
            }

            return sum2;
        }
    }
}