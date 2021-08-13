using System;

namespace Burkardt.Function
{
    public static class Zeta
    {
        public static double zeta_m1(double p, double tol)

            //****************************************************************************80
            //
            //  Purpose:
            //  
            //    ZETA_M1 estimates the Riemann Zeta function minus 1.
            //
            //  Discussion:
            //
            //    This function includes the Euler-McLaurin correction.
            //
            //    ZETA_M1 ( P ) = ZETA ( P ) - 1
            //
            //    ZETA(P) has the form 1 + small terms.  Computing ZETA(P)-1
            //    allows for greater accuracy in the small terms.
            //
            //  Definition:
            //
            //    For 1 < P, the Riemann Zeta function is defined as:
            //
            //      ZETA ( P ) = Sum ( 1 <= N < Infinity ) 1 / N^P
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Thompson,
            //    Atlas for Computing Mathematical Functions,
            //    Wiley, 1997,
            //    ISBN: 0471181714,
            //    LC: QA331 T385
            //
            //  Parameters:
            //
            //    Input, double P, the Math.Power to which the integers are raised.
            //    P must be greater than 1.
            //
            //    Input, double TOL, the requested relative tolerance.
            //
            //    Output, double ZETA_M1, an approximation to the Riemann
            //    Zeta function minus 1.
            //
        {
            double base_;
            int k;
            int n;
            double negp;
            double nsterm;
            double value;

            if (p <= 1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("ZETA_M1 - Fatal error!");
                Console.WriteLine("  Exponent P <= 1.0.");
            }

            nsterm = p * (p + 1.0) * (p + 2.0) * (p + 3.0) * (p + 4.0)
                     / 30240.0;

            base_ = nsterm * Math.Pow(2.0, p) / tol;

            n = (int)Math.Pow(base_, 1.0 / (p + 5.0));
            if (n < 10)
            {
                n = 10;
            }

            negp = -p;
            value = 0.0;
            for (k = 2; k < n; k++)
            {
                base_ = (double)(k);
                value = value + Math.Pow(base_, negp);
            }

            /*
            Euler-McLaurin correction.
            */
            base_ = (double)(n);

            value = value + Math.Pow(base_, negp)
                * (0.5 + (double)(n) / (p - 1.0)
                       + p * (1.0 -
                              (p + 1.0) * (p + 2.0) / (double)(60 * n * n))
                       / (double)(12 * n)
                       + nsterm / Math.Pow(base_, p + 5.0));

            return value;
        }

        public static double zeta_naive(double p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZETA_NAIVE estimates the Riemann Zeta function.
            //
            //  Definition:
            //
            //    For 1 < P, the Riemann Zeta function is defined as:
            //
            //      ZETA ( P ) = Sum ( 1 <= N < +oo ) 1 / N^P
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input, double P, the Math.Power to which the integers are raised.
            //    P must be greater than 1.  
            //
            //    Output, double ZETA_NAIVE, an approximation to the Riemann
            //    Zeta function.
            //
        {
            int n;
            double value;
            double value_old;

            if (p <= 1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("ZETA_NAIVE - Fatal error!");
                Console.WriteLine("  Exponent P <= 1.0.");
                return (1);
            }

            value = 0.0;
            n = 0;

            for (;;)
            {
                n = n + 1;
                value_old = value;
                value = value + 1.0 / Math.Pow((double)n, p);

                if (value <= value_old || 1000 <= n)
                {
                    break;
                }

            }

            return value;
        }

    }
}