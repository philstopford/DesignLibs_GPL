using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_power(double r, int p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_POWER computes an integer Math.Power of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the base.
        //
        //    Input, int P, the Math.Power, which may be negative.
        //
        //    Output, double R8_POWER, the value of R^P.
        //
    {
        double value = 0;
        switch (p)
        {
            //
            //  Special case.  R^0 = 1.
            //
            case 0:
                value = 1.0;
                break;
            //
            default:
            {
                switch (r)
                {
                    case 0.0 when 0 < p:
                        value = 0.0;
                        break;
                    case 0.0:
                        value = Math.Pow(r, p);
                        break;
                    default:
                    {
                        value = p switch
                        {
                            >= 1 => Math.Pow(r, p),
                            _ => Math.Pow(r, p)
                        };

                        break;
                    }
                }

                break;
            }
        }

        return value;
    }

    public static double r8_power_fast(double r, int p, ref int mults)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_POWER_FAST computes the P-th Math.Power of R, for real R and integer P.
        //
        //  Discussion:
        //
        //    Obviously, R^P can be computed using P-1 multiplications.
        //
        //    However, R^P can also be computed using at most 2*LOG2(P) multiplications.
        //    To do the calculation this way, let N = LOG2(P).
        //    Compute A, A^2, A^4, ..., A^N by N-1 successive squarings.
        //    Start the value of R^P at A, and each time that there is a 1 in
        //    the binary expansion of P, multiply by the current result of the squarings.
        //
        //    This algorithm is not optimal.  For small exponents, and for special
        //    cases, the result can be computed even more quickly.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the base.
        //
        //    Input, int P, the Math.Power, which may be negative.
        //
        //    Output, int &MULTS, the number of multiplications and divisions.
        //
        //    Output, double R8_POWER_FAST, the value of R^P.
        //
    {
        double value = 0;

        mults = 0;
        switch (r)
        {
            //
            //  Special bases.
            //
            case 1.0:
                value = 1.0;
                return value;
            case -1.0:
            {
                value = (p % 2) switch
                {
                    1 => -1.0,
                    _ => 1.0
                };

                return value;
            }
            case 0.0 when p <= 0:
                Console.WriteLine("");
                Console.WriteLine("R8_POWER_FAST - Fatal error!");
                Console.WriteLine("  Base is zero, and exponent is negative.");
                return 1;
            case 0.0:
                value = 0.0;
                return value;
        }

        switch (p)
        {
            //
            //  Special Math.Powers.
            //
            case -1:
                value = 1.0 / r;
                mults += 1;
                return value;
            case 0:
                value = 1.0;
                return value;
            case 1:
                value = r;
                return value;
        }

        //
        //  Some work to do.
        //
        int p_mag = Math.Abs(p);
        int p_sign = i4_sign(p);

        value = 1.0;
        double r2 = r;

        while (0 < p_mag)
        {
            switch (p_mag % 2)
            {
                case 1:
                    value *= r2;
                    mults += 1;
                    break;
            }

            p_mag /= 2;
            r2 *= r2;
            mults += 1;
        }

        switch (p_sign)
        {
            case -1:
                value = 1.0 / value;
                mults += 1;
                break;
        }

        return value;
    }

}