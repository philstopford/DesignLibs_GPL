using System;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

public static class Zernike
{
    public static double zernike_poly(int m, int n, double rho)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZERNIKE_POLY evaluates a Zernike polynomial at RHO.
        //
        //  Discussion:
        //
        //    This routine uses the facts that:
        //
        //    *) R^M_N = 0 if M < 0, or N < 0, or N < M.
        //    *) R^M_M = RHO^M
        //    *) R^M_N = 0 if mod ( N - M, 2 ) = 1.
        //
        //    and the recursion:
        //
        //    R^M_(N+2) = A * [ ( B * RHO^2 - C ) * R^M_N - D * R^M_(N-2) ]
        //
        //    where
        //
        //    A = ( N + 2 ) / ( ( N + 2 )^2 - M^2 )
        //    B = 4 * ( N + 1 )
        //    C = ( N + M )^2 / N + ( N - M + 2 )^2 / ( N + 2 )
        //    D = ( N^2 - M^2 ) / N
        //
        //    I wish I could clean up the recursion in the code, but for
        //    now, I have to treat the case M = 0 specially.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein,
        //    CRC Concise Encyclopedia of Mathematics,
        //    CRC Press, 2002,
        //    Second edition,
        //    ISBN: 1584883472,
        //    LC: QA5.W45.
        //
        //  Parameters:
        //
        //    Input, int M, the upper index.
        //
        //    Input, int N, the lower index.
        //
        //    Input, double RHO, the radial coordinate.
        //
        //    Output, double ZERNIKE_POLY, the value of the Zernike
        //    polynomial R^M_N at the point RHO.
        //
    {
        double a;
        double b;
        double c;
        double d;
        int nn;
        double z;
        double zm2;
        double zp2;
        switch (m)
        {
            //
            //  Do checks.
            //
            case < 0:
                z = 0.0;
                return z;
        }

        switch (n)
        {
            case < 0:
                z = 0.0;
                return z;
        }

        if (n < m)
        {
            z = 0.0;
            return z;
        }

        switch ((n - m) % 2)
        {
            case 1:
                z = 0.0;
                return z;
        }

        zm2 = 0.0;
        z = Math.Pow(rho, m);

        switch (m)
        {
            case 0 when n == 0:
                return z;
            case 0:
            {
                zm2 = z;
                z = 2.0 * rho * rho - 1.0;

                for (nn = m + 2; nn <= n - 2; nn += 2)
                {
                    a = (nn + 2)
                        / (double)((nn + 2) * (nn + 2) - m * m);

                    b = 4 * (nn + 1);

                    c = (nn + m) * (nn + m) / (double)nn
                        + (nn - m + 2) * (nn - m + 2)
                        / (double)(nn + 2);

                    d = (nn * nn - m * m) / (double)nn;

                    zp2 = a * ((b * rho * rho - c) * z - d * zm2);
                    zm2 = z;
                    z = zp2;
                }

                break;
            }
            default:
            {
                for (nn = m; nn <= n - 2; nn += 2)
                {
                    a = (nn + 2)
                        / (double)((nn + 2) * (nn + 2) - m * m);

                    b = 4 * (nn + 1);

                    c = (nn + m) * (nn + m) / (double)nn
                        + (nn - m + 2) * (nn - m + 2)
                        / (double)(nn + 2);

                    d = (nn * nn - m * m) / (double)nn;

                    zp2 = a * ((b * rho * rho - c) * z - d * zm2);
                    zm2 = z;
                    z = zp2;
                }

                break;
            }
        }

        return z;
    }

    public static double[] zernike_poly_coef(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZERNIKE_POLY_COEF: coefficients of a Zernike polynomial.
        //
        //  Discussion:
        //
        //    With our coefficients stored in C(0:N), the
        //    radial function R^M_N(RHO) is given by
        //
        //      R^M_N(RHO) = C(0)
        //                 + C(1) * RHO
        //                 + C(2) * RHO^2
        //                 + ...
        //                 + C(N) * RHO^N
        //
        //    and the odd and even Zernike polynomials are
        //
        //      Z^M_N(RHO,PHI,odd)  = R^M_N(RHO) * sin(PHI)
        //      Z^M_N(RHO,PHI,even) = R^M_N(RHO) * cos(PHI)
        //
        //    The first few "interesting" values of R are:
        //
        //    R^0_0 = 1
        //
        //    R^1_1 = RHO
        //
        //    R^0_2 = 2 * RHO^2 - 1
        //    R^2_2 =     RHO^2
        //
        //    R^1_3 = 3 * RHO^3 - 2 * RHO
        //    R^3_3 =     RHO^3
        //
        //    R^0_4 = 6 * RHO^4 - 6 * RHO^2 + 1
        //    R^2_4 = 4 * RHO^4 - 3 * RHO^2
        //    R^4_4 =     RHO^4
        //
        //    R^1_5 = 10 * RHO^5 - 12 * RHO^3 + 3 * RHO
        //    R^3_5 =  5 * RHO^5 -  4 * RHO^3
        //    R^5_5 =      RHO^5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein,
        //    CRC Concise Encyclopedia of Mathematics,
        //    CRC Press, 2002,
        //    Second edition,
        //    ISBN: 1584883472,
        //    LC: QA5.W45.
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the polynomial.
        //    Normally, 0 <= M <= N and 0 <= N.
        //
        //    Output, double ZERNIKE_POLY_COEF[N+1], the coefficients of the polynomial.
        //
    {
        double[] c;
        int l;
        int nm_minus;
        int nm_plus;

        c = new double[n + 1];

        typeMethods.r8vec_zero(n + 1, ref c);

        switch (n)
        {
            case < 0:
                return c;
        }

        switch (m)
        {
            case < 0:
                return c;
        }

        if (n < m)
        {
            return c;
        }

        switch ((n - m) % 2)
        {
            case 1:
                return c;
        }

        nm_plus = (m + n) / 2;
        nm_minus = (n - m) / 2;

        c[n] = typeMethods.r8_choose(n, nm_plus);

        for (l = 0; l <= nm_minus - 1; l++)
        {
            c[n - 2 * l - 2] = -(double)((nm_plus - l) * (nm_minus - l))
                * c[n - 2 * l] / ((n - l) * (l + 1));

        }

        return c;
    }

}