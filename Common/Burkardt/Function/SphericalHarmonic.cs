using System;

namespace Burkardt.Function;

public static class SphericalHarmonic
{
    public static void spherical_harmonic(int l, int m, double theta, double phi,
            ref double[] c, ref double[] s )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERICAL_HARMONIC evaluates spherical harmonic functions.
        //
        //  Discussion:
        //
        //    The spherical harmonic function Y(L,M,THETA,PHI,X) is the
        //    angular part of the solution to Laplace's equation in spherical
        //    coordinates.
        //
        //    Y(L,M,THETA,PHI,X) is related to the associated Legendre
        //    function as follows:
        //
        //      Y(L,M,THETA,PHI,X) = FACTOR * P(L,M,cos(THETA)) * exp ( i * M * PHI )
        //
        //    Here, FACTOR is a normalization factor:
        //
        //      FACTOR = sqrt ( ( 2 * L + 1 ) * ( L - M )! / ( 4 * PI * ( L + M )! ) )
        //
        //    In Mathematica, a spherical harmonic function can be evaluated by
        //
        //      SphericalHarmonicY [ l, m, theta, phi ]
        //
        //    Note that notational tradition in physics requires that THETA
        //    and PHI represent the reverse of what they would normally mean
        //    in mathematical notation; that is, THETA goes up and down, and
        //    PHI goes around.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Eric Weisstein,
        //    CRC Concise Encyclopedia of Mathematics,
        //    CRC Press, 2002,
        //    Second edition,
        //    ISBN: 1584883472,
        //    LC: QA5.W45.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input, int L, the first index of the spherical harmonic function.
        //    Normally, 0 <= L.
        //
        //    Input, int M, the second index of the spherical harmonic function.
        //    Normally, -L <= M <= L.
        //
        //    Input, double THETA, the polar angle, for which
        //    0 <= THETA <= PI.
        //
        //    Input, double PHI, the longitudinal angle, for which
        //    0 <= PHI <= 2*PI.
        //
        //    Output, double C[L+1], S[L+1], the real and imaginary
        //    parts of the functions Y(L,0:L,THETA,PHI).
        //
    {
        int i;

        int m_abs = Math.Abs(m);

        double[] plm = new double[l + 1];

        PolynomialNS.Legendre.legendre_associated_normalized(l, m_abs, Math.Cos(theta), ref plm);

        double angle = m * phi;

        switch (m)
        {
            case >= 0:
            {
                for (i = 0; i <= l; i++)
                {
                    c[i] = plm[i] * Math.Cos(angle);
                    s[i] = plm[i] * Math.Sin(angle);
                }

                break;
            }
            default:
            {
                for (i = 0; i <= l; i++)
                {
                    c[i] = -plm[i] * Math.Cos(angle);
                    s[i] = -plm[i] * Math.Sin(angle);
                }

                break;
            }
        }
    }

}