using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double psi(double xx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
        //
        //  Discussion:
        //
        //    The main computation involves evaluation of rational Chebyshev
        //    approximations.  PSI was written at Argonne National Laboratory
        //    for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2021
        //
        //  Author:
        //
        //    Barry Brown, James Lovato, Kathy Russell.
        //
        //  Reference:
        //
        //    William Cody, Strecok and Thacher,
        //    Chebyshev Approximations for the Psi Function,
        //    Mathematics of Computation,
        //    Volume 27, 1973, pages 123-127.
        //
        //  Parameters:
        //
        //    Input, double *XX, the argument of the psi function.
        //
        //    Output, double PSI, the value of the psi function.  PSI
        //    is assigned the value 0 when the psi function is undefined.
        //
    {
        const double dx0 = 1.461632144968362341262659542325721325e0;
        const double piov4 = .785398163397448e0;
        double[] p1 =  {
                .895385022981970e-02,.477762828042627e+01,.142441585084029e+03,
                .118645200713425e+04,.363351846806499e+04,.413810161269013e+04,
                .130560269827897e+04
            }
            ;
        double[] p2 =  {
                -.212940445131011e+01,-.701677227766759e+01,-.448616543918019e+01,
                -.648157123766197e+00
            }
            ;
        double[] q1 =  {
                .448452573429826e+02,.520752771467162e+03,.221000799247830e+04,
                .364127349079381e+04,.190831076596300e+04,.691091682714533e-05
            }
            ;
        double[] q2 =  {
                .322703493791143e+02,.892920700481861e+02,.546117738103215e+02,
                .777788548522962e+01
            }
            ;
        const int K1 = 3;
        const int K2 = 1;
        double psi, den, upper, w;
        int i;
        //
        //  MACHINE DEPENDENT CONSTANTS:
        //
        //  XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
        //  WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
        //  AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
        //  ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
        //  PSI MAY BE REPRESENTED AS ALOG(X).
        //
        //  XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
        //  MAY BE REPRESENTED BY 1/X.
        //
        double xmax1 = ipmpar(K1);
        xmax1 = Math.Min(xmax1, 1.0e0 / dpmpar(K2));
        double xsmall = 1e-9;
        double x = xx;
        double aug = 0.0e0;
        switch (x)
        {
            case >= 0.5e0:
                goto S50;
        }

        //
        //  X < 0.5,  USE REFLECTION FORMULA
        //  PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
        //
        if (Math.Abs(x) > xsmall)
        {
            goto S10;
        }

        switch (x)
        {
            case 0.0e0:
                goto S100;
        }

        //
        //     0 < ABS(X) <= XSMALL.  USE 1/X AS A SUBSTITUTE
        //     FOR  PI*COTAN(PI*X)
        //
        aug = -(1.0e0 / x);
        goto S40;
        S10:
        //
        //  REDUCTION OF ARGUMENT FOR COTAN
        //
        w = -x;
        double sgn = piov4;
        switch (w)
        {
            case > 0.0e0:
                goto S20;
        }

        w = -w;
        sgn = -sgn;
        S20:
        //
        //  MAKE AN ERROR EXIT IF X <= -XMAX1
        //
        if (w >= xmax1)
        {
            goto S100;
        }

        int nq = (int)Math.Truncate(w);
        w -= nq;
        nq = (int)Math.Truncate(w * 4.0e0);
        w = 4.0e0 * (w - nq * .25e0);
        //
        //  W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
        //  ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
        //  QUADRANT AND DETERMINE SIGN
        //
        int n = nq / 2;
        if (n + n != nq)
        {
            w = 1.0e0 - w;
        }

        double z = piov4 * w;
        int m = n / 2;
        if (m + m != n)
        {
            sgn = -sgn;
        }

        //
        //  DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
        //
        n = (nq + 1) / 2;
        m = n / 2;
        m += m;
        if (m != n)
        {
            goto S30;
        }

        switch (z)
        {
            //
            //  CHECK FOR SINGULARITY
            //
            case 0.0e0:
                goto S100;
        }

        //
        //  USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
        //  SIN/COS AS A SUBSTITUTE FOR TAN
        //
        aug = sgn * (Math.Cos(z) / Math.Sin(z) * 4.0e0);
        goto S40;
        S30:
        aug = sgn * (Math.Sin(z) / Math.Cos(z) * 4.0e0);
        S40:
        x = 1.0e0 - x;
        S50:
        switch (x)
        {
            case > 3.0e0:
                goto S70;
        }

        //
        //  0.5 <= X <= 3.0
        //
        den = x;
        upper = p1[0] * x;
        for (i = 1; i <= 5; i++)
        {
            den = (den + q1[i - 1]) * x;
            upper = (upper + p1[i + 1 - 1]) * x;
        }

        den = (upper + p1[6]) / (den + q1[5]);
        double xmx0 = x - dx0;
        psi = den * xmx0 + aug;
        return psi;
        S70:
        //
        //  IF X .GE. XMAX1, PSI = LN(X)
        //
        if (x >= xmax1)
        {
            goto S90;
        }

        //
        //  3.0 < X < XMAX1
        //
        w = 1.0e0 / (x * x);
        den = w;
        upper = p2[0] * w;
        for (i = 1; i <= 3; i++)
        {
            den = (den + q2[i - 1]) * w;
            upper = (upper + p2[i + 1 - 1]) * w;
        }

        aug = upper / (den + q2[3]) - 0.5e0 / x + aug;
        S90:
        psi = aug + Math.Log(x);
        return psi;
        S100:
        //
        //  ERROR RETURN
        //
        psi = 0.0e0;
        return psi;
    }

    public static void psi_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PSI_VALUES returns some values of the Psi or Digamma function.
        //
        //  Discussion:
        //
        //    PSI(X) = d LN ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
        //
        //    PSI(1) = - Euler's constant.
        //
        //    PSI(X+1) = PSI(X) + 1 / X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz and Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    US Department of Commerce, 1964.
        //
        //  Parameters:
        //
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 11;

        double[] fx_vec =  {
                -0.5772156649E+00, -0.4237549404E+00, -0.2890398966E+00,
                -0.1691908889E+00, -0.0613845446E+00, -0.0364899740E+00,
                0.1260474528E+00, 0.2085478749E+00, 0.2849914333E+00,
                0.3561841612E+00, 0.4227843351E+00
            }
            ;
        double[] x_vec =  {
                1.0E+00, 1.1E+00, 1.2E+00,
                1.3E+00, 1.4E+00, 1.5E+00,
                1.6E+00, 1.7E+00, 1.8E+00,
                1.9E+00, 2.0E+00
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            x = 0.0E+00;
            fx = 0.0E+00;
        }
        else
        {
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}