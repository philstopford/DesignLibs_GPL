using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {

        public static double dpmpar(int i)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DPMPAR provides machine constants for double precision arithmetic.
            //
            //  Discussion:
            //
            //     DPMPAR PROVIDES THE double PRECISION MACHINE CONSTANTS FOR
            //     THE COMPUTER BEING USED.   It is assumed that tHE
            //     double PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
            //     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX.
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
            //    ALFRED H. MORRIS, JR
            //
            //  Input:
            //
            //    int I, a value of 1, 2 or 3.
            //
            //  Output:
            //
            //    double DPMPAR, the value requested by I.
            //    DPMPAR(1) = B^(1 - M), THE MACHINE PRECISION,
            //    DPMPAR(2) = B^(EMIN - 1), THE SMALLEST MAGNITUDE,
            //    DPMPAR(3) = B^EMAX*(1 - B^(-M)), THE LARGEST MAGNITUDE.
        {
            double b;
            double binv;
            double bm1;
            int emax;
            int emin;
            int ibeta;
            int K1 = 4;
            int K2 = 8;
            int K3 = 9;
            int K4 = 10;
            int m;
            double one;
            double value;
            double w;
            double z;

            if (i == 1)
            {
                b = ipmpar(K1);
                m = ipmpar(K2);
                value = Math.Pow(b, (double) (1 - m));
            }
            else if (i == 2)
            {
                b = ipmpar(K1);
                emin = ipmpar(K3);
                one = 1.0;
                binv = one / b;
                w = Math.Pow(b, (double) (emin + 2));
                value = w * binv * binv * binv;
            }
            else
            {
                ibeta = ipmpar(K1);
                m = ipmpar(K2);
                emax = ipmpar(K4);
                b = ibeta;
                bm1 = ibeta - 1;
                one = 1.0;
                z = Math.Pow(b, (double) (m - 1));
                w = ((z - one) * b + bm1) / (b * z);
                z = Math.Pow(b, (double) (emax - 2));
                value = w * z * b * b;
            }

            return value;
        }

        public static void erf_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ERF_VALUES returns some values of the ERF or "error" function.
            //
            //  Definition:
            //
            //    ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
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
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.0000000000E+00, 0.1124629160E+00, 0.2227025892E+00, 0.3286267595E+00,
                0.4283923550E+00, 0.5204998778E+00, 0.6038560908E+00, 0.6778011938E+00,
                0.7421009647E+00, 0.7969082124E+00, 0.8427007929E+00, 0.8802050696E+00,
                0.9103139782E+00, 0.9340079449E+00, 0.9522851198E+00, 0.9661051465E+00,
                0.9763483833E+00, 0.9837904586E+00, 0.9890905016E+00, 0.9927904292E+00,
                0.9953222650E+00
            };
            double[] x_vec =
            {
                0.0E+00, 0.1E+00, 0.2E+00, 0.3E+00,
                0.4E+00, 0.5E+00, 0.6E+00, 0.7E+00,
                0.8E+00, 0.9E+00, 1.0E+00, 1.1E+00,
                1.2E+00, 1.3E+00, 1.4E+00, 1.5E+00,
                1.6E+00, 1.7E+00, 1.8E+00, 1.9E+00,
                2.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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

        public static double error_f(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ERROR_F evaluates the error function ERF.
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
            //  Parameters:
            //
            //    Input, double *X, the argument.
            //
            //    Output, double ERROR_F, the value of the error function at X.
            //
        {
            double[] a =
            {
                .771058495001320e-04, -.133733772997339e-02, .323076579225834e-01,
                .479137145607681e-01, .128379167095513e+00
            };
            double ax;
            double[] b =
            {
                .301048631703895e-02, .538971687740286e-01, .375795757275549e+00
            };
            double bot;
            double c = .564189583547756e0;
            double erf1;
            double[] p =
            {
                -1.36864857382717e-07, 5.64195517478974e-01, 7.21175825088309e+00,
                4.31622272220567e+01, 1.52989285046940e+02, 3.39320816734344e+02,
                4.51918953711873e+02, 3.00459261020162e+02
            };
            double[] q =
            {
                1.00000000000000e+00, 1.27827273196294e+01, 7.70001529352295e+01,
                2.77585444743988e+02, 6.38980264465631e+02, 9.31354094850610e+02,
                7.90950925327898e+02, 3.00459260956983e+02
            };
            double[] r =
            {
                2.10144126479064e+00, 2.62370141675169e+01, 2.13688200555087e+01,
                4.65807828718470e+00, 2.82094791773523e-01
            };
            double[] s =
            {
                9.41537750555460e+01, 1.87114811799590e+02, 9.90191814623914e+01,
                1.80124575948747e+01
            };
            double t;
            double top;
            double x2;

            ax = Math.Abs(x);

            if (ax <= 0.5e0)
            {
                t = x * x;
                top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0e0;
                bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0e0;
                erf1 = x * (top / bot);
                return erf1;
            }

            if (ax > 4.0e0) goto S20;
            top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[
                7];
            bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[
                7];
            erf1 = 0.5e0 + (0.5e0 - Math.Exp(-(x * x)) * top / bot);
            if (x < 0.0e0) erf1 = -erf1;
            return erf1;

            S20:
            if (ax >= 5.8e0) goto S30;
            x2 = x * x;
            t = 1.0e0 / x2;
            top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
            bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0e0;
            erf1 = (c - top / (x2 * bot)) / ax;
            erf1 = 0.5e0 + (0.5e0 - Math.Exp(-x2) * erf1);
            if (x < 0.0e0) erf1 = -erf1;
            return erf1;

            S30:
            erf1 = fifdsign(1.0e0, x);

            return erf1;
        }

        public static double error_fc(int ind, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ERROR_FC evaluates the complementary error function ERFC.
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
            //  Parameters:
            //
            //    Input, int *IND, chooses the scaling.
            //    If IND is nonzero, then the value returned has been multiplied by
            //    EXP(X*X).
            //
            //    Input, double *X, the argument of the function.
            //
            //    Output, double ERROR_FC, the value of the complementary
            //    error function.
            //
        {
            double c = .564189583547756e0;
            double[] a =
            {
                .771058495001320e-04, -.133733772997339e-02, .323076579225834e-01,
                .479137145607681e-01, .128379167095513e+00
            };
            double[] b =
            {
                .301048631703895e-02, .538971687740286e-01, .375795757275549e+00
            };
            double[] p =
            {
                -1.36864857382717e-07, 5.64195517478974e-01, 7.21175825088309e+00,
                4.31622272220567e+01, 1.52989285046940e+02, 3.39320816734344e+02,
                4.51918953711873e+02, 3.00459261020162e+02
            };
            double[] q =
            {
                1.00000000000000e+00, 1.27827273196294e+01, 7.70001529352295e+01,
                2.77585444743988e+02, 6.38980264465631e+02, 9.31354094850610e+02,
                7.90950925327898e+02, 3.00459260956983e+02
            };
            double[] r =
            {
                2.10144126479064e+00, 2.62370141675169e+01, 2.13688200555087e+01,
                4.65807828718470e+00, 2.82094791773523e-01
            };
            double[] s =
            {
                9.41537750555460e+01, 1.87114811799590e+02, 9.90191814623914e+01,
                1.80124575948747e+01
            };
            int K1 = 1;
            double erfc1, ax, bot, e, t, top, w;

            //
            //  ABS(X) <= 0.5
            //
            ax = Math.Abs(x);
            if (ax > 0.5e0) goto S10;
            t = x * x;
            top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0e0;
            bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0e0;
            erfc1 = 0.5e0 + (0.5e0 - x * (top / bot));
            if (ind != 0) erfc1 = Math.Exp(t) * erfc1;
            return erfc1;
            S10:
            //
            //  0.5 < ABS(X) <= 4
            //
            if (ax > 4.0e0) goto S20;
            top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[
                7];
            bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[
                7];
            erfc1 = top / bot;
            goto S40;
            S20:
            //
            //  4 < ABS(X)
            //
            if (x <= -5.6e0) goto S60;
            if (ind != 0) goto S30;
            if (x > 100.0e0) goto S70;
            if (x * x > -exparg(K1)) goto S70;
            S30:
            t = Math.Pow(1.0e0 / x, 2.0);
            top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
            bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0e0;
            erfc1 = (c - t * top / bot) / ax;
            S40:
            //
            //  FINAL ASSEMBLY
            //
            if (ind == 0) goto S50;
            if (x < 0.0e0) erfc1 = 2.0e0 * Math.Exp(x * x) - erfc1;
            return erfc1;
            S50:
            w = x * x;
            t = w;
            e = w - t;
            erfc1 = (0.5e0 + (0.5e0 - e)) * Math.Exp(-t) * erfc1;
            if (x < 0.0e0) erfc1 = 2.0e0 - erfc1;
            return erfc1;
            S60:
            //
            //  LIMIT VALUE FOR LARGE NEGATIVE X
            //
            erfc1 = 2.0e0;
            if (ind != 0) erfc1 = 2.0e0 * Math.Exp(x * x);
            return erfc1;
            S70:
            //
            //  LIMIT VALUE FOR LARGE POSITIVE X WHEN IND = 0
            //
            erfc1 = 0.0e0;
            return erfc1;
        }

        public static double esum(int mu, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ESUM evaluates exp ( MU + X ).
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
            //  Parameters:
            //
            //    Input, int *MU, part of the argument.
            //
            //    Input, double *X, part of the argument.
            //
            //    Output, double ESUM, the value of exp ( MU + X ).
            //
        {
            double esum, w;

            if (x > 0.0e0) goto S10;
            if (mu < 0) goto S20;
            w = (double) mu + x;
            if (w > 0.0e0) goto S20;
            esum = Math.Exp(w);
            return esum;
            S10:
            if (mu > 0) goto S20;
            w = (double) mu + x;
            if (w < 0.0e0) goto S20;
            esum = Math.Exp(w);
            return esum;
            S20:
            w = mu;
            esum = Math.Exp(w) * Math.Exp(x);
            return esum;
        }

        public static double eval_pol(double[] a, int n, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EVAL_POL evaluates a polynomial at X.
            //
            //  Discussion:
            //
            //    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X^N
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
            //  Input:
            //
            //    double A(0:N), coefficients of the polynomial.
            //
            //    int *N, length of A.
            //
            //    double *X, the point at which the polynomial
            //    is to be evaluated.
            //
            //  Output:
            //
            //    double EVAL_POL, the value of the polynomial at X.
            //
        {
            double devlpl, term;
            int i;

            term = a[n - 1];
            for (i = n - 1 - 1; i >= 0; i--)
            {
                term = a[i] + term * x;
            }

            devlpl = term;
            return devlpl;
        }

        public static double exparg(int l)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXPARG returns the largest or smallest legal argument for EXP.
            //
            //  Discussion:
            //
            //    Only an approximate limit for the argument of EXP is desired.
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
            //  Input:
            //
            //    int *L, indicates which limit is desired.
            //    If L = 0, then the largest positive argument for EXP is desired.
            //    Otherwise, the largest negative argument for EXP for which the
            //    result is nonzero is desired.
            //
            //  Output:
            //
            //    double EXPARG, the desired value.
            //
        {
            int b;
            double exparg;
            int K1 = 4;
            int K2 = 9;
            int K3 = 10;
            double lnb;
            int m;

            b = ipmpar(K1);
            if (b != 2) goto S10;
            lnb = .69314718055995e0;
            goto S40;
            S10:
            if (b != 8) goto S20;
            lnb = 2.0794415416798e0;
            goto S40;
            S20:
            if (b != 16) goto S30;
            lnb = 2.7725887222398e0;
            goto S40;
            S30:
            lnb = Math.Log((double) b);
            S40:
            if (l == 0) goto S50;
            m = ipmpar(K2) - 1;
            exparg = 0.99999e0 * ((double) m * lnb);
            return exparg;
            S50:
            m = ipmpar(K3);
            exparg = 0.99999e0 * ((double) m * lnb);
            return exparg;
        }

        public static double fifdsign(double mag, double sign)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
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
            //  Parameters:
            //
            //  mag     -     magnitude
            //  sign    -     sign to be transfered
            //
        {
            if (mag < 0) mag = -mag;
            if (sign < 0) mag = -mag;
            return mag;

        }

        public static double fpser(double a, double b, double x, double eps)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FPSER evaluates IX(A,B)(X) for very small B.
            //
            //  Discussion:
            //
            //    This routine is appropriate for use when
            //
            //      B < min ( EPS, EPS * A )
            //
            //    and
            //
            //      X <= 0.5.
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
            //  Parameters:
            //
            //    Input, double *A, *B, parameters of the function.
            //
            //    Input, double *X, the point at which the function is to
            //    be evaluated.
            //
            //    Input, double *EPS, a tolerance.
            //
            //    Output, double FPSER, the value of IX(A,B)(X).
            //
        {
            double an;
            double c;
            double fpser;
            int K1 = 1;
            double s;
            double t;
            double tol;

            fpser = 1.0e0;
            if (a <= 1e-3 * eps) goto S10;
            fpser = 0.0e0;
            t = a * Math.Log(x);
            if (t < exparg(K1)) return fpser;
            fpser = Math.Exp(t);
            S10:
            //
            //                NOTE THAT 1/B(A,B) = B
            //
            fpser = b / a * fpser;
            tol = eps / a;
            an = a + 1.0e0;
            t = x;
            s = t / an;
            S20:
            an = an + 1.0e0;
            t = x * t;
            c = t / an;
            s = s + c;
            if (Math.Abs(c) > tol) goto S20;
            fpser = fpser * (1.0e0 + a * s);

            return fpser;
        }

        public static double gam1(double a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5D+00 <= A <= 1.5
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
            //  Parameters:
            //
            //    Input, double *A, forms the argument of the Gamma function.
            //
            //    Output, double GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
            //
        {
            double s1 = .273076135303957e+00;
            double s2 = .559398236957378e-01;
            double[] p =
            {
                .577215664901533e+00, -.409078193005776e+00, -.230975380857675e+00,
                .597275330452234e-01, .766968181649490e-02, -.514889771323592e-02,
                .589597428611429e-03
            };
            double[] q =
            {
                .100000000000000e+01, .427569613095214e+00, .158451672430138e+00,
                .261132021441447e-01, .423244297896961e-02
            };
            double[] r =
            {
                -.422784335098468e+00, -.771330383816272e+00, -.244757765222226e+00,
                .118378989872749e+00, .930357293360349e-03, -.118290993445146e-01,
                .223047661158249e-02, .266505979058923e-03, -.132674909766242e-03
            };
            double gam1, bot, d, t, top, w, T1;

            t = a;
            d = a - 0.5e0;
            if (d > 0.0e0) t = d - 0.5e0;
            T1 = t;
            if (T1 < 0) goto S40;
            else if (T1 == 0) goto S10;
            else goto S20;
            S10:
            gam1 = 0.0e0;
            return gam1;
            S20:
            top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]) * t + p[0];
            bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.0e0;
            w = top / bot;
            if (d > 0.0e0) goto S30;
            gam1 = a * w;
            return gam1;
            S30:
            gam1 = t / a * (w - 0.5e0 - 0.5e0);
            return gam1;
            S40:
            top = (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]) * t + r[3]) * t + r[2]) * t + r[1]) * t +
                  r[0];
            bot = (s2 * t + s1) * t + 1.0e0;
            w = top / bot;
            if (d > 0.0e0) goto S50;
            gam1 = a * (w + 0.5e0 + 0.5e0);
            return gam1;
            S50:
            gam1 = t * w / a;
            return gam1;
        }

        public static double gsumln(double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GSUMLN evaluates the function ln(Gamma(A + B)).
            //
            //  Discussion:
            //
            //    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
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
            //  Parameters:
            //
            //    Input, double *A, *B, values whose sum is the argument of
            //    the Gamma function.
            //
            //    Output, double GSUMLN, the value of ln(Gamma(A+B)).
            //
        {
            double gsumln;
            double T1;
            double T2;
            double x;

            x = a + b - 2e0;
            if (x > 0.25e0) goto S10;
            T1 = 1.0e0 + x;
            gsumln = gamma_ln1(T1);
            return gsumln;
            S10:
            if (x > 1.25e0) goto S20;
            gsumln = gamma_ln1(x) + alnrel(x);
            return gsumln;
            S20:
            T2 = x - 1.0e0;
            gsumln = gamma_ln1(T2) + Math.Log(x * (1.0e0 + x));
            return gsumln;
        }

        public static int ipmpar(int i)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    IPMPAR returns integer machine constants.
            //
            //  Discussion:
            //
            //    Input arguments 1 through 3 are queries about integer arithmetic.
            //    We assume integers are represented in the N-digit, base-A form
            //
            //      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
            //
            //    where 0 <= X(0:N-1) < A.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //    Then:
            //
            //      IPMPAR(1) = A, the base of integer arithmetic;
            //      IPMPAR(2) = N, the number of base A digits;
            //      IPMPAR(3) = A^N - 1, the largest magnitude.
            //
            //    It is assumed that the single and double precision floating
            //    point arithmetics have the same base, say B, and that the
            //    nonzero numbers are represented in the form
            //
            //      sign * (B^E) * (X(1)/B + ... + X(M)/B^M)
            //
            //    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
            //    EMIN <= E <= EMAX.
            //
            //    Input argument 4 is a query about the base of real arithmetic:
            //
            //      IPMPAR(4) = B, the base of single and double precision arithmetic.
            //
            //    Input arguments 5 through 7 are queries about single precision
            //    floating point arithmetic:
            //
            //     IPMPAR(5) = M, the number of base B digits for single precision.
            //     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
            //     IPMPAR(7) = EMAX, the largest exponent E for single precision.
            //
            //    Input arguments 8 through 10 are queries about double precision
            //    floating point arithmetic:
            //
            //     IPMPAR(8) = M, the number of base B digits for double precision.
            //     IPMPAR(9) = EMIN, the smallest exponent E for double precision.
            //     IPMPAR(10) = EMAX, the largest exponent E for double precision.
            //
            //  Author:
            //
            //    Barry Brown, James Lovato, Kathy Russell.
            //
            //  Reference:
            //
            //    Phyllis Fox, Andrew Hall, and Norman Schryer,
            //    Algorithm 528,
            //    Framework for a Portable FORTRAN Subroutine Library,
            //    ACM Transactions on Mathematical Software,
            //    Volume 4, 1978, pages 176-188.
            //
            //  Parameters:
            //
            //    Input, int *I, the index of the desired constant.
            //
            //    Output, int IPMPAR, the value of the desired constant.
            //
        {
            int[] imach = new int [11];
            int ipmpar;
            //     MACHINE CONSTANTS FOR AMDAHL MACHINES.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 16;
            //   imach[5] = 6;
            //   imach[6] = -64;
            //   imach[7] = 63;
            //   imach[8] = 14;
            //   imach[9] = -64;
            //   imach[10] = 63;
            //
            //     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
            //       PC 7300, AND AT&T 6300.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -125;
            //   imach[7] = 128;
            //   imach[8] = 53;
            //   imach[9] = -1021;
            //   imach[10] = 1024;
            //
            //     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
            //
            //   imach[1] = 2;
            //   imach[2] = 33;
            //   imach[3] = 8589934591;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -256;
            //   imach[7] = 255;
            //   imach[8] = 60;
            //   imach[9] = -256;
            //   imach[10] = 255;
            //
            //     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
            //
            //   imach[1] = 2;
            //   imach[2] = 39;
            //   imach[3] = 549755813887;
            //   imach[4] = 8;
            //   imach[5] = 13;
            //   imach[6] = -50;
            //   imach[7] = 76;
            //   imach[8] = 26;
            //   imach[9] = -50;
            //   imach[10] = 76;
            //
            //     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
            //
            //   imach[1] = 2;
            //   imach[2] = 39;
            //   imach[3] = 549755813887;
            //   imach[4] = 8;
            //   imach[5] = 13;
            //   imach[6] = -50;
            //   imach[7] = 76;
            //   imach[8] = 26;
            //   imach[9] = -32754;
            //   imach[10] = 32780;
            //
            //     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
            //       60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
            //       ARITHMETIC (NOS OPERATING SYSTEM).
            //
            //   imach[1] = 2;
            //   imach[2] = 48;
            //   imach[3] = 281474976710655;
            //   imach[4] = 2;
            //   imach[5] = 48;
            //   imach[6] = -974;
            //   imach[7] = 1070;
            //   imach[8] = 95;
            //   imach[9] = -926;
            //   imach[10] = 1070;
            //
            //     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
            //       ARITHMETIC (NOS/VE OPERATING SYSTEM).
            //
            //   imach[1] = 2;
            //   imach[2] = 63;
            //   imach[3] = 9223372036854775807;
            //   imach[4] = 2;
            //   imach[5] = 48;
            //   imach[6] = -4096;
            //   imach[7] = 4095;
            //   imach[8] = 96;
            //   imach[9] = -4096;
            //   imach[10] = 4095;
            //
            //     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
            //
            //   imach[1] = 2;
            //   imach[2] = 63;
            //   imach[3] = 9223372036854775807;
            //   imach[4] = 2;
            //   imach[5] = 47;
            //   imach[6] = -8189;
            //   imach[7] = 8190;
            //   imach[8] = 94;
            //   imach[9] = -8099;
            //   imach[10] = 8190;
            //
            //     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
            //
            //   imach[1] = 2;
            //   imach[2] = 15;
            //   imach[3] = 32767;
            //   imach[4] = 16;
            //   imach[5] = 6;
            //   imach[6] = -64;
            //   imach[7] = 63;
            //   imach[8] = 14;
            //   imach[9] = -64;
            //   imach[10] = 63;
            //
            //     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
            //       THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
            //       5/7/9 AND THE SEL SYSTEMS 85/86.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 16;
            //   imach[5] = 6;
            //   imach[6] = -64;
            //   imach[7] = 63;
            //   imach[8] = 14;
            //   imach[9] = -64;
            //   imach[10] = 63;
            //
            //     MACHINE CONSTANTS FOR THE IBM PC.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -125;
            //   imach[7] = 128;
            //   imach[8] = 53;
            //   imach[9] = -1021;
            //   imach[10] = 1024;
            //
            //     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
            //       MACFORTRAN II.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -125;
            //   imach[7] = 128;
            //   imach[8] = 53;
            //   imach[9] = -1021;
            //   imach[10] = 1024;
            //
            //     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -127;
            //   imach[7] = 127;
            //   imach[8] = 56;
            //   imach[9] = -127;
            //   imach[10] = 127;
            //
            //     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
            //
            //   imach[1] = 2;
            //   imach[2] = 35;
            //   imach[3] = 34359738367;
            //   imach[4] = 2;
            //   imach[5] = 27;
            //   imach[6] = -128;
            //   imach[7] = 127;
            //   imach[8] = 54;
            //   imach[9] = -101;
            //   imach[10] = 127;
            //
            //     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
            //
            //   imach[1] = 2;
            //   imach[2] = 35;
            //   imach[3] = 34359738367;
            //   imach[4] = 2;
            //   imach[5] = 27;
            //   imach[6] = -128;
            //   imach[7] = 127;
            //   imach[8] = 62;
            //   imach[9] = -128;
            //   imach[10] = 127;
            //
            //     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
            //       32-BIT INTEGER ARITHMETIC.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -127;
            //   imach[7] = 127;
            //   imach[8] = 56;
            //   imach[9] = -127;
            //   imach[10] = 127;
            //
            //     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -125;
            //   imach[7] = 128;
            //   imach[8] = 53;
            //   imach[9] = -1021;
            //   imach[10] = 1024;
            //
            //     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
            //       SERIES (MIPS R3000 PROCESSOR).
            //
            //   imach[1] = 2;
            //   imach[2] = 31;
            //   imach[3] = 2147483647;
            //   imach[4] = 2;
            //   imach[5] = 24;
            //   imach[6] = -125;
            //   imach[7] = 128;
            //   imach[8] = 53;
            //   imach[9] = -1021;
            //   imach[10] = 1024;
            //
            //     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
            //       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
            //       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

            imach[1] = 2;
            imach[2] = 31;
            imach[3] = 2147483647;
            imach[4] = 2;
            imach[5] = 24;
            imach[6] = -125;
            imach[7] = 128;
            imach[8] = 53;
            imach[9] = -1021;
            imach[10] = 1024;

            ipmpar = imach[i];
            return ipmpar;
        }

        public static double dexpm1(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DEXPM1 evaluates the function EXP(X) - 1.
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
            //    Armido DiDinato and Alfred Morris,
            //    Algorithm 708:
            //    Significant Digit Computation of the Incomplete Beta Function Ratios,
            //    ACM Transactions on Mathematical Software,
            //    Volume 18, 1993, pages 360-373.
            //
            //  Parameters:
            //
            //    Input, double *X, the value at which exp(X)-1 is desired.
            //
            //    Output, double DEXPM1, the value of exp(X)-1.
            //
        {
            double p1 = 0.914041914819518e-09;
            double p2 = 0.238082361044469e-01;
            double q1 = -0.499999999085958e+00;
            double q2 = 0.107141568980644e+00;
            double q3 = -0.119041179760821e-01;
            double q4 = 0.595130811860248e-03;
            double dexpm1;
            double w;

            if (Math.Abs(x) <= 0.15e0)
            {
                dexpm1 = x * (((
                                   p2 * x
                                   + p1) * x
                               + 1.0e0)
                              / ((((
                                       q4 * x
                                       + q3) * x
                                   + q2) * x
                                  + q1) * x
                                 + 1.0e0));
            }
            else if (x <= 0.0e0)
            {
                w = Math.Exp(x);
                dexpm1 = w - 0.5e0 - 0.5e0;
            }
            else
            {
                w = Math.Exp(x);
                dexpm1 = w * (0.5e0 + (0.5e0 - 1.0e0 / w));
            }

            return dexpm1;
        }
        //****************************************************************************80

        public static double dinvnr(double p, double q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DINVNR computes the inverse of the normal distribution.
            //
            //  Discussion:
            //
            //    Returns X such that CUMNOR(X)  =   P,  i.e., the  integral from -
            //    infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
            //
            //    The rational function on page 95 of Kennedy and Gentle is used as a start
            //    value for the Newton method of finding roots.
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
            //    Kennedy and Gentle,
            //    Statistical Computing,
            //    Marcel Dekker, NY, 1980,
            //    QA276.4  K46
            //
            //  Parameters:
            //
            //    Input, double *P, *Q, the probability, and the complementary
            //    probability.
            //
            //    Output, double DINVNR, the argument X for which the
            //    Normal CDF has the value P.
            //
        {
            double maxit = 100;
            double eps = (1.0e-13);
            double r2pi = 0.3989422804014326e0;
            double nhalf = (-0.5e0);

            double dennor(double x)
            {
                return (r2pi * Math.Exp(nhalf * (x) * (x)));
            }

            double ccum = 0;
            double cum = 0;
            double dinvnr;
            double dx;
            int i;
            double pp;
            bool qporq = false;
            double strtx;
            double xcur;
            //
            //  FIND MINIMUM OF P AND Q
            //
            qporq = p <= q;
            if (!qporq) goto S10;
            pp = p;
            goto S20;
            S10:
            pp = q;
            S20:
            //
            //  INITIALIZATION STEP
            //
            strtx = stvaln(pp);
            xcur = strtx;
            //
            //  NEWTON INTERATIONS
            //
            for (i = 1; i <= maxit; i++)
            {
                cumnor(xcur, ref cum, ref ccum);
                dx = (cum - pp) / dennor(xcur);
                xcur = xcur - dx;
                if (Math.Abs(dx / xcur) < eps) goto S40;
            }

            dinvnr = strtx;
            //
            //  IF WE GET HERE, NEWTON HAS FAILED
            //
            if (!qporq) dinvnr = -dinvnr;
            return dinvnr;
            S40:
            //
            //  IF WE GET HERE, NEWTON HAS SUCCEDED
            //
            dinvnr = xcur;
            if (!qporq) dinvnr = -dinvnr;
            return dinvnr;
        }

        public static void dinvr(ref E0000Data data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DINVR bounds the zero of the function and invokes DZROR.
            //
            //  Discussion:
            //
            //    This routine seeks to find bounds on a root of the function and
            //    invokes ZROR to perform the zero finding.  STINVR must have been
            //    called before this routine in order to set its parameters.
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
            //    J C P Bus and T J Dekker,
            //    Two Efficient Algorithms with Guaranteed Convergence for
            //    Finding a Zero of a Function,
            //    ACM Transactions on Mathematical Software,
            //    Volume 1, Number 4, pages 330-345, 1975.
            //
            //  Parameters:
            //
            //    Input/output, integer STATUS.  At the beginning of a zero finding
            //    problem, STATUS should be set to 0 and INVR invoked.  The value
            //    of parameters other than X will be ignored on this call.
            //    If INVR needs the function to be evaluated, it will set STATUS to 1
            //    and return.  The value of the function should be set in FX and INVR
            //    again called without changing any of its other parameters.
            //    If INVR finishes without error, it returns with STATUS 0, and X an
            //    approximate root of F(X).
            //    If INVR cannot bound the function, it returns a negative STATUS and
            //    sets QLEFT and QHI.
            //
            //    Output, double precision X, the value at which F(X) is to be evaluated.
            //
            //    Input, double precision FX, the value of F(X) calculated by the user
            //    on the previous call, when INVR returned with STATUS = 1.
            //
            //    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
            //    case, QLEFT is TRUE if the stepping search terminated unsucessfully
            //    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
            //
            //    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
            //    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
            //    if F(X) < Y.
            //
        {
            E0000 ( 0, ref data );
        }

        public class E0000Data
        {
            public int ientry { get; set; }
            public int status { get; set; }
            public double x { get; set; }
            public double fx { get; set; }
            public bool qleft  { get; set; }
            public bool qhi  { get; set; }
            public double zabsst  { get; set; }
            public double zabsto  { get; set; }
            public double zbig  { get; set; }
            public double zsmall  { get; set; }
            public double zrelst { get; set; }
            public double zrelto { get; set; }
            public double zstpmu  { get; set; }
        }
        
        public static double dlanor(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DLANOR evaluates the logarithm of the asymptotic Normal CDF.
            //
            //  Discussion:
            //
            //    This routine computes the logarithm of the cumulative normal distribution
            //    from abs ( x ) to infinity for  5 <= abs ( X ).
            //
            //    The relative error at X = 5 is about 0.5E-5.
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
            //    Milton Abramowitz and Irene Stegun,
            //    Handbook of Mathematical Functions
            //    1966, Formula 26.2.12.
            //
            //  Parameters:
            //
            //    Input, double *X, the value at which the Normal CDF is to be
            //    evaluated.  It is assumed that 5 <= abs ( X ).
            //
            //    Output, double DLANOR, the logarithm of the asymptotic
            //    Normal CDF.
            //
        {
            double approx;
            double[] coef =
            {
                -1.0e0, 3.0e0, -15.0e0, 105.0e0, -945.0e0, 10395.0e0, -135135.0e0, 2027025.0e0,
                -34459425.0e0, 654729075.0e0, -13749310575e0, 316234143225.0e0
            };
            double correc;
            double dlanor;
            double dlsqpi = 0.91893853320467274177e0;
            int K1 = 12;
            double T2;
            double xx;
            double xx2;

            xx = Math.Abs(x);
            if (xx < 5.0e0)
            {
                throw new Exception("DLANOR: Argument too small.");
            }

            approx = -dlsqpi - 0.5e0 * xx * xx - Math.Log(xx);
            xx2 = xx * xx;
            T2 = 1.0e0 / xx2;
            correc = eval_pol(coef, K1, T2) / xx2;
            correc = alnrel(correc);
            dlanor = approx + correc;

            return dlanor;
        }

        public static void dstinv(ref E0000Data edata)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSTINV seeks a value X such that F(X) = Y.
            //
            //  Discussion:
            //
            //    DSTINV is the double precision set inverse finder.
            //    It uses reverse communication.
            //
            //    Given a monotone function F and a value Y, it finds X
            //    such that F(X) = Y.
            //
            //    This routine sets quantities needed by INVR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  More Precise Description of INVR -
            //
            //     F must be a monotone function, the results of QMFINV are
            //     otherwise undefined.  QINCR must be .TRUE. if F is non-
            //     decreasing and .FALSE. if F is non-increasing.
            //     QMFINV will return .TRUE. if and only if F(SMALL) and
            //     F(BIG) bracket Y, i. e.,
            //          QINCR is .TRUE. and F(SMALL)<=Y<=F(BIG) or
            //          QINCR is .FALSE. and F(BIG)<=Y<=F(SMALL)
            //     if QMFINV returns .TRUE., then the X returned satisfies
            //     the following condition.  let
            //               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
            //     then if QINCR is .TRUE.,
            //          F(X-TOL(X)) <= Y <= F(X+TOL(X))
            //     and if QINCR is .FALSE.
            //          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
            //
            //                              Method
            //     Compares F(X) with Y for the input value of X then uses QINCR
            //     to determine whether to step left or right to bound the
            //     desired x.  the initial step size is
            //          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
            //     Iteratively steps right or left until it bounds X.
            //     At each step which doesn't bound X, the step size is doubled.
            //     The routine is careful never to step beyond SMALL or BIG.  If
            //     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
            //     after setting QLEFT and QHI.
            //     If X is successfully bounded then Algorithm R of the paper
            //     'Two Efficient Algorithms with Guaranteed Convergence for
            //     Finding a Zero of a Function' by J. C. P. Bus and
            //     T. J. Dekker in ACM Transactions on Mathematical
            //     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
            //     to find the zero of the function F(X)-Y. This is routine
            //     QRZERO.
            //
            //  Parameters:
            //
            //    double SMALL --> The left endpoint of the interval to be
            //          searched for a solution.
            //
            //    double BIG --> The right endpoint of the interval to be
            //          searched for a solution.
            //
            //    double ABSSTP, RELSTP --> The initial step size in the search
            //          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
            //
            //    double STPMUL --> When a step doesn't bound the zero, the step
            //                size is multiplied by STPMUL and another step
            //                taken.  A popular value is 2.0
            //
            //    double ABSTOL, RELTOL --> Two numbers that determine the accuracy
            //          of the solution.  See function for a precise definition.
            //
        {
            E0000(1,ref edata);
        }

        public static double dstrem(double z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
            //
            //  Discussion:
            //
            //    This routine returns
            //
            //      ln ( Gamma ( Z ) ) - Sterling ( Z )
            //
            //    where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
            //
            //    Sterling(Z) = ln ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
            //
            //    If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
            //    with values calculated using Maple.
            //
            //    Otherwise, the difference is computed explicitly.
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
            //  Parameters:
            //
            //    Input, double *Z, the value at which the Sterling
            //    remainder is to be calculated.  Z must be positive.
            //
            //    Output, double DSTREM, the Sterling remainder.
            //
        {
            double hln2pi = 0.91893853320467274178e0;
            int ncoef = 10;

            double[] coef =
                {
                    0.0e0, 0.0833333333333333333333333333333e0,
                    -0.00277777777777777777777777777778e0, 0.000793650793650793650793650793651e0,
                    -0.000595238095238095238095238095238e0,
                    0.000841750841750841750841750841751e0, -0.00191752691752691752691752691753e0,
                    0.00641025641025641025641025641026e0, -0.0295506535947712418300653594771e0,
                    0.179644372368830573164938490016e0
                }
                ;
            int K1 = 10;
            double dstrem;
            double sterl;
            double T2;
            //
            //    For information, here are the next 11 coefficients of the
            //    remainder term in Sterling's formula
            //            -1.39243221690590111642743221691
            //            13.4028640441683919944789510007
            //            -156.848284626002017306365132452
            //            2193.10333333333333333333333333
            //            -36108.7712537249893571732652192
            //            691472.268851313067108395250776
            //            -0.152382215394074161922833649589D8
            //            0.382900751391414141414141414141D9
            //            -0.108822660357843910890151491655D11
            //            0.347320283765002252252252252252D12
            //            -0.123696021422692744542517103493D14
            //
            if (z <= 0.0e0)
            {
                throw new Exception("Zero or negative argument in DSTREM");
            }

            if (!(z > 6.0e0)) goto S10;
            T2 = 1.0e0 / Math.Pow(z, 2.0);
            dstrem = eval_pol(coef, K1, T2) * z;
            goto S20;
            S10:
            sterl = hln2pi + (z - 0.5e0) * Math.Log(z) - z;
            dstrem = gamma_log(z) - sterl;
            S20:
            return dstrem;
        }

        public static void dstzr(ref E0001Data edata)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSTXR sets quantities needed by the zero finder.
            //
            //  Discussion:
            //
            //     Double precision SeT ZeRo finder - Reverse communication version
            //                              Function
            //     Sets quantities needed by ZROR.  The function of ZROR
            //     and the quantities set is given here.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Concise Description
            //
            //    Given a function F
            //     find XLO such that F(XLO) = 0.
            //          More Precise Description -
            //     Input condition. F is a double function of a single
            //     double argument and XLO and XHI are such that
            //          F(XLO)*F(XHI)  <=  0.0
            //     If the input condition is met, QRZERO returns .TRUE.
            //     and output values of XLO and XHI satisfy the following
            //          F(XLO)*F(XHI)  <= 0.
            //          ABS(F(XLO)  <= ABS(F(XHI)
            //          ABS(XLO-XHI)  <= TOL(X)
            //     where
            //          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
            //     If this algorithm does not find XLO and XHI satisfying
            //     these conditions then QRZERO returns .FALSE.  This
            //     implies that the input condition was not met.
            //
            //  Parameters:
            //
            //     XLO --> The left endpoint of the interval to be
            //           searched for a solution.
            //                    XLO is DOUBLE PRECISION
            //     XHI --> The right endpoint of the interval to be
            //           for a solution.
            //                    XHI is DOUBLE PRECISION
            //     ABSTOL, RELTOL --> Two numbers that determine the accuracy
            //                      of the solution.  See function for a
            //                      precise definition.
            //                    ABSTOL is DOUBLE PRECISION
            //                    RELTOL is DOUBLE PRECISION
            //
            //                              Method
            //     Algorithm R of the paper 'Two Efficient Algorithms with
            //     Guaranteed Convergence for Finding a Zero of a Function'
            //     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
            //     Mathematical Software, Volume 1, no. 4 page 330
            //     (Dec. '75) is employed to find the zero of F(X)-Y.
            //
        {
            E0001(1,ref edata);
        }

        public class E0001Data
        {
            public int IENTRY { get; set; }
            public int status  { get; set; }
            public double x { get; set; }
            public double fx  { get; set; }
            public double xlo { get; set; }
            public double xhi { get; set; }
            public bool qleft  { get; set; }
            public bool qhi  { get; set; }
            public double zabstl { get; set; }
            public double zreltl  { get; set; }
            public double zxhi  { get; set; }
            public double zxlo  { get; set; }
        }

        public static void dzror(ref E0000Data edata, ref E0001Data e1data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DZROR seeks the zero of a function using reverse communication.
            //
            //  Discussion:
            //
            //    Performs the zero finding.  STZROR must have been called before
            //    this routine in order to set its parameters.
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
            //  Parameters:
            //
            //    int STATUS <--> At the beginning of a zero finding problem, STATUS
            //                 should be set to 0 and ZROR invoked.  (The value
            //                 of other parameters will be ignored on this call.)
            //
            //                 When ZROR needs the function evaluated, it will set
            //                 STATUS to 1 and return.  The value of the function
            //                 should be set in FX and ZROR again called without
            //                 changing any of its other parameters.
            //
            //                 When ZROR has finished without error, it will return
            //                 with STATUS 0.  In that case (XLO,XHI) bound the answe
            //
            //                 If ZROR finds an error (which implies that F(XLO)-Y an
            //                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
            //                 this case, XLO and XHI are undefined.
            //
            //    double X <-- The value of X at which F(X) is to be evaluated.
            //
            //    double FX --> The value of F(X) calculated when ZROR returns with
            //            STATUS = 1.
            //
            //    double XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
            //             inverval in X containing the solution below.
            //
            //    double XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
            //             inverval in X containing the solution above.
            //
            //    bool QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
            //                at XLO.  If it is .FALSE. the search terminated
            //                unsucessfully at XHI.
            //
            //    bool QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
            //              search and .FALSE. if F(X) < Y at the
            //              termination of the search.
            //
        {
            e1data.status = edata.status;
            e1data.qleft = edata.qleft;
            e1data.qhi = edata.qhi;
            
            e1data.x = edata.x;
            e1data.fx = edata.fx;
            
            E0001(0, ref e1data);
        }

        public static double rcomp(double a, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RCOMP evaluates exp(-X) * X**A / Gamma(A).
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
            //  Input:
            //
            //    double *A, *X, arguments of the quantity to be computed.
            //
            //  Output:
            //
            //    double RCOMP, the value of exp(-X) * X^A / Gamma(A).
            //
            //  Local:
            //
            //    double RT2PIN = 1/SQRT(2*PI)
            //
        {
            double rt2pin = .398942280401433e0;
            double rcomp, t, t1, u;
            rcomp = 0.0e0;
            if (a >= 20.0e0) goto S20;
            t = a * Math.Log(x) - x;
            if (a >= 1.0e0) goto S10;
            rcomp = a * Math.Exp(t) * (1.0e0 + gam1(a));
            return rcomp;
            S10:
            rcomp = Math.Exp(t) / gamma_x(a);
            return rcomp;
            S20:
            u = x / a;
            if (u == 0.0e0) return rcomp;
            t = Math.Pow(1.0e0 / a, 2.0);
            t1 = (((0.75e0 * t - 1.0e0) * t + 3.5e0) * t - 105.0e0) / (a * 1260.0e0);
            t1 = t1 - (a * rlog(u));
            rcomp = rt2pin * Math.Sqrt(a) * Math.Exp(t1);
            return rcomp;
        }

        public static double rexp(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    REXP evaluates the function EXP(X) - 1.
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
            //  Input:
            //
            //    double *X, the argument of the function.
            //
            //  Output:
            //
            //    double REXP, the value of EXP(X)-1.
            //
        {
            double p1 = .914041914819518e-09;
            double p2 = .238082361044469e-01;
            double q1 = -.499999999085958e+00;
            double q2 = .107141568980644e+00;
            double q3 = -.119041179760821e-01;
            double q4 = .595130811860248e-03;
            double rexp, w;

            if (Math.Abs(x) > 0.15e0) goto S10;
            rexp = x * (((p2 * x + p1) * x + 1.0e0) / ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.0e0));
            return rexp;
            S10:
            w = Math.Exp(x);
            if (x > 0.0e0) goto S20;
            rexp = w - 0.5e0 - 0.5e0;
            return rexp;
            S20:
            rexp = w * (0.5e0 + (0.5e0 - 1.0e0 / w));
            return rexp;
        }

        public static double rlog(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RLOG computes  X - 1 - LN(X).
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
            //  Parameters:
            //
            //    Input, double *X, the argument of the function.
            //
            //    Output, double RLOG, the value of the function.
            //
        {
            double a = .566749439387324e-01;
            double b = .456512608815524e-01;
            double p0 = .333333333333333e+00;
            double p1 = -.224696413112536e+00;
            double p2 = .620886815375787e-02;
            double q1 = -.127408923933623e+01;
            double q2 = .354508718369557e+00;
            double rlog, r, t, u, w, w1;

            if (x < 0.61e0 || x > 1.57e0) goto S40;
            if (x < 0.82e0) goto S10;
            if (x > 1.18e0) goto S20;
            //
            //  ARGUMENT REDUCTION
            //
            u = x - 0.5e0 - 0.5e0;
            w1 = 0.0e0;
            goto S30;
            S10:
            u = x - 0.7e0;
            u /= 0.7e0;
            w1 = a - u * 0.3e0;
            goto S30;
            S20:
            u = 0.75e0 * x - 1e0;
            w1 = b + u / 3.0e0;
            S30:
            //
            //  SERIES EXPANSION
            //
            r = u / (u + 2.0e0);
            t = r * r;
            w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0e0);
            rlog = 2.0e0 * t * (1.0e0 / (1.0e0 - r) - r * w) + w1;
            return rlog;
            S40:
            r = x - 0.5e0 - 0.5e0;
            rlog = r - Math.Log(x);
            return rlog;
        }

        public static double rlog1(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RLOG1 evaluates the function X - ln ( 1 + X ).
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
            //  Parameters:
            //
            //    Input, double *X, the argument.
            //
            //    Output, double RLOG1, the value of X - ln ( 1 + X ).
            //
        {
            double a = .566749439387324e-01;
            double b = .456512608815524e-01;
            double h;
            double p0 = .333333333333333e+00;
            double p1 = -.224696413112536e+00;
            double p2 = .620886815375787e-02;
            double q1 = -.127408923933623e+01;
            double q2 = .354508718369557e+00;
            double rlog1, r, t, w, w1;

            if (x < -0.39e0 || x > 0.57e0) goto S40;
            if (x < -0.18e0) goto S10;
            if (x > 0.18e0) goto S20;
            //
            //  ARGUMENT REDUCTION
            //
            h = x;
            w1 = 0.0e0;
            goto S30;
            S10:
            h = x + 0.3e0;
            h /= 0.7e0;
            w1 = a - h * 0.3e0;
            goto S30;
            S20:
            h = 0.75e0 * x - 0.25e0;
            w1 = b + h / 3.0e0;
            S30:
            //
            //  SERIES EXPANSION
            //
            r = h / (h + 2.0e0);
            t = r * r;
            w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0e0);
            rlog1 = 2.0e0 * t * (1.0e0 / (1.0e0 - r) - r * w) + w1;
            return rlog1;
            S40:
            w = x + 0.5e0 + 0.5e0;
            rlog1 = x - Math.Log(w);
            return rlog1;
        }
        
        public static double stvaln ( double p )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STVALN provides starting values for the inverse of the normal distribution.
            //
            //  Discussion:
            //
            //    The routine returns X such that
            //      P = CUMNOR(X),
            //    that is,
            //      P = Integral from -infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU.
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
            //    Kennedy and Gentle,
            //    Statistical Computing,
            //    Marcel Dekker, NY, 1980, page 95,
            //    QA276.4  K46
            //
            //  Parameters:
            //
            //    Input, double *P, the probability whose normal deviate
            //    is sought.
            //
            //    Output, double STVALN, the normal deviate whose probability
            //    is P.
            //
        {
            double[] xden = {
                0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
                0.38560700634e-2
            };
            double[] xnum = {
                -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
                -0.453642210148e-4
            };
            int K1 = 5;
            double sign;
            double stvaln,y,z;

            if(!(p <= 0.5e0)) goto S10;
            sign = -1.0e0;
            z = p;
            goto S20;
            S10:
            sign = 1.0e0;
            z = 1.0e0-p;
            S20:
            y = Math.Sqrt(-(2.0e0*Math.Log(z)));
            stvaln = y+ eval_pol ( xnum, K1, y ) / eval_pol ( xden, K1, y );
            stvaln = sign*stvaln;
            return stvaln;
        }


    }
}