using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
       public static double beta(double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA evaluates the beta function.
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
            //  Input:
            //
            //    double A, B, the arguments of the beta function.
            //
            //  Output:
            //
            //    double BETA, the value of the beta function.
            //
        {
            double value;

            value = Math.Exp(beta_log(a, b));

            return value;
        }

        public static double beta_asym(double a, double b, double lambda, double eps)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
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
            //    Input, double *A, *B, the parameters of the function.
            //    A and B should be nonnegative.  It is assumed that both A and B
            //    are greater than or equal to 15.
            //
            //    Input, double *LAMBDA, the value of ( A + B ) * Y - B.
            //    It is assumed that 0 <= LAMBDA.
            //
            //    Input, double *EPS, the tolerance.
            //
            //  Local:
            //
            //    int NUM,  THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
            //    ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
            //    THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
            //
            //    E0 = 2/SQRT(PI)
            //
            //    E1 = 2^(-3/2)
            //
        {
            double[] a0 = new double[21];
            double[] b0 = new double[21];
            double bsum;
            double[] c = new double[21];
            double[] d = new double[21];
            double dsum;
            double e0 = 1.12837916709551e0;
            double e1 = 0.353553390593274e0;
            double f;
            double h;
            double h2;
            double hn;
            double j0;
            double j1;
            int K3 = 1;
            int num = 20;
            double r;
            double r0;
            double r1;
            double value;

            double s,
                sum,
                t,
                t0,
                t1,
                u,
                w,
                w0,
                z,
                z0,
                z2,
                zn,
                znm1;
            int i, im1, imj, j, m, mm1, mmj, n, np1;
            double T1, T2;

            value = 0.0e0;
            if (a >= b) goto S10;
            h = a / b;
            r0 = 1.0e0 / (1.0e0 + h);
            r1 = (b - a) / b;
            w0 = 1.0e0 / Math.Sqrt(a * (1.0e0 + h));
            goto S20;
            S10:
            h = b / a;
            r0 = 1.0e0 / (1.0e0 + h);
            r1 = (b - a) / a;
            w0 = 1.0e0 / Math.Sqrt(b * (1.0e0 + h));
            S20:
            T1 = -(lambda / a);
            T2 = lambda / b;
            f = a * rlog1(T1) + b * rlog1(T2);
            t = Math.Exp(-f);
            if (t == 0.0e0) return value;
            z0 = Math.Sqrt(f);
            z = 0.5e0 * (z0 / e1);
            z2 = f + f;
            a0[0] = 2.0e0 / 3.0e0 * r1;
            c[0] = -(0.5e0 * a0[0]);
            d[0] = -c[0];
            j0 = 0.5e0 / e0 * error_fc(K3, z0);
            j1 = e1;
            sum = j0 + d[0] * w0 * j1;
            s = 1.0e0;
            h2 = h * h;
            hn = 1.0e0;
            w = w0;
            znm1 = z;
            zn = z2;
            for (n = 2; n <= num; n += 2)
            {
                hn = h2 * hn;
                a0[n - 1] = 2.0e0 * r0 * (1.0e0 + h * hn) / ((double) n + 2.0e0);
                np1 = n + 1;
                s = s + hn;
                a0[np1 - 1] = 2.0e0 * r1 * s / ((double) n + 3.0e0);
                for (i = n; i <= np1; i++)
                {
                    r = -(0.5e0 * ((double) i + 1.0e0));
                    b0[0] = r * a0[0];
                    for (m = 2; m <= i; m++)
                    {
                        bsum = 0.0e0;
                        mm1 = m - 1;
                        for (j = 1; j <= mm1; j++)
                        {
                            mmj = m - j;
                            bsum = bsum + (((double) j * r - (double) mmj) * a0[j - 1] * b0[mmj - 1]);
                        }

                        b0[m - 1] = r * a0[m - 1] + bsum / (double) m;
                    }

                    c[i - 1] = b0[i - 1] / ((double) i + 1.0e0);
                    dsum = 0.0e0;
                    im1 = i - 1;
                    for (j = 1; j <= im1; j++)
                    {
                        imj = i - j;
                        dsum = dsum + (d[imj - 1] * c[j - 1]);
                    }

                    d[i - 1] = -(dsum + c[i - 1]);
                }

                j0 = e1 * znm1 + ((double) n - 1.0e0) * j0;
                j1 = e1 * zn + (double) n * j1;
                znm1 = z2 * znm1;
                zn = z2 * zn;
                w = w0 * w;
                t0 = d[n - 1] * w * j0;
                w = w0 * w;
                t1 = d[np1 - 1] * w * j1;
                sum = sum + (t0 + t1);
                if (Math.Abs(t0) + Math.Abs(t1) <= eps * sum) goto S80;
            }

            S80:
            u = Math.Exp(-bcorr(a, b));
            value = e0 * t * u * sum;
            return value;
        }

        public static double beta_frac(double a, double b, double x, double y, double lambda,
                double eps)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
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
            //    double *A, *B, the parameters of the function.
            //    A and B should be nonnegative.  It is assumed that both A and
            //    B are greater than 1.
            //
            //    double *X, *Y.  X is the argument of the
            //    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
            //
            //    double *LAMBDA, the value of ( A + B ) * Y - B.
            //
            //    double *EPS, a tolerance.
            //
            //  Output:
            //
            //    double BETA_FRAC, the value of the continued
            //    fraction approximation for IX(A,B).
            //
        {
            double alpha;
            double an;
            double anp1;
            double beta;
            double bfrac;
            double bn;
            double bnp1;
            double c;
            double c0;
            double c1;
            double e;
            double n;
            double p;
            double r;
            double r0;
            double s;
            double t;
            double w;
            double yp1;

            bfrac = beta_rcomp(a, b, x, y);

            if (bfrac == 0.0e0)
            {
                return bfrac;
            }

            c = 1.0e0 + lambda;
            c0 = b / a;
            c1 = 1.0e0 + 1.0e0 / a;
            yp1 = y + 1.0e0;
            n = 0.0e0;
            p = 1.0e0;
            s = a + 1.0e0;
            an = 0.0e0;
            bn = anp1 = 1.0e0;
            bnp1 = c / c1;
            r = c1 / c;
            //
            //  CONTINUED FRACTION CALCULATION
            //
            S10:
            n = n + 1.0e0;
            t = n / a;
            w = n * (b - n) * x;
            e = a / s;
            alpha = p * (p + c0) * e * e * (w * x);
            e = (1.0e0 + t) / (c1 + t + t);
            beta = n + w / s + e * (c + n * yp1);
            p = 1.0e0 + t;
            s = s + 2.0e0;
            //
            //  UPDATE AN, BN, ANP1, AND BNP1
            //
            t = alpha * an + beta * anp1;
            an = anp1;
            anp1 = t;
            t = alpha * bn + beta * bnp1;
            bn = bnp1;
            bnp1 = t;
            r0 = r;
            r = anp1 / bnp1;

            if (Math.Abs(r - r0) <= (eps) * r)
            {
                goto S20;
            }

            //
            //  RESCALE AN, BN, ANP1, AND BNP1
            //
            an /= bnp1;
            bn /= bnp1;
            anp1 = r;
            bnp1 = 1.0e0;
            goto S10;
            //
            //  TERMINATION
            //
            S20:
            bfrac = bfrac * r;
            return bfrac;
        }

        public static void beta_grat(double a, double b, double x, double y, ref double w,
                double eps, ref int ierr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
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
            //    Input, double *A, *B, the parameters of the function.
            //    A and B should be nonnegative.  It is assumed that 15 <= A
            //    and B <= 1, and that B is less than A.
            //
            //    Input, double *X, *Y.  X is the argument of the
            //    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
            //
            //    Input/output, double *W, a quantity to which the
            //    result of the computation is to be added on output.
            //
            //    Input, double *EPS, a tolerance.
            //
            //    Output, int *IERR, an error flag, which is 0 if no error
            //    was detected.
            //
        {
            double bm1;
            double bp2n;
            double[] c = new double[30];
            double cn;
            double coef;
            double[] d = new double[30];
            double dj;
            int i;
            double j;
            double l;
            double lnx;
            int n;
            double n2;
            int nm1;
            double nu;
            double p;
            double q;
            double r;
            double s;
            double sum;
            double t;
            double T1;
            double t2;
            double u;
            double v;
            double z;

            bm1 = b - 0.5e0 - 0.5e0;
            nu = a + 0.5e0 * bm1;
            if (y > 0.375e0) goto S10;
            T1 = -y;
            lnx = alnrel(T1);
            goto S20;
            S10:
            lnx = Math.Log(x);
            S20:
            z = -(nu * lnx);
            if (b * z == 0.0e0) goto S70;
            //
            //  COMPUTATION OF THE EXPANSION
            //  SET R = EXP(-Z)*Z^B/GAMMA(B)
            //
            r = b * (1.0e0 + gam1(b)) * Math.Exp(b * Math.Log(z));
            r = r * (Math.Exp(a * lnx) * Math.Exp(0.5e0 * bm1 * lnx));
            u = algdiv(b, a) + b * Math.Log(nu);
            u = r * Math.Exp(-u);
            if (u == 0.0e0) goto S70;
            gamma_rat1(b, z, r, p, q, eps);
            v = 0.25e0 * Math.Pow(1.0e0 / nu, 2.0);
            t2 = 0.25e0 * lnx * lnx;
            l = w / u;
            j = q / r;
            sum = j;
            t = cn = 1.0e0;
            n2 = 0.0e0;
            for (n = 1; n <= 30; n++)
            {
                bp2n = b + n2;
                j = (bp2n * (bp2n + 1.0e0) * j + (z + bp2n + 1.0e0) * t) * v;
                n2 = n2 + 2.0e0;
                t = t * t2;
                cn /= (n2 * (n2 + 1.0e0));
                c[n - 1] = cn;
                s = 0.0e0;
                if (n == 1) goto S40;
                nm1 = n - 1;
                coef = b - (double) n;
                for (i = 1; i <= nm1; i++)
                {
                    s = s + (coef * c[i - 1] * d[n - i - 1]);
                    coef = coef + b;
                }

                S40:
                d[n - 1] = bm1 * cn + s / (double) n;
                dj = d[n - 1] * j;
                sum = sum + dj;
                if (sum <= 0.0e0) goto S70;
                if (Math.Abs(dj) <= eps * (sum + l)) goto S60;
            }

            S60:
            //
            //  ADD THE RESULTS TO W
            //
            ierr = 0;
            w = w + (u * sum);
            return;
            S70:
            //
            //  THE EXPANSION CANNOT BE COMPUTED
            //
            ierr = 1;
        }

        public static void beta_inc(double a, double b, double x, double y, ref double w,
                ref double w1, ref int ierr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_INC evaluates the incomplete beta function IX(A,B).
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
            //    Alfred H Morris, Jr.
            //
            //  Parameters:
            //
            //    Input, double *A, *B, the parameters of the function.
            //    A and B should be nonnegative.
            //
            //    Input, double *X, *Y.  X is the argument of the
            //    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
            //
            //    Output, double *W, *W1, the values of IX(A,B) and
            //    1-IX(A,B).
            //
            //    Output, int *IERR, the error flag.
            //    0, no error was detected.
            //    1, A or B is negative;
            //    2, A = B = 0;
            //    3, X < 0 or 1 < X;
            //    4, Y < 0 or 1 < Y;
            //    5, X + Y /= 1;
            //    6, X = A = 0;
            //    7, Y = B = 0.
            //
        {
            double a0;
            double b0;
            int K1 = 1;
            double eps;
            double lambda;
            double t;
            double x0;
            double y0;
            double z;

            int ierr1, ind, n;
            double T2, T3, T4, T5;
            //
            //  EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
            //  NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
            //
            eps = dpmpar(K1);
            w = w1 = 0.0e0;
            if (a < 0.0e0 || b < 0.0e0) goto S270;
            if (a == 0.0e0 && b == 0.0e0) goto S280;
            if (x < 0.0e0 || x > 1.0e0) goto S290;
            if (y < 0.0e0 || y > 1.0e0) goto S300;
            z = x + y - 0.5e0 - 0.5e0;
            if (Math.Abs(z) > 3.0e0 * eps) goto S310;
            ierr = 0;
            if (x == 0.0e0) goto S210;
            if (y == 0.0e0) goto S230;
            if (a == 0.0e0) goto S240;
            if (b == 0.0e0) goto S220;
            eps = Math.Max(eps, 1e-15);
            if (Math.Max(a, b) < 1e-3 * eps) goto S260;
            ind = 0;
            a0 = a;
            b0 = b;
            x0 = x;
            y0 = y;
            if (Math.Min(a0, b0) > 1.0e0) goto S40;
            //
            //  PROCEDURE FOR A0 <= 1 OR B0 <= 1
            //
            if (x <= 0.5e0) goto S10;
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
            S10:
            if (b0 < Math.Min(eps, eps * a0)) goto S90;
            if (a0 < Math.Min(eps, eps * b0) && b0 * x0 <= 1.0e0) goto S100;
            if (Math.Max(a0, b0) > 1.0e0) goto S20;
            if (a0 >= Math.Min(0.2e0, b0)) goto S110;
            if (Math.Pow(x0, a0) <= 0.9e0) goto S110;
            if (x0 >= 0.3e0) goto S120;
            n = 20;
            goto S140;
            S20:
            if (b0 <= 1.0e0) goto S110;
            if (x0 >= 0.3e0) goto S120;
            if (x0 >= 0.1e0) goto S30;
            if (Math.Pow(x0 * b0, a0) <= 0.7e0) goto S110;
            S30:
            if (b0 > 15.0e0) goto S150;
            n = 20;
            goto S140;
            S40:
            //
            //  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
            //
            if (a > b) goto S50;
            lambda = a - (a + b) * x;
            goto S60;
            S50:
            lambda = (a + b) * y - b;
            S60:
            if (lambda >= 0.0e0) goto S70;
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
            lambda = Math.Abs(lambda);
            S70:
            if (b0 < 40.0e0 && b0 * x0 <= 0.7e0) goto S110;
            if (b0 < 40.0e0) goto S160;
            if (a0 > b0) goto S80;
            if (a0 <= 100.0e0) goto S130;
            if (lambda > 0.03e0 * a0) goto S130;
            goto S200;
            S80:
            if (b0 <= 100.0e0) goto S130;
            if (lambda > 0.03e0 * b0) goto S130;
            goto S200;
            S90:
            //
            //  EVALUATION OF THE APPROPRIATE ALGORITHM
            //
            w = fpser(a0, b0, x0, eps);
            w1 = 0.5e0 + (0.5e0 - w);
            goto S250;
            S100:
            w1 = apser(a0, b0, x0, eps);
            w = 0.5e0 + (0.5e0 - w1);
            goto S250;
            S110:
            w = beta_pser(a0, b0, x0, eps);
            w1 = 0.5e0 + (0.5e0 - w);
            goto S250;
            S120:
            w1 = beta_pser(b0, a0, y0, eps);
            w = 0.5e0 + (0.5e0 - w1);
            goto S250;
            S130:
            T2 = 15.0e0 * eps;
            w = beta_frac(a0, b0, x0, y0, lambda, T2);
            w1 = 0.5e0 + (0.5e0 - w);
            goto S250;
            S140:
            w1 = beta_up(b0, a0, y0, x0, n, eps);
            b0 = b0 + (double) n;
            S150:
            T3 = 15.0e0 * eps;
            beta_grat(b0, a0, y0, x0, w1, T3, ierr1);
            w = 0.5e0 + (0.5e0 - w1);
            goto S250;
            S160:
            n = (int) b0;
            b0 = b0 - (double) n;
            if (b0 != 0.0e0) goto S170;
            n = n - 1;
            b0 = 1.0e0;
            S170:
            w = beta_up(b0, a0, y0, x0, n, eps);
            if (x0 > 0.7e0) goto S180;
            w = w + beta_pser(a0, b0, x0, eps);
            w1 = 0.5e0 + (0.5e0 - w);
            goto S250;
            S180:
            if (a0 > 15.0e0) goto S190;
            n = 20;
            w = w + beta_up(a0, b0, x0, y0, n, eps);
            a0 = a0 + (double) n;
            S190:
            T4 = 15.0e0 * eps;
            beta_grat(a0, b0, x0, y0, w, T4, ierr1);
            w1 = 0.5e0 + (0.5e0 - w);
            goto S250;
            S200:
            T5 = 100.0e0 * eps;
            w = beta_asym(a0, b0, lambda, T5);
            w1 = 0.5e0 + (0.5e0 - w);
            goto S250;
            S210:
            //
            //  TERMINATION OF THE PROCEDURE
            //
            if (a == 0.0e0) goto S320;
            S220:
            w = 0.0e0;
            w1 = 1.0e0;
            return;
            S230:
            if (b == 0.0e0) goto S330;
            S240:
            w = 1.0e0;
            w1 = 0.0e0;
            return;
            S250:
            if (ind == 0) return;
            t = w;
            w = w1;
            w1 = t;
            return;
            S260:
            //
            //  PROCEDURE FOR A AND B < 1.E-3 * EPS
            //
            w = b / (a + b);
            w1 = a / (a + b);
            return;
            S270:
            //
            //  ERROR RETURN
            //
            ierr = 1;
            return;
            S280:
            ierr = 2;
            return;
            S290:
            ierr = 3;
            return;
            S300:
            ierr = 4;
            return;
            S310:
            ierr = 5;
            return;
            S320:
            ierr = 6;
            return;
            S330:
            ierr = 7;
            return;
        }

        public static void beta_inc_values(ref int n_data, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_INC_VALUES returns some values of the incomplete Beta function.
            //
            //  Discussion:
            //
            //    The incomplete Beta function may be written
            //
            //      BETA_INC(A,B,X) = Integral (0 to X) T^(A-1) * (1-T)^(B-1) dT
            //                      / Integral (0 to 1) T^(A-1) * (1-T)^(B-1) dT
            //
            //    Thus,
            //
            //      BETA_INC(A,B,0.0) = 0.0
            //      BETA_INC(A,B,1.0) = 1.0
            //
            //    Note that in Mathematica, the expressions:
            //
            //      BETA[A,B]   = Integral (0 to 1) T^(A-1) * (1-T)^(B-1) dT
            //      BETA[X,A,B] = Integral (0 to X) T^(A-1) * (1-T)^(B-1) dT
            //
            //    and thus, to evaluate the incomplete Beta function requires:
            //
            //      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
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
            //    Karl Pearson,
            //    Tables of the Incomplete Beta Function,
            //    Cambridge University Press, 1968.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, double *A, *B, the parameters of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 30;

            double[] a_vec =  {
                0.5E+00, 0.5E+00, 0.5E+00, 1.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00,
                2.0E+00, 2.0E+00, 2.0E+00, 2.0E+00,
                2.0E+00, 2.0E+00, 2.0E+00, 2.0E+00,
                2.0E+00, 5.5E+00, 10.0E+00, 10.0E+00,
                10.0E+00, 10.0E+00, 20.0E+00, 20.0E+00,
                20.0E+00, 20.0E+00, 20.0E+00, 30.0E+00,
                30.0E+00, 40.0E+00
            }
            ;
            double[] b_vec =  {
                0.5E+00, 0.5E+00, 0.5E+00, 0.5E+00,
                0.5E+00, 0.5E+00, 0.5E+00, 1.0E+00,
                2.0E+00, 2.0E+00, 2.0E+00, 2.0E+00,
                2.0E+00, 2.0E+00, 2.0E+00, 2.0E+00,
                2.0E+00, 5.0E+00, 0.5E+00, 5.0E+00,
                5.0E+00, 10.0E+00, 5.0E+00, 10.0E+00,
                10.0E+00, 20.0E+00, 20.0E+00, 10.0E+00,
                10.0E+00, 20.0E+00
            }
            ;
            double[] fx_vec =  {
                0.0637686E+00, 0.2048328E+00, 1.0000000E+00, 0.0E+00,
                0.0050126E+00, 0.0513167E+00, 0.2928932E+00, 0.5000000E+00,
                0.028E+00, 0.104E+00, 0.216E+00, 0.352E+00,
                0.500E+00, 0.648E+00, 0.784E+00, 0.896E+00,
                0.972E+00, 0.4361909E+00, 0.1516409E+00, 0.0897827E+00,
                1.0000000E+00, 0.5000000E+00, 0.4598773E+00, 0.2146816E+00,
                0.9507365E+00, 0.5000000E+00, 0.8979414E+00, 0.2241297E+00,
                0.7586405E+00, 0.7001783E+00
            }
            ;
            double[] x_vec =  {
                0.01E+00, 0.10E+00, 1.00E+00, 0.0E+00,
                0.01E+00, 0.10E+00, 0.50E+00, 0.50E+00,
                0.1E+00, 0.2E+00, 0.3E+00, 0.4E+00,
                0.5E+00, 0.6E+00, 0.7E+00, 0.8E+00,
                0.9E+00, 0.50E+00, 0.90E+00, 0.50E+00,
                1.00E+00, 0.50E+00, 0.80E+00, 0.60E+00,
                0.80E+00, 0.50E+00, 0.60E+00, 0.70E+00,
                0.80E+00, 0.70E+00
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0.0E+00;
                b = 0.0E+00;
                x = 0.0E+00;
                fx = 0.0E+00;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static double beta_log(double a0, double b0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_LOG evaluates the logarithm of the beta function.
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
            //    Input, double *A0, *B0, the parameters of the function.
            //    A0 and B0 should be nonnegative.
            //
            //    Output, double *BETA_LOG, the value of the logarithm
            //    of the Beta function.
            //
        {
            double a;
            double b;
            double c;
            double e = .918938533204673e0;
            double h;
            int i;
            int n;
            double T1;
            double u;
            double v;
            double value;
            double w;
            double z;

            a = Math.Min(a0, b0);
            b = Math.Max(a0, b0);
            if (a >= 8.0e0) goto S100;
            if (a >= 1.0e0) goto S20;
            //
            //  PROCEDURE WHEN A < 1
            //
            if (b >= 8.0e0) goto S10;
            T1 = a + b;
            value = gamma_log(a) + (gamma_log(b) - gamma_log(T1));
            return value;
            S10:
            value = gamma_log(a) + algdiv(a, b);
            return value;
            S20:
            //
            //  PROCEDURE WHEN 1 <= A < 8
            //
            if (a > 2.0e0) goto S40;
            if (b > 2.0e0) goto S30;
            value = gamma_log(a) + gamma_log(b) - gsumln(a, b);
            return value;
            S30:
            w = 0.0e0;
            if (b < 8.0e0) goto S60;
            value = gamma_log(a) + algdiv(a, b);
            return value;
            S40:
            //
            //  REDUCTION OF A WHEN B <= 1000
            //
            if (b > 1000.0e0) goto S80;
            n = (int) (a - 1.0e0);
            w = 1.0e0;
            for (i = 1; i <= n; i++)
            {
                a = a - 1.0e0;
                h = a / b;
                w = w * (h / (1.0e0 + h));
            }

            w = Math.Log(w);
            if (b < 8.0e0) goto S60;
            value = w + gamma_log(a) + algdiv(a, b);
            return value;
            S60:
            //
            //  REDUCTION OF B WHEN B < 8
            //
            n = (int) (b - 1.0e0);
            z = 1.0e0;
            for (i = 1; i <= n; i++)
            {
                b = b - 1.0e0;
                z = z * (b / (a + b));
            }

            value = w + Math.Log(z) + (gamma_log(a) + (gamma_log(b) - gsumln(a, b)));
            return value;
            S80:
            //
            //  REDUCTION OF A WHEN 1000 < B.
            //
            n = (int) (a - 1.0e0);
            w = 1.0e0;
            for (i = 1; i <= n; i++)
            {
                a = a - 1.0e0;
                w = w * (a / (1.0e0 + a / b));
            }

            value = Math.Log(w) - (double) n * Math.Log(b) + (gamma_log(a) + algdiv(a, b));
            return value;
            S100:
            //
            //  PROCEDURE WHEN 8 <= A.
            //
            w = bcorr(a, b);
            h = a / b;
            c = h / (1.0e0 + h);
            u = -((a - 0.5e0) * Math.Log(c));
            v = b * alnrel(h);
            if (u <= v) goto S110;
            value = -(0.5e0 * Math.Log(b)) + e + w - v - u;
            return value;
            S110:
            value = -(0.5e0 * Math.Log(b)) + e + w - u - v;
            return value;
        }

        public static double beta_pser(double a, double b, double x, double eps)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
            //
            //  Discussion:
            //
            //    BETA_PSER is used when B <= 1 or B*X <= 0.7.
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
            //    double *A, *B, the parameters.
            //
            //    double *X, the point where the function
            //    is to be evaluated.
            //
            //    double *EPS, the tolerance.
            //
            //  Output:
            //
            //    Output, double BETA_PSER, the approximate value of IX(A,B)(X).
            //
        {
            double a0;
            double apb;
            double b0;
            double bpser;
            double c;
            int i;
            int m;
            double n;
            double sum;
            double t;
            double tol;
            double u;
            double w;
            double z;

            bpser = 0.0e0;
            if (x == 0.0e0)
            {
                return bpser;
            }

            //
            //  COMPUTE THE FACTOR X^A/(A*BETA(A,B))
            //
            a0 = Math.Min(a, b);
            if (a0 < 1.0e0) goto S10;
            z = a * Math.Log(x) - beta_log(a, b);
            bpser = Math.Exp(z) / a;
            goto S100;
            S10:
            b0 = Math.Max(a, b);
            if (b0 >= 8.0e0) goto S90;
            if (b0 > 1.0e0) goto S40;
            //
            //  PROCEDURE FOR A0 < 1 AND B0 <= 1
            //
            bpser = Math.Pow(x, a);
            if (bpser == 0.0e0) return bpser;
            apb = a + b;
            if (apb > 1.0e0) goto S20;
            z = 1.0e0 + gam1(apb);
            goto S30;
            S20:
            u = a + b - 1e0;
            z = (1.0e0 + gam1(u)) / apb;
            S30:
            c = (1.0e0 + gam1(a)) * (1.0e0 + gam1(b)) / z;
            bpser = bpser * (c * (b / apb));
            goto S100;
            S40:
            //
            //  PROCEDURE FOR A0 < 1 AND 1 < B0 < 8
            //
            u = gamma_ln1(a0);
            m = (int) (b0 - 1.0e0);
            if (m < 1) goto S60;
            c = 1.0e0;
            for (i = 1; i <= m; i++)
            {
                b0 = b0 - 1.0e0;
                c = c * (b0 / (a0 + b0));
            }

            u = Math.Log(c) + u;
            S60:
            z = a * Math.Log(x) - u;
            b0 = b0 - 1.0e0;
            apb = a0 + b0;
            if (apb > 1.0e0) goto S70;
            t = 1.0e0 + gam1(apb);
            goto S80;
            S70:
            u = a0 + b0 - 1e0;
            t = (1.0e0 + gam1(u)) / apb;
            S80:
            bpser = Math.Exp(z) * (a0 / a) * (1.0e0 + gam1(b0)) / t;
            goto S100;
            S90:
            //
            //  PROCEDURE FOR A0 < 1 AND B0 .GE. 8
            //
            u = gamma_ln1(a0) + algdiv(a0, b0);
            z = a * Math.Log(x) - u;
            bpser = a0 / a * Math.Exp(z);
            S100:
            if (bpser == 0.0e0 || a <= 0.1e0 * eps) return bpser;
            //
            //  COMPUTE THE SERIES
            //
            sum = n = 0.0e0;
            c = 1.0e0;
            tol = eps / a;
            S110:
            n = n + 1.0e0;
            c = c * ((0.5e0 + (0.5e0 - b / n)) * x);
            w = c / (a + n);
            sum = sum + w;
            if (Math.Abs(w) > tol) goto S110;
            bpser = bpser * (1.0e0 + a * sum);

            return bpser;
        }

        public static double beta_rcomp(double a, double b, double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_RCOMP evaluates X^A * Y^B / Beta(A,B).
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
            //    double *A, *B, the parameters of the Beta function.
            //    A and B should be nonnegative.
            //
            //    double *X, *Y, define the numerator of the fraction.
            //
            //  Output:
            //
            //    double BETA_RCOMP, the value of X^A * Y^B / Beta(A,B).
            //
            //  Local:
            //
            //    double CONST = 1/SQRT(2*PI)
            //
        {
            double a0;
            double apb;
            double b0;
            double brcomp;
            double c;
            double Const = 0.398942280401433e0;
            double e;
            double h;
            int i;
            double lambda;
            double lnx;
            double lny;
            int n;
            double t;
            double T1;
            double T2;
            double u;
            double v;
            double x0;
            double y0;
            double z;

            brcomp = 0.0e0;
            if (x == 0.0e0 || y == 0.0e0) return brcomp;
            a0 = Math.Min(a, b);
            if (a0 >= 8.0e0) goto S130;
            if (x > 0.375e0) goto S10;
            lnx = Math.Log(x);
            T1 = -x;
            lny = alnrel(T1);
            goto S30;
            S10:
            if (y > 0.375e0) goto S20;
            T2 = -y;
            lnx = alnrel(T2);
            lny = Math.Log(y);
            goto S30;
            S20:
            lnx = Math.Log(x);
            lny = Math.Log(y);
            S30:
            z = a * lnx + b * lny;
            if (a0 < 1.0e0) goto S40;
            z = z - beta_log(a, b);
            brcomp = Math.Exp(z);
            return brcomp;
            S40:
            //
            //  PROCEDURE FOR A < 1 OR B < 1
            //
            b0 = Math.Max(a, b);
            if (b0 >= 8.0e0) goto S120;
            if (b0 > 1.0e0) goto S70;
            //
            //  ALGORITHM FOR B0 <= 1
            //
            brcomp = Math.Exp(z);
            if (brcomp == 0.0e0) return brcomp;
            apb = a + b;
            if (apb > 1.0e0) goto S50;
            z = 1.0e0 + gam1(apb);
            goto S60;
            S50:
            u = a + b - 1e0;
            z = (1.0e0 + gam1(u)) / apb;
            S60:
            c = (1.0e0 + gam1(a)) * (1.0e0 + gam1(b)) / z;
            brcomp = brcomp * (a0 * c) / (1.0e0 + a0 / b0);
            return brcomp;
            S70:
            //
            //  ALGORITHM FOR 1 < B0 < 8
            //
            u = gamma_ln1(a0);
            n = (int) (b0 - 1.0e0);
            if (n < 1) goto S90;
            c = 1.0e0;
            for (i = 1; i <= n; i++)
            {
                b0 = b0 - 1.0e0;
                c = c * (b0 / (a0 + b0));
            }

            u = Math.Log(c) + u;
            S90:
            z = z - u;
            b0 = b0 - 1.0e0;
            apb = a0 + b0;
            if (apb > 1.0e0) goto S100;
            t = 1.0e0 + gam1(apb);
            goto S110;
            S100:
            u = a0 + b0 - 1e0;
            t = (1.0e0 + gam1(u)) / apb;
            S110:
            brcomp = a0 * Math.Exp(z) * (1.0e0 + gam1(b0)) / t;
            return brcomp;
            S120:
            //
            //  ALGORITHM FOR 8 <= B0
            //
            u = gamma_ln1(a0) + algdiv(a0, b0);
            brcomp = a0 * Math.Exp(z - u);
            return brcomp;
            S130:
            //
            //  PROCEDURE FOR A .GE. 8 AND B .GE. 8
            //
            if (a > b) goto S140;
            h = a / b;
            x0 = h / (1.0e0 + h);
            y0 = 1.0e0 / (1.0e0 + h);
            lambda = a - (a + b) * x;
            goto S150;
            S140:
            h = b / a;
            x0 = 1.0e0 / (1.0e0 + h);
            y0 = h / (1.0e0 + h);
            lambda = (a + b) * y - b;
            S150:
            e = -(lambda / a);
            if (Math.Abs(e) > 0.6e0) goto S160;
            u = rlog1(e);
            goto S170;
            S160:
            u = e - Math.Log(x / x0);
            S170:
            e = lambda / b;
            if (Math.Abs(e) > 0.6e0) goto S180;
            v = rlog1(e);
            goto S190;
            S180:
            v = e - Math.Log(y / y0);
            S190:
            z = Math.Exp(-(a * u + b * v));
            brcomp = Const * Math.Sqrt(b * x0) * z * Math.Exp(-bcorr(a, b));
            return brcomp;
        }

        public static double beta_rcomp1(int mu, double a, double b, double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_RCOMP1 evaluates exp(MU) * X^A * Y^B / Beta(A,B).
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
            //    Input, int MU, ?
            //
            //    Input, double A, B, the parameters of the Beta function.
            //    A and B should be nonnegative.
            //
            //    Input, double X, Y, ?
            //
            //    Output, double BETA_RCOMP1, the value of
            //    exp(MU) * X^A * Y^B / Beta(A,B).
            //
            //  Local:
            //
            //    double CONST = 1/SQRT(2*PI)
            //
        {
            double a0;
            double apb;
            double b0;
            double brcmp1;
            double c;
            double Const = 0.398942280401433e0;
            double e;
            double h;
            int i;
            double lambda;
            double lnx;
            double lny;
            int n;
            double t;
            double T1;
            double T2;
            double T3;
            double T4;
            double u;
            double v;
            double x0;
            double y0;
            double z;

            a0 = Math.Min(a, b);
            if (a0 >= 8.0e0) goto S130;
            if (x > 0.375e0) goto S10;
            lnx = Math.Log(x);
            T1 = -x;
            lny = alnrel(T1);
            goto S30;
            S10:
            if (y > 0.375e0) goto S20;
            T2 = -y;
            lnx = alnrel(T2);
            lny = Math.Log(y);
            goto S30;
            S20:
            lnx = Math.Log(x);
            lny = Math.Log(y);
            S30:
            z = a * lnx + b * lny;
            if (a0 < 1.0e0) goto S40;
            z = z - beta_log(a, b);
            brcmp1 = esum(mu, z);
            return brcmp1;
            S40:
            //
            //  PROCEDURE FOR A < 1 OR B < 1
            //
            b0 = Math.Max(a, b);
            if (b0 >= 8.0e0) goto S120;
            if (b0 > 1.0e0) goto S70;
            //
            //  ALGORITHM FOR B0 <= 1
            //
            brcmp1 = esum(mu, z);
            if (brcmp1 == 0.0e0) return brcmp1;
            apb = a + b;
            if (apb > 1.0e0) goto S50;
            z = 1.0e0 + gam1(apb);
            goto S60;
            S50:
            u = a + b - 1e0;
            z = (1.0e0 + gam1(u)) / apb;
            S60:
            c = (1.0e0 + gam1(a)) * (1.0e0 + gam1(b)) / z;
            brcmp1 = brcmp1 * (a0 * c) / (1.0e0 + a0 / b0);
            return brcmp1;
            S70:
            //
            //  ALGORITHM FOR 1 < B0 < 8
            //
            u = gamma_ln1(a0);
            n = (int) (b0 - 1.0e0);
            if (n < 1) goto S90;
            c = 1.0e0;
            for (i = 1; i <= n; i++)
            {
                b0 = b0 - 1.0e0;
                c = c * (b0 / (a0 + b0));
            }

            u = Math.Log(c) + u;
            S90:
            z = z - u;
            b0 = b0 - 1.0e0;
            apb = a0 + b0;
            if (apb > 1.0e0) goto S100;
            t = 1.0e0 + gam1(apb);
            goto S110;
            S100:
            u = a0 + b0 - 1e0;
            t = (1.0e0 + gam1(u)) / apb;
            S110:
            brcmp1 = a0 * esum(mu, z) * (1.0e0 + gam1(b0)) / t;
            return brcmp1;
            S120:
            //
            //  ALGORITHM FOR 8 <= B0
            //
            u = gamma_ln1(a0) + algdiv(a0, b0);
            T3 = z - u;
            brcmp1 = a0 * esum(mu, T3);
            return brcmp1;
            S130:
            //
            // PROCEDURE FOR A .GE. 8 AND B .GE. 8
            //
            if (a > b) goto S140;
            h = a / b;
            x0 = h / (1.0e0 + h);
            y0 = 1.0e0 / (1.0e0 + h);
            lambda = a - (a + b) * x;
            goto S150;
            S140:
            h = b / a;
            x0 = 1.0e0 / (1.0e0 + h);
            y0 = h / (1.0e0 + h);
            lambda = (a + b) * y - b;
            S150:
            e = -(lambda / a);
            if (Math.Abs(e) > 0.6e0) goto S160;
            u = rlog1(e);
            goto S170;
            S160:
            u = e - Math.Log(x / x0);
            S170:
            e = lambda / b;
            if (Math.Abs(e) > 0.6e0) goto S180;
            v = rlog1(e);
            goto S190;
            S180:
            v = e - Math.Log(y / y0);
            S190:
            T4 = -(a * u + b * v);
            z = esum(mu, T4);
            brcmp1 = Const * Math.Sqrt(b * x0) * z * Math.Exp(-bcorr(a, b));
            return brcmp1;
        }

        public static double beta_up(double a, double b, double x, double y, int n,
                double eps)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
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
            //    Input, double *A, *B, the parameters of the function.
            //    A and B should be nonnegative.
            //
            //    Input, double *X, *Y, ?
            //
            //    Input, int *N, ?
            //
            //    Input, double *EPS, the tolerance.
            //
            //    Output, double BETA_UP, the value of IX(A,B) - IX(A+N,B).
            //
        {
            double ap1;
            double apb;
            double bup;
            double d;
            int i;
            int k;
            int K1 = 1;
            int K2 = 0;
            int kp1;
            double l;
            int mu;
            int nm1;
            double r;
            double t;
            double w;
            //
            //  OBTAIN THE SCALING FACTOR EXP(-MU) AND
            //  EXP(MU)*(X^A*Y^B/BETA(A,B))/A
            //
            apb = a + b;
            ap1 = a + 1.0e0;
            mu = 0;
            d = 1.0e0;
            if (n == 1 || a < 1.0e0) goto S10;
            if (apb < 1.1e0 * ap1) goto S10;
            mu = (int) Math.Abs(exparg(K1));
            k = (int) exparg(K2);
            if (k < mu) mu = k;
            t = mu;
            d = Math.Exp(-t);
            S10:
            bup = beta_rcomp1(mu, a, b, x, y) / a;
            if (n == 1 || bup == 0.0e0) return bup;
            nm1 = n - 1;
            w = d;
            //
            //  LET K BE THE INDEX OF THE MAXIMUM TERM
            //
            k = 0;
            if (b <= 1.0e0) goto S50;
            if (y > 1e-4) goto S20;
            k = nm1;
            goto S30;
            S20:
            r = (b - 1.0e0) * x / y - a;
            if (r < 1.0e0) goto S50;
            t = (double) nm1;
            k = nm1;
            if (r < t) k = (int) r;
            S30:
            //
            //  ADD THE INCREASING TERMS OF THE SERIES
            //
            for (i = 1; i <= k; i++)
            {
                l = i - 1;
                d = (apb + l) / (ap1 + l) * x * d;
                w = w + d;
            }

            if (k == nm1) goto S70;
            S50:
            //
            //  ADD THE REMAINING TERMS OF THE SERIES
            //
            kp1 = k + 1;
            for (i = kp1; i <= nm1; i++)
            {
                l = i - 1;
                d = (apb + l) / (ap1 + l) * x * d;
                w = w + d;
                if (d <= eps * w) goto S70;
            }

            S70:
            //
            //  TERMINATE THE PROCEDURE
            //
            bup = bup * w;
            return bup;
        }


    }
}