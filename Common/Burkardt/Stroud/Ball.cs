using System;
using Burkardt.Quadrature;

namespace Burkardt.Stroud
{
    public static class Ball
    {
        
        public static double ball_f1_nd(int setting, Func<int, int, double[], double> func, int n, double[] center,
                double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_F1_ND approximates an integral inside a ball in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
            //
            //  Discussion:
            //
            //    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( int n, double x[] ), the name of the user supplied
            //    function which evaluates F at the N-vector X.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double CENTER[N], the center of the ball.
            //
            //    Input, double R, the radius of the ball.
            //
            //    Output, double BALL_F1_ND, the approximate integral of the function.
            //
        {
            int i;
            int ihi;
            int itemp;
            int j;
            int k;
            int khi;
            int ktemp;
            double quad;
            double result;
            double t;
            double temp;
            double u;
            double u2;
            double v;
            double volume;
            double w;
            double[] x;
            double y;

            if (r == 0.0)
            {
                result = 0.0;
                return result;
            }

            x = new double[n];

            u2 = (1.0 - 2.0 * Math.Sqrt(1.0 / (double)(n + 4)))
                 / (double)(n + 2);
            u = Math.Sqrt(u2);
            for (i = 0; i < n; i++)
            {
                x[i] = center[i] - r * u;
            }

            w = 1.0 / (double)((n + 1) * (int)Math.Pow(2, n));

            quad = 0.0;
            ihi = (int)Math.Pow(2, n);

            for (i = 0; i < ihi; i++)
            {
                itemp = i;
                for (j = 0; j < n; j++)
                {
                    u = (center[j] - x[j]) / r;

                    if ((itemp % 2) == 1)
                    {
                        x[j] = center[j] - Math.Abs(x[j] - center[j]);
                    }
                    else
                    {
                        x[j] = center[j] + Math.Abs(x[j] - center[j]);
                    }

                    itemp = itemp / 2;
                }

                quad = quad + w * func(setting, n, x);
            }

            temp = Math.Sqrt((double)(n + 4));

            t = Math.Sqrt(2.0 * (double)(n + 1) / (double)(n + 2))
                / ((double)(n) * temp);
            ;

            y = (1.0 + 2.0 / ((double)(n) * temp)) / (double)(n + 2);
            v = Math.Sqrt(y - t);
            u = Math.Sqrt(y + (double)(n - 1) * t);

            khi = (int)Math.Pow(2, n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = center[j] - r * v;
                }

                x[i] = center[i] - r * u;

                for (k = 0; k < khi; k++)
                {
                    ktemp = k;
                    for (j = 0; j < n; j++)
                    {
                        if ((ktemp % 2) == 1)
                        {
                            x[j] = center[j] - Math.Abs(x[j] - center[j]);
                        }
                        else
                        {
                            x[j] = center[j] + Math.Abs(x[j] - center[j]);
                        }

                        ktemp = ktemp / 2;
                    }
                    
                    quad = quad + w * func(setting, n, x);
                }

                x[i] = center[i] - r * v;
            }

            volume = ball_volume_nd(n, r);
            result = quad * volume;

            return result;
        }

        public static double ball_f3_nd(int setting, Func<int, int, double[], double> func, int n, double[] center,
                double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_F3_ND approximates an integral inside a ball in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
            //
            //  Discussion:
            //
            //    A 2**(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( int n, double x[] ), the name of the user supplied
            //    function which evaluates F at the N-vector X.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double CENTER[N], the center of the ball.
            //
            //    Input, double R, the radius of the ball.
            //
            //    Output, double BALL_F3_ND, the approximate integral of the function.
            //
        {
            int i;
            int j;
            int jhi;
            int jtemp;
            int k;
            double quad;
            double result;
            double ri;
            double s;
            double volume;
            double weight;
            double[] x;

            if (r == 0.0)
            {
                result = 0.0;
                return result;
            }

            x = new double[n];

            quad = 0.0;
            //
            //  The first point is the center of the ball.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = center[i];
            }

            weight = 4.0 / (double)((int)Math.Pow(n + 2, 2));
            quad = quad + weight * func(setting, n, x);

            s = 1.0 / Math.Sqrt((double)(n + 4));

            for (i = 0; i < n; i++)
            {
                ri = Math.Sqrt((double)(i + 3) / (double)(n + 4));
                //
                //  Set up the first point, with (I) zeroes, RI, and then N-I-1 S's.
                //
                for (j = 0; j < n; j++)
                {
                    if (j < i)
                    {
                        x[j] = center[j];
                    }
                    else if (j == i)
                    {
                        x[j] = center[j] + r * ri;
                    }
                    else
                    {
                        x[j] = center[j] + r * s;
                    }
                }

                weight = Math.Pow(2.0, i + 1 - n) * (double)(n + 4)
                         / (double)((i + 2) * (i + 3) * (n + 2));
                //
                //  Now go through all sign permutations of the basic point.
                //
                jhi = (int)Math.Pow(2, n - i);

                for (j = 0; j < jhi; j++)
                {
                    jtemp = j;

                    for (k = i; k < n; k++)
                    {
                        if ((jtemp % 2) == 1)
                        {
                            x[k] = center[k] - Math.Abs(x[k] - center[k]);
                        }
                        else
                        {
                            x[k] = center[k] + Math.Abs(x[k] - center[k]);
                        }

                        jtemp = jtemp / 2;
                    }

                    quad = quad + weight * func(setting, n, x);
                }
            }

            volume = ball_volume_nd(n, r);
            result = quad * volume;

            return result;
        }

        public static double ball_monomial_nd(int n, int[] p, double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_MONOMIAL_ND integrates a monomial on a ball in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) <= R * R
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Gerald Folland,
            //    How to Integrate a Polynomial Over a Sphere,
            //    American Mathematical Monthly,
            //    Volume 108, Number 5, May 2001, pages 446-448.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, int P[N], the exponents of X(1) through X(N) in the monomial.
            //    The exponents must be nonnegative.
            //
            //    Input, double R, the radius of the ball.
            //
            //    Output, double BALL_MONOMIAL_ND, the integral of
            //    X1**P(1)*X2**P(2)*...*XN**P(N) over the ball.
            //
        {
            int i;
            double power;
            double value;

            power = (double)(n);
            for (i = 0; i < n; i++)
            {
                power = power + (double)p[i];
            }

            value = Sphere.sphere_unit_monomial_nd(n, p) * Math.Pow(r, power) / power;

            return value;
        }

        public static double ball_unit_07_3d(int setting, Func<int, double, double, double, double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_UNIT_07_3D approximates an integral inside the unit ball in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z <= 1.
            //
            //  Discussion:
            //
            //    A 64 point 7-th degree formula is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y, double z ), the name of 
            //    the user supplied function which evaluates F(X,Y,Z).
            //
            //    Output, double BALL_UNIT_07_3D, the approximate integral of the function.
            //
        {
            int order = 4;

            double angle;
            int i;
            int j;
            int k;
            
            double quad;
            double result;
            double volume;
            double w;
            double[] weight1 =
            {
                0.19455533421780251826,
                0.13877799911553081506,
                0.13877799911553081506,
                0.19455533421780251826
            };
            double[] weight2 = new double[4];
            double[] weight3 = new double[4];
            double x;
            double[] xtab1 =
            {
                -0.906179845938663992797626878299,
                -0.538469310105683091036314420700,
                0.538469310105683091036314420700,
                0.906179845938663992797626878299
            };
            double[] xtab2 = new double[4];
            double[] xtab3 = new double[4];
            double y;
            double z;
            //
            //  Set XTAB2 and WEIGHT2.
            //
            for (j = 0; j < order; j++)
            {
                angle = Math.PI * (double)(2 * j - 1) / (double)(2 * order);
                xtab2[j] = Math.Cos(angle);
            }

            for (j = 0; j < order; j++)
            {
                weight2[j] = 1.0;
            }

            //
            //  Set XTAB3 and WEIGHT3 for the interval [-1,1].
            //
            LegendreQuadrature.legendre_set(order, ref xtab3, ref weight3);

            w = 3.0 / 16.0;

            quad = 0.0;

            for (i = 0; i < order; i++)
            {
                for (j = 0; j < order; j++)
                {
                    for (k = 0; k < order; k++)
                    {
                        x = xtab1[i] * Math.Sqrt(1.0 - xtab2[j] * xtab2[j])
                                     * Math.Sqrt(1.0 - xtab3[k] * xtab3[k]);
                        y = xtab1[i] * xtab2[j] * Math.Sqrt(1.0 - xtab3[k] * xtab3[k]);
                        z = xtab1[i] * xtab3[k];

                        quad = quad + w * weight1[i] * weight2[j] * weight3[k]
                            * func(setting, x, y, z);
                    }
                }
            }

            volume = ball_unit_volume_3d();
            result = quad * volume;

            return result;
        }

        public static double ball_unit_14_3d(int setting, Func<int, double, double, double, double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_UNIT_14_3D approximates an integral inside the unit ball in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z <= 1.
            //
            //  Discussion:
            //
            //    A 288 point 14-th degree formula is used, Stroud number S3:14-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y, double z ), the name of 
            //    the user supplied function which evaluates F(X,Y,Z).
            //
            //    Output, double RESULT, the approximate integral of the function.
            //
        {
            int i;
            int j;
            int k;
            int l;
            int m;
            int n;
            double quad;
            double[] r =
            {
                0.968160240, 0.836031107, 0.613371433, 0.324253423
            };
            double result;
            double temp;
            double volume;
            double w1;
            double w2;
            double[] weight =
            {
                0.076181268, 0.126263673, 0.098048133, 0.032840260
            };
            double x;
            double[] xtab =
            {
                -0.151108275, 0.315838353, 0.346307112,
                -0.101808787, -0.409228403
            };
            double y;
            double[] ytab =
            {
                0.155240600, 0.257049387, 0.666277790,
                0.817386065, 0.501547712
            };
            double z;
            double[] ztab =
            {
                0.976251323, 0.913330032, 0.660412970,
                0.567022920, 0.762221757
            };

            quad = 0.0;

            for (m = 0; m < 4; m++)
            {
                w1 = 125.0 * weight[m] / 3360.0;
                x = 0.525731112 * r[m];
                y = 0.850650808 * r[m];
                z = 0.0;


                for (j = 0; j < 2; j++)
                {
                    x = -x;
                    for (k = 0; k < 2; k++)
                    {
                        y = -y;
                        for (l = 0; l < 3; l++)
                        {
                            temp = z;
                            z = y;
                            y = x;
                            x = temp;
                            quad = quad + w1 * func(setting, x, y, z);
                        }
                    }
                }

                w2 = 143.0 * weight[m] / 3360.0;

                for (n = 0; n < 5; n++)
                {
                    x = xtab[n] * r[m];
                    y = ytab[n] * r[m];
                    z = ztab[n] * r[m];

                    for (i = 0; i < 3; i++)
                    {
                        temp = x;
                        x = z;
                        z = -y;
                        y = -temp;

                        for (j = 0; j < 3; j++)
                        {
                            temp = z;
                            z = y;
                            y = x;
                            x = temp;

                            quad = quad + w2 * func(setting, x, y, z);
                        }

                        y = -y;
                        z = -z;
                        quad = quad + w2 * func(setting, x, y, z);
                    }
                }
            }

            volume = ball_unit_volume_3d();
            result = quad * volume;

            return result;
        }

        public static double ball_unit_15_3d(int setting, Func<int, double, double, double, double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_UNIT_15_3D approximates an integral inside the unit ball in 3D.
            //
            //  Integration region:
            //
            //    X * X + Y * Y + Z * Z <= 1.
            //
            //  Discussion:
            //
            //    A 512 point 15-th degree formula is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 October 2000
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y, double z ), the name of the 
            //    user supplied function which evaluates F(X,Y,Z).
            //
            //    Output, double BALL_UNIT_15_3D, the approximate integral of the function.
            //
        {
            double cj;
            double ck;
            int i;
            int j;
            int k;
            int order1 = 4;
            int order2 = 8;
            
            double quad;
            double result;
            double sj;
            double sk;
            double volume;
            double w;
            double[] weight1 =
            {
                0.0328402599, 0.0980481327, 0.1262636728, 0.0761812678
            };
            double[] weight2;
            double x;
            double[] xtab1 =
            {
                0.3242534234, 0.6133714327, 0.8360311073, 0.9681602395
            };
            double[] xtab2;
            double y;
            double z;

            xtab2 = new double[order2];
            weight2 = new double[order2];

            LegendreQuadrature.legendre_set(order2, ref xtab2, ref weight2);

            w = 3.0 / 32.0;

            quad = 0.0;

            for (i = 0; i < order1; i++)
            {
                for (j = 0; j < order2; j++)
                {
                    sj = xtab2[j];
                    cj = Math.Sqrt(1.0 - sj * sj);

                    for (k = 1; k <= 16; k++)
                    {
                        sk = Math.Sin((double)(k) * Math.PI / 8.0);
                        ck = Math.Cos((double)(k) * Math.PI / 8.0);
                        x = xtab1[i] * cj * ck;
                        y = xtab1[i] * cj * sk;
                        z = xtab1[i] * sj;
                        quad = quad + w * weight1[i] * weight2[j] * func(setting, x, y, z);
                    }
                }
            }

            volume = ball_unit_volume_3d();
            result = quad * volume;

            return result;
        }

        public static double ball_unit_f1_nd(int setting, Func<int, int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_UNIT_F1_ND approximates an integral inside the unit ball in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) <= 1.
            //
            //  Discussion:
            //
            //    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( int n, double x[] ), the name of the 
            //    user supplied function which evaluates F at the N-vector X.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double BALL_UNIT_F1_ND, the approximate integral of the function.
            //
        {
            int i;
            int ihi;
            int itemp;
            int j;
            int k;
            int khi;
            int ktemp;
            double quad;
            double result;
            double t;
            double temp;
            double u;
            double u2;
            double v;
            double volume;
            double w;
            double[] x;
            double y;

            x = new double[n];

            u2 = (1.0 - 2.0 * Math.Sqrt(1.0 / (double)(n + 4)))
                 / (double)(n + 2);
            u = Math.Sqrt(u2);
            for (i = 0; i < n; i++)
            {
                x[i] = -u;
            }

            w = 1.0 / (double)((n + 1) * (int)Math.Pow(2, n));

            quad = 0.0;
            ihi = (int)Math.Pow(2, n);

            for (i = 0; i < ihi; i++)
            {
                itemp = i;

                for (j = 0; j < n; j++)
                {
                    if ((itemp % 2) == 1)
                    {
                        x[j] = -Math.Abs(x[j]);
                    }
                    else
                    {
                        x[j] = Math.Abs(x[j]);
                    }

                    itemp = itemp / 2;
                }

                quad = quad + w * func(setting, n, x);
            }

            temp = Math.Sqrt((double)(n + 4));

            t = Math.Sqrt(2.0 * (double)(n + 1) / (double)(n + 2))
                / ((double)(n) * temp);

            y = (1.0 + 2.0 / ((double)(n) * temp))
                / (double)(n + 2);
            v = Math.Sqrt(y - t);
            u = Math.Sqrt(y + (double)(n - 1) * t);

            khi = (int)Math.Pow(2, n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = -v;
                }

                x[i] = -u;

                for (k = 0; k < khi; k++)
                {
                    ktemp = k;

                    for (j = 0; j < n; j++)
                    {
                        if ((ktemp % 2) == 1)
                        {
                            x[j] = -Math.Abs(x[j]);
                        }
                        else
                        {
                            x[j] = Math.Abs(x[j]);
                        }

                        ktemp = ktemp / 2;
                    }

                    quad = quad + w * func(setting, n, x);
                }

                x[i] = -v;
            }

            volume = ball_unit_volume_nd(n);
            result = quad * volume;

            return result;
        }

        public static double ball_unit_f3_nd(int setting, Func<int, int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_UNIT_F3_ND approximates an integral inside the unit ball in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) <= 1.
            //
            //  Discussion:
            //
            //    A 2^(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( int n, double x[] ), the name of the 
            //    user supplied function which evaluates F at the N-vector X.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double BALL_UNIT_F3_ND, the approximate integral of the function.
            //
        {
            int i;
            int j;
            int jtemp;
            int k;
            double quad;
            double result;
            double ri;
            double s;
            double volume;
            double weight;
            double[] x;

            quad = 0.0;
            //
            //  The first point is the center of the ball.
            //
            x = new double[n];
            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            weight = 4.0 / (double)((n + 2) * (n + 2));
            quad = quad + weight * func(setting, n, x);

            s = 1.0 / Math.Sqrt((double)(n + 4));

            for (i = 0; i < n; i++)
            {
                ri = Math.Sqrt((double)(i + 3) / (double)(n + 4));
                //
                //  Set up the first point, with (I-1) zeroes, RI, and then N-I S's.
                //
                for (j = 0; j < n; j++)
                {
                    if (j < i)
                    {
                        x[j] = 0.0;
                    }
                    else if (j == i)
                    {
                        x[j] = ri;
                    }
                    else
                    {
                        x[j] = s;
                    }
                }

                weight = Math.Pow(2.0, i + 1 - n) * (double)(n + 4)
                         / (double)((i + 2) * (i + 3) * (n + 2));
                //
                //  Now go through all sign permutations of the basic point.
                //
                for (j = 0; j < (int)Math.Pow(2, n - i); j++)
                {
                    jtemp = j;

                    for (k = i; k < n; k++)
                    {
                        if ((jtemp % 2) == 1)
                        {
                            x[k] = -Math.Abs(x[k]);
                        }
                        else
                        {
                            x[k] = Math.Abs(x[k]);
                        }

                        jtemp = jtemp / 2;
                    }

                    quad = quad + weight * func(setting, n, x);
                }
            }

            volume = ball_unit_volume_nd(n);
            result = quad * volume;

            return result;
        }

        public static double ball_unit_volume_3d()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_UNIT_VOLUME_3D computes the volume of the unit ball in 3D.
            //
            //  Integration region:
            //
            //    X * X + Y * Y + Z * Z <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double BALL_UNIT_VOLUME_3D, the volume of the ball.
            //
        {
            
            double value;

            value = (4.0 / 3.0) * Math.PI;

            return value;
        }

        public static double ball_unit_volume_nd(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_UNIT_VOLUME_ND computes the volume of the unit ball in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) <= 1.
            //
            //  Discussion:
            //
            //    N  Volume
            //
            //    2             PI
            //    3  (4/3)    * PI
            //    4  (1/2)    * PI^2
            //    5  (8/15)   * PI^2
            //    6  (1/6)    * PI^3
            //    7  (16/105) * PI^3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double BALL_UNIT_VOLUME_ND, the volume of the ball.
            //
        {
            int i;
            int m;
            
            double volume;

            if ((n % 2) == 0)
            {
                m = n / 2;
                volume = Math.Pow(Math.PI, m);
                for (i = 1; i <= m; i++)
                {
                    volume = volume / (double)(i);
                }
            }
            else
            {
                m = (n - 1) / 2;
                volume = Math.Pow(Math.PI, m) * (int)Math.Pow(2, n);
                for (i = m + 1; i <= 2 * m + 1; i++)
                {
                    volume = volume / (double)(i);
                }
            }

            return volume;
        }

        public static double ball_volume_3d(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_VOLUME_3D computes the volume of a ball in 3D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y + Z*Z <= R * R
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the ball.
            //
            //    Output, double BALL_VOLUME_3D, the volume of the ball.
            //
        {
            
            double volume;

            volume = (4.0 / 3.0) * Math.PI * r * r * r;

            return volume;
        }

        public static double ball_volume_nd(int n, double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL_VOLUME_ND computes the volume of a ball in ND.
            //
            //  Integration region:
            //
            //    sum ( X(1:N)^2 ) <= R * R
            //
            //  Discussion:
            //
            //    N  Volume
            //
            //    2             PI   * R^2
            //    3  (4/3)    * PI   * R^3
            //    4  (1/2)    * PI^2 * R^4
            //    5  (8/15)   * PI^2 * R^5
            //    6  (1/6)    * PI^3 * R^6
            //    7  (16/105) * PI^3 * R^7
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double R, the radius of the ball.
            //
            //    Output, double BALL_VOLUME_ND, the volume of the ball.
            //
        {
            double volume;

            volume = ball_unit_volume_nd(n) * Math.Pow(r, n);

            return volume;
        }

    }
}