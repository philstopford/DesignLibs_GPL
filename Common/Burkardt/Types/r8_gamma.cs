using System;
using Burkardt.Probability;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_gamma(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA evaluates Gamma(X) for a real argument.
        //
        //  Discussion:
        //
        //    This routine calculates the gamma function for a real argument X.
        //
        //    Computation is based on an algorithm outlined in reference 1.
        //    The program uses rational functions that approximate the gamma
        //    function to at least 20 significant decimal digits.  Coefficients
        //    for the approximation over the interval (1,2) are unpublished.
        //    Those for the approximation for 12 <= X are from reference 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody, Laura Stoltz.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody,
        //    An Overview of Software Development for Special Functions,
        //    in Numerical Analysis Dundee, 1975,
        //    edited by GA Watson,
        //    Lecture Notes in Mathematics 506,
        //    Springer, 1976.
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
        //    Charles Mesztenyi, John Rice, Henry Thatcher,
        //    Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968,
        //    LC: QA297.C64.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the function.
        //
        //    Output, double R8_GAMMA, the value of the function.
        //
    {
        //
        //  Coefficients for minimax approximation over (12, INF).
        //
        double[] c =
            {
                -1.910444077728E-03,
                8.4171387781295E-04,
                -5.952379913043012E-04,
                7.93650793500350248E-04,
                -2.777777777777681622553E-03,
                8.333333333333333331554247E-02,
                5.7083835261E-03
            }
            ;
        const double eps = 2.22E-16;
        const double half = 0.5;
        int i;
        const double one = 1.0;
        double[] p =
            {
                -1.71618513886549492533811E+00,
                2.47656508055759199108314E+01,
                -3.79804256470945635097577E+02,
                6.29331155312818442661052E+02,
                8.66966202790413211295064E+02,
                -3.14512729688483675254357E+04,
                -3.61444134186911729807069E+04,
                6.64561438202405440627855E+04
            }
            ;
        double[] q =
            {
                -3.08402300119738975254353E+01,
                3.15350626979604161529144E+02,
                -1.01515636749021914166146E+03,
                -3.10777167157231109440444E+03,
                2.25381184209801510330112E+04,
                4.75584627752788110767815E+03,
                -1.34659959864969306392456E+05,
                -1.15132259675553483497211E+05
            }
            ;
        double res;
        const double sqrtpi = 0.9189385332046727417803297;
        const double twelve = 12.0;
        const double two = 2.0;
        double value = 0;
        const double xbig = 171.624;
        const double xinf = 1.79E+308;
        const double xminin = 2.23E-308;
        double y1;
        const double zero = 0.0;

        bool parity = false;
        double fact = one;
        int n = 0;
        double y = x;
        //
        //  Argument is negative.
        //
        if (y <= zero)
        {
            y = -x;
            y1 = (int) y;
            res = y - y1;

            if (Math.Abs(res - zero) > double.Epsilon)
            {
                if (Math.Abs(y1 - (int) (y1 * half) * two) > double.Epsilon)
                {
                    parity = true;
                }

                fact = -Math.PI / Math.Sin(Math.PI * res);
                y += one;
            }
            else
            {
                res = xinf;
                value = res;
                return value;
            }
        }

        switch (y)
        {
            //
            //  Argument is positive.
            //
            //
            //  Argument < EPS.
            //
            case < eps and >= xminin:
                res = one / y;
                break;
            case < eps:
                res = xinf;
                value = res;
                return value;
            case < twelve:
            {
                y1 = y;
                //
                //  0.0 < argument < 1.0.
                //
                double z;
                if (y < one)
                {
                    z = y;
                    y += one;
                }
                //
                //  1.0 < argument < 12.0.
                //  Reduce argument if necessary.
                //
                else
                {
                    n = (int) y - 1;
                    y -= n;
                    z = y - one;
                }

                //
                //  Evaluate approximation for 1.0 < argument < 2.0.
                //
                double xnum = zero;
                double xden = one;
                for (i = 0; i < 8; i++)
                {
                    xnum = (xnum + p[i]) * z;
                    xden = xden * z + q[i];
                }

                res = xnum / xden + one;
                //
                //  Adjust result for case  0.0 < argument < 1.0.
                //
                if (y1 < y)
                {
                    res /= y1;
                }
                //
                //  Adjust result for case 2.0 < argument < 12.0.
                //
                else if (y < y1)
                {
                    for (i = 1; i <= n; i++)
                    {
                        res *= y;
                        y += one;
                    }
                }

                break;
            }
            //
            //  Evaluate for 12.0 <= argument.
            //
            case <= xbig:
            {
                double ysq = y * y;
                double sum = c[6];
                for (i = 0; i < 6; i++)
                {
                    sum = sum / ysq + c[i];
                }

                sum = sum / y - y + sqrtpi;
                sum += (y - half) * Math.Log(y);
                res = Math.Exp(sum);
                break;
            }
            default:
                res = xinf;
                value = res;
                return value;
        }

        res = parity switch
        {
            //
            //  Final adjustments and return.
            //
            true => -res,
            _ => res
        };

        if (Math.Abs(fact - one) > double.Epsilon)
        {
            res = fact / res;
        }

        value = res;

        return value;
    }

    public static double r8_gamma_inc(double p, double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_INC computes the incomplete Gamma function.
        //
        //  Discussion:
        //
        //    GAMMA_INC(P,       0) = 0,
        //    GAMMA_INC(P,Infinity) = 1.
        //
        //    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T^(P-1) EXP(-T) DT / GAMMA(P).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by B L Shea.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    B L Shea,
        //    Chi-squared and Incomplete Gamma Integral,
        //    Algorithm AS239,
        //    Applied Statistics,
        //    Volume 37, Number 3, 1988, pages 466-473.
        //
        //  Parameters:
        //
        //    Input, double P, the exponent parameter.
        //    0.0 < P.
        //
        //    Input, double X, the integral limit parameter.
        //    If X is less than or equal to 0, the value is returned as 0.
        //
        //    Output, double R8_GAMMA_INC, the value of the function.
        //
    {
        double a;
        double arg;
        double c;
        const double exp_arg_min = -88.0;
        const double overflow = 1.0E+37;
        const double plimit = 1000.0;
        double pn1;
        const double tol = 1.0E-07;
        const double xbig = 1.0E+08;

        double value = 0;

        switch (p)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("R8_GAMMA_INC - Fatal error!");
                Console.WriteLine("  Parameter P <= 0.");
                return 1;
        }

        switch (x)
        {
            case <= 0.0:
                value = 0.0;
                return value;
        }

        //
        //  Use a normal approximation if PLIMIT < P.
        //
        if (plimit < p)
        {
            pn1 = 3.0 * Math.Sqrt(p) * (Math.Pow(x / p, 1.0 / 3.0)
                + 1.0 / (9.0 * p) - 1.0);
            value = Normal.normal_01_cdf(pn1);
            return value;
        }

        //
        //  Is X extremely large compared to P?
        //
        if (xbig < x)
        {
            value = 1.0;
            return value;
        }

        //
        //  Use Pearson's series expansion.
        //  (P is not large enough to force overflow in the log of Gamma.
        //
        if (x <= 1.0 || x < p)
        {
            arg = p * Math.Log(x) - x - Helpers.LogGamma(p + 1.0);
            c = 1.0;
            value = 1.0;
            a = p;

            for (;;)
            {
                a += 1.0;
                c = c * x / a;
                value += c;

                if (c <= tol)
                {
                    break;
                }
            }

            arg += Math.Log(value);

            value = exp_arg_min <= arg ? Math.Exp(arg) : 0.0;
        }
        //
        //  Use a continued fraction expansion.
        //
        else
        {
            arg = p * Math.Log(x) - x - Helpers.LogGamma(p);
            a = 1.0 - p;
            double b = a + x + 1.0;
            c = 0.0;
            pn1 = 1.0;
            double pn2 = x;
            double pn3 = x + 1.0;
            double pn4 = x * b;
            value = pn3 / pn4;

            for (;;)
            {
                a += 1.0;
                b += 2.0;
                c += 1.0;
                double pn5 = b * pn3 - a * c * pn1;
                double pn6 = b * pn4 - a * c * pn2;

                switch (Math.Abs(pn6))
                {
                    case > 0.0:
                    {
                        double rn = pn5 / pn6;

                        if (Math.Abs(value - rn) <= Math.Min(tol, tol * rn))
                        {
                            arg += Math.Log(value);

                            if (exp_arg_min <= arg)
                            {
                                value = 1.0 - Math.Exp(arg);
                            }
                            else
                            {
                                value = 1.0;
                            }

                            return value;
                        }

                        value = rn;
                        break;
                    }
                }

                pn1 = pn3;
                pn2 = pn4;
                pn3 = pn5;
                pn4 = pn6;
                //
                //  Rescale terms in continued fraction if terms are large.
                //
                if (!(overflow <= Math.Abs(pn5)))
                {
                    continue;
                }

                pn1 /= overflow;
                pn2 /= overflow;
                pn3 /= overflow;
                pn4 /= overflow;
            }
        }

        return value;
    }

    public static double r8_gamma_log(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_LOG evaluates the logarithm of the gamma function.
        //
        //  Discussion:
        //
        //    This routine calculates the LOG(GAMMA) function for a positive real
        //    argument X.  Computation is based on an algorithm outlined in
        //    references 1 and 2.  The program uses rational functions that
        //    theoretically approximate LOG(GAMMA) to at least 18 significant
        //    decimal digits.  The approximation for X > 12 is from reference
        //    3, while approximations for X < 12.0 are similar to those in
        //    reference 1, but are unpublished.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody, Laura Stoltz.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody, Kenneth Hillstrom,
        //    Chebyshev Approximations for the Natural Logarithm of the
        //    Gamma Function,
        //    Mathematics of Computation,
        //    Volume 21, Number 98, April 1967, pages 198-203.
        //
        //    Kenneth Hillstrom,
        //    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
        //    May 1969.
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
        //    Charles Mesztenyi, John Rice, Henry Thatcher,
        //    Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968,
        //    LC: QA297.C64.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the function.
        //
        //    Output, double R8_GAMMA_LOG, the value of the function.
        //
    {
        double[] c =
        {
            -1.910444077728E-03,
            8.4171387781295E-04,
            -5.952379913043012E-04,
            7.93650793500350248E-04,
            -2.777777777777681622553E-03,
            8.333333333333333331554247E-02,
            5.7083835261E-03
        };
        const double d1 = -5.772156649015328605195174E-01;
        const double d2 = 4.227843350984671393993777E-01;
        const double d4 = 1.791759469228055000094023;
        const double frtbig = 2.25E+76;
        double[] p1 =
        {
            4.945235359296727046734888,
            2.018112620856775083915565E+02,
            2.290838373831346393026739E+03,
            1.131967205903380828685045E+04,
            2.855724635671635335736389E+04,
            3.848496228443793359990269E+04,
            2.637748787624195437963534E+04,
            7.225813979700288197698961E+03
        };
        double[] p2 =
        {
            4.974607845568932035012064,
            5.424138599891070494101986E+02,
            1.550693864978364947665077E+04,
            1.847932904445632425417223E+05,
            1.088204769468828767498470E+06,
            3.338152967987029735917223E+06,
            5.106661678927352456275255E+06,
            3.074109054850539556250927E+06
        };
        double[] p4 =
        {
            1.474502166059939948905062E+04,
            2.426813369486704502836312E+06,
            1.214755574045093227939592E+08,
            2.663432449630976949898078E+09,
            2.940378956634553899906876E+10,
            1.702665737765398868392998E+11,
            4.926125793377430887588120E+11,
            5.606251856223951465078242E+11
        };
        double[] q1 =
        {
            6.748212550303777196073036E+01,
            1.113332393857199323513008E+03,
            7.738757056935398733233834E+03,
            2.763987074403340708898585E+04,
            5.499310206226157329794414E+04,
            6.161122180066002127833352E+04,
            3.635127591501940507276287E+04,
            8.785536302431013170870835E+03
        };
        double[] q2 =
        {
            1.830328399370592604055942E+02,
            7.765049321445005871323047E+03,
            1.331903827966074194402448E+05,
            1.136705821321969608938755E+06,
            5.267964117437946917577538E+06,
            1.346701454311101692290052E+07,
            1.782736530353274213975932E+07,
            9.533095591844353613395747E+06
        };
        double[] q4 =
        {
            2.690530175870899333379843E+03,
            6.393885654300092398984238E+05,
            4.135599930241388052042842E+07,
            1.120872109616147941376570E+09,
            1.488613728678813811542398E+10,
            1.016803586272438228077304E+11,
            3.417476345507377132798597E+11,
            4.463158187419713286462081E+11
        };
        double res;
        const double sqrtpi = 0.9189385332046727417803297;
        const double xbig = 2.55E+305;
        const double xinf = 1.79E+308;

        double y = x;

        switch (y)
        {
            case > 0.0 and <= xbig:
            {
                if (y <= r8_epsilon())
                {
                    res = -Math.Log(y);
                }
                //
                //  EPS < X <= 1.5.
                //
                else
                {
                    double corr;
                    double xden;
                    double xm2;
                    double xnum;
                    int i;
                    switch (y)
                    {
                        case <= 1.5:
                        {
                            double xm1;
                            switch (y)
                            {
                                case < 0.6796875:
                                    corr = -Math.Log(y);
                                    xm1 = y;
                                    break;
                                default:
                                    corr = 0.0;
                                    xm1 = y - 0.5 - 0.5;
                                    break;
                            }

                            switch (y)
                            {
                                case <= 0.5:
                                case >= 0.6796875:
                                {
                                    xden = 1.0;
                                    xnum = 0.0;
                                    for (i = 0; i < 8; i++)
                                    {
                                        xnum = xnum * xm1 + p1[i];
                                        xden = xden * xm1 + q1[i];
                                    }

                                    res = corr + xm1 * (d1 + xm1 * (xnum / xden));
                                    break;
                                }
                                default:
                                {
                                    xm2 = y - 0.5 - 0.5;
                                    xden = 1.0;
                                    xnum = 0.0;
                                    for (i = 0; i < 8; i++)
                                    {
                                        xnum = xnum * xm2 + p2[i];
                                        xden = xden * xm2 + q2[i];
                                    }

                                    res = corr + xm2 * (d2 + xm2 * (xnum / xden));
                                    break;
                                }
                            }

                            break;
                        }
                        //
                        //  1.5 < X <= 4.0.
                        //
                        case <= 4.0:
                        {
                            xm2 = y - 2.0;
                            xden = 1.0;
                            xnum = 0.0;
                            for (i = 0; i < 8; i++)
                            {
                                xnum = xnum * xm2 + p2[i];
                                xden = xden * xm2 + q2[i];
                            }

                            res = xm2 * (d2 + xm2 * (xnum / xden));
                            break;
                        }
                        //
                        //  4.0 < X <= 12.0.
                        //
                        case <= 12.0:
                        {
                            double xm4 = y - 4.0;
                            xden = -1.0;
                            xnum = 0.0;
                            for (i = 0; i < 8; i++)
                            {
                                xnum = xnum * xm4 + p4[i];
                                xden = xden * xm4 + q4[i];
                            }

                            res = d4 + xm4 * (xnum / xden);
                            break;
                        }
                        //
                        default:
                        {
                            res = 0.0;

                            if (y <= frtbig)
                            {
                                res = c[6];
                                double ysq = y * y;
                                for (i = 0; i < 6; i++)
                                {
                                    res = res / ysq + c[i];
                                }
                            }

                            res /= y;
                            corr = Math.Log(y);
                            res = res + sqrtpi - 0.5 * corr;
                            res += y * (corr - 1.0);
                            break;
                        }
                    }
                }

                break;
            }
            //
            default:
                res = xinf;
                break;
        }

        //
        //  Final adjustments and return.
        //
        return res;
    }
        
    public static double r8_gamma_log_int(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_LOG_INT computes the logarithm of Gamma of an integer N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument of the logarithm of the Gamma function.
        //    0 < N.
        //
        //    Output, double R8_GAMMA_LOG_INT, the logarithm of
        //    the Gamma function of N.
        //
    {
        switch (n)
        {
            case <= 0:
                Console.WriteLine(" ");
                Console.WriteLine("R8_GAMMA_LOG_INT - Fatal error!");
                Console.WriteLine("  Illegal input value of N = " + n + "");
                Console.WriteLine("  But N must be strictly positive.");
                return 1;
            default:
            {
                double value = Helpers.LogGamma(n);

                return value;
            }
        }
    }
}