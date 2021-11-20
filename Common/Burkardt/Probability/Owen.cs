using System;

namespace Burkardt.Probability;

public static class Owen
{
    public static void owen_values(ref int n_data, ref double h, ref double a, ref double t )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OWEN_VALUES returns some values of Owen's T function.
        //
        //  Discussion:
        //
        //    Owen's T function is useful for computation of the bivariate normal
        //    distribution and the distribution of a skewed normal distribution.
        //
        //    Although it was originally formulated in terms of the bivariate
        //    normal function, the function can be defined more directly as
        //
        //      T(H,A) = 1 / ( 2 * Math.PI ) *
        //        Integral ( 0 <= X <= A ) e^(H^2*(1+X^2)/2) / (1+X^2) dX
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Mike Patefield, David Tandy,
        //    Fast and Accurate Calculation of Owen's T Function,
        //    Journal of Statistical Software,
        //    Volume 5, Number 5, 2000, pages 1-25.
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &H, a parameter.
        //
        //    Output, double &A, the upper limit of the integral.
        //
        //    Output, double &T, the value of the function.
        //
    {
        const int N_MAX = 28;

        double[] a_vec =
        {
            0.2500000000000000E+00,
            0.4375000000000000E+00,
            0.9687500000000000E+00,
            0.0625000000000000E+00,
            0.5000000000000000E+00,
            0.9999975000000000E+00,
            0.5000000000000000E+00,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.5000000000000000E+00,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.5000000000000000E+00,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.5000000000000000E+00,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.5000000000000000E+00,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.1000000000000000E+02,
            0.1000000000000000E+03
        };

        double[] h_vec =
        {
            0.0625000000000000E+00,
            6.5000000000000000E+00,
            7.0000000000000000E+00,
            4.7812500000000000E+00,
            2.0000000000000000E+00,
            1.0000000000000000E+00,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.2500000000000000E+00,
            0.2500000000000000E+00,
            0.2500000000000000E+00,
            0.2500000000000000E+00,
            0.1250000000000000E+00,
            0.1250000000000000E+00,
            0.1250000000000000E+00,
            0.1250000000000000E+00,
            0.7812500000000000E-02,
            0.7812500000000000E-02,
            0.7812500000000000E-02,
            0.7812500000000000E-02,
            0.7812500000000000E-02,
            0.7812500000000000E-02
        };

        double[] t_vec =
        {
            3.8911930234701366E-02,
            2.0005773048508315E-11,
            6.3990627193898685E-13,
            1.0632974804687463E-07,
            8.6250779855215071E-03,
            6.6741808978228592E-02,
            0.4306469112078537E-01,
            0.6674188216570097E-01,
            0.7846818699308410E-01,
            0.7929950474887259E-01,
            0.6448860284750376E-01,
            0.1066710629614485E+00,
            0.1415806036539784E+00,
            0.1510840430760184E+00,
            0.7134663382271778E-01,
            0.1201285306350883E+00,
            0.1666128410939293E+00,
            0.1847501847929859E+00,
            0.7317273327500385E-01,
            0.1237630544953746E+00,
            0.1737438887583106E+00,
            0.1951190307092811E+00,
            0.7378938035365546E-01,
            0.1249951430754052E+00,
            0.1761984774738108E+00,
            0.1987772386442824E+00,
            0.2340886964802671E+00,
            0.2479460829231492E+00
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            h = 0.0;
            a = 0.0;
            t = 0.0;
        }
        else
        {
            h = h_vec[n_data - 1];
            a = a_vec[n_data - 1];
            t = t_vec[n_data - 1];
        }
    }

    public static double tfn(double h, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TFN calculates the T function of Owen.
        //
        //  Discussion:
        //
        //    Owen's T function is useful for computation of the bivariate normal
        //    distribution and the distribution of a skewed normal distribution.
        //
        //    Although it was originally formulated in terms of the bivariate
        //    normal function, the function can be defined more directly as
        //
        //      T(H,A) = 1 / ( 2 * PI ) *
        //        Integral ( 0 <= X <= A ) e^( -H^2 * (1+X^2) / 2 ) / (1+X^2) dX
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 December 2004
        //
        //  Author:
        //
        //    FORTRAN77 original version by J C Young, C E Minder.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    D B Owen,
        //    Tables for computing the bivariate normal distribution,
        //    Annals of Mathematical Statistics,
        //    Volume 27, pages 1075-1090, 1956.
        //
        //    J C Young, C E Minder,
        //    Algorithm AS 76,
        //    An Algorithm Useful in Calculating Non-Central T and
        //      Bivariate Normal Distributions,
        //    Applied Statistics,
        //    Volume 23, Number 3, 1974, pages 455-457.
        //
        //  Parameters:
        //
        //    Input, double H, A, the arguments of the T function.
        //
        //    Output, double TFN, the value of the T function.
        //
    {
        const int NGAUSS = 10;

        const double two_pi_inverse = 0.1591549430918953;
        const double tv1 = 1.0E-35;
        const double tv2 = 15.0;
        const double tv3 = 15.0;
        const double tv4 = 1.0E-05;
        double value;
        double[] weight =  {
            0.666713443086881375935688098933E-01,
            0.149451349150580593145776339658,
            0.219086362515982043995534934228,
            0.269266719309996355091226921569,
            0.295524224714752870173892994651,
            0.295524224714752870173892994651,
            0.269266719309996355091226921569,
            0.219086362515982043995534934228,
            0.149451349150580593145776339658,
            0.666713443086881375935688098933E-01
        };
        double[] xtab =  {
            -0.973906528517171720077964012084,
            -0.865063366688984510732096688423,
            -0.679409568299024406234327365115,
            -0.433395394129247190799265943166,
            -0.148874338981631210884826001130,
            0.148874338981631210884826001130,
            0.433395394129247190799265943166,
            0.679409568299024406234327365115,
            0.865063366688984510732096688423,
            0.973906528517171720077964012084
        };
        switch (Math.Abs(h))
        {
            //
            //  Test for H near zero.
            //
            case < tv1:
                value = Math.Atan(a) * two_pi_inverse;
                break;
            //
            //  Test for large values of abs(H).
            //
            case > tv2:
                value = 0.0;
                break;
            //
            default:
            {
                if (Math.Abs(a) < tv1)
                {
                    value = 0.0;
                }
                //
                //  Test whether abs(A) is so large that it must be truncated.
                //  If so, the truncated value is H2.
                //
                else
                {
                    double hs = -0.5 * h * h;
                    double h2 = a;
                    double as_ = a * a;
                    //
                    //  Computation of truncation point by Newton iteration.
                    //
                    double rt;
                    if (tv3 <= Math.Log(1.0 + as_) - hs * as_)
                    {
                        double h1 = 0.5 * a;
                        as_ = 0.25 * as_;

                        for (;;)
                        {
                            rt =  as_ +1.0;
                            h2 = h1 + (hs * as_ +tv3 - Math.Log(rt) )
                                / (2.0 * h1 * (1.0 / rt - hs));
                            as_ = h2* h2;

                            if (Math.Abs(h2 - h1) < tv4)
                            {
                                break;
                            }

                            h1 = h2;
                        }
                    }

                    //
                    //  Gaussian quadrature on the interval [0,H2].
                    //
                    rt = 0.0;
                    int i;
                    for (i = 0; i < NGAUSS; i++)
                    {
                        double x = 0.5 * h2 * (xtab[i] + 1.0);
                        rt += weight[i] * Math.Exp(hs * (1.0 + x * x)) / (1.0 + x * x);
                    }

                    value = rt * 0.5 * h2 * two_pi_inverse;
                }

                break;
            }
        }

        return value;

    }

}