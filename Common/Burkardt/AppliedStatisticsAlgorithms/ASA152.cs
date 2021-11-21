using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static double chyper(bool point, int kk, int ll, int mm, int nn, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHYPER computes point or cumulative hypergeometric probabilities.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Richard Lund.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    PR Freeman,
        //    Algorithm AS 59:
        //    Hypergeometric Probabilities,
        //    Applied Statistics,
        //    Volume 22, Number 1, 1973, pages 130-133.
        //
        //    Richard Lund,
        //    Algorithm AS 152:
        //    Cumulative hypergeometric probabilities,
        //    Applied Statistics,
        //    Volume 29, Number 2, 1980, pages 221-223.
        //
        //    BL Shea,
        //    Remark AS R77:
        //    A Remark on Algorithm AS 152: Cumulative hypergeometric probabilities,
        //    Applied Statistics,
        //    Volume 38, Number 1, 1989, pages 199-204.
        //
        //  Parameters:
        //
        //    Input, bool POINT, is TRUE if the point probability is desired,
        //    and FALSE if the cumulative probability is desired.
        //
        //    Input, int KK, the sample size.
        //    0 <= KK <= MM.
        //
        //    Input, int LL, the number of successes in the sample.
        //    0 <= LL <= KK.
        //
        //    Input, int MM, the population size that was sampled.
        //    0 <= MM.
        //
        //    Input, int NN, the number of "successes" in the population.
        //    0 <= NN <= MM.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no error occurred.
        //    nonzero, an error occurred.
        //
        //    Output, double CHYPER, the PDF (point probability) of
        //    exactly LL successes out of KK samples, or the CDF (cumulative
        //    probability) of up to LL successes out of KK samples.
        //
    {
        const double elimit = -88.0;
        const int mbig = 600;
        const int mvbig = 1000;
        const double rootpi = 2.506628274631001;
        const double scale = 1.0E+35;

        ifault = 0;

        int k = kk + 1;
        int l = ll + 1;
        int m = mm + 1;
        int n = nn + 1;

        bool dir = true;
        //
        //  Check arguments are within permitted limits.
        //
        double value = 0.0;

        if (n < 1 || m < n || k < 1 || m < k)
        {
            ifault = 1;
            return value;
        }

        if (l < 1 || m - n < k - l)
        {
            ifault = 2;
            return value;
        }

        value = point switch
        {
            false => 1.0,
            _ => value
        };

        if (n < l || k < l)
        {
            ifault = 2;
            return value;
        }

        ifault = 0;
        value = 1.0;

        if (k == 1 || k == m || n == 1 || n == m)
        {
            return value;
        }

        switch (point)
        {
            case false when ll == Math.Min(kk, nn):
                return value;
        }

        double p = nn / (double) (mm - nn);

        if (16.0 * Math.Max(p, 1.0 / p)
            < Math.Min(kk, mm - kk) &&
            mvbig < mm && -100.0 < elimit)
        {
            //
            //  Use a normal approximation.
            //
            double mean = kk * nn / (double) mm;

            double sig = Math.Sqrt(mean * ((mm - nn) / (double) mm)
                                        * ((mm - kk) / (double) (mm - 1)));

            switch (point)
            {
                case true:
                {
                    double arg = -0.5 * Math.Pow((ll - mean) / sig, 2);
                    if (elimit <= arg)
                    {
                        value = Math.Exp(arg) / (sig * rootpi);
                    }
                    else
                    {
                        value = 0.0;
                    }

                    break;
                }
                default:
                    value = alnorm((ll + 0.5 - mean) / sig, false);
                    break;
            }
        }
        else
        {
            //
            //  Calculate exact hypergeometric probabilities.
            //  Interchange K and N if this saves calculations.
            //
            int i;
            if (Math.Min(n - 1, m - n) < Math.Min(k - 1, m - k))
            {
                i = k;
                k = n;
                n = i;
            }

            if (m - k < k - 1)
            {
                dir = !dir;
                l = n - l + 1;
                k = m - k + 1;
            }

            int j;
            if (mbig < mm)
            {
                //
                //  Take logarithms of factorials.
                //
                p = alnfac(nn)
                    - alnfac(mm)
                    + alnfac(mm - kk)
                    + alnfac(kk)
                    + alnfac(mm - nn)
                    - alnfac(ll)
                    - alnfac(nn - ll)
                    - alnfac(kk - ll)
                    - alnfac(mm - nn - kk + ll);

                value = elimit <= p ? Math.Exp(p) : 0.0;
            }
            else
            {
                //
                //  Use Freeman/Lund algorithm.
                //
                for (i = 1; i <= l - 1; i++)
                {
                    value = value * ((k - i) * (n - i))
                            / ((l - i) * (m - i));
                }

                if (l != k)
                {
                    j = m - n + l;
                    for (i = l; i <= k - 1; i++)
                    {
                        value = value * (j - i) / (m - i);
                    }
                }
            }

            switch (point)
            {
                case true:
                    return value;
            }

            switch (value)
            {
                case 0.0:
                {
                    //
                    //  We must recompute the point probability since it has underflowed.
                    //
                    if (mm <= mbig)
                    {
                        p = alnfac(nn)
                            - alnfac(mm)
                            + alnfac(kk)
                            + alnfac(mm - nn)
                            - alnfac(ll)
                            - alnfac(nn - ll)
                            - alnfac(kk - ll)
                            - alnfac(mm - nn - kk + ll)
                            + alnfac(mm - kk);
                    }

                    p += Math.Log(scale);

                    if (p < elimit)
                    {
                        ifault = 3;
                        if ((nn * kk + nn + kk + 1)
                            / (double) (mm + 2) < ll)
                        {
                            value = 1.0;
                        }

                        return value;
                    }

                    p = Math.Exp(p);
                    break;
                }
                //
                default:
                    p = value * scale;
                    break;
            }

            double pt = 0.0;
            int nl = n - l;
            int kl = k - l;
            int mnkl = m - n - kl + 1;

            if (l <= kl)
            {
                for (i = 1; i <= l - 1; i++)
                {
                    p = p * ((l - i) * (mnkl - i)) /
                        ((nl + i) * (kl + i));
                    pt += p;
                }
            }
            else
            {
                dir = !dir;
                for (j = 0; j <= kl - 1; j++)
                {
                    p = p * ((nl - j) * (kl - j))
                        / ((l + j) * (mnkl + j));
                    pt += p;
                }
            }

            ifault = p switch
            {
                0.0 => 3,
                _ => ifault
            };

            switch (dir)
            {
                case true:
                    value += pt / scale;
                    break;
                default:
                    value = 1.0 - pt / scale;
                    break;
            }
        }

        return value;
    }

    public static void hypergeometric_cdf_values(ref int n_data, ref int sam, ref int suc, ref int pop,
            ref int n, ref double fx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
        //
        //  Discussion:
        //
        //    CDF(X)(A,B) is the probability of at most X successes in A trials,
        //    given that the probability of success on a single trial is B.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`DiscreteDistributions`]
        //      dist = HypergeometricDistribution [ sam, suc, pop ]
        //      CDF [ dist, n ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2004
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
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //    Daniel Zwillinger,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition, CRC Press, 1996, pages 651-652.
        //
        //  Parameters:
        //
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int *SAM, int *SUC, int *POP, the sample size, 
        //    success size, and population parameters of the function.
        //
        //    Output, int *N, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 16;

        double[] fx_vec =  {
                0.6001858177500578E-01,
                0.2615284665839845E+00,
                0.6695237889132748E+00,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.5332595856827856E+00,
                0.1819495964117640E+00,
                0.4448047017527730E-01,
                0.9999991751316731E+00,
                0.9926860896560750E+00,
                0.8410799901444538E+00,
                0.3459800113391901E+00,
                0.0000000000000000E+00,
                0.2088888139634505E-02,
                0.3876752992448843E+00,
                0.9135215248834896E+00
            }
            ;

        int[] n_vec =  {
                7, 8, 9, 10,
                6, 6, 6, 6,
                6, 6, 6, 6,
                0, 0, 0, 0
            }
            ;

        int[] pop_vec =  {
                100, 100, 100, 100,
                100, 100, 100, 100,
                100, 100, 100, 100,
                90, 200, 1000, 10000
            }
            ;

        int[] sam_vec =  {
                10, 10, 10, 10,
                6, 7, 8, 9,
                10, 10, 10, 10,
                10, 10, 10, 10
            }
            ;

        int[] suc_vec =  {
                90, 90, 90, 90,
                90, 90, 90, 90,
                10, 30, 50, 70,
                90, 90, 90, 90
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
            sam = 0;
            suc = 0;
            pop = 0;
            n = 0;
            fx = 0.0;
        }
        else
        {
            sam = sam_vec[n_data - 1];
            suc = suc_vec[n_data - 1];
            pop = pop_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static void hypergeometric_pdf_values(ref int n_data, ref int sam, ref int suc, ref int pop,
            ref int n, ref double fx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
        //
        //  Discussion:
        //
        //    CDF(X)(A,B) is the probability of X successes in A trials,
        //    given that the probability of success on a single trial is B.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      dist = HypergeometricDistribution [ sam, suc, pop ]
        //      PDF [ dist, n ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 January 2008
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
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //    Daniel Zwillinger,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition, CRC Press, 1996, pages 651-652.
        //
        //  Parameters:
        //
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int *SAM, int *SUC, int *POP, the sample size, 
        //    success size, and population parameters of the function.
        //
        //    Output, int *N, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 16;

        double[] fx_vec =  {
                0.05179370533242827E+00,
                0.2015098848089788E+00,
                0.4079953223292903E+00,
                0.3304762110867252E+00,
                0.5223047493549780E+00,
                0.3889503452643453E+00,
                0.1505614239732950E+00,
                0.03927689321042477E+00,
                0.00003099828465518108E+00,
                0.03145116093938197E+00,
                0.2114132170316862E+00,
                0.2075776621999210E+00,
                0.0000000000000000E+00,
                0.002088888139634505E+00,
                0.3876752992448843E+00,
                0.9135215248834896E+00
            }
            ;

        int[] n_vec =  {
                7, 8, 9, 10,
                6, 6, 6, 6,
                6, 6, 6, 6,
                0, 0, 0, 0
            }
            ;

        int[] pop_vec =  {
                100, 100, 100, 100,
                100, 100, 100, 100,
                100, 100, 100, 100,
                90, 200, 1000, 10000
            }
            ;

        int[] sam_vec =  {
                10, 10, 10, 10,
                6, 7, 8, 9,
                10, 10, 10, 10,
                10, 10, 10, 10
            }
            ;

        int[] suc_vec =  {
                90, 90, 90, 90,
                90, 90, 90, 90,
                10, 30, 50, 70,
                90, 90, 90, 90
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
            sam = 0;
            suc = 0;
            pop = 0;
            n = 0;
            fx = 0.0;
        }
        else
        {
            sam = sam_vec[n_data - 1];
            suc = suc_vec[n_data - 1];
            pop = pop_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}