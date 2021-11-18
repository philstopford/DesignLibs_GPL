using System;

namespace Burkardt.PDFLib;

public static partial class PDF
{

    public static double i4_binomial_pdf(int n, double p, int k)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BINOMIAL_PDF evaluates the binomial PDF.
        //
        //  Discussion:
        //
        //    pdf(n,p,k) = C(n,k) p^k (1-p)^(n-k)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2018
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the number of binomial trials.
        //    0 < N.
        //
        //    Input, double P, the probability of a success in one trial.
        //
        //    Input, int K, the number of successes.
        //
        //    Output, double I4_BINOMIAL_PDF, the probability of K successes
        //    in N trials with a per-trial success probability of P.
        //
    {
        double value = 0;

        switch (k)
        {
            case < 0:
                value = 0.0;
                break;
            default:
            {
                if (k <= n)
                {
                    value = r8_choose(n, k) * Math.Pow(p, k) * Math.Pow(1.0 - p, n - k);
                }
                else
                {
                    value = 0.0;
                }

                break;
            }
        }

        return value;
    }

    public static int i4_binomial_sample(int n, double pp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BINOMIAL_SAMPLE generates a binomial random deviate.
        //
        //  Discussion:
        //
        //    This procedure generates a single random deviate from a binomial
        //    distribution whose number of trials is N and whose
        //    probability of an event in each trial is P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Barry Brown, James Lovato.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Voratas Kachitvichyanukul, Bruce Schmeiser,
        //    Binomial Random Variate Generation,
        //    Communications of the ACM,
        //    Volume 31, Number 2, February 1988, pages 216-222.
        //
        //  Parameters:
        //
        //    Input, int N, the number of binomial trials, from which a
        //    random deviate will be generated.
        //    0 < N.
        //
        //    Input, double PP, the probability of an event in each trial of
        //    the binomial distribution from which a random deviate is to be generated.
        //    0.0 < PP < 1.0.
        //
        //    Output, int I4_BINOMIAL_SAMPLE, a random deviate from the
        //    distribution.
        //
    {
        double al;
        double alv;
        double amaxp;
        double c;
        double f;
        double f1;
        double f2;
        double ffm;
        double fm;
        double g;
        int i;
        int ix = 0;
        int ix1;
        int k;
        int m;
        int mp;
        double p;
        double p1;
        double p2;
        double p3;
        double p4;
        double q;
        double qn;
        double r;
        double t;
        double u;
        double v;
        int value;
        double w;
        double w2;
        double x;
        double x1;
        double x2;
        double xl;
        double xll;
        double xlr;
        double xm;
        double xnp;
        double xnpq;
        double xr;
        double ynorm;
        double z;
        double z2;

        switch (pp)
        {
            case <= 0.0:
            case >= 1.0:
                Console.WriteLine("");
                Console.WriteLine("I4_BINOMIAL_SAMPLE - Fatal error!");
                Console.WriteLine("  PP is out of range.");
                return 1;
        }

        p = Math.Min(pp, 1.0 - pp);
        q = 1.0 - p;
        xnp = n * p;

        switch (xnp)
        {
            case < 30.0:
            {
                qn = Math.Pow(q, n);
                r = p / q;
                g = r * (n + 1);

                for (;;)
                {
                    ix = 0;
                    f = qn;
                    u = r8_uniform_01_sample();

                    for (;;)
                    {
                        if (u < f)
                        {
                            ix = pp switch
                            {
                                > 0.5 => n - ix,
                                _ => ix
                            };

                            value = ix;
                            return value;
                        }

                        if (110 < ix)
                        {
                            break;
                        }

                        u -= f;
                        ix += 1;
                        f *= (g / ix - r);
                    }
                }
            }
        }

        ffm = xnp + p;
        m = (int) ffm;
        fm = m;
        xnpq = xnp * q;
        p1 = (int) (2.195 * Math.Sqrt(xnpq) - 4.6 * q) + 0.5;
        xm = fm + 0.5;
        xl = xm - p1;
        xr = xm + p1;
        c = 0.134 + 20.5 / (15.3 + fm);
        al = (ffm - xl) / (ffm - xl * p);
        xll = al * (1.0 + 0.5 * al);
        al = (xr - ffm) / (xr * q);
        xlr = al * (1.0 + 0.5 * al);
        p2 = p1 * (1.0 + c + c);
        p3 = p2 + c / xll;
        p4 = p3 + c / xlr;
        //
        //  Generate a variate.
        //
        for (;;)
        {
            u = r8_uniform_01_sample() * p4;
            v = r8_uniform_01_sample();
            //
            //  Triangle
            //
            if (u < p1)
            {
                ix = pp switch
                {
                    > 0.5 => n - ix,
                    _ => (int) (xm - p1 * v + u)
                };

                value = ix;
                return value;
            }

            //
            //  Parallelogram
            //
            if (u <= p2)
            {
                x = xl + (u - p1) / c;
                v = v * c + 1.0 - Math.Abs(xm - x) / p1;

                switch (v)
                {
                    case <= 0.0:
                    case > 1.0:
                        continue;
                    default:
                        ix = (int) x;
                        break;
                }
            }
            else if (u <= p3)
            {
                ix = (int) (xl + Math.Log(v) / xll);
                switch (ix)
                {
                    case < 0:
                        continue;
                    default:
                        v = v * (u - p2) * xll;
                        break;
                }
            }
            else
            {
                ix = (int) (xr - Math.Log(v) / xlr);
                if (n < ix)
                {
                    continue;
                }

                v = v * (u - p3) * xlr;
            }

            k = Math.Abs(ix - m);

            if (k <= 20 || xnpq / 2.0 - 1.0 <= k)
            {
                f = 1.0;
                r = p / q;
                g = (n + 1) * r;

                if (m < ix)
                {
                    mp = m + 1;
                    for (i = mp; i <= ix; i++)
                    {
                        f *= (g / i - r);
                    }
                }
                else if (ix < m)
                {
                    ix1 = ix + 1;
                    for (i = ix1; i <= m; i++)
                    {
                        f /= (g / i - r);
                    }
                }

                if (v <= f)
                {
                    ix = pp switch
                    {
                        > 0.5 => n - ix,
                        _ => ix
                    };

                    value = ix;
                    return value;
                }
            }
            else
            {
                amaxp = k / xnpq * ((k * (k / 3.0
                                          + 0.625) + 0.1666666666666) / xnpq + 0.5);
                ynorm = -(double) (k * k) / (2.0 * xnpq);
                alv = Math.Log(v);

                if (alv < ynorm - amaxp)
                {
                    ix = pp switch
                    {
                        > 0.5 => n - ix,
                        _ => ix
                    };

                    value = ix;
                    return value;
                }

                if (ynorm + amaxp < alv)
                {
                    continue;
                }

                x1 = ix + 1;
                f1 = fm + 1.0;
                z = n + 1 - fm;
                w = n - ix + 1;
                z2 = z * z;
                x2 = x1 * x1;
                f2 = f1 * f1;
                w2 = w * w;

                t = xm * Math.Log(f1 / x1) + (n - m + 0.5) * Math.Log(z / w)
                                           + (ix - m) * Math.Log(w * p / (x1 * q))
                                           + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0
                                               / f2) / f2) / f2) / f2) / f1 / 166320.0
                                           + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0
                                               / z2) / z2) / z2) / z2) / z / 166320.0
                                           + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0
                                               / x2) / x2) / x2) / x2) / x1 / 166320.0
                                           + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0
                                               / w2) / w2) / w2) / w2) / w / 166320.0;

                if (alv <= t)
                {
                    ix = pp switch
                    {
                        > 0.5 => n - ix,
                        _ => ix
                    };

                    value = ix;
                    return value;
                }
            }
        }
    }

    public static int i4_uniform_sample(int a, int b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_UNIFORM_SAMPLE returns a scaled pseudorandom I4 between A and B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, B, the limits of the interval.
        //
        //    Output, int I4_UNIFORM_AB, a number between A and B.
        //
    {
        int a2;
        int b2;
        double u;
        int value;
        //
        //  We prefer A < B.
        //
        a2 = Math.Min(a, b);
        b2 = Math.Max(a, b);

        u = r8_uniform_01_sample();
        //
        //  Scale to [A2-0.5,B2+0.5].
        //
        u = (1.0 - u) * (a2 - 0.5)
            + u * (b2 + 0.5);
        //
        //  Round.
        //
        value = (int) Math.Round(u);
        //
        //  Enforce limits.
        //
        if (value < a2)
        {
            value = a2;
        }

        if (b2 < value)
        {
            value = b2;
        }

        return value;
    }

    public static double i4vec_multinomial_pdf(int n, double[] p, int m, int[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_MULTINOMIAL_PDF evaluates the multinomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 June 2013
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the number of trials.
        //
        //    Input, double P[M], the probability of each outcome
        //    on any single trial.
        //
        //    Input, int M, the number of possible outcomes
        //    of a single trial.
        //
        //    Input, int X[M], the results of N trials,
        //    with X(I) the number of times outcome I occurred.
        //
        //    Output, double I4VEC_MULTINOMIAL_PDF, the probability
        //    density function evaluated at X.
        //
    {
        int bot;
        int c;
        int i;
        int j;
        double pdf;
        int top;
        //
        //  The combinatorial coefficient is an integer.
        //
        c = 1;
        top = n;
        for (i = 0; i < m; i++)
        {
            bot = 1;
            for (j = 0; j < x[i]; j++)
            {
                c = c * top / bot;
                top -= 1;
                bot += 1;
            }
        }

        pdf = c;
        for (i = 0; i < m; i++)
        {
            pdf *= Math.Pow(p[i], x[i]);
        }

        return pdf;
    }

    public static int[] i4vec_multinomial_sample(int n, double[] p, int ncat)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_MULTINOMIAL_SAMPLE generates a multinomial random deviate.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Barry Brown, James Lovato.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Luc Devroye,
        //    Non-Uniform Random Variate Generation,
        //    Springer, 1986,
        //    ISBN: 0387963057,
        //    LC: QA274.D48.
        //
        //  Parameters:
        //
        //    Input, int N, the number of events, which will be
        //    classified into one of the NCAT categories.
        //
        //    Input, double P[NCAT-1].  P(I) is the probability that an event
        //    will be classified into category I.  Thus, each P(I) must be between 
        //    0.0 and 1.0.  Only the first NCAT-1 values of P must be defined since 
        //    P(NCAT) would be 1.0 minus the sum of the first NCAT-1 P's.
        //
        //    Input, int NCAT, the number of categories.
        //
        //    Output, int I4VEC_MULTINOMIAL_SAMPLE[NCAT], a random observation from 
        //    the multinomial distribution.  All IX(i) will be nonnegative and their 
        //    sum will be N.
        //
    {
        int i;
        int icat;
        int[] ix;
        int ntot;
        double prob;
        double ptot;

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("I4VEC_MULTINOMIAL_SAMPLE - Fatal error!");
                Console.WriteLine("  N < 0");
                return new int[1];
        }

        switch (ncat)
        {
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("I4VEC_MULTINOMIAL_SAMPLE - Fatal error!");
                Console.WriteLine("  NCAT <= 1");
                return new int[1];
        }

        for (i = 0; i < ncat - 1; i++)
        {
            switch (p[i])
            {
                case < 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("I4VEC_MULTINOMIAL_SAMPLE - Fatal error!");
                    Console.WriteLine("  Some P(i) < 0.");
                    return new int[1];
                case > 1.0:
                    Console.WriteLine("");
                    Console.WriteLine("I4VEC_MULTINOMIAL_SAMPLE - Fatal error!");
                    Console.WriteLine("  Some 1 < P(i).");
                    return new int[1];
            }
        }

        ptot = 0.0;
        for (i = 0; i < ncat - 1; i++)
        {
            ptot += p[i];
        }

        switch (ptot)
        {
            case > 0.99999:
                Console.WriteLine("");
                Console.WriteLine("I4VEC_MULTINOMIAL_SAMPLE - Fatal error!");
                Console.WriteLine("  1.0 < Sum of P().");
                return new int[1];
        }

        //
        //  Initialize variables.
        //
        ntot = n;
        ptot = 1.0;

        ix = new int[ncat];
        for (i = 0; i < ncat; i++)
        {
            ix[i] = 0;
        }

        //
        //  Generate the observation.
        //
        for (icat = 0; icat < ncat - 1; icat++)
        {
            prob = p[icat] / ptot;
            ix[icat] = i4_binomial_sample(ntot, prob);
            ntot -= ix[icat];
            switch (ntot)
            {
                case <= 0:
                    return ix;
                default:
                    ptot -= p[icat];
                    break;
            }
        }

        ix[ncat - 1] = ntot;

        return ix;
    }


}