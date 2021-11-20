using System;

namespace Burkardt.Probability;

public static class Ribesl
{
    public static int ribesl(double x, double alpha, int nb, int ize, ref double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RIBESL calculates I Bessel function with non-integer orders.
        //
        //  Discussion:
        //
        //    This routine calculates Bessel functions I SUB(N+ALPHA) (X)
        //    for non-negative argument X, and non-negative order N+ALPHA,
        //    with or without exponential scaling.
        //
        //    This program is based on a program written by David
        //    Sookne that computes values of the Bessel functions J or
        //    I of real argument and integer order.  Modifications include
        //    the restriction of the computation to the I Bessel function
        //    of non-negative real argument, the extension of the computation
        //    to arbitrary positive order, the inclusion of optional
        //    exponential scaling, and the elimination of most underflow.
        //
        //    In case of an error, NCALC will not equal NB, and not all I's are
        //    calculated to the desired accuracy.
        //
        //    If NCALC < 0:  An argument is out of range. For example,
        //    NB <= 0, IZE is not 1 or 2, or IZE = 1 and EXPARG <= ABS(X)
        //    In this case, the B-vector is not calculated, and NCALC is
        //    set to MIN(NB,0)-1 so that NCALC /= NB.
        //
        //    If 0 < NCALC < NB, then not all requested function values could
        //    be calculated accurately.  This usually occurs because NB is
        //    much larger than ABS(X).  In this case, B(N) is calculated
        //    to the desired accuracy for N <= NCALC, but precision
        //    is lost for NCALC < N <= NB.  If B(N) does not vanish
        //    for NCALC < N (because it is too small to be represented),
        //    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
        //    significant figures of B(N) can be trusted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 March 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Frank Olver, David Sookne,
        //    A Note on Backward Recurrence Algorithms,
        //    Mathematics of Computation,
        //    Volume 26, 1972, pages 941-947.
        //
        //    David Sookne,
        //    Bessel Functions of Real Argument and Integer Order,
        //    NBS Journal of Research B,
        //    Volume 77B, 1973, pages 125-132.
        //
        //    William Cody,
        //    Algorithm 597:
        //    Sequence of Modified Bessel Functions of the First Kind,
        //    ACM Transactions of Mathematical Software,
        //    Volume 9, Number 2, June 1983, pages 242-245.
        //
        //  Parameters:
        //
        //    Input, double X, the argument for which the functions
        //    are to be calculated.
        //
        //    Input, double ALPHA,the fractional part of the order
        //    for which the functions are to be calculated.
        //    0 <= ALPHA < 1.0.
        //
        //    Input, int NB, the number of functions to be calculated.
        //    The first function calculated is of order ALPHA, and the
        //    last is of order (NB - 1 + ALPHA).  1 <= NB.
        //
        //    Input, int IZE, scaling option.
        //    1, unscaled I's are to calculated,
        //    2, exponentially scaled I's are to be calculated.
        //
        //    Output, double B[NB], the values of the functions
        //    I(ALPHA,X) through I(NB-1+ALPHA,X), with scaling if requested.
        //
        //    Output, int RIBESL, the value of NCALC, the error indicator.
        //    If NCALC = NB, then all the requested values were calculated
        //    to the desired accuracy.
        //
        //  Local Parameeters:
        //
        //    BETA, the radix for the floating-point system.
        //
        //    MINEXP, smallest representable power of BETA.
        //
        //    MAXEXP, smallest power of BETA that overflows
        //
        //    IT, number of bits in the mantissa of a working precision variable.
        //
        //    NSIG, decimal significance desired.  Should be set to
        //    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
        //    in decreased accuracy while setting NSIG higher will
        //    increase CPU time without increasing accuracy.  The
        //    truncation error is limited to a relative error of
        //    T=.5*10^(-NSIG).
        //
        //    ENTEN, 10.0^K, where K is the largest integer such that
        //    ENTEN is machine-representable in working precision
        //
        //    ENSIG, 10.0^NSIG
        //
        //    RTNSIG, 10.0^(-K) for the smallest integer K such that
        //    NSIG/4 <= K.
        //
        //    ENMTEN, smallest ABS(X) such that X/4 does not underflow
        //
        //    XLARGE, upper limit on the magnitude of X when IZE=2.  Bear
        //    in mind that if ABS(X)=N, then at least N iterations
        //    of the backward recursion will be executed.  The value
        //    of 10.0^4 is used on every machine.
        //
        //    EXPARG, largest working precision argument that the library
        //    EXP routine can handle and upper limit on the
        //    magnitude of X when IZE=1; approximately log(BETA^MAXEXP).
        //
        //    Approximate values for some important machines are:
        //
        //                        beta       minexp      maxexp       it
        //
        //  CRAY-1        (S.P.)    2        -8193        8191        48
        //  Cyber 180/855
        //    under NOS   (S.P.)    2         -975        1070        48
        //  IEEE (IBM/XT,
        //    SUN, etc.)  (S.P.)    2         -126         128        24
        //  IEEE (IBM/XT,
        //    SUN, etc.)  (D.P.)    2        -1022        1024        53
        //  IBM 3033      (D.P.)   16          -65          63        14
        //  VAX           (S.P.)    2         -128         127        24
        //  VAX D-Format  (D.P.)    2         -128         127        56
        //  VAX G-Format  (D.P.)    2        -1024        1023        53
        //
        //
        //                        NSIG       ENTEN       ENSIG      RTNSIG
        //
        // CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
        // Cyber 180/855
        //   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
        // IEEE (IBM/XT,
        //   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
        // IEEE (IBM/XT,
        //   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
        // IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
        // VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
        // VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
        // VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
        //
        //
        //                         ENMTEN      XLARGE   EXPARG
        //
        // CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677
        // Cyber 180/855
        //   under NOS   (S.P.)   1.25E-293    1.0E+4     741
        // IEEE (IBM/XT,
        //   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88
        // IEEE (IBM/XT,
        //   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
        // IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
        // VAX           (S.P.)   1.17E-38     1.0E+4      88
        // VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
        // VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
        //
    {
        const double constant = 1.585;
        double empal;
        const double enmten = 8.9E-308;
        const double ensig = 1.0E+16;
        const double enten = 1.0E+308;
        const double exparg = 709.0;
        const double half = 0.5;
        int n;
        int ncalc;
        const int nsig = 16;
        const double one = 1.0;
        const double rtnsig = 1.0E-04;
        double tempa;
        double tempb;
        double tempc;
        double tover;
        const double two = 2.0;
        const double xlarge = 1.0E+04;
        const double zero = 0.0;
        switch (nb)
        {
            //
            //  Check for X, NB, OR IZE out of range.
            //
            case <= 0:
                ncalc = Math.Min(nb, 0) - 1;
                return ncalc;
        }

        switch (x)
        {
            case < 0.0:
                ncalc = Math.Min(nb, 0) - 1;
                return ncalc;
        }

        switch (alpha)
        {
            case < 0.0:
                ncalc = Math.Min(nb, 0) - 1;
                return ncalc;
            case >= 1.0:
                ncalc = Math.Min(nb, 0) - 1;
                return ncalc;
        }

        switch (ize)
        {
            case 1 when exparg < x:
                ncalc = Math.Min(nb, 0) - 1;
                return ncalc;
            case 2 when xlarge < x:
                ncalc = Math.Min(nb, 0) - 1;
                return ncalc;
        }

        //
        //  Use 2-term ascending series for small X.
        //
        ncalc = nb;
        int magx = (int) x;
        //
        //  Initialize the forward sweep, the P-sequence of Olver.
        //
        if (rtnsig <= x)
        {
            int nbmx = nb - magx;
            n = magx + 1;
            double en = n + n + (alpha + alpha);
            double plast = one;
            double p = en / x;
            //
            //  Calculate general significance test.
            //
            double test = ensig + ensig;

            if (5 * nsig < 2 * magx)
            {
                test = Math.Sqrt(test * p);
            }
            else
            {
                test /= Math.Pow(constant, magx);
            }

            //
            //  Calculate P-sequence until N = NB-1.  Check for possible overflow.
            //
            bool flag = false;

            int l;
            int nend;
            double pold;
            switch (nbmx)
            {
                case >= 3:
                {
                    tover = enten / ensig;
                    int nstart = magx + 2;
                    nend = nb - 1;

                    int k;
                    for (k = nstart; k <= nend; k++)
                    {
                        n = k;
                        en += two;
                        pold = plast;
                        plast = p;
                        p = en * plast / x + pold;
                        //
                        //  To avoid overflow, divide P-sequence by TOVER.  Calculate
                        //  P-sequence until 1 < ABS(P).
                        //
                        if (!(tover < p))
                        {
                            continue;
                        }

                        tover = enten;
                        p /= tover;
                        plast /= tover;
                        double psave = p;
                        double psavel = plast;
                        nstart = n + 1;

                        for (;;)
                        {
                            n += 1;
                            en += two;
                            pold = plast;
                            plast = p;
                            p = en * plast / x + pold;

                            if (1.0 < p)
                            {
                                break;
                            }
                        }

                        tempb = en / x;
                        //
                        //  Calculate backward test, and find NCALC, the highest N
                        //  such that the test is passed.
                        //
                        test = pold * plast / ensig;
                        test *= half - half / (tempb * tempb);
                        p = plast * tover;
                        n -= 1;
                        en -= two;
                        nend = Math.Min(nb, n);

                        ncalc = nend + 1;

                        for (l = nstart; l <= nend; l++)
                        {
                            pold = psavel;
                            psavel = psave;
                            psave = en * psavel / x + pold;

                            if (!(test < psave * psavel))
                            {
                                continue;
                            }

                            ncalc = l;
                            break;
                        }

                        ncalc -= 1;
                        flag = true;
                        break;
                    }

                    switch (flag)
                    {
                        case false:
                            n = nend;
                            en = n + n + (alpha + alpha);
                            //
                            //  Calculate special significance test for 2 < NBMX.
                            //
                            test = Math.Max(test, Math.Sqrt(plast * ensig) * Math.Sqrt(p + p));
                            break;
                    }

                    break;
                }
            }

            switch (flag)
            {
                //
                //  Calculate P-sequence until significance test passed.
                //
                case false:
                {
                    for (;;)
                    {
                        n += 1;
                        en += two;
                        pold = plast;
                        plast = p;
                        p = en * plast / x + pold;

                        if (test <= p)
                        {
                            break;
                        }
                    }

                    break;
                }
            }

            //
            //  Initialize the backward recursion and the normalization sum.
            //
            n += 1;
            en += two;
            tempb = zero;
            tempa = one / p;
            double em = n - one;
            empal = em + alpha;
            double emp2al = em - one + (alpha + alpha);
            double total = tempa * empal * emp2al / em;
            nend = n - nb;
            switch (nend)
            {
                //
                //  N < NB, so store B(N) and set higher orders to zero.
                //
                case < 0:
                {
                    b[n - 1] = tempa;
                    nend = -nend;

                    for (l = 1; l <= nend; l++)
                    {
                        b[n + l - 1] = zero;
                    }

                    nend = n - 2;
                    switch (nend)
                    {
                        //
                        //  Calculate via difference equation and store B(N), until N = 2.
                        //
                        case > 0:
                        {
                            for (l = 1; l <= nend; l++)
                            {
                                n -= 1;
                                en -= two;
                                b[n - 1] = en * b[n] / x + b[n + 1];
                                em -= one;
                                emp2al -= one;
                                emp2al = n switch
                                {
                                    2 => one,
                                    _ => emp2al
                                };

                                empal -= one;
                                total = (total + b[n - 1] * empal) * emp2al / em;
                            }

                            break;
                        }
                    }

                    //
                    //  Calculate B(1).
                    //
                    b[0] = two * empal * b[1] / x + b[2];

                    total = total + total + b[0];
                    break;
                }
                //
                default:
                {
                    switch (nend)
                    {
                        case > 0:
                        {
                            for (l = 1; l <= nend; l++)
                            {
                                n -= 1;
                                en -= two;
                                tempc = tempb;
                                tempb = tempa;
                                tempa = en * tempb / x + tempc;
                                em -= one;
                                emp2al -= one;

                                if (n == 1)
                                {
                                    break;
                                }

                                emp2al = n switch
                                {
                                    2 => one,
                                    _ => emp2al
                                };

                                empal -= one;
                                total = (total + tempa * empal) * emp2al / em;
                            }

                            break;
                        }
                    }

                    //
                    //  Store B(NB).
                    //
                    b[n - 1] = tempa;

                    switch (nb)
                    {
                        case <= 1:
                            total = total + total + tempa;
                            break;
                        //
                        default:
                        {
                            n -= 1;
                            en -= two;
                            b[n - 1] = en * tempa / x + tempb;

                            switch (n)
                            {
                                case > 1:
                                {
                                    em -= one;
                                    emp2al -= one;

                                    emp2al = n switch
                                    {
                                        2 => one,
                                        _ => emp2al
                                    };

                                    empal -= one;
                                    total = (total + b[n - 1] * empal) * emp2al / em;

                                    nend = n - 2;
                                    switch (nend)
                                    {
                                        //
                                        //  Calculate via difference equation and store B(N), until N = 2.
                                        //
                                        case > 0:
                                        {
                                            for (l = 1; l <= nend; l++)
                                            {
                                                n -= 1;
                                                en -= two;
                                                b[n - 1] = en * b[n] / x + b[n + 1];
                                                em -= one;
                                                emp2al -= one;
                                                emp2al = n switch
                                                {
                                                    2 => one,
                                                    _ => emp2al
                                                };

                                                empal -= one;
                                                total = (total + b[n - 1] * empal) * emp2al / em;
                                            }

                                            break;
                                        }
                                    }

                                    //
                                    //  Calculate B(1).
                                    //
                                    b[0] = two * empal * b[1] / x + b[2];
                                    break;
                                }
                            }

                            total = total + total + b[0];
                            break;
                        }
                    }

                    break;
                }
            }

            //
            //  Normalize.  Divide all B(N) by TOTAL.
            //
            if (Math.Abs(alpha - zero) > double.Epsilon)
            {
                total = total * Helpers.Gamma(one + alpha) * Math.Pow(x * half, -alpha);
            }

            switch (ize)
            {
                case 1:
                    total *= Math.Exp(-x);
                    break;
            }

            tempa = enmten;

            switch (total)
            {
                case > 1.0:
                    tempa *= total;
                    break;
            }

            for (n = 1; n <= nb; n++)
            {
                if (b[n - 1] < tempa)
                {
                    b[n - 1] = zero;
                }

                b[n - 1] /= total;
            }

            return ncalc;
        }
        //
        //  Two-term ascending series for small X.
        //

        tempa = one;
        empal = one + alpha;
        double halfx = zero;

        if (enmten < x)
        {
            halfx = half * x;
        }

        if (Math.Abs(alpha - zero) > double.Epsilon)
        {
            tempa = Math.Pow(halfx, alpha) / Helpers.Gamma(empal);
        }

        switch (ize)
        {
            case 2:
                tempa *= Math.Exp(-x);
                break;
        }

        tempb = zero;

        if (one < x + one)
        {
            tempb = halfx * halfx;
        }

        b[0] = tempa + tempa * tempb / empal;

        if (Math.Abs(x - zero) > double.Epsilon && Math.Abs(b[0] - zero) <= double.Epsilon)
        {
            ncalc = 0;
        }

        switch (nb)
        {
            case > 1 when Math.Abs(x - zero) <= double.Epsilon:
            {
                int i;
                for (i = 1; i < nb; i++)
                {
                    b[i] = zero;
                }

                break;
            }
            //
            //  Calculate higher-order functions.
            //
            case > 1:
            {
                tempc = halfx;
                tover = (enmten + enmten) / x;

                if (Math.Abs(tempb - zero) > double.Epsilon)
                {
                    tover = enmten / tempb;
                }

                for (n = 2; n <= nb; n++)
                {
                    tempa /= empal;
                    empal += one;
                    tempa *= tempc;

                    if (tempa <= tover * empal)
                    {
                        tempa = zero;
                    }

                    b[n - 1] = tempa + tempa * tempb / empal;

                    if (Math.Abs(b[n - 1] - zero) <= double.Epsilon && n < ncalc)
                    {
                        ncalc = n - 1;
                    }
                }

                break;
            }
        }

        return ncalc;
    }
}