using System;
using System.Diagnostics;
using Burkardt.Types;

namespace Burkardt.Chebyshev;

public static class ChebyshevSeries
{
    public static double[] dfrnt(double[] xx, int npl)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DFRNT determines the derivative of a Chebyshev series.
        //
        //  Discussion:
        //
        //    This routine computes the Chebyshev series of the derivative of a 
        //    function whose Chebyshev series is given.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double XX[NPL], the given Chebyshev series.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Output, double DFRNT[NPL], the Chebyshev series for the
        //    derivative.
        //
    {
        double dl;
        double dn;
        int k;
        int l;
        double[] x2;
        double xxl;
        double xxn;

        x2 = new double[npl];
        dn = npl - 1;
        xxn = xx[npl - 2];
        x2[npl - 2] = 2.0 * xx[npl - 1] * dn;
        x2[npl - 1] = 0.0;

        for (k = 3; k <= npl; k++)
        {
            l = npl - k + 1;
            dl = l;
            xxl = xx[l - 1];
            x2[l - 1] = x2[l + 1] + 2.0 * xxn * dl;
            xxn = xxl;
        }

        return x2;
    }

    public static double echeb(double x, double[] coef, int npl)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ECHEB evaluates a Chebyshev series at a point.
        //
        //  Discussion:
        //
        //    This routine evaluates a Chebyshev series at a point in [-1,+1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NPL], the Chebyshev series.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Output, double ECHEB, the value of the Chebyshev series at X.
        //
    {
        double br;
        double brp2 = 0;
        double brpp;
        double fx;
        int j;
        int k;

        br = 0.0;
        brpp = 0.0;

        for (k = 1; k <= npl; k++)
        {
            j = npl - k + 1;
            brp2 = brpp;
            brpp = br;
            br = 2.0 * x * brpp - brp2 + coef[j - 1];
        }

        fx = 0.5 * (br - brp2);
        return fx;
    }

    public static double edcheb(double x, double[] coef, int npl)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EDCHEB evaluates the derivative of a Chebyshev series at a point.
        //
        //  Discussion:
        //
        //    This routine evaluates the derivative of a Chebyshev series 
        //    at a point in [-1,+1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NPL], the Chebyshev series.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Output, double EDCHEB, the value of the derivative of the
        //    Chebyshev series at X.
        //
    {
        double bf = 0;
        double bj = 0;
        double bjp2;
        double bjpl;
        double dj;
        double fx;
        int j;
        int k;
        int n;
        double xj;
        double xjp2;
        double xjpl;

        xjp2 = 0.0;
        xjpl = 0.0;
        bjp2 = 0.0;
        bjpl = 0.0;
        n = npl - 1;

        for (k = 1; k <= n; k++)
        {
            j = npl - k;
            dj = j;
            xj = 2.0 * coef[j] * dj + xjp2;
            bj = 2.0 * x * bjpl - bjp2 + xj;
            bf = bjp2;
            bjp2 = bjpl;
            bjpl = bj;
            xjp2 = xjpl;
            xjpl = xj;
        }

        fx = 0.5 * (bj - bf);
        return fx;
    }

    public static double[] invert(double[] x, double[] xx, int npl, int net)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERT computes the inverse Chebyshev series.
        //
        //  Discussion:
        //
        //    This routine computes the inverse of a Chebyshev series, starting with
        //    an initial approximation XX. 
        //
        //    The routine uses the Euler method and computes all powers EPS^K 
        //    up to K=2^(NET+1), where EPS = 1 - X * ( XX inverse ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double X[NPL], the Chebyshev series.
        //
        //    Input, double XX[NPL], an initial approximation for the
        //    inverse Chebyshev series.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Input, int NET, the number of iterations to take.
        //
        //    Output, double INVERT[NPL], the estimated Chebyshev
        //    series of the inverse function.
        //
    {
        int k;
        double s;
        double[] w2;
        double[] ww;
        double[] xnvse;

        ww = mltply_new(x, xx, npl);

        s = -1.0;
        typeMethods.r8vec_scale(s, npl, ref ww);
        ww[0] += 2.0;

        w2 = mltply_new(ww, ww, npl);
        ww[0] = 2.0 * ww[0];

        xnvse = new double[npl];

        for (k = 1; k <= net; k++)
        {
            mltply(ww, w2, npl, xnvse);

            typeMethods.r8vec_add(npl, xnvse, ref ww);

            mltply(w2, w2, npl, xnvse);

            typeMethods.r8vec_copy(npl, xnvse, ref w2);
        }

        mltply(ww, xx, npl, xnvse);

        return xnvse;
    }

    public static void mltply(double[] xx, double[] x2, int npl, double[] x3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MLTPLY_NEW multiplies two Chebyshev series.
        //
        //  Discussion:
        //
        //    This routine multiplies two given Chebyshev series, XX and X2,
        //    to produce an output Chebyshev series, X3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double XX[NPL], the first Chebyshev series.
        //
        //    Input, double X2[NPL], the second Chebyshev series.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Output, double X3[NPL], the Chebyshev series of the
        //    product.
        //
    {
        double ex;
        int k;
        int l;
        int m;
        int mm;

        for (k = 1; k <= npl; k++)
        {
            ex = 0.0;
            mm = npl - k + 1;
            for (m = 1; m <= mm; m++)
            {
                l = m + k - 1;
                ex = ex + xx[m - 1] * x2[l - 1] + xx[l - 1] * x2[m - 1];
            }

            x3[k - 1] = 0.5 * ex;
        }

        x3[0] -= 0.5 * xx[0] * x2[0];

        for (k = 3; k <= npl; k++)
        {
            ex = 0.0;
            mm = k - 1;
            for (m = 2; m <= mm; m++)
            {
                l = k - m + 1;
                ex += xx[m - 1] * x2[l - 1];
            }

            x3[k - 1] = 0.5 * ex + x3[k - 1];
        }
    }

    public static double[] mltply_new ( double[] xx, double[] x2, int npl )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MLTPLY_NEW multiplies two Chebyshev series.
        //
        //  Discussion:
        //
        //    This routine multiplies two given Chebyshev series, XX and X2,
        //    to produce an output Chebyshev series, X3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double XX[NPL], the first Chebyshev series.
        //
        //    Input, double X2[NPL], the second Chebyshev series.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Output, double MLTPLY_NEW[NPL], the Chebyshev series of the
        //    product.
        //
    {
        double[] x3;

        x3 = typeMethods.r8vec_zero_new ( npl );

        mltply ( xx, x2, npl, x3 );

        return x3;
    }

    public static double[] ntgrt(double[] xx, int npl)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NTGRT determines the integral of a Chebyshev series.
        //
        //  Discussion:
        //
        //    This routine computes the Chebyshev series for the integral of a 
        //    function whose Chebyshev series is given.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double XX[NPL], the Chebyshev series.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Output, double NTGRT[NPL], the Chebyshev series for the
        //    integral of the function.
        //
    {
        double dk;
        int k;
        int n;
        double term;
        double[] x2;
        double xpr;

        x2 = new double[npl];

        xpr = xx[0];
        x2[0] = 0.0;
        n = npl - 1;

        for (k = 2; k <= n; k++)
        {
            dk = k - 1;
            term = (xpr - xx[k]) / (2.0 * dk);
            xpr = xx[k - 1];
            x2[k - 1] = term;
        }

        dk = n;
        x2[npl - 1] = xpr / (2.0 * dk);

        return x2;
    }

    public static double[] binom ( double[] x, double[] xx, int npl, int m, int nt )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOM: binomial expansion series for the (-1/M) power of a Chebyshev series.
        //
        //  Discussion:
        //
        //    This routine uses a certain number of terms of the binomial expansion 
        //    series to estimate the (-1/M) power of a given Chebyshev series. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Broucke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 4, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double X[NPL], the given Chebyshev series.
        //
        //    Input, double XX[NPL], an initial estimate for
        //    the Chebyshev series for the input function raised to the (-1/M) power.
        //
        //    Input, int NPL, the number of terms in the 
        //    Chebyshev series.
        //
        //    Input, int M, defines the exponent, (-1/M).
        //    0 < M.
        //
        //    Input, int NT, the number of terms of the binomial
        //    series to be used.
        //
        //    Output, double BINOM[NPL], the estimated Chebyshev series
        //    for the input function raised to the (-1/M) power.
        //
    {
        double alfa;
        double coef;
        double dkm2;
        double dkmm;
        double dm;
        int k;
        double[] w2;
        double[] w3;
        double[] ww;
        double[] xa;

        dm = m;
        alfa = -1.0 / dm;

        ww = typeMethods.r8vec_copy_new(npl, x);

        w2 = new double[npl];

        for (k = 1; k <= m; k++)
        {
            mltply(ww, xx, npl, w2);
            typeMethods.r8vec_copy(npl, w2, ref ww);
        }

        ww[0] -= 2.0;

        xa = typeMethods.r8vec_zero_new(npl);
        xa[0] = 2.0;

        w3 = typeMethods.r8vec_copy_new(npl, xa);

        for (k = 2; k <= nt; k++)
        {
            dkmm = k - 1;
            dkm2 = k - 2;
            coef = (alfa - dkm2) / dkmm;

            mltply(w3, ww, npl, w2);

            typeMethods.r8vec_copy(npl, w2, ref w3);
            typeMethods.r8vec_scale(coef, npl, ref w3);
            typeMethods.r8vec_add(npl, w3, ref xa);
        }

        mltply(xa, xx, npl, w2);

        typeMethods.r8vec_copy(npl, w2, ref xa);

        return xa;
    }

    public static double echebser0(double x, double[] coef, int nc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ECHEBSER0 evaluates a Chebyshev series.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of 
        //    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
        //    gives an example for treating the first derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ECHEBSER0, the value of the Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];

        x2 = 2.0 * x;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;
        }

        value = 0.5 * (b0 - b2);

        return value;
    }

    public static double echebser1(double x, double[] coef, int nc, ref double y1 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ECHEBSER1 evaluates a Chebyshev series and first derivative.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of 
        //    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
        //    gives an example for treating the first derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ECHEBSER1, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];
        c0 = coef[nc - 1];

        x2 = 2.0 * x;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;

            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + x2 * c1;
                    break;
            }
        }

        y1 = c0 - c2;

        value = 0.5 * (b0 - b2);

        return value;
    }

    public static double echebser2(double x, double[] coef, int nc, ref double y1, ref double y2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ECHEBSER2 evaluates a Chebyshev series and two derivatives.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of 
        //    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
        //    gives an example for treating the first derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ECHEBSER2, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
        //    Output, double &Y2, the value of the 2nd derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        double d0;
        double d1 = 0.0;
        double d2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];
        c0 = coef[nc - 1];
        d0 = coef[nc - 1];

        x2 = 2.0 * x;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;

            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + x2 * c1;
                    break;
            }

            switch (i)
            {
                case > 1:
                    d2 = d1;
                    d1 = d0;
                    d0 = c0 - d2 + x2 * d1;
                    break;
            }
        }

        y2 = (d0 - d2) * 4.0;
        y1 = c0 - c2;

        value = 0.5 * (b0 - b2);

        return value;
    }

    public static double echebser3(double x, double[] coef, int nc, ref double y1, ref double y2,
            ref  double y3 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ECHEBSER3 evaluates a Chebyshev series and three derivatives.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of 
        //    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
        //    gives an example for treating the first derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ECHEBSER3, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
        //    Output, double &Y2, the value of the 2nd derivative of the
        //    Chebyshev series at X.
        //
        //    Output, double &Y3, the value of the 3rd derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        double d0;
        double d1 = 0.0;
        double d2 = 0.0;
        double e0;
        double e1 = 0.0;
        double e2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];
        c0 = coef[nc - 1];
        d0 = coef[nc - 1];
        e0 = coef[nc - 1];

        x2 = 2.0 * x;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;

            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + x2 * c1;
                    break;
            }

            switch (i)
            {
                case > 1:
                    d2 = d1;
                    d1 = d0;
                    d0 = c0 - d2 + x2 * d1;
                    break;
            }

            switch (i)
            {
                case > 2:
                    e2 = e1;
                    e1 = e0;
                    e0 = d0 - e2 + x2 * e1;
                    break;
            }
        }

        y3 = (e0 - e2) * 24.0;
        y2 = (d0 - d2) * 4.0;
        y1 = c0 - c2;

        value = 0.5 * (b0 - b2);

        return value;
    }

    public static double echebser4(double x, double[] coef, int nc, ref double y1, ref double y2,
            ref double y3, ref double y4 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ECHEBSER4 evaluates a Chebyshev series and four derivatives.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of 
        //    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
        //    gives an example for treating the first derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 April 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ECHEBSER3, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, &Y2, &Y3, &Y4, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        double d0;
        double d1 = 0.0;
        double d2 = 0.0;
        double e0;
        double e1 = 0.0;
        double e2 = 0.0;
        double f0;
        double f1 = 0.0;
        double f2 = 0.0;
        int i;
        double y0;

        b0 = coef[nc - 1];
        c0 = coef[nc - 1];
        d0 = coef[nc - 1];
        e0 = coef[nc - 1];
        f0 = coef[nc - 1];

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + 2.0 * x * b1;

            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + 2.0 * x * c1;
                    break;
            }

            switch (i)
            {
                case > 1:
                    d2 = d1;
                    d1 = d0;
                    d0 = c0 - d2 + 2.0 * x * d1;
                    break;
            }

            switch (i)
            {
                case > 2:
                    e2 = e1;
                    e1 = e0;
                    e0 = d0 - e2 + 2.0 * x * e1;
                    break;
            }

            switch (i)
            {
                case > 3:
                    f2 = f1;
                    f1 = f0;
                    f0 = e0 - f2 + 2.0 * x * f1;
                    break;
            }
        }

        y0 = (b0 - b2) / 2.0;
        y1 = c0 - c2;
        y2 = (d0 - d2) * 2.0 * 2.0;
        y3 = (e0 - e2) * 6.0 * 4.0;
        y4 = (f0 - f2) * 24.0 * 8.0;

        return y0;
    }

    public static double evenchebser0(double x, double[] coef, int nc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EVENCHEBSER0 evaluates an even Chebyshev series.
        //
        //  Discussion:
        //
        //    This function implements Clenshaw's modification of his
        //    algorithm for even series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double EVENCHEBSER0, the value of the Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];

        x2 = 4.0 * x * x - 2.0;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;
        }

        value = 0.5 * (b0 - b2);

        return value;
    }

    public static double evenchebser1(double x, double[] coef, int nc, ref double y1 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EVENCHEBSER1 evaluates an even Chebyshev series and first derivative.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of 
        //    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
        //    gives an example for treating the first derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double EVENCHEBSER1, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];
        c0 = coef[nc - 1];

        x2 = 4.0 * x * x - 2.0;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;
            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + x2 * c1;
                    break;
            }
        }

        y1 = (c0 - c2) * 4.0 * x;

        value = 0.5 * (b0 - b2);

        return value;
    }

    public static double evenchebser2(double x, double[] coef, int nc, ref double y1, ref double y2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EVENCHEBSER2 evaluates an even Chebyshev series and first two derivatives.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of
        //    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
        //    gives an example for treating the first derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double EVENCHEBSER2, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
        //    Output, double &Y2, the value of the 2nd derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        double d0;
        double d1 = 0.0;
        double d2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];
        c0 = coef[nc - 1];
        d0 = coef[nc - 1];

        x2 = 4.0 * x * x - 2.0;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;
            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + x2 * c1;
                    break;
            }

            switch (i)
            {
                case > 1:
                    d2 = d1;
                    d1 = d0;
                    d0 = c0 - d2 + x2 * d1;
                    break;
            }
        }

        y2 = (d0 - d2) * 64.0 * x * x + (c0 - c2) * 4.0;
        y1 = (c0 - c2) * 4.0 * x;

        value = 0.5 * (b0 - b2);

        return value;
    }

    public static double oddchebser0(double x, double[] coef, int nc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ODDCHEBSER0 evaluates an odd Chebyshev series.
        //
        //  Discussion:
        //
        //    This function implements Clenshaw's modification of  his algorithm
        //    for odd series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ODDCHEBSER0, the value of the Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        int i;
        double value = 0;
        double x2;

        b0 = coef[nc - 1];

        x2 = 4.0 * x * x - 2.0;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = coef[i] - b2 + x2 * b1;
        }

        value = (b0 - b1) * x;

        return value;
    }

    public static double oddchebser1(double x, double[] coef, int nc, ref double y1 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ODDCHEBSER1 evaluates an odd Chebyshev series and the first derivative.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of
        //    Clenshaw's algorithm. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ODDCHEBSER1, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        double coefi;
        int i;
        double value = 0;
        double x2;

        coefi = 2.0 * coef[nc - 1];
        b0 = coefi;
        c0 = coefi;

        x2 = 4.0 * x * x - 2.0;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            coefi = 2.0 * coef[i] - coefi;
            b0 = coefi - b2 + x2 * b1;
            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + x2 * c1;
                    break;
            }
        }

        y1 = (c0 - c2) * 4.0 * x * x + (b0 - b2) * 0.5;
        value = (b0 - b2) * 0.5 * x;

        return value;
    }

    public static double oddchebser2(double x, double[] coef, int nc, ref double y1, ref double y2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ODDCHEBSER2 evaluates an odd Chebyshev series and first two derivatives.
        //
        //  Discussion:
        //
        //    This function implements a modification and extension of
        //    Clenshaw's algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2014
        //
        //  Author:
        //
        //    Manfred Zimmer
        //
        //  Reference:
        //
        //    Charles Clenshaw,
        //    Mathematical Tables, Volume 5,
        //    Chebyshev series for mathematical functions,
        //    London, 1962.
        //
        //    Gerhard Maess,
        //    Vorlesungen ueber Numerische Mathematik II, Analysis,
        //    Berlin, Akademie_Verlag, 1984-1988,
        //    ISBN: 978-3764318840,
        //    LC: QA297.M325.��
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //    -1 <= X <= +1.
        //
        //    Input, double COEF[NC], the Chebyshev series.
        //
        //    Input, int NC, the number of terms in the Chebyshev series.
        //    0 < NC.
        //
        //    Output, double ODDCHEBSER2, the value of the Chebyshev series at X.
        //
        //    Output, double &Y1, the value of the 1st derivative of the
        //    Chebyshev series at X.
        //
        //    Output, double &Y2, the value of the 2nd derivative of the
        //    Chebyshev series at X.
        //
    {
        Debug.Assert
        (
            0 < nc &&
            Math.Abs(x) <= 1.0
        );

        double b0;
        double b1 = 0.0;
        double b2 = 0.0;
        double c0;
        double c1 = 0.0;
        double c2 = 0.0;
        double d0;
        double d1 = 0.0;
        double d2 = 0.0;
        double coefi;
        int i;
        double value = 0;
        double x2;

        coefi = 2.0 * coef[nc - 1];
        b0 = coefi;
        c0 = coefi;
        d0 = coefi;

        x2 = 4.0 * x * x - 2.0;

        for (i = nc - 2; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            coefi = 2.0 * coef[i] - coefi;
            b0 = coefi - b2 + x2 * b1;
            switch (i)
            {
                case > 0:
                    c2 = c1;
                    c1 = c0;
                    c0 = b0 - c2 + x2 * c1;
                    break;
            }

            switch (i)
            {
                case > 1:
                    d2 = d1;
                    d1 = d0;
                    d0 = c0 - d2 + x2 * d1;
                    break;
            }
        }

        value = (b0 - b2) * 0.5 * x;

        x2 = x * x;
        y1 = (c0 - c2) * 4.0 * x2 + (b0 - b2) * 0.5;
        y2 = ((d0 - d2) * 64.0 * x2 + (c0 - c2) * 12.0) * x;

        return value;
    }
}