using System;
using System.Globalization;
using Burkardt.Types;
using Burkardt.Weight;

namespace Burkardt.Quadrature;

public static class CHKQFS
{
    public static void chkqfs(double[] t, double[] wts, int[] mlt, int nt, int nwts, int[] ndx,
            int key, ref double[] w, int mop, int mex, int kind, double alpha, double beta,
            int lo )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHKQFS checks the polynomial accuracy of a quadrature formula.
        //
        //  Discussion:
        //
        //    This routine will optionally print weights, and results of a moments test.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 February 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, double T[NT], the knots.
        //
        //    Input, double WTS[NWTS], the weights.
        //
        //    Input, int MLT[NT], the multiplicity of the knots.
        //
        //    Input, int NT, the number of knots.
        //
        //    Input, int NWTS, the number of weights.
        //
        //    Input, int NDX[NT], used to index the array WTS.  
        //    If KEY = 1, then NDX need not be preset.  For more details see the 
        //    comments in CAWIQ.
        //
        //    Input, int KEY, indicates the structure of the WTS
        //    array.  It will normally be set to 1.  This will cause the weights to be 
        //    packed sequentially in array WTS.  For more details see the comments 
        //    in CAWIQ.
        //
        //    Input/output, double W[MEX], the moments array.
        //    This is input only if KIND = 0.
        //
        //    Input, int MOP, the expected order of precision of the
        //    quadrature formula.
        //
        //    Input, int MEX, the number of moments to be tested.
        //    MEX must be at least 1.  Set MEX = 1 and LO < 0 for no moment check.
        //
        //    Input, int KIND, the rule.
        //    0, unknown weight function (the user must set the first MEX moments in
        //       array W in this case.)
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        //    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
        //    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
        //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        //    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
        //    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        //
        //    Input, double ALPHA, the value of Alpha, if needed.
        //
        //    Input, double BETA, the value of Beta, if needed.
        //
        //    Input, int LO, selects the action to carry out.
        //     > 0, print weights and moment tests.
        //     = 0, print nothing.  Dompute moment test.
        //     < 0, print weights only.  Do not compute moment tests.
        //
        //  Local Parameters:
        //
        //    Local, double E[MEX], ER[MEX], the absolute and relative 
        //    errors of the quadrature formula applied to (X-DEL)^n.
        //
        //    Local, double QM[MEX], the values of the quadrature formula
        //    applied to (X-DEL)^N.
        //
    {
        double[] e;
        double ek;
        double emn;
        double emx;
        double erest;
        double ern;
        double erx;
        double[] er;
        int i;
        int j;
        int jl;
        int k;
        int kindp;
        int kjl;
        int l;
        int m;
        int mx;
        double px;
        double tmp;
        double tmpx;
        double prec;
        double[] qm;
        //
        //  KIND may be set to -1 to allow printing of moments only.
        //
        //  This feature is only used internally, by CHKQF.
        //
        kindp = Math.Max(0, kind);

        if (lo != 0 && kind != -1)
        {
            if (kindp != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  Interpolatory quadrature formula");
                Console.WriteLine("");
                Console.WriteLine("  Type  Interval       Weight function               Name");
                Console.WriteLine("");
                switch (kindp)
                {
                    case 1:
                        Console.WriteLine("    1    (-1,+1)            1.0                    Legendre");
                        break;
                    case 2:
                        Console.WriteLine("    2    (-1,+1)    ((b-x)*(x-a))^(-0.5)          Chebyshev Type 1");
                        break;
                    case 3:
                        Console.WriteLine("    3    (-1,+1)    ((b-x)*(x-a))^alpha           Gegenbauer");
                        break;
                    case 4:
                        Console.WriteLine("    4    (-1,+1)  (b-x)^alpha*(x-a)^beta          Jacobi");
                        break;
                    case 5:
                        Console.WriteLine("    5   (a,+oo)   (x-a)^alpha*exp(-b*(x-a))      Gen Laguerre");
                        break;
                    case 6:
                        Console.WriteLine("    6  (-oo,+oo) |x-a|^alpha*exp(-b*(x-a)^2)  Gen Hermite");
                        break;
                    case 7:
                        Console.WriteLine("    7    (-1,+1)    |x-(a+b)/2.0|^alpha        Exponential");
                        break;
                    case 8:
                        Console.WriteLine("    8   (0,+oo)    (x-a)^alpha*(x+b)^beta         Rational");
                        break;
                    case 9:
                        Console.WriteLine("    9    (-1,+1)    ((b-x)*(x-a))^(+0.5)          Chebyshev Type 2");
                        break;
                }

                switch (kindp)
                {
                    case >= 3 and <= 8:
                        Console.WriteLine("                  alpha      " + alpha + "");
                        break;
                }

                switch (kindp)
                {
                    case 4:
                    case 8:
                        Console.WriteLine("                  beta       " + beta + "");
                        break;
                }

            }

            if (kind != -1)
            {
                prec = typeMethods.r8_epsilon();
                Console.WriteLine("");
                Console.WriteLine("  Machine precision = " + prec + "");
            }

            Console.WriteLine("");
            Console.WriteLine("           Knots               Mult                Weights");
            Console.WriteLine("");

            for (i = 1; i <= nt; i++)
            {
                k = Math.Abs(ndx[i - 1]);
                if (k != 0)
                {
                    Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                      + t[i - 1].ToString("0.#################").PadLeft(26)
                                      + mlt[i - 1].ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                      + wts[k - 1].ToString("0.#################").PadLeft(26) + "");
                    for (j = k + 1; j <= k + mlt[i - 1] - 1; j++)
                    {
                        Console.WriteLine("                                  "
                                          + wts[j - 1].ToString("0.#################").PadLeft(26) + "");
                    }
                }
            }
        }

        switch (lo)
        {
            case < 0:
                return;
        }

        //
        //  Compute the moments in W.
        //
        if (kindp != 0)
        {
            w = WM.wm(mex, kindp, alpha, beta);
        }

        e = new double[mex];
        er = new double[mex];
        qm = new double[mex];

        for (j = 0; j < mex; j++)
        {
            qm[j] = 0.0;
        }

        erest = 0.0;

        for (k = 1; k <= nt; k++)
        {
            tmp = 1.0;
            l = Math.Abs(ndx[k - 1]);
            switch (l)
            {
                case 0:
                    continue;
            }

            erest += Math.Abs(wts[l - 1]);
            for (j = 1; j <= mex; j++)
            {
                qm[j - 1] += tmp * wts[l - 1];
                tmpx = tmp;
                px = 1.0;
                for (jl = 2; jl <= Math.Min(mlt[k - 1], mex - j + 1); jl++)
                {
                    kjl = j + jl - 1;
                    tmpx *= (kjl - 1);
                    qm[kjl - 1] += tmpx * wts[l + jl - 2] / px;
                    switch (key)
                    {
                        case <= 0:
                            px *= jl;
                            break;
                    }
                }

                tmp *= t[k - 1];
            }

        }

        for (j = 0; j < mex; j++)
        {
            e[j] = w[j] - qm[j];
            er[j] = e[j] / (Math.Abs(w[j]) + 1.0);
        }

        //
        //  For some strange weight functions W(1) may vanish.
        //
        erest /= (Math.Abs(w[0]) + 1.0);

        switch (lo)
        {
            case > 0:
            {
                m = mop + 1;
                mx = Math.Min(mop, mex);

                emx = Math.Abs(e[0]);
                emn = emx;
                erx = Math.Abs(er[0]);
                ern = erx;
                for (k = 1; k < mx; k++)
                {
                    emx = Math.Max(Math.Abs(e[k]), emx);
                    emn = Math.Min(Math.Abs(e[k]), emn);
                    erx = Math.Max(Math.Abs(er[k]), erx);
                    ern = Math.Min(Math.Abs(er[k]), ern);
                }

                Console.WriteLine("");
                Console.WriteLine("  Comparison of moments");
                Console.WriteLine("");
                Console.WriteLine("  Order of precision " + mop + "");
                Console.WriteLine("  Errors :    Absolute    Relative");
                Console.WriteLine("  ---------+-------------------------");
                Console.WriteLine("  Minimum :" + emn.ToString("0.###").PadLeft(12)
                                                + "  " + ern.ToString("0.###").PadLeft(12) + "");
                Console.WriteLine("  Maximum :" + emx.ToString("0.###").PadLeft(12)
                                                + "  " + erx.ToString("0.###").PadLeft(12) + "");
                Console.WriteLine("");
                Console.WriteLine("  Weights ratio       "
                                  + erest.ToString("0.###").PadLeft(12) + "");

                if (m <= mex)
                {
                    ek = e[m - 1];
                    for (j = 1; j <= mop; j++)
                    {
                        ek /= j;
                    }

                    Console.WriteLine("  Error in " + mop + "th power "
                                      + e[m - 1].ToString("0.###").PadLeft(12) + "");
                    Console.WriteLine("  Error constant      "
                                      + ek.ToString("0.###").PadLeft(12) + "");
                }

                Console.WriteLine("");
                Console.WriteLine("  Moments:");
                Console.WriteLine("");
                Console.WriteLine("            True             from QF            Error      Relative");
                Console.WriteLine("");
                for (j = 1; j <= mx; j++)
                {
                    Console.WriteLine(j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                      + w[j - 1].ToString("0.##########").PadLeft(19)
                                      + qm[j - 1].ToString("0.##########").PadLeft(19)
                                      + e[j - 1].ToString("0.###").PadLeft(12)
                                      + er[j - 1].ToString("0.###").PadLeft(12) + "");
                }

                Console.WriteLine("");
                for (j = m; j <= mex; j++)
                {
                    Console.WriteLine(j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                      + w[j - 1].ToString("0.##########").PadLeft(19)
                                      + qm[j - 1].ToString("0.##########").PadLeft(19)
                                      + e[j - 1].ToString("0.###").PadLeft(12)
                                      + er[j - 1].ToString("0.###").PadLeft(12) + "");
                }

                break;
            }
        }
    }
}