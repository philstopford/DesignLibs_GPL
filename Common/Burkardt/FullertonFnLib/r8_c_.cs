using System;

namespace Burkardt.FullertonFnLib;

public static partial class FullertonLib
{
    public class r8CBRTData
    {
        public int niter;
        public r8PakData pakdata = new();
    }
    public static double r8_cbrt(ref r8CBRTData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CBRT computes the cube root of an R8.
        //
        //  Discussion:
        //
        //    The approximation is a generalized Chebyshev series converted
        //    to polynomial form.  The approximation is nearly best in the 
        //    sense of relative error with 4.085 digits accuracy.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the number whose square root is desired.
        //
        //    Output, double R8_CBRT, the cube root of X.
        //
    {
        double[] cbrt2 = {
                0.62996052494743658238360530363911,
                0.79370052598409973737585281963615,
                1.0,
                1.25992104989487316476721060727823,
                1.58740105196819947475170563927231
            }
            ;
        int n = 0;
        double y = 0;

        data.niter = data.niter switch
        {
            0 => (int) (1.443 * Math.Log(-0.106 * Math.Log(0.1 * r8_mach(3))) + 1.0),
            _ => data.niter
        };

        double value = 0.0;

        if (x != 0.0)
        {
            r8_upak(Math.Abs(x), ref y, ref n);
            int ixpnt = n / 3;
            int irem = n - 3 * ixpnt + 3;

            value = 0.439581 + y * (
                0.928549 + y * (
                    -0.512653 + y *
                    0.144586));

            int iter;
            for (iter = 1; iter <= data.niter; iter++)
            {
                double vsq = value * value;
                value += (y - value * vsq) / (3.0 * vsq);
            }

            value = x switch
            {
                < 0.0 => -Math.Abs(value),
                _ => +Math.Abs(value)
            };

            value = r8_pak( ref data.pakdata, cbrt2[irem - 1] * value, ixpnt);
        }

        return value;
    }

    public static double r8_chi(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHI evaluates the hyperbolic cosine integral of an R8 argument.
        //
        //  Discussion:
        //
        //    The hyperbolic cosine integral is defined by
        //
        //      CHI(X) = gamma + log ( x ) 
        //        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
        //
        //    where gamma is Euler's constant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_CHI, the hyperbolic cosine integral
        //    evaluated at X.
        //
    {
        double value = 0;

        r8E1Data data = new();

        value = 0.5 * (r8_ei(x) - r8_e1(ref data, x));

        return value;
    }

    public class r8CHUData
    {
        public double eps;
        public r8CHUScaledData scaledData = new();
        public r8PochData pochData = new();
        public r8Poch1Data poch1Data = new();
        public r8ExprelData expreldata = new();
        public r8GamrData gamrdata = new();

    }
    public static double r8_chu(ref r8CHUData data, ref r8GammaData gdata, double a, double b, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHU evaluates the confluent hypergeometric function of R8 arguments.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters.
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_CHU, the function value.
        //
    {
        int i;
        int m;

        double sum;
        double t;
        double value;
        double xi;
        double xi1;

        data.eps = data.eps switch
        {
            0.0 => r8_mach(3),
            _ => data.eps
        };

        switch (x)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_CHU - Fatal error!");
                Console.WriteLine("  X < 0.");
                return 1;
            case 0.0 when 1.0 <= b:
                Console.WriteLine("");
                Console.WriteLine("R8_CHU - Fatal error!");
                Console.WriteLine("  X = 0 and 1 <= B.");
                return 1;
            case 0.0:
                value = r8_gamma(ref gdata, 1.0 - b) / r8_gamma(ref gdata, 1.0 + a - b);
                return value;
        }

        if (r8_max(Math.Abs(a), 1.0)
            * r8_max(Math.Abs(1.0 + a - b), 1.0) < 0.99 * Math.Abs(x))
        {
            value = r8_power(x, -a) * r8_chu_scaled(ref data.scaledData, a, b, x);
            return value;
        }

        double aintb = b switch
        {
            //
            //  The ascending series will be used, because the descending rational
            //  approximation (which is based on the asymptotic series) is unstable.
            //
            >= 0.0 => r8_aint(b + 0.5),
            _ => r8_aint(b - 0.5)
        };

        double beps = b - aintb;
        int n = (int) aintb;
        double alnx = Math.Log(x);
        double xtoeps = Math.Exp(-beps * alnx);
        switch (n)
        {
            //
            //  Evaluate the finite sum.
            //
            //  Consider the case b < 1.0 first.
            //
            case < 1:
            {
                sum = 1.0;
                t = 1.0;
                m = -n;
                for (i = 1; i <= m; i++)
                {
                    xi1 = i - 1;
                    t = t * (a + xi1) * x / ((b + xi1) * (xi1 + 1.0));
                    sum += t;
                }

                sum = r8_poch(ref data.pochData, ref gdata, 1.0 + a - b, -a) * sum;
                break;
            }
            //
            default:
            {
                sum = 0.0;
                m = n - 2;

                switch (m)
                {
                    case >= 0:
                    {
                        t = 1.0;
                        sum = 1.0;

                        for (i = 1; i <= m; i++)
                        {
                            xi = i;
                            t = t * (a - b + xi) * x / ((1.0 - b + xi) * xi);
                            sum += t;
                        }

                        sum = r8_gamma(ref gdata, b - 1.0) * r8_gamr(ref data.gamrdata, ref gdata, a)
                                                           * r8_power(x, 1 - n) * xtoeps * sum;
                        break;
                    }
                }

                break;
            }
        }

        int istrt = n switch
        {
            //
            //  Next evaluate the infinite sum.
            //
            < 1 => 1 - n,
            _ => 0
        };

        xi = istrt;

        double factor = r8_mop(n) * r8_gamr(ref data.gamrdata, ref gdata, 1.0 + a - b) * r8_power(x, xi);

        if (beps != 0.0)
        {
            factor = factor * beps * Math.PI / Math.Sin(beps * Math.PI);
        }

        double pochai = r8_poch(ref data.pochData, ref gdata, a, xi);
        double gamri1 = r8_gamr(ref data.gamrdata, ref gdata, xi + 1.0);
        double gamrni = r8_gamr(ref data.gamrdata, ref gdata, aintb + xi);
        double b0 = factor * r8_poch(ref data.pochData, ref gdata, a, xi - beps)
                           * gamrni * r8_gamr(ref data.gamrdata, ref gdata, xi + 1.0 - beps);
        switch (Math.Abs(xtoeps - 1.0))
        {
            //
            //  x^(-beps) is close to 1.0, so we must be careful in evaluating the
            //  differences.
            //
            case <= 0.5:
            {
                double pch1ai = r8_poch1(ref data.poch1Data, ref gdata, a + xi, -beps);
                double pch1i = r8_poch1(ref data.poch1Data, ref gdata, xi + 1.0 - beps, beps);
                double c0 = factor * pochai * gamrni * gamri1 * (
                    -r8_poch1(ref data.poch1Data, ref gdata, b + xi, -beps) + pch1ai
                    - pch1i + beps * pch1ai * pch1i);
                //
                //  xeps1 = (1.0 - x^(-beps))/beps = (x^(-beps) - 1.0)/(-beps)
                //
                double xeps1 = alnx * r8_exprel( ref data.expreldata,-beps * alnx);

                value = sum + c0 + xeps1 * b0;
                double xn = n;

                for (i = 1; i <= 1000; i++)
                {
                    xi = istrt + i;
                    xi1 = istrt + i - 1;
                    b0 = (a + xi1 - beps) * b0 * x
                         / ((xn + xi1) * (xi - beps));
                    c0 = (a + xi1) * c0 * x / ((b + xi1) * xi)
                         - ((a - 1.0) * (xn + 2.0 * xi - 1.0)
                            + xi * (xi - beps)) * b0
                         / (xi * (b + xi1) * (a + xi1 - beps));
                    t = c0 + xeps1 * b0;
                    value += t;
                    if (Math.Abs(t) < data.eps * Math.Abs(value))
                    {
                        return value;
                    }
                }

                Console.WriteLine("");
                Console.WriteLine("R8_CHU - Fatal error!");
                Console.WriteLine("  No convergence in 1000 terms.");
                return 1;
            }
        }

        //
        //  x^(-beps) is very different from 1.0, so the straightforward
        //  formulation is stable.
        //
        double a0 = factor * pochai * r8_gamr(ref data.gamrdata, ref gdata, b + xi) * gamri1 / beps;
        b0 = xtoeps * b0 / beps;

        value = sum + a0 - b0;

        for (i = 1; i <= 1000; i++)
        {
            xi = istrt + i;
            xi1 = istrt + i - 1;
            a0 = (a + xi1) * a0 * x / ((b + xi1) * xi);
            b0 = (a + xi1 - beps) * b0 * x
                 / ((aintb + xi1) * (xi - beps));
            t = a0 - b0;
            value += t;
            if (Math.Abs(t) < data.eps * Math.Abs(value))
            {
                return value;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("R8_CHU - Fatal error!");
        Console.WriteLine("  No convergence in 1000 terms.");
        return 1;
    }

    public class r8CHUScaledData
    {
        public double eps;

    }
    public static double r8_chu_scaled(ref r8CHUScaledData data, double a, double b, double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHU_SCALED: scaled confluent hypergeometric function of R8 arguments.
        //
        //  Discussion:
        //
        //    Evaluate, for large z, z^a * u(a,b,z)  where U is the logarithmic
        //    confluent hypergeometric function.  A rational approximation due to
        //    Y L Luke is used.  When U is not in the asymptotic region, that is, when A
        //    or B is large compared with Z, considerable significance loss occurs.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters.
        //
        //    Input, double Z, the argument.
        //
        //    Output, double R8_CHU_SCALED, the function value.
        //
    {
        double[] aa = new double[4];
        double[] bb = new double[4];
        int i;
        //static double sqeps = 0.0;

        data.eps = data.eps switch
        {
            0.0 => 4.0 * r8_mach(4),
            _ => data.eps
        };

        double bp = 1.0 + a - b;
        double ab = a * bp;
        double ct2 = 2.0 * (z - ab);
        double sab = a + bp;

        bb[0] = 1.0;
        aa[0] = 1.0;

        double ct3 = sab + 1.0 + ab;
        bb[1] = 1.0 + 2.0 * z / ct3;
        aa[1] = 1.0 + ct2 / ct3;

        double anbn = ct3 + sab + 3.0;
        double ct1 = 1.0 + 2.0 * z / anbn;
        bb[2] = 1.0 + 6.0 * ct1 * z / ct3;
        aa[2] = 1.0 + 6.0 * ab / anbn + 3.0 * ct1 * ct2 / ct3;

        for (i = 4; i <= 300; i++)
        {
            double x2i1 = 2 * i - 3;
            ct1 = x2i1 / (x2i1 - 2.0);
            anbn = anbn + x2i1 + sab;
            ct2 = (x2i1 - 1.0) / anbn;
            double c2 = x2i1 * ct2 - 1.0;
            double d1z = x2i1 * 2.0 * z / anbn;

            ct3 = sab * ct2;
            double g1 = d1z + ct1 * (c2 + ct3);
            double g2 = d1z - c2;
            double g3 = ct1 * (1.0 - ct3 - 2.0 * ct2);

            bb[3] = g1 * bb[2] + g2 * bb[1] + g3 * bb[0];
            aa[3] = g1 * aa[2] + g2 * aa[1] + g3 * aa[0];

            double value = aa[3] / bb[3];

            if (Math.Abs(value - aa[0] / bb[0]) < data.eps * Math.Abs(value))
            {
                return value;
            }

            int j;
            for (j = 0; j < 3; j++)
            {
                aa[j] = aa[j + 1];
                bb[j] = bb[j + 1];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("R8_CHU_SCALED - Fatal error!");
        Console.WriteLine("  No convergence after 300 terms.");
        return 1;
    }

    public class r8CIData
    {
        public int nci;
        public double xsml;
        public r8SifgData sifgdata = new();
    }
    public static double r8_ci( ref r8CIData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CI evaluates the cosine integral Ci of an R8 argument.
        //
        //  Discussion:
        //
        //    The cosine integral is defined by
        //
        //      CI(X) = - integral ( X <= T < Infinity ) ( cos ( T ) ) / T  dT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_CI, the cosine integral Ci evaluated at X.
        //
    {
        double[] cics = {
                -0.34004281856055363156281076633129873,
                -1.03302166401177456807159271040163751,
                0.19388222659917082876715874606081709,
                -0.01918260436019865893946346270175301,
                0.00110789252584784967184098099266118,
                -0.00004157234558247208803840231814601,
                0.00000109278524300228715295578966285,
                -0.00000002123285954183465219601280329,
                0.00000000031733482164348544865129873,
                -0.00000000000376141547987683699381798,
                0.00000000000003622653488483964336956,
                -0.00000000000000028911528493651852433,
                0.00000000000000000194327860676494420,
                -0.00000000000000000001115183182650184,
                0.00000000000000000000005527858887706,
                -0.00000000000000000000000023907013943,
                0.00000000000000000000000000091001612,
                -0.00000000000000000000000000000307233,
                0.00000000000000000000000000000000926
            }
            ;
        double f = 0;
        double g = 0;
        double value;
        double y;

        switch (data.nci)
        {
            case 0:
                data.nci = r8_inits(cics, 19, 0.1 * r8_mach(3));
                data.xsml = Math.Sqrt(r8_mach(3));
                break;
        }

        switch (x)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_CI - Fatal error!");
                Console.WriteLine("  X <= 0.0.");
                return 1;
        }

        if (x <= data.xsml)
        {
            y = -1.0;
            value = Math.Log(x) - 0.5 + r8_csevl(y, cics, data.nci);
        }
        else
        {
            switch (x)
            {
                case <= 4.0:
                    y = (x * x - 8.0) * 0.125;
                    value = Math.Log(x) - 0.5 + r8_csevl(y, cics, data.nci);
                    break;
                default:
                    r8_sifg(ref data.sifgdata, x, ref f, ref g);
                    double sinx = Math.Sin(x);
                    value = f * sinx - g * Math.Cos(x);
                    break;
            }
        }

        return value;
    }

    public class r8CINData
    {
        public int ncin;
        public double xmin;
        public r8SifgData sifgdata = new();
    }
        
    public static double r8_cin( ref r8CINData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CIN evaluates the alternate cosine integral Cin of an R8 argument.
        //
        //  Discussion:
        //
        //    CIN(X) = gamma + log(X) 
        //      + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_CIN, the cosine integral Cin evaluated at X.
        //
    {
        double[] cincs = {
                0.37074501750909688741654801228564992,
                -0.05893574896364446831956864397363697,
                0.00538189642113569124048745326203340,
                -0.00029860052841962135319594906563410,
                0.00001095572575321620077031054467306,
                -0.00000028405454877346630491727187731,
                0.00000000546973994875384912457861806,
                -0.00000000008124187461318157083277452,
                0.00000000000095868593117706609013181,
                -0.00000000000000920266004392351031377,
                0.00000000000000007325887999017895024,
                -0.00000000000000000049143726675842909,
                0.00000000000000000000281577746753902,
                -0.00000000000000000000001393986788501,
                0.00000000000000000000000006022485646,
                -0.00000000000000000000000000022904717,
                0.00000000000000000000000000000077273,
                -0.00000000000000000000000000000000233
            }
            ;
        const double eul = 0.57721566490153286060651209008240;
        double f = 0;
        double g = 0;
        double value;

        switch (data.ncin)
        {
            case 0:
                data.ncin = r8_inits(cincs, 18, 0.1 * r8_mach(3));
                data.xmin = Math.Sqrt(r8_mach(1));
                break;
        }

        double absx = Math.Abs(x);

        if (absx <= data.xmin)
        {
            value = 0.0;
        }
        else
        {
            switch (absx)
            {
                case <= 4.0:
                    value = r8_csevl((x * x - 8.0) * 0.125, cincs, data.ncin) * x * x;
                    break;
                default:
                    r8_sifg(ref data.sifgdata, absx, ref f, ref g);
                    double sinx = Math.Sin(absx);
                    value = -f * sinx + g * Math.Cos(absx) + Math.Log(absx) + eul;
                    break;
            }
        }

        return value;
    }

    public class r8CINHData
    {
        public int ncinh;
        public double xmin;
        public double xsml;

    }
    public static double r8_cinh( ref r8CINHData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CINH: alternate hyperbolic cosine integral Cinh of an R8 argument.
        //
        //  Discussion:
        //
        //    Cinh ( x ) = Integral ( 0 <= t <= x ) ( cosh ( t ) - 1 ) dt / t
        //
        //    The original text of this program had a mistake:
        //      y = x * x / 9.0 - 1.0
        //    has been corrected to
        //      y = x * x / 4.5 - 1.0
        //    JVB, 27 March 2010
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_CINH, the hyperbolic cosine integral Cinh
        //    evaluated at X.
        //
    {
        double[] cinhcs = {
                0.1093291636520734431407425199795917,
                0.0573928847550379676445323429825108,
                0.0028095756978830353416404208940774,
                0.0000828780840721356655731765069792,
                0.0000016278596173914185577726018815,
                0.0000000227809519255856619859083591,
                0.0000000002384484842463059257284002,
                0.0000000000019360829780781957471028,
                0.0000000000000125453698328172559683,
                0.0000000000000000663637449497262300,
                0.0000000000000000002919639263594744,
                0.0000000000000000000010849123956107,
                0.0000000000000000000000034499080805,
                0.0000000000000000000000000094936664,
                0.0000000000000000000000000000228291,
                0.0000000000000000000000000000000484
            }
            ;
        const double eul = 0.57721566490153286060651209008240;
        double value;

        switch (data.ncinh)
        {
            case 0:
                data.ncinh = r8_inits(cinhcs, 16, 0.1 * r8_mach(3));
                data.xsml = Math.Sqrt(r8_mach(3));
                data.xmin = 2.0 * Math.Sqrt(r8_mach(1));
                break;
        }

        double absx = Math.Abs(x);

        switch (x)
        {
            case 0.0:
                value = 0.0;
                break;
            default:
            {
                if (absx <= data.xmin)
                {
                    value = 0.0;
                }
                else
                {
                    double y;
                    if (x <= data.xsml)
                    {
                        y = -1.0;
                        value = x * x * (0.25 + r8_csevl(y, cinhcs, data.ncinh));
                    }
                    else
                    {
                        switch (x)
                        {
                            case <= 3.0:
                                y = x * x / 4.5 - 1.0;
                                value = x * x * (0.25 + r8_csevl(y, cinhcs, data.ncinh));
                                break;
                            default:
                                value = r8_chi(absx) - eul - Math.Log(absx);
                                break;
                        }
                    }
                }

                break;
            }
        }

        return value;
    }

    public class r8CosData
    {
        public int ntsn;
        public double xmax;
        public double xsml;
        public double xwarn;
        public r8SqrtData sqrtdata = new();
    }
    public static double r8_cos( ref r8CosData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COS evaluates the cosine of an R8 argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_COS, the cosine of X.
        //
    {
        const double pi2 = 1.57079632679489661923132169163975;
        const double pi2rec = 0.63661977236758134307553505349006;
        const double pihi = 3.140625;
        const double pilo = 9.6765358979323846264338327950288E-04;
        const double pirec = 0.31830988618379067153776752674503;
        double[] sincs = {
                -0.374991154955873175839919279977323464,
                -0.181603155237250201863830316158004754,
                0.005804709274598633559427341722857921,
                -0.000086954311779340757113212316353178,
                0.000000754370148088851481006839927030,
                -0.000000004267129665055961107126829906,
                0.000000000016980422945488168181824792,
                -0.000000000000050120578889961870929524,
                0.000000000000000114101026680010675628,
                -0.000000000000000000206437504424783134,
                0.000000000000000000000303969595918706,
                -0.000000000000000000000000371357734157,
                0.000000000000000000000000000382486123,
                -0.000000000000000000000000000000336623,
                0.000000000000000000000000000000000256
            }
            ;
        double value;

        switch (data.ntsn)
        {
            case 0:
                data.ntsn = r8_inits(sincs, 15, 0.1 * r8_mach(3));
                data.xsml = r8_sqrt( ref data.sqrtdata, 2.0 * r8_mach(3));
                data.xmax = 1.0 / r8_mach(4);
                data.xwarn = r8_sqrt(ref data.sqrtdata, data.xmax);
                break;
        }

        double absx = Math.Abs(x);
        double y = absx + pi2;

        if (data.xmax < y)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_COS - Warning!");
            Console.WriteLine("  No precision because |X| is big.");
            value = 0.0;
            return value;
        }

        if (data.xwarn < y)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_COS - Warning!");
            Console.WriteLine("  Answer < half precision because |X| is big.");
        }

        value = 1.0;

        if (absx < data.xsml)
        {
            return value;
        }

        double xn = (int) (y * pirec + 0.5);
        int n2 = (int) (r8_mod(xn, 2.0) + 0.5);
        xn -= 0.5;
        double f = absx - xn * pihi - xn * pilo;

        xn = 2.0 * (f * pi2rec) * (f * pi2rec) - 1.0;
        value = f + f * r8_csevl(xn, sincs, data.ntsn);

        if (n2 != 0)
        {
            value = -value;
        }

        value = value switch
        {
            < -1.0 => -1.0,
            > 1.0 => 1.0,
            _ => value
        };

        return value;
    }

    public static double r8_cos_deg(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COS_DEG evaluates the cosine of an R8 argument in degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument in degrees.
        //
        //    Output, double R8_COS_DEG, the cosine of X.
        //
    {
        const double raddeg = 0.017453292519943295769236907684886;

        double value = Math.Cos(raddeg * x);

        switch (x % 90.0)
        {
            case 0.0:
            {
                int n = (int) (Math.Abs(x) / 90.0 + 0.5);
                n %= 2;

                switch (n)
                {
                    case 1:
                        value = 0.0;
                        break;
                    default:
                    {
                        value = value switch
                        {
                            < 0.0 => -1.0,
                            _ => +1.0
                        };

                        break;
                    }
                }

                break;
            }
        }

        return value;
    }

    public class r8CoshData
    {
        public double ymax;
    }
        
    public static double r8_cosh( ref r8CoshData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COSH evaluates the hyperbolic cosine of an R8 argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_COSH, the hyperbolic cosine of X.
        //
    {
        double value;

        data.ymax = data.ymax switch
        {
            0.0 => 1.0 / Math.Sqrt(r8_mach(3)),
            _ => data.ymax
        };

        double y = Math.Exp(Math.Abs(x));

        if (y < data.ymax)
        {
            value = 0.5 * (y + 1.0 / y);
        }
        else
        {
            value = 0.5 * y;
        }

        return value;
    }


    public class r8CotData
    {
        public int nterms;
        public double sqeps;
        public double xmax;
        public double xmin;
        public double xsml;

    }
    public static double r8_cot( ref r8CotData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COT evaluates the cotangent of an R8 argument.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Wayne Fullerton.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Wayne Fullerton,
        //    Portable Special Function Routines,
        //    in Portability of Numerical Software,
        //    edited by Wayne Cowell,
        //    Lecture Notes in Computer Science, Volume 57,
        //    Springer 1977,
        //    ISBN: 978-3-540-08446-4,
        //    LC: QA297.W65.
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_COT, the cotangent of X.
        //
    {
        double[] cotcs = {
                +0.240259160982956302509553617744970,
                -0.165330316015002278454746025255758E-01,
                -0.429983919317240189356476228239895E-04,
                -0.159283223327541046023490851122445E-06,
                -0.619109313512934872588620579343187E-09,
                -0.243019741507264604331702590579575E-11,
                -0.956093675880008098427062083100000E-14,
                -0.376353798194580580416291539706666E-16,
                -0.148166574646746578852176794666666E-18,
                -0.583335658903666579477984000000000E-21,
                -0.229662646964645773928533333333333E-23,
                -0.904197057307483326719999999999999E-26,
                -0.355988551920600064000000000000000E-28,
                -0.140155139824298666666666666666666E-30,
                -0.551800436872533333333333333333333E-33
            }
            ;
        const double pi2rec = 0.011619772367581343075535053490057;
        double value;

        switch (data.nterms)
        {
            case 0:
                data.nterms = r8_inits(cotcs, 15, 0.1 * r8_mach(3));
                data.xmax = 1.0 / r8_mach(4);
                data.xsml = Math.Sqrt(3.0 * r8_mach(3));
                data.xmin = Math.Exp(r8_max(Math.Log(r8_mach(1)),
                    -Math.Log(r8_mach(2))) + 0.01);
                data.sqeps = Math.Sqrt(r8_mach(4));
                break;
        }

        double y = Math.Abs(x);

        if (y < data.xmin)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_COT - Fatal error!");
            Console.WriteLine("  |X| is too small.");
            return 1;
        }

        if (data.xmax < y)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_COT - Fatal error!");
            Console.WriteLine("  |X| is too big.");
            return 1;
        }

        //
        //  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
        //  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
        //  = aint(.625*y) + aint(z) + rem(z)
        //
        double ainty = r8_aint(y);
        double yrem = y - ainty;
        double prodbg = 0.625 * ainty;
        ainty = r8_aint(prodbg);
        y = prodbg - ainty + 0.625 * yrem + y * pi2rec;
        double ainty2 = r8_aint(y);
        ainty += ainty2;
        y -= ainty2;

        int ifn = (int) (ainty % 2.0);
        y = ifn switch
        {
            1 => 1.0 - y,
            _ => y
        };

        switch (Math.Abs(x))
        {
            case > 0.5 when y < Math.Abs(x) * data.sqeps:
                Console.WriteLine("");
                Console.WriteLine("R8_COT - Warning!");
                Console.WriteLine("  Answer less than half precision.");
                Console.WriteLine("  |X| too big, or X nearly a nonzero multiple of pi.");
                return 1;
        }

        switch (y)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_COT - Fatal error!");
                Console.WriteLine("  X is a multiple of pi.");
                return 1;
        }

        if (y <= data.xsml)
        {
            value = 1.0 / y;
        }
        else
        {
            switch (y)
            {
                case <= 0.25:
                    value = (0.5 + r8_csevl(32.0 * y * y - 1.0, cotcs, data.nterms)) / y;
                    break;
                case <= 0.5:
                    value = (0.5 + r8_csevl(8.0 * y * y - 1.0,
                        cotcs, data.nterms)) / (0.5 * y);

                    value = (value * value - 1.0) * 0.5 / value;
                    break;
                default:
                    value = (0.5 + r8_csevl(2.0 * y * y - 1.0,
                        cotcs, data.nterms)) / (0.25 * y);
                    value = (value * value - 1.0) * 0.5 / value;
                    value = (value * value - 1.0) * 0.5 / value;
                    break;
            }
        }

        value = ifn switch
        {
            1 => -value,
            _ => x switch
            {
                < 0.0 => -Math.Abs(value),
                _ => +Math.Abs(value)
            }
        };

        return value;
    }

    public static double r8_csevl(double x, double[] a, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CSEVL evaluates a Chebyshev series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2011
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    Volume 16, Number 4, April 1973, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Input, double CS[N], the Chebyshev coefficients.
        //
        //    Input, int N, the number of Chebyshev coefficients.
        //
        //    Output, double R8_CSEVL, the Chebyshev series evaluated at X.
        //
    {
        double b2 = 0;
        int i;

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8_CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms <= 0.");
                return 1;
            case > 1000:
                Console.WriteLine("");
                Console.WriteLine("R8_CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms greater than 1000.");
                return 1;
        }

        switch (x)
        {
            case < -1.1:
            case > 1.1:
                Console.WriteLine("");
                Console.WriteLine("R8_CSEVL - Fatal error!");
                Console.WriteLine("  X outside (-1,+1).");
                return 1;
        }

        double twox = 2.0 * x;
        double b1 = 0.0;
        double b0 = 0.0;

        for (i = n - 1; 0 <= i; i--)
        {
            b2 = b1;
            b1 = b0;
            b0 = twox * b1 - b2 + a[i];
        }

        double value = 0.5 * (b0 - b2);

        return value;
    }
}