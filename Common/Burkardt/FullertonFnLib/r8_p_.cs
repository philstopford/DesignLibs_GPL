using System;
using Burkardt.Types;

namespace Burkardt.FullertonFnLib;

public static partial class FullertonLib
{
    public class r8PakData
    {
        public int nmax;
        public int nmin;

    }
    public static double r8_pak(ref r8PakData data, double y, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PAK packs a base 2 exponent into an R8.
        //
        //  Discussion:
        //
        //    This routine is almost the inverse of R8_UPAK.  It is not exactly 
        //    the inverse, because abs(x) need not be between 0.5 and 1.0.  
        //    If both R8_PAK and 2.0^n were known to be in range, we could compute
        //    R8_PAK = x * 2.0^n .
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
        //  Parameters:
        //
        //    Input, double Y, the mantissa.
        //
        //    Input, int N, the exponent.
        //
        //    Output, double R8_PAK, the packed value.
        //
    {
        const double aln210 = 3.321928094887362347870319429489;
        int ny = 0;
        double value = 0;

        switch (data.nmin)
        {
            case 0:
            {
                double aln2b = 1.0;
                if (i4_mach(10) != 2)
                {
                    aln2b = r8_mach(5) * aln210;
                }

                data.nmin = (int)(aln2b * i4_mach(15));
                data.nmax = (int)(aln2b * i4_mach(16));
                break;
            }
        }

        r8_upak(y, ref value, ref ny);

        int nsum = n + ny;

        if (nsum < data.nmin)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_PAK - Warning!");
            Console.WriteLine("  Packed number underflows.");
            value = 0.0;
            return value;
        }

        if (data.nmax < nsum)
        {
            Console.WriteLine("");
            Console.WriteLine("R8_PAK - Fatal error!");
            Console.WriteLine("  Packed number overflows.");
            return 1;
        }

        while (nsum < 0)
        {
            value = 0.5 * value;
            nsum += 1;
        }

        while (0 < nsum)
        {
            value = 2.0 * value;
            nsum -= 1;
        }

        return value;
    }

    public class r8PochData
    {
        public double eps;

        public r8FacData facdata = new();

        public r8LnrelData lnreldata = new();
        public r8LgmcData lgmcdata = new();
        public r8GamrData gamrdata = new();
        public r8LgamsData lgamsdata = new();
    }
    public static double r8_poch(ref r8PochData data, ref r8GammaData gdata, double a, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_POCH evaluates Pochhammer's function of R8 arguments.
        //
        //  Discussion:
        //
        //    POCH ( A, X ) = Gamma ( A + X ) / Gamma ( A ).
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
        //    Input, double A, X, the arguments.
        //
        //    Output, double R8_POCH, the Pochhammer function of A and X.
        //
    {
        double alnga = 0;
        double alngax = 0;
        int n;
            
        double sgnga = 0;
        double sgngax = 0;
        //static double sqeps = 0.0;
        double value = 0;

        data.eps = data.eps switch
        {
            0.0 => r8_mach(4),
            _ => data.eps
        };

        double ax = a + x;

        switch (ax)
        {
            case <= 0.0 when Math.Abs(r8_aint(ax) - ax) <= typeMethods.r8_epsilon():
            {
                if (0.0 < a || Math.Abs(r8_aint(a) - a) > typeMethods.r8_epsilon())
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_POCH - Fatal error!");
                    Console.WriteLine("  A + X is nonpositive int,");
                    Console.WriteLine("  but A is not.");
                    return 1;
                }

                switch (x)
                {
                    //
                    //  We know here that both A+X and A are non-positive integers.
                    //
                    case 0.0:
                        value = 1.0;
                        break;
                    default:
                    {
                        if (-20.0 < r8_min(a + x, a))
                        {
                            n = (int) x;
                            int ia = (int) a;
                            value = r8_mop(n) * r8_fac(ref data.facdata, -ia) / r8_fac(ref data.facdata, -ia - n);
                        }
                        else
                        {
                            n = (int) x;
                            value = r8_mop(n) * Math.Exp((a - 0.5)
                                                         * r8_lnrel(ref data.lnreldata, x / (a - 1.0))
                                                         + x * Math.Log(-a + 1.0 - x) - x
                                                         + r8_lgmc(ref data.lgmcdata, -a + 1.0)
                                                         - r8_lgmc(ref data.lgmcdata, -a - x + 1.0));
                        }

                        break;
                    }
                }

                return value;
            }
        }

        switch (a)
        {
            //
            //  A + X is not zero or a negative integer.
            //
            case <= 0.0 when Math.Abs(r8_aint(a) - a) <= typeMethods.r8_epsilon():
                value = 0.0;
                return value;
        }

        n = (int)Math.Abs(x);
        //
        //  x is a small non-positive integer, presummably a common case.
        //
        if (Math.Abs(n - x) <= typeMethods.r8_epsilon() && n <= 20)
        {
            value = 1.0;
            int i;
            for (i = 1; i <= n; i++)
            {
                value *= a + (i - 1);
            }

            return value;
        }

        double absax = Math.Abs(a + x);
        double absa = Math.Abs(a);

        if (r8_max(absax, absa) <= 20.0)
        {
            value = r8_gamma(ref gdata, a + x) * r8_gamr(ref data.gamrdata, ref gdata, a);
            return value;
        }

        if (0.5 * absa < Math.Abs(x))
        {
            r8_lgams(ref data.lgamsdata, ref gdata, a + x, ref alngax, ref sgngax);
            r8_lgams(ref data.lgamsdata, ref gdata, a, ref alnga, ref sgnga);
            value = sgngax * sgnga * Math.Exp(alngax - alnga);
            return value;
        }

        double b = a switch
        {
            //
            //  abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
            //  a+x and a must have the same sign.  for negative a, we use
            //  gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
            //  sin(pi*a)/sin(pi*(a+x))
            //
            < 0.0 => -a - x + 1.0,
            _ => a
        };

        value = Math.Exp((b - 0.5) * r8_lnrel(ref data.lnreldata, x / b)
            + x * Math.Log(b + x) - x + r8_lgmc(ref data.lgmcdata, b + x) - r8_lgmc(ref data.lgmcdata, b));

        if (0.0 <= a || value == 0.0)
        {
            return value;
        }

        double cospix = Math.Cos(Math.PI * x);
        double sinpix = Math.Sin(Math.PI * x);
        double cospia = Math.Cos(Math.PI * a);
        double sinpia = Math.Sin(Math.PI * a);

        double errpch = Math.Abs(x) * (1.0 + Math.Log(b));
        double den = cospix + cospia * sinpix / sinpia;
        double err = (Math.Abs(x) * (Math.Abs(sinpix)
                                     + Math.Abs(cospia * cospix / sinpia))
                      + Math.Abs(a * sinpix) / sinpia / sinpia) * Math.PI;
        err = errpch + err / Math.Abs(den);

        value /= den;

        return value;
    }

    public class r8Poch1Data
    {
        public double alneps;
        public double sqtbig;

        public r8PsiData psidata = new();
        public r8CotData cotdata = new();
        public r8ExprelData expreldata = new();
        public r8PochData pochdata = new();
    }
    public static double r8_poch1(ref r8Poch1Data data, ref r8GammaData gdata, double a, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_POCH1 evaluates a quantity related to Pochhammer's symbol.
        //
        //  Discussion:
        //
        //    Evaluate a generalization of Pochhammer's symbol for special
        //    situations that require especially accurate values when x is small in
        //      poch1(a,x) = (poch(a,x)-1)/x
        //                 = (gamma(a+x)/gamma(a) - 1.0)/x .
        //    This specification is particularly suited for stably computing
        //    expressions such as
        //      (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
        //           = poch1(a,x) - poch1(b,x)
        //    Note that poch1(a,0.0) = psi(a)
        //
        //    When abs(x) is so small that substantial cancellation will occur if
        //    the straightforward formula is used, we  use an expansion due
        //    to fields and discussed by y. l. luke, the special functions and their
        //    approximations, vol. 1, academic press, 1969, page 34.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2011
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
        //    Input, double A, the parameter.
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double R8_POCH1, the value of the function.
        //
    {
        double[] bern = {
                +0.833333333333333333333333333333333E-01,
                -0.138888888888888888888888888888888E-02,
                +0.330687830687830687830687830687830E-04,
                -0.826719576719576719576719576719576E-06,
                +0.208767569878680989792100903212014E-07,
                -0.528419013868749318484768220217955E-09,
                +0.133825365306846788328269809751291E-10,
                -0.338968029632258286683019539124944E-12,
                +0.858606205627784456413590545042562E-14,
                -0.217486869855806187304151642386591E-15,
                +0.550900282836022951520265260890225E-17,
                -0.139544646858125233407076862640635E-18,
                +0.353470703962946747169322997780379E-20,
                -0.895351742703754685040261131811274E-22,
                +0.226795245233768306031095073886816E-23,
                -0.574472439520264523834847971943400E-24,
                +0.145517247561486490186626486727132E-26,
                -0.368599494066531017818178247990866E-28,
                +0.933673425709504467203255515278562E-30,
                -0.236502241570062993455963519636983E-31
            }
            ;
        double[] gbern = new double[21];
        int ii;

        double value;

        switch (data.sqtbig)
        {
            case 0.0:
                data.sqtbig = 1.0 / Math.Sqrt(24.0 * r8_mach(1));
                data.alneps = Math.Log(r8_mach(3));
                break;
        }

        switch (x)
        {
            case 0.0:
                value = r8_psi(ref data.psidata, a);
                return value;
        }

        double absx = Math.Abs(x);
        double absa = Math.Abs(a);

        if (0.1 * absa < absx || 0.1 < absx * Math.Log(r8_max(absa, 2.0)))
        {
            value = r8_poch(ref data.pochdata, ref gdata, a, x);
            value = (value - 1.0) / x;
            return value;
        }

        double bp = a switch
        {
            < -0.5 => 1.0 - a - x,
            _ => a
        };

        int incr = bp switch
        {
            < 10.0 => (int) r8_aint(11.0 - bp),
            _ => 0
        };

        double b = bp + incr;

        double var = b + 0.5 * (x - 1.0);
        double alnvar = Math.Log(var);
        double q = x * alnvar;
        double poly1 = 0.0;

        if (var < data.sqtbig)
        {
            double var2 = 1.0 / var / var;

            double rho = 0.5 * (x + 1.0);
            gbern[0] = 1.0;
            gbern[1] = -rho / 12.0;
            double term = var2;
            poly1 = gbern[1] * term;

            int nterms = (int) (-0.5 * data.alneps / alnvar + 1.0);

            switch (nterms)
            {
                case > 20:
                    Console.WriteLine("");
                    Console.WriteLine("R8_POCH1 - Fatal error!");
                    Console.WriteLine(" 20 < NTERMS.");
                    return 1;
            }

            int k;
            for (k = 2; k <= nterms; k++)
            {
                double gbk = 0.0;
                int j;
                for (j = 1; j <= k; j++)
                {
                    int ndx = k - j + 1;
                    gbk += bern[ndx - 1] * gbern[j - 1];
                }

                gbern[k] = -rho * gbk / k;
                term = term * (2 * k - 2 - x)
                            * (2 * k - 1 - x) * var2;
                poly1 += gbern[k] * term;
            }
        }

        poly1 = (x - 1.0) * poly1;
        value = r8_exprel(ref data.expreldata, q) * (alnvar + q * poly1) + poly1;
        //
        //  we have r8_poch1(b,x), but bp is small, so we use backwards recursion
        //  to obtain r8_poch1(bp,x).
        //
        for (ii = 1; ii <= incr; ii++)
        {
            int i = incr - ii;
            double binv = 1.0 / (bp + i);
            value = (value - binv) / (1.0 + x * binv);
        }

        if (Math.Abs(bp - a) <= typeMethods.r8_epsilon())
        {
            return value;
        }

        //
        //  we have r8_poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
        //  formula to obtain r8_poch1(a,x).
        //
        double sinpxx = Math.Sin(Math.PI * x) / x;
        double sinpx2 = Math.Sin(0.5 * Math.PI * x);
        double trig = sinpxx * r8_cot(ref data.cotdata, Math.PI * b) - 2.0 * sinpx2 * (sinpx2 / x);

        value = trig + (1.0 + x * trig) * value;

        return value;
    }

    public static double r8_power(double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_POWER evaluates A^B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 August 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the base.
        //
        //    Input, double B, the exponent.
        //
        //    Output, double R8_POWER, the value of A^B.
        //
    {
        double value = 0;

        value = Math.Pow(a, b);

        return value;
    }

    public class r8PsiData
    {
        public double dxrel;
        public int ntapsi;
        public int ntpsi;
        public double xbig;

        public r8CotData cotdata = new();
    }
        
    public static double r8_psi( ref r8PsiData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PSI evaluates the psi function of an R8 argument.
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
        //    Output, double R8_PSI, the psi function of X.
        //
    {
        double[] apsics = {
                -0.832710791069290760174456932269E-03,
                -0.416251842192739352821627121990E-03,
                +0.103431560978741291174463193961E-06,
                -0.121468184135904152987299556365E-09,
                +0.311369431998356155521240278178E-12,
                -0.136461337193177041776516100945E-14,
                +0.902051751315416565130837974000E-17,
                -0.831542997421591464829933635466E-19,
                +0.101224257073907254188479482666E-20,
                -0.156270249435622507620478933333E-22,
                +0.296542716808903896133226666666E-24,
                -0.674686886765702163741866666666E-26,
                +0.180345311697189904213333333333E-27,
                -0.556901618245983607466666666666E-29,
                +0.195867922607736251733333333333E-30,
                -0.775195892523335680000000000000E-32
            }
            ;

        double[] psics = {
                -0.38057080835217921520437677667039E-01,
                +0.49141539302938712748204699654277,
                -0.56815747821244730242892064734081E-01,
                +0.83578212259143131362775650747862E-02,
                -0.13332328579943425998079274172393E-02,
                +0.22031328706930824892872397979521E-03,
                -0.37040238178456883592889086949229E-04,
                +0.62837936548549898933651418717690E-05,
                -0.10712639085061849855283541747074E-05,
                +0.18312839465484165805731589810378E-06,
                -0.31353509361808509869005779796885E-07,
                +0.53728087762007766260471919143615E-08,
                -0.92116814159784275717880632624730E-09,
                +0.15798126521481822782252884032823E-09,
                -0.27098646132380443065440589409707E-10,
                +0.46487228599096834872947319529549E-11,
                -0.79752725638303689726504797772737E-12,
                +0.13682723857476992249251053892838E-12,
                -0.23475156060658972717320677980719E-13,
                +0.40276307155603541107907925006281E-14,
                -0.69102518531179037846547422974771E-15,
                +0.11856047138863349552929139525768E-15,
                -0.20341689616261559308154210484223E-16,
                +0.34900749686463043850374232932351E-17,
                -0.59880146934976711003011081393493E-18,
                +0.10273801628080588258398005712213E-18,
                -0.17627049424561071368359260105386E-19,
                +0.30243228018156920457454035490133E-20,
                -0.51889168302092313774286088874666E-21,
                +0.89027730345845713905005887487999E-22,
                -0.15274742899426728392894971904000E-22,
                +0.26207314798962083136358318079999E-23,
                -0.44964642738220696772598388053333E-24,
                +0.77147129596345107028919364266666E-25,
                -0.13236354761887702968102638933333E-25,
                +0.22709994362408300091277311999999E-26,
                -0.38964190215374115954491391999999E-27,
                +0.66851981388855302310679893333333E-28,
                -0.11469986654920864872529919999999E-28,
                +0.19679385886541405920515413333333E-29,
                -0.33764488189750979801907200000000E-30,
                +0.57930703193214159246677333333333E-31
            }
            ;
        double value = 0;

        switch (data.ntpsi)
        {
            case 0:
                data.ntpsi = r8_inits(psics, 42, 0.1 * r8_mach(3));
                data.ntapsi = r8_inits(apsics, 16, 0.1 * r8_mach(3));
                data.xbig = 1.0 / Math.Sqrt(r8_mach(3));
                data.dxrel = Math.Sqrt(r8_mach(4));
                break;
        }

        double y = Math.Abs(x);

        switch (y)
        {
            case < 10.0:
            {
                int n = (int) x;
                switch (x)
                {
                    case < 0.0:
                        n -= 1;
                        break;
                }

                y = x - n;
                n -= 1;
                value = r8_csevl(2.0 * y - 1.0, psics, data.ntpsi);

                int i;
                switch (n)
                {
                    case 0:
                        return value;
                    case < 0:
                    {
                        n = -n;

                        switch (x)
                        {
                            case 0.0:
                                Console.WriteLine("");
                                Console.WriteLine("R8_PSI - Fatal error!");
                                Console.WriteLine("  X is zero.");
                                return 1;
                            case < 0.0 when x + (n - 2) == 0.0:
                                Console.WriteLine("");
                                Console.WriteLine("R8_PSI - Fatal error!");
                                Console.WriteLine("  X is a negative int.");
                                return 1;
                            case < -0.5 when Math.Abs((x - r8_aint(x - 0.5)) / x) < data.dxrel:
                                Console.WriteLine("");
                                Console.WriteLine("R8_PSI - Warning!");
                                Console.WriteLine("  Answer is less than half precision");
                                Console.WriteLine("  because X is near a negative int.");
                                break;
                        }

                        for (i = 1; i <= n; i++)
                        {
                            value -= 1.0 / (x + (i - 1));
                        }

                        break;
                    }
                    case > 0:
                    {
                        for (i = 1; i <= n; i++)
                        {
                            value += 1.0 / (y + i);
                        }

                        break;
                    }
                }

                break;
            }
            default:
            {
                double aux = y < data.xbig ? r8_csevl(8.0 / y / y - 1.0, apsics, data.ntapsi) : 0.0;

                value = x switch
                {
                    < 0.0 => Math.Log(Math.Abs(x)) - 0.5 / x + aux - Math.PI * r8_cot(ref data.cotdata, Math.PI * x),
                    > 0.0 => Math.Log(x) - 0.5 / x + aux,
                    _ => value
                };

                break;
            }
        }

        return value;
    }
}