﻿using System;

namespace Burkardt.FullertonFnLib
{
    public static partial class FullertonLib
    {
        public static double r8_pak(double y, int n)

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
            double aln2b;
            int nmax = 0;
            int nmin = 0;
            int nsum;
            int ny = 0;
            double value = 0;

            if (nmin == 0)
            {
                aln2b = 1.0;
                if (i4_mach(10) != 2)
                {
                    aln2b = r8_mach(5) * aln210;
                }

                nmin = (int)(aln2b * (double) (i4_mach(15)));
                nmax = (int)(aln2b * (double) (i4_mach(16)));
            }

            r8_upak(y, ref value, ref ny);

            nsum = n + ny;

            if (nsum < nmin)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_PAK - Warning!");
                Console.WriteLine("  Packed number underflows.");
                value = 0.0;
                return value;
            }

            if (nmax < nsum)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_PAK - Fatal error!");
                Console.WriteLine("  Packed number overflows.");
                return (1);
            }

            while (nsum < 0)
            {
                value = 0.5 * value;
                nsum = nsum + 1;
            }

            while (0 < nsum)
            {
                value = 2.0 * value;
                nsum = nsum - 1;
            }

            return value;
        }

        public static double r8_poch(double a, double x)

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
            double absa;
            double absax;
            double alnga = 0;
            double alngax = 0;
            double ax;
            double b;
            double cospia;
            double cospix;
            double den;
            double eps = 0.0;
            double err;
            double errpch;
            int i;
            int ia;
            int n;
            const double pi = 3.141592653589793238462643383279503;
            double sgnga = 0;
            double sgngax = 0;
            double sinpia;
            double sinpix;
            //static double sqeps = 0.0;
            double value;

            if (eps == 0.0)
            {
                eps = r8_mach(4);
                //  sqeps = sqrt ( eps );
            }

            ax = a + x;

            if (ax <= 0.0 && r8_aint(ax) == ax)
            {
                if (0.0 < a || r8_aint(a) != a)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_POCH - Fatal error!");
                    Console.WriteLine("  A + X is nonpositive int,");
                    Console.WriteLine("  but A is not.");
                    return (1);
                }

                //
                //  We know here that both A+X and A are non-positive integers.
                //
                if (x == 0.0)
                {
                    value = 1.0;
                }
                else if (-20.0 < r8_min(a + x, a))
                {
                    n = (int) (x);
                    ia = (int) (a);
                    value = r8_mop(n) * r8_fac(-ia) / r8_fac(-ia - n);
                }
                else
                {
                    n = (int) (x);
                    value = r8_mop(n) * Math.Exp((a - 0.5)
                                            * r8_lnrel(x / (a - 1.0))
                                            + x * Math.Log(-a + 1.0 - x) - x
                                            + r8_lgmc(-a + 1.0)
                                            - r8_lgmc(-a - x + 1.0));
                }

                return value;
            }

            //
            //  A + X is not zero or a negative integer.
            //
            if (a <= 0.0 && r8_aint(a) == a)
            {
                value = 0.0;
                return value;
            }

            n = (int)Math.Abs(x);
            //
            //  x is a small non-positive integer, presummably a common case.
            //
            if ((double) (n) == x && n <= 20)
            {
                value = 1.0;
                for (i = 1; i <= n; i++)
                {
                    value = value * (a + (double) (i - 1));
                }

                return value;
            }

            absax = Math.Abs(a + x);
            absa = Math.Abs(a);

            if (r8_max(absax, absa) <= 20.0)
            {
                value = r8_gamma(a + x) * r8_gamr(a);
                return value;
            }

            if (0.5 * absa < Math.Abs(x))
            {
                r8_lgams(a + x, ref alngax, ref sgngax);
                r8_lgams(a, ref alnga, ref sgnga);
                value = sgngax * sgnga * Math.Exp(alngax - alnga);
                return value;
            }

            //
            //  abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
            //  a+x and a must have the same sign.  for negative a, we use
            //  gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
            //  sin(pi*a)/sin(pi*(a+x))
            //
            if (a < 0.0)
            {
                b = -a - x + 1.0;
            }
            else
            {
                b = a;
            }

            value = Math.Exp((b - 0.5) * r8_lnrel(x / b)
                + x * Math.Log(b + x) - x + r8_lgmc(b + x) - r8_lgmc(b));

            if (0.0 <= a || value == 0.0)
            {
                return value;
            }

            cospix = Math.Cos(pi * x);
            sinpix = Math.Sin(pi * x);
            cospia = Math.Cos(pi * a);
            sinpia = Math.Sin(pi * a);

            errpch = Math.Abs(x) * (1.0 + Math.Log(b));
            den = cospix + cospia * sinpix / sinpia;
            err = (Math.Abs(x) * (Math.Abs(sinpix)
                                  + Math.Abs(cospia * cospix / sinpia))
                   + Math.Abs(a * sinpix) / sinpia / sinpia) * pi;
            err = errpch + err / Math.Abs(den);

            value = value / den;

            return value;
        }

        public static double r8_poch1(double a, double x)

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
            double absa;
            double absx;
            double alneps = 0.0;
            double alnvar;
            double b;
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
            double binv;
            double bp;
            double[] gbern = new double[21];
            double gbk;
            int i;
            int ii;
            int incr;
            int j;
            int k;
            int ndx;
            int nterms;
            const double pi = 3.141592653589793238462643383279503;
            double poly1;
            double q;
            double rho;
            double sinpxx;
            double sinpx2;
            double sqtbig = 0.0;
            double term;
            double trig;
            double value;
            double var;
            double var2;

            if (sqtbig == 0.0)
            {
                sqtbig = 1.0 / Math.Sqrt(24.0 * r8_mach(1));
                alneps = Math.Log(r8_mach(3));
            }

            if (x == 0.0)
            {
                value = r8_psi(a);
                return value;
            }

            absx = Math.Abs(x);
            absa = Math.Abs(a);

            if (0.1 * absa < absx || 0.1 < absx * Math.Log(r8_max(absa, 2.0)))
            {
                value = r8_poch(a, x);
                value = (value - 1.0) / x;
                return value;
            }

            if (a < -0.5)
            {
                bp = 1.0 - a - x;
            }
            else
            {
                bp = a;
            }

            if (bp < 10.0)
            {
                incr = (int)r8_aint(11.0 - bp);
            }
            else
            {
                incr = 0;
            }

            b = bp + (double) (incr);

            var = b + 0.5 * (x - 1.0);
            alnvar = Math.Log(var);
            q = x * alnvar;
            poly1 = 0.0;

            if (var < sqtbig)
            {
                var2 = 1.0 / var / var;

                rho = 0.5 * (x + 1.0);
                gbern[0] = 1.0;
                gbern[1] = -rho / 12.0;
                term = var2;
                poly1 = gbern[1] * term;

                nterms = (int) (-0.5 * alneps / alnvar + 1.0);

                if (20 < nterms)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_POCH1 - Fatal error!");
                    Console.WriteLine(" 20 < NTERMS.");
                    return (1);
                }

                for (k = 2; k <= nterms; k++)
                {
                    gbk = 0.0;
                    for (j = 1; j <= k; j++)
                    {
                        ndx = k - j + 1;
                        gbk = gbk + bern[ndx - 1] * gbern[j - 1];
                    }

                    gbern[k] = -rho * gbk / (double) (k);
                    term = term * ((double) (2 * k - 2) - x)
                                * ((double) (2 * k - 1) - x) * var2;
                    poly1 = poly1 + gbern[k] * term;
                }
            }

            poly1 = (x - 1.0) * poly1;
            value = r8_exprel(q) * (alnvar + q * poly1) + poly1;
            //
            //  we have r8_poch1(b,x), but bp is small, so we use backwards recursion
            //  to obtain r8_poch1(bp,x).
            //
            for (ii = 1; ii <= incr; ii++)
            {
                i = incr - ii;
                binv = 1.0 / (bp + (double) (i));
                value = (value - binv) / (1.0 + x * binv);
            }

            if (bp == a)
            {
                return value;
            }

            //
            //  we have r8_poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
            //  formula to obtain r8_poch1(a,x).
            //
            sinpxx = Math.Sin(pi * x) / x;
            sinpx2 = Math.Sin(0.5 * pi * x);
            trig = sinpxx * r8_cot(pi * b) - 2.0 * sinpx2 * (sinpx2 / x);

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
            double value;

            value = Math.Pow(a, b);

            return value;
        }

        public static double r8_psi(double x)

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
            double aux;
            double dxrel = 0.0;
            int i;
            int n;
            int ntapsi = 0;
            int ntpsi = 0;
            const double pi = 3.14159265358979323846264338327950;
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
            double xbig = 0.0;
            double y;

            if (ntpsi == 0)
            {
                ntpsi = r8_inits(psics, 42, 0.1 * r8_mach(3));
                ntapsi = r8_inits(apsics, 16, 0.1 * r8_mach(3));
                xbig = 1.0 / Math.Sqrt(r8_mach(3));
                dxrel = Math.Sqrt(r8_mach(4));
            }

            y = Math.Abs(x);

            if (y < 10.0)
            {
                n = (int) (x);
                if (x < 0.0)
                {
                    n = n - 1;
                }

                y = x - (double) (n);
                n = n - 1;
                value = r8_csevl(2.0 * y - 1.0, psics, ntpsi);

                if (n == 0)
                {
                    return value;
                }
                else if (n < 0)
                {
                    n = -n;

                    if (x == 0.0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_PSI - Fatal error!");
                        Console.WriteLine("  X is zero.");
                        return (1);
                    }

                    if (x < 0.0 && x + (double) (n - 2) == 0.0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_PSI - Fatal error!");
                        Console.WriteLine("  X is a negative int.");
                        return (1);
                    }

                    if (x < -0.5 && Math.Abs((x - r8_aint(x - 0.5)) / x) < dxrel)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_PSI - Warning!");
                        Console.WriteLine("  Answer is less than half precision");
                        Console.WriteLine("  because X is near a negative int.");
                    }

                    for (i = 1; i <= n; i++)
                    {
                        value = value - 1.0 / (x + (double) (i - 1));
                    }
                }
                else if (0 < n)
                {
                    for (i = 1; i <= n; i++)
                    {
                        value = value + 1.0 / (y + (double) (i));
                    }

                }
            }
            else
            {
                if (y < xbig)
                {
                    aux = r8_csevl(8.0 / y / y - 1.0, apsics, ntapsi);
                }
                else
                {
                    aux = 0.0;
                }

                if (x < 0.0)
                {
                    value = Math.Log(Math.Abs(x)) - 0.5 / x + aux
                            - pi * r8_cot(pi * x);
                }
                else if (0.0 < x)
                {
                    value = Math.Log(x) - 0.5 / x + aux;
                }
            }

            return value;
        }
    }
}