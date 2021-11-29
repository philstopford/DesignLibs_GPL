using System;

namespace Burkardt.FullertonFnLib;

public static partial class FullertonLib
{
    public class r8GAMIData
    {
        public r8GamitData gamitdata = new();
        public r8LngamData lngamdata = new();
    }
    public static double r8_gami(ref r8GAMIData data, ref r8GammaData gdata, double a, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMI evaluates the incomplete gamma function for an R8 argument.
        //
        //  Discussion:
        //
        //    GAMI = Integral ( 0 <= T <= X ) exp ( - t ) * t^( a - 1 )  dt
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
        //    Input, double A, the parameter.
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_GAMI, the value of the incomplete gamma function.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GAMI - Fatal error!");
                Console.WriteLine("  A <= 0.");
                return 1;
            default:
                double value = 0;
                switch (x)
                {
                    case < 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("R8_GAMI - Fatal error!");
                        Console.WriteLine("  X < 0.");
                        return 1;
                    case 0.0:
                        value = 0.0;
                        break;
                    default:
                        double factor = Math.Exp(r8_lngam(ref data.lngamdata, ref gdata, a) + a * Math.Log(x));
                        value = factor * r8_gamit(ref data.gamitdata, ref gdata, a, x);
                        break;
                }

                return value;
        }
    }

    public class r8GamicData
    {
        public double alneps;
        public double eps;
        public r8GMICData gmicdata = new();
        public r8GMITData gmitdata = new();
        public r8LngamData lngamdata = new();
        public r8LgicData lgicdata = new();
        public r8LgitData lgitdata = new();
        public r8LgamsData lgamsdata = new();
    }
        
    public static double r8_gamic( ref r8GamicData data, ref r8GammaData gdata, double a, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMIC evaluates the complementary incomplete gamma function.
        //
        //  Discussion:
        //
        //    GAMIC = integral ( x <= t < oo ) exp(-t) * t^(a-1) dt
        //
        //    GAMIC is evaluated for arbitrary real values of A and non-negative
        //    values X (even though GAMIC is defined for X < 0.0), except that
        //    for X = 0 and A <= 0.0, GAMIC is undefined.
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
        //    Walter Gautschi,
        //    A Computational Procedure for Incomplete Gamma Functions,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 4, December 1979, pages 466-481.
        //
        //  Parameters:
        //
        //    Input, double A, the parameter.
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double R8_GAMIC, the value of the incomplete 
        //    gamma function.
        //
    {
        double algap1 = 0;
        double alngs = 0;
        //static double bot = 0.0;
        //int ma;
        double sgng;
        double sgngam = 0;
        double sgngs = 0;
        //static double sqeps = 0.0;
        double t;
        double value;

        switch (data.eps)
        {
            case 0.0:
                data.eps = 0.5 * r8_mach(3);
                //  sqeps = sqrt ( r8_mach ( 4 ) );
                data.alneps = -Math.Log(r8_mach(3));
                //  bot = log ( r8_mach ( 1 ) );
                break;
        }

        switch (x)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GAMIC - Fatal error!");
                Console.WriteLine("  X < 0.");
                return 1;
            case 0.0 when a <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GAMIC - Fatal error!");
                Console.WriteLine("  X = 0 and A <= 0.");
                return 1;
            case 0.0:
                value = Math.Exp(r8_lngam(ref data.lngamdata, ref gdata, a + 1.0) - Math.Log(a));

                return value;
        }

        double alx = Math.Log(x);
        double sga = a switch
        {
            < 0.0 => -1.0,
            _ => +1.0
        };

        double ainta = r8_aint(a + 0.5 * sga);
        double aeps = a - ainta;

        int izero = 0;

        switch (x)
        {
            case < 1.0:
            {
                switch (a)
                {
                    case <= 0.5 when Math.Abs(aeps) <= 0.001:
                    {
                        double e = -ainta switch
                        {
                            <= 1.0 => 2.0,
                            _ => 2.0 * (-ainta + 2.0) / (ainta * ainta - 1.0)
                        };

                        e -= alx * r8_power(x, -0.001);

                        if (e * Math.Abs(aeps) <= data.eps)
                        {
                            value = r8_gmic(ref data.gmicdata, ref gdata, a, x, alx);
                            return value;
                        }

                        break;
                    }
                }

                r8_lgams(ref data.lgamsdata, ref gdata, a + 1.0, ref algap1, ref sgngam);
                double gstar = r8_gmit(ref data.gmitdata, ref gdata, a, x, algap1, sgngam, alx);

                switch (gstar)
                {
                    case 0.0:
                        izero = 1;
                        break;
                    default:
                        alngs = Math.Log(Math.Abs(gstar));
                        sgngs = r8_sign(gstar);
                        break;
                }

                break;
            }
            default:
            {
                if (a < x)
                {
                    value = Math.Exp(r8_lgic(ref data.lgicdata, a, x, alx));
                    return value;
                }

                sgngam = 1.0;
                algap1 = r8_lngam(ref data.lngamdata, ref gdata, a + 1.0);
                sgngs = 1.0;
                alngs = r8_lgit(ref data.lgitdata, a, x, algap1);
                break;
            }
        }

        double h = 1.0;

        if (izero != 1)
        {
            t = a * alx + alngs;

            if (data.alneps < t)
            {
                sgng = -sgngs * sga * sgngam;
                t = t + algap1 - Math.Log(Math.Abs(a));
                value = sgng * Math.Exp(t);
                return value;
            }

            if (-data.alneps < t)
            {
                h = 1.0 - sgngs * Math.Exp(t);
            }
        }

        sgng = r8_sign(h) * sga * sgngam;
        t = Math.Log(Math.Abs(h)) + algap1 - Math.Log(Math.Abs(a));
        value = sgng * Math.Exp(t);

        return value;
    }

    public class r8GamitData
    {
        public double alneps;
        public r8GMITData gmitdata = new();
        public r8LngamData lngamdata = new();
        public r8LgicData lgicdata = new();
        public r8LgitData lgitdata = new();
        public r8LgamsData lgamsdata = new();
        public r8GamrData gamrdata = new();

    }
    public static double r8_gamit( ref r8GamitData data, ref r8GammaData gdata, double a, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMIT evaluates Tricomi's incomplete gamma function for an R8 argument.
        //
        //  Discussion:
        //
        //      GAMIT = x^(-a) / gamma(a) 
        //        * Integral ( 0 <= t <= x ) exp(-t) * t^(a-1) dt
        //
        //    with analytic continuation for a <= 0.0.  Gamma(x) is the complete
        //    gamma function of X.  GAMIT is evaluated for arbitrary real values of
        //    A and for non-negative values of X (even though GAMIT is defined for
        //    X < 0.0).
        //
        //    A slight deterioration of 2 or 3 digits accuracy will occur when
        //    gamit is very large or very small in absolute value, because log-
        //    arithmic variables are used.  Also, if the parameter A is very close
        //    to a negative integer (but not a negative integer), there is a loss
        //    of accuracy.
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
        //    Walter Gautschi,
        //    A Computational Procedure for Incomplete Gamma Functions,
        //    ACM Transactions on Mathematical Software,
        //    Volume 5, Number 4, December 1979, pages 466-481.
        //
        //  Parameters:
        //
        //    Input, double A, the parameter.
        //
        //    Input, double X, the argument.
        //
        //    Output, double R8_GAMIT, the function value.
        //
    {
        double algap1 = 0;
        double alx;
        //static double bot = 0.0;
        double sgngam = 0;
        //static double sqeps = 0.0;
        double t;
        double value;

        data.alneps = data.alneps switch
        {
            0.0 => -Math.Log(r8_mach(3)),
            _ => data.alneps
        };

        switch (x)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GAMIT - Fatal error!");
                Console.WriteLine("  X is negative.");
                return 1;
            case 0.0:
                alx = 0.0;
                break;
            default:
                alx = Math.Log(x);
                break;
        }

        double sga = a switch
        {
            < 0.0 => -1.0,
            _ => +1.0
        };

        double ainta = r8_aint(a + 0.5 * sga);
        double aeps = a - ainta;

        switch (x)
        {
            case 0.0:
            {
                if (0.0 < ainta || aeps != 0.0)
                {
                    value = r8_gamr(ref data.gamrdata, ref gdata, a + 1.0);
                }
                else
                {
                    value = 0.0;
                }

                return value;
            }
            case <= 1.0:
            {
                if (-0.5 <= a || aeps != 0.0)
                {
                    r8_lgams(ref data.lgamsdata, ref gdata, a + 1.0, ref algap1, ref sgngam);
                }

                value = r8_gmit(ref data.gmitdata, ref gdata, a, x, algap1, sgngam, alx);
                return value;
            }
        }

        if (x <= a)
        {
            t = r8_lgit(ref data.lgitdata, a, x, r8_lngam(ref data.lngamdata, ref gdata, a + 1.0));
            value = Math.Exp(t);
            return value;
        }

        double alng = r8_lgic(ref data.lgicdata, a, x, alx);
        //
        //  Evaluate in terms of log (r8_gamic (a, x))
        //
        double h = 1.0;

        if (aeps != 0.0 || 0.0 < ainta)
        {
            r8_lgams(ref data.lgamsdata, ref gdata, a + 1.0, ref algap1, ref sgngam);
            t = Math.Log(Math.Abs(a)) + alng - algap1;

            if (data.alneps < t)
            {
                t -= a * alx;
                value = -sga * sgngam * Math.Exp(t);
                return value;
            }

            if (-data.alneps < t)
            {
                h = 1.0 - sga * sgngam * Math.Exp(t);
            }
        }

        t = -a * alx + Math.Log(Math.Abs(h));

        value = h switch
        {
            < 0.0 => -Math.Exp(t),
            _ => +Math.Exp(t)
        };

        return value;
    }

    public static void r8_gaml(ref double xmin, ref double xmax )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAML evaluates bounds for an R8 argument of the gamma function.
        //
        //  Discussion:
        //
        //    This function calculates the minimum and maximum legal bounds 
        //    for X in the evaluation of GAMMA ( X ).
        //
        //    XMIN and XMAX are not the only bounds, but they are the only 
        //    non-trivial ones to calculate.
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
        //    Output, double &XMIN, &XMAX, the bounds.
        //
    {
        int i;

        double alnsml = Math.Log(r8_mach(1));
        xmin = -alnsml;

        for (i = 1; i <= 10; i++)
        {
            double xold = xmin;
            double xln = Math.Log(xmin);
            xmin -= xmin * ((xmin + 0.5) * xln - xmin
                                               - 0.2258 + alnsml) / (xmin * xln + 0.5);

            switch (Math.Abs(xmin - xold))
            {
                case < 0.005:
                {
                    xmin = -xmin + 0.01;

                    double alnbig = Math.Log(r8_mach(2));
                    xmax = alnbig;

                    int j;
                    for (j = 1; j <= 10; j++)
                    {
                        xold = xmax;
                        xln = Math.Log(xmax);
                        xmax -= xmax * ((xmax - 0.5) * xln - xmax
                            + 0.9189 - alnbig) / (xmax * xln - 0.5);

                        switch (Math.Abs(xmax - xold))
                        {
                            case < 0.005:
                                xmax -= 0.01;
                                xmin = r8_max(xmin, -xmax + 1.0);
                                return;
                        }
                    }

                    Console.WriteLine("");
                    Console.WriteLine("R8_GAML - Fatal error!");
                    Console.WriteLine("  Unable to find XMAX.");
                    return;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("R8_GAML - Fatal error!");
        Console.WriteLine("  Unable to find XMIN.");
    }

    public class r8GammaData
    {
        public double dxrel;
        public double xmax;
        public double xmin;
        public double xsml;
        public int ngcs;
        public r8LgmcData lgmcdata = new();
    }

    public static double r8_gamma(ref r8GammaData data, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA evaluates the gamma function of an R8 argument.
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
        //    Output, double R8_GAMMA, the gamma function of X.
        //
    {
        double[] gcs = {
                +0.8571195590989331421920062399942E-02,
                +0.4415381324841006757191315771652E-02,
                +0.5685043681599363378632664588789E-01,
                -0.4219835396418560501012500186624E-02,
                +0.1326808181212460220584006796352E-02,
                -0.1893024529798880432523947023886E-03,
                +0.3606925327441245256578082217225E-04,
                -0.6056761904460864218485548290365E-05,
                +0.1055829546302283344731823509093E-05,
                -0.1811967365542384048291855891166E-06,
                +0.3117724964715322277790254593169E-07,
                -0.5354219639019687140874081024347E-08,
                +0.9193275519859588946887786825940E-09,
                -0.1577941280288339761767423273953E-09,
                +0.2707980622934954543266540433089E-10,
                -0.4646818653825730144081661058933E-11,
                +0.7973350192007419656460767175359E-12,
                -0.1368078209830916025799499172309E-12,
                +0.2347319486563800657233471771688E-13,
                -0.4027432614949066932766570534699E-14,
                +0.6910051747372100912138336975257E-15,
                -0.1185584500221992907052387126192E-15,
                +0.2034148542496373955201026051932E-16,
                -0.3490054341717405849274012949108E-17,
                +0.5987993856485305567135051066026E-18,
                -0.1027378057872228074490069778431E-18,
                +0.1762702816060529824942759660748E-19,
                -0.3024320653735306260958772112042E-20,
                +0.5188914660218397839717833550506E-21,
                -0.8902770842456576692449251601066E-22,
                +0.1527474068493342602274596891306E-22,
                -0.2620731256187362900257328332799E-23,
                +0.4496464047830538670331046570666E-24,
                -0.7714712731336877911703901525333E-25,
                +0.1323635453126044036486572714666E-25,
                -0.2270999412942928816702313813333E-26,
                +0.3896418998003991449320816639999E-27,
                -0.6685198115125953327792127999999E-28,
                +0.1146998663140024384347613866666E-28,
                -0.1967938586345134677295103999999E-29,
                +0.3376448816585338090334890666666E-30,
                -0.5793070335782135784625493333333E-31
            }
            ;

        const double sq2pil = 0.91893853320467274178032973640562;
        double value;

        switch (data.ngcs)
        {
            case 0:
                data.ngcs = r8_inits(gcs, 42, 0.1 * r8_mach(3));
                r8_gaml(ref data.xmin, ref data.xmax);
                data.xsml = Math.Exp(r8_max(Math.Log(r8_mach(1)),
                    -Math.Log(r8_mach(2))) + 0.01);
                data.dxrel = Math.Sqrt(r8_mach(4));
                break;
        }

        double y = Math.Abs(x);

        switch (y)
        {
            case <= 10.0:
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
                value = 0.9375 + r8_csevl(2.0 * y - 1.0, gcs, data.ngcs);

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
                                Console.WriteLine("R8_GAMMA - Fatal error!");
                                Console.WriteLine("  X is 0.");
                                return 1;
                            case < 0.0 when x + (n - 2) == 0.0:
                                Console.WriteLine("");
                                Console.WriteLine("R8_GAMMA - Fatal error!");
                                Console.WriteLine("  X is a negative int.");
                                return 1;
                            case < -0.5 when Math.Abs((x - r8_aint(x - 0.5)) / x) < data.dxrel:
                                Console.WriteLine("");
                                Console.WriteLine("R8_GAMMA - Warning!");
                                Console.WriteLine("  X too near a negative int,");
                                Console.WriteLine("  answer is half precision.");
                                break;
                        }

                        if (y < data.xsml)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("R8_GAMMA - Fatal error!");
                            Console.WriteLine("  X is so close to zero that Gamma overflows.");
                            return 1;
                        }

                        for (i = 1; i <= n; i++)
                        {
                            value /= x + (i - 1);
                        }

                        break;
                    }
                    default:
                    {
                        switch (n)
                        {
                            default:
                            {
                                for (i = 1; i <= n; i++)
                                {
                                    value = (y + i) * value;
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
            default:
            {
                if (data.xmax < x)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Fatal error!");
                    Console.WriteLine("  X so big that Gamma overflows.");
                    return 1;
                }

                //
                //  Underflow.
                //
                if (x < data.xmin)
                {
                    value = 0.0;
                    return value;
                }

                value = Math.Exp((y - 0.5) * Math.Log(y) - y + sq2pil + r8_lgmc(ref data.lgmcdata, y));

                switch (x)
                {
                    case > 0.0:
                        return value;
                }

                if (Math.Abs((x - r8_aint(x - 0.5)) / x) < data.dxrel)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Warning!");
                    Console.WriteLine("  X too near a negative int,");
                    Console.WriteLine("  answer is half precision.");
                }

                double sinpiy = Math.Sin(Math.PI * y);

                switch (sinpiy)
                {
                    case 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("R8_GAMMA - Fatal error!");
                        Console.WriteLine("  X is a negative int.");
                        return 1;
                }

                value = -Math.PI / (y * sinpiy * value);
                break;
            }
        }

        return value;
    }

    public class r8GamrData
    {
        public r8LgamsData lgamsdata = new();
    }
    public static double r8_gamr(ref r8GamrData data, ref r8GammaData gdata, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMR evaluates the reciprocal gamma function of an R8 argument.
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
        //    Output, double R8_GAMR, the value of the reciprocal gamma
        //    function at X.
        //
    {
        double alngx = 0;
        double value = 0;
        double sgngx = 0;

        switch (x)
        {
            case <= 0.0 when Math.Abs(r8_aint(x) - x) <= double.Epsilon:
                value = 0.0;
                break;
            default:
            {
                switch (Math.Abs(x))
                {
                    case <= 10.0:
                        value = 1.0 / r8_gamma(ref gdata, x);
                        break;
                    default:
                        r8_lgams(ref data.lgamsdata, ref gdata, x, ref alngx, ref sgngx);
                        value = sgngx * Math.Exp(-alngx);
                        break;
                }

                break;
            }
        }

        return value;
    }

    public class r8GMICData
    {
        public double bot;
        public double eps;
        public r8LngamData lngamdata = new();
    }
        
    public static double r8_gmic( ref r8GMICData data, ref r8GammaData gdata, double a, double x, double alx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GMIC: complementary incomplete gamma, small X, A near negative int.
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
        //    Input, double A, the parameter.
        //
        //    Input, double X, the argument.
        //
        //    Input, double ALX, the logarithm of X.
        //
        //    Output, double R8_GMIC, the complementary incomplete 
        //    gamma function.
        //
    {
        const double euler = 0.57721566490153286060651209008240;
        int k;

        switch (data.eps)
        {
            case 0.0:
                data.eps = 0.5 * r8_mach(3);
                data.bot = Math.Log(r8_mach(1));
                break;
        }

        switch (a)
        {
            case > 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GMIC - Fatal error!");
                Console.WriteLine("  A must be near a negative int.");
                return 1;
        }

        switch (x)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GMIC - Fatal error!");
                Console.WriteLine("  X <= 0.");
                return 1;
        }

        int m = -(int) (a - 0.5);
        double fm = m;

        double te = 1.0;
        double t = 1.0;
        double s = t;
        bool converged = false;

        for (k = 1; k <= 200; k++)
        {
            double fkp1 = k + 1;
            te = -x * te / (fm + fkp1);
            t = te / fkp1;
            s += t;
            if (!(Math.Abs(t) < data.eps * s))
            {
                continue;
            }

            converged = true;
            break;
        }

        switch (converged)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("R8_GMIC - Fatal error!");
                Console.WriteLine("  No convergence after 200 iterations.");
                return 1;
        }

        double value = -alx - euler + x * s / (fm + 1.0);

        switch (m)
        {
            case 0:
                return value;
            case 1:
                value = -value - 1.0 + 1.0 / x;
                return value;
        }

        te = fm;
        t = 1.0;
        s = t;
        int mm1 = m - 1;

        for (k = 1; k <= mm1; k++)
        {
            double fk = k;
            te = -x * te / fk;
            t = te / (fm - fk);
            s += t;
            if (Math.Abs(t) < data.eps * Math.Abs(s))
            {
                break;
            }
        }

        for (k = 1; k <= m; k++)
        {
            value += 1.0 / k;
        }

        double sgng = (m % 2) switch
        {
            1 => -1.0,
            _ => +1.0
        };

        double alng = Math.Log(value) - r8_lngam(ref data.lngamdata, ref gdata, fm + 1.0);

        if (data.bot < alng)
        {
            value = sgng * Math.Exp(alng);
        }
        else
        {
            value = 0.0;
        }

        if (s != 0.0)
        {
            value += r8_sign(s) * Math.Exp(-fm * alx + Math.Log(Math.Abs(s) / fm));
        }

        return value;
    }

    public class r8GMITData
    {
        public double bot;
        public double eps;
        public r8LngamData lngamdata = new();
    }
    public static double r8_gmit( ref r8GMITData data, ref r8GammaData gdata, double a, double x, double algap1, double sgngam, double alx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GMIT: Tricomi's incomplete gamma function for small X.
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
        //    Input, double A, the parameter.
        //
        //    Input, double X, the argument.
        //
        //    Input, double ALGAP1, the logarithm of Gamma ( A + 1 ).
        //
        //    Input, double SGNGAM, the sign of Gamma ( A + 1 ).
        //
        //    Input, double ALX, the logarithm of X.
        //
        //    Output, double R8_GMIT, the Tricomi incomplete gamma function.
        //
    {
        double algs;
        int k;
        double value;

        switch (data.eps)
        {
            case 0.0:
                data.eps = 0.5 * r8_mach(3);
                data.bot = Math.Log(r8_mach(1));
                break;
        }

        switch (x)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8_GMIT - Fatal error!");
                Console.WriteLine("  X <= 0.");
                return 1;
        }

        int ma = a switch
        {
            < 0.0 => (int) (a - 0.5),
            _ => (int) (a + 0.5)
        };

        double aeps = a - ma;

        double ae = a switch
        {
            < -0.5 => aeps,
            _ => a
        };

        double t = 1.0;
        double te = ae;
        double s = t;
        bool converged = false;

        for (k = 1; k <= 200; k++)
        {
            double fk = k;
            te = -x * te / fk;
            t = te / (ae + fk);
            s += t;
            if (!(Math.Abs(t) < data.eps * Math.Abs(s)))
            {
                continue;
            }

            converged = true;
            break;
        }

        switch (converged)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("R8_GMIT - Fatal error!");
                Console.WriteLine("  No convergence in 200 iterations.");
                return 1;
        }

        switch (a)
        {
            case >= -0.5:
                algs = -algap1 + Math.Log(s);
                value = Math.Exp(algs);
                return value;
        }

        algs = -r8_lngam(ref data.lngamdata, ref gdata, 1.0 + aeps) + Math.Log(s);
        s = 1.0;
        int m = -ma - 1;
        t = 1.0;

        for (k = 1; k <= m; k++)
        {
            t = x * t / (aeps - (m + 1 - k));
            s += t;
            if (Math.Abs(t) < data.eps * Math.Abs(s))
            {
                break;
            }
        }

        value = 0.0;
        algs = -(double) ma * Math.Log(x) + algs;

        if (s == 0.0 || aeps == 0.0)
        {
            value = Math.Exp(algs);
            return value;
        }

        double sgng2 = sgngam * r8_sign(s);
        double alg2 = -x - algap1 + Math.Log(Math.Abs(s));

        if (data.bot < alg2)
        {
            value = sgng2 * Math.Exp(alg2);
        }

        if (data.bot < algs)
        {
            value += Math.Exp(algs);
        }

        return value;
    }
}