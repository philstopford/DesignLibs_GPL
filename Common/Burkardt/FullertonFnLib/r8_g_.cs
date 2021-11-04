using System;

namespace Burkardt.FullertonFnLib
{
    public static partial class FullertonLib
    {
        public class r8GAMIData
        {
            public r8GamitData gamitdata = new r8GamitData();
            public r8LngamData lngamdata = new r8LngamData();
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
            double factor;
            double value;

            if (a <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GAMI - Fatal error!");
                Console.WriteLine("  A <= 0.");
                return (1);
            }

            if (x < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GAMI - Fatal error!");
                Console.WriteLine("  X < 0.");
                return (1);
            }
            else if (x == 0.0)
            {
                value = 0.0;
            }
            else
            {
                factor = Math.Exp(r8_lngam(ref data.lngamdata, ref gdata, a) + a * Math.Log(x));
                value = factor * r8_gamit(ref data.gamitdata, ref gdata, a, x);
            }

            return value;
        }

        public class r8GamicData
        {
            public double alneps = 0.0;
            public double eps = 0.0;
            public r8GMICData gmicdata = new r8GMICData();
            public r8GMITData gmitdata = new r8GMITData();
            public r8LngamData lngamdata = new r8LngamData();
            public r8LgicData lgicdata = new r8LgicData();
            public r8LgitData lgitdata = new r8LgitData();
            public r8LgamsData lgamsdata = new r8LgamsData();
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
            double aeps;
            double ainta;
            double algap1 = 0;
            double alngs = 0;
            double alx;
            //static double bot = 0.0;
            double e;
            double gstar;
            double h;
            int izero;
            //int ma;
            double sga;
            double sgng;
            double sgngam = 0;
            double sgngs = 0;
            //static double sqeps = 0.0;
            double t;
            double value;

            if (data.eps == 0.0)
            {
                data.eps = 0.5 * r8_mach(3);
                //  sqeps = sqrt ( r8_mach ( 4 ) );
                data.alneps = -Math.Log(r8_mach(3));
                //  bot = log ( r8_mach ( 1 ) );
            }

            if (x < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GAMIC - Fatal error!");
                Console.WriteLine("  X < 0.");
                return (1);
            }

            if (x == 0.0)
            {
                if (a <= 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMIC - Fatal error!");
                    Console.WriteLine("  X = 0 and A <= 0.");
                    return (1);
                }

                value = Math.Exp(r8_lngam(ref data.lngamdata, ref gdata, a + 1.0) - Math.Log(a));

                return value;
            }

            alx = Math.Log(x);
            if (a < 0.0)
            {
                sga = -1.0;
            }
            else
            {
                sga = +1.0;
            }

            ainta = r8_aint(a + 0.5 * sga);
            aeps = a - ainta;

            izero = 0;

            if (x < 1.0)
            {
                if (a <= 0.5 && Math.Abs(aeps) <= 0.001)
                {
                    if (-ainta <= 1.0)
                    {
                        e = 2.0;
                    }
                    else
                    {
                        e = 2.0 * (-ainta + 2.0) / (ainta * ainta - 1.0);
                    }

                    e = e - alx * r8_power(x, -0.001);

                    if (e * Math.Abs(aeps) <= data.eps)
                    {
                        value = r8_gmic(ref data.gmicdata, ref gdata, a, x, alx);
                        return value;
                    }
                }

                r8_lgams(ref data.lgamsdata, ref gdata, a + 1.0, ref algap1, ref sgngam);
                gstar = r8_gmit(ref data.gmitdata, ref gdata, a, x, algap1, sgngam, alx);

                if (gstar == 0.0)
                {
                    izero = 1;
                }
                else
                {
                    alngs = Math.Log(Math.Abs(gstar));
                    sgngs = r8_sign(gstar);
                }
            }
            else
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
            }

            h = 1.0;

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
            public double alneps = 0.0;
            public r8GMITData gmitdata = new r8GMITData();
            public r8LngamData lngamdata = new r8LngamData();
            public r8LgicData lgicdata = new r8LgicData();
            public r8LgitData lgitdata = new r8LgitData();
            public r8LgamsData lgamsdata = new r8LgamsData();
            public r8GamrData gamrdata = new r8GamrData();

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
            double aeps;
            double ainta;
            double algap1 = 0;
            double alng;
            double alx;
            //static double bot = 0.0;
            double h;
            double sga;
            double sgngam = 0;
            //static double sqeps = 0.0;
            double t;
            double value;

            if (data.alneps == 0.0)
            {
                data.alneps = -Math.Log(r8_mach(3));
                //  sqeps = sqrt ( r8_mach ( 4 ) );
                //  bot = log ( r8_mach ( 1 ) );
            }

            if (x < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GAMIT - Fatal error!");
                Console.WriteLine("  X is negative.");
                return (1);
            }
            else if (x == 0.0)
            {
                alx = 0.0;
            }
            else
            {
                alx = Math.Log(x);
            }

            if (a < 0.0)
            {
                sga = -1.0;
            }
            else
            {
                sga = +1.0;
            }

            ainta = r8_aint(a + 0.5 * sga);
            aeps = a - ainta;

            if (x == 0.0)
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

            if (x <= 1.0)
            {
                if (-0.5 <= a || aeps != 0.0)
                {
                    r8_lgams(ref data.lgamsdata, ref gdata, a + 1.0, ref algap1, ref sgngam);
                }

                value = r8_gmit(ref data.gmitdata, ref gdata, a, x, algap1, sgngam, alx);
                return value;
            }

            if (x <= a)
            {
                t = r8_lgit(ref data.lgitdata, a, x, r8_lngam(ref data.lngamdata, ref gdata, a + 1.0));
                value = Math.Exp(t);
                return value;
            }

            alng = r8_lgic(ref data.lgicdata, a, x, alx);
            //
            //  Evaluate in terms of log (r8_gamic (a, x))
            //
            h = 1.0;

            if (aeps != 0.0 || 0.0 < ainta)
            {
                r8_lgams(ref data.lgamsdata, ref gdata, a + 1.0, ref algap1, ref sgngam);
                t = Math.Log(Math.Abs(a)) + alng - algap1;

                if (data.alneps < t)
                {
                    t = t - a * alx;
                    value = -sga * sgngam * Math.Exp(t);
                    return value;
                }

                if (-data.alneps < t)
                {
                    h = 1.0 - sga * sgngam * Math.Exp(t);
                }
            }

            t = -a * alx + Math.Log(Math.Abs(h));

            if (h < 0.0)
            {
                value = -Math.Exp(t);
            }
            else
            {
                value = +Math.Exp(t);
            }

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
            double alnbig;
            double alnsml;
            int i;
            int j;
            double xln;
            double xold;

            alnsml = Math.Log(r8_mach(1));
            xmin = -alnsml;

            for (i = 1; i <= 10; i++)
            {
                xold = xmin;
                xln = Math.Log(xmin);
                xmin = xmin - xmin * ((xmin + 0.5) * xln - xmin
                                                         - 0.2258 + alnsml) / (xmin * xln + 0.5);

                if (Math.Abs(xmin - xold) < 0.005)
                {
                    xmin = -xmin + 0.01;

                    alnbig = Math.Log(r8_mach(2));
                    xmax = alnbig;

                    for (j = 1; j <= 10; j++)
                    {
                        xold = xmax;
                        xln = Math.Log(xmax);
                        xmax = xmax - xmax * ((xmax - 0.5) * xln - xmax
                            + 0.9189 - alnbig) / (xmax * xln - 0.5);

                        if (Math.Abs(xmax - xold) < 0.005)
                        {
                            xmax = xmax - 0.01;
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

            Console.WriteLine("");
            Console.WriteLine("R8_GAML - Fatal error!");
            Console.WriteLine("  Unable to find XMIN.");
            return;
        }

        public class r8GammaData
        {
            public double dxrel = 0.0;
            public double xmax = 0.0;
            public double xmin = 0.0;
            public double xsml = 0.0;
            public int ngcs = 0;
            public r8LgmcData lgmcdata = new r8LgmcData();
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
            int i;
            int n;
            
            double sinpiy;
            const double sq2pil = 0.91893853320467274178032973640562;
            double value;
            double y;

            if (data.ngcs == 0)
            {
                data.ngcs = r8_inits(gcs, 42, 0.1 * r8_mach(3));
                r8_gaml(ref data.xmin, ref data.xmax);
                data.xsml = Math.Exp(r8_max(Math.Log(r8_mach(1)),
                    -Math.Log(r8_mach(2))) + 0.01);
                data.dxrel = Math.Sqrt(r8_mach(4));
            }

            y = Math.Abs(x);

            if (y <= 10.0)
            {
                n = (int) (x);
                if (x < 0.0)
                {
                    n = n - 1;
                }

                y = x - (double) (n);
                n = n - 1;
                value = 0.9375 + r8_csevl(2.0 * y - 1.0, gcs, data.ngcs);

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
                        Console.WriteLine("R8_GAMMA - Fatal error!");
                        Console.WriteLine("  X is 0.");
                        return (1);
                    }

                    if (x < 0.0 && x + (double) (n - 2) == 0.0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_GAMMA - Fatal error!");
                        Console.WriteLine("  X is a negative int.");
                        return (1);
                    }

                    if (x < -0.5 && Math.Abs((x - r8_aint(x - 0.5)) / x) < data.dxrel)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_GAMMA - Warning!");
                        Console.WriteLine("  X too near a negative int,");
                        Console.WriteLine("  answer is half precision.");
                    }

                    if (y < data.xsml)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_GAMMA - Fatal error!");
                        Console.WriteLine("  X is so close to zero that Gamma overflows.");
                        return (1);
                    }

                    for (i = 1; i <= n; i++)
                    {
                        value = value / (x + (double) (i - 1));
                    }

                }
                else if (n == 0)
                {
                }
                else
                {
                    for (i = 1; i <= n; i++)
                    {
                        value = (y + (double) (i)) * value;
                    }
                }
            }
            else
            {
                if (data.xmax < x)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Fatal error!");
                    Console.WriteLine("  X so big that Gamma overflows.");
                    return (1);
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

                if (0.0 < x)
                {
                    return value;
                }

                if (Math.Abs((x - r8_aint(x - 0.5)) / x) < data.dxrel)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Warning!");
                    Console.WriteLine("  X too near a negative int,");
                    Console.WriteLine("  answer is half precision.");
                }

                sinpiy = Math.Sin(Math.PI * y);

                if (sinpiy == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Fatal error!");
                    Console.WriteLine("  X is a negative int.");
                    return (1);
                }

                value = -Math.PI / (y * sinpiy * value);
            }

            return value;
        }

        public class r8GamrData
        {
            public r8LgamsData lgamsdata = new r8LgamsData();
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
            double value;
            double sgngx = 0;

            if (x <= 0.0 && r8_aint(x) == x)
            {
                value = 0.0;
            }
            else if (Math.Abs(x) <= 10.0)
            {
                value = 1.0 / r8_gamma(ref gdata, x);
            }
            else
            {
                r8_lgams(ref data.lgamsdata, ref gdata, x, ref alngx, ref sgngx);
                value = sgngx * Math.Exp(-alngx);
            }

            return value;
        }

        public class r8GMICData
        {
            public double bot = 0.0;
            public double eps = 0.0;
            public r8LngamData lngamdata = new r8LngamData();
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
            double alng;
            bool converged;
            const double euler = 0.57721566490153286060651209008240;
            double fk;
            double fkp1;
            double fm;
            int k;
            int m;
            int mm1;
            double s;
            double sgng;
            double t;
            double te;
            double value;

            if (data.eps == 0.0)
            {
                data.eps = 0.5 * r8_mach(3);
                data.bot = Math.Log(r8_mach(1));
            }

            if (0.0 < a)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GMIC - Fatal error!");
                Console.WriteLine("  A must be near a negative int.");
                return (1);
            }

            if (x <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GMIC - Fatal error!");
                Console.WriteLine("  X <= 0.");
                return (1);
            }

            m = -(int) (a - 0.5);
            fm = (double) (m);

            te = 1.0;
            t = 1.0;
            s = t;
            converged = false;

            for (k = 1; k <= 200; k++)
            {
                fkp1 = (double) (k + 1);
                te = -x * te / (fm + fkp1);
                t = te / fkp1;
                s = s + t;
                if (Math.Abs(t) < data.eps * s)
                {
                    converged = true;
                    break;
                }
            }

            if (!converged)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GMIC - Fatal error!");
                Console.WriteLine("  No convergence after 200 iterations.");
                return (1);
            }

            value = -alx - euler + x * s / (fm + 1.0);

            if (m == 0)
            {
                return value;
            }
            else if (m == 1)
            {
                value = -value - 1.0 + 1.0 / x;
                return value;
            }

            te = fm;
            t = 1.0;
            s = t;
            mm1 = m - 1;

            for (k = 1; k <= mm1; k++)
            {
                fk = (double) (k);
                te = -x * te / fk;
                t = te / (fm - fk);
                s = s + t;
                if (Math.Abs(t) < data.eps * Math.Abs(s))
                {
                    break;
                }
            }

            for (k = 1; k <= m; k++)
            {
                value = value + 1.0 / (double) (k);
            }

            if ((m % 2) == 1)
            {
                sgng = -1.0;
            }
            else
            {
                sgng = +1.0;
            }

            alng = Math.Log(value) - r8_lngam(ref data.lngamdata, ref gdata, fm + 1.0);

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
                value = value
                        + r8_sign(s) * Math.Exp(-fm * alx + Math.Log(Math.Abs(s) / fm));
            }

            return value;
        }

        public class r8GMITData
        {
            public double bot = 0.0;
            public double eps = 0.0;
            public r8LngamData lngamdata = new r8LngamData();
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
            double ae;
            double aeps;
            double alg2;
            double algs;
            bool converged;
            double fk;
            int k;
            int m;
            int ma;
            double s;
            double sgng2;
            double t;
            double te;
            double value;

            if (data.eps == 0.0)
            {
                data.eps = 0.5 * r8_mach(3);
                data.bot = Math.Log(r8_mach(1));
            }

            if (x <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GMIT - Fatal error!");
                Console.WriteLine("  X <= 0.");
                return (1);
            }

            if (a < 0.0)
            {
                ma = (int) (a - 0.5);
            }
            else
            {
                ma = (int) (a + 0.5);
            }

            aeps = a - (double) (ma);

            if (a < -0.5)
            {
                ae = aeps;
            }
            else
            {
                ae = a;
            }

            t = 1.0;
            te = ae;
            s = t;
            converged = false;

            for (k = 1; k <= 200; k++)
            {
                fk = (double) (k);
                te = -x * te / fk;
                t = te / (ae + fk);
                s = s + t;
                if (Math.Abs(t) < data.eps * Math.Abs(s))
                {
                    converged = true;
                    break;
                }
            }

            if (!converged)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_GMIT - Fatal error!");
                Console.WriteLine("  No convergence in 200 iterations.");
                return (1);
            }

            if (-0.5 <= a)
            {
                algs = -algap1 + Math.Log(s);
                value = Math.Exp(algs);
                return value;
            }

            algs = -r8_lngam(ref data.lngamdata, ref gdata, 1.0 + aeps) + Math.Log(s);
            s = 1.0;
            m = -ma - 1;
            t = 1.0;

            for (k = 1; k <= m; k++)
            {
                t = x * t / (aeps - (double) (m + 1 - k));
                s = s + t;
                if (Math.Abs(t) < data.eps * Math.Abs(s))
                {
                    break;
                }
            }

            value = 0.0;
            algs = -(double) (ma) * Math.Log(x) + algs;

            if (s == 0.0 || aeps == 0.0)
            {
                value = Math.Exp(algs);
                return value;
            }

            sgng2 = sgngam * r8_sign(s);
            alg2 = -x - algap1 + Math.Log(Math.Abs(s));

            if (data.bot < alg2)
            {
                value = sgng2 * Math.Exp(alg2);
            }

            if (data.bot < algs)
            {
                value = value + Math.Exp(algs);
            }

            return value;
        }
    }
}