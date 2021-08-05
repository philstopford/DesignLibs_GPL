using System;

namespace Burkardt.FullertonFnLib
{
    public static partial class FullertonLib
    {
        public class r8LBetaData
        {
            public r8LgmcData lgmcdata = new r8LgmcData();
            public r8LngamData lngamdata = new r8LngamData();
            public r8LnrelData lnreldata = new r8LnrelData();
        }
        public static double r8_lbeta(ref r8LBetaData data, ref r8GammaData gdata, double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LBETA evaluates the logarithm of the beta function of R8 arguments.
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
            //    Input, double A, B, the arguments.
            //
            //    Output, double R8_LBETA, the logarithm of the beta function of A
            //    and B.
            //
        {
            double corr;
            double p;
            double q;
            const double sq2pil = 0.91893853320467274178032973640562;
            double value;

            p = r8_min(a, b);
            q = r8_max(a, b);

            if (p <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LBETA - Fatal error!");
                Console.WriteLine("  Both arguments must be greater than 0.");
                return (1);
            }
            else if (p < 10.0 && q <= 10.0)
            {
                value = Math.Log(r8_gamma(ref gdata, p) * (r8_gamma(ref gdata, q) / r8_gamma(ref gdata, p + q)));
            }
            else if (p < 10.0)
            {
                corr = r8_lgmc(ref data.lgmcdata, q) - r8_lgmc(ref data.lgmcdata, p + q);

                value = r8_lngam(ref data.lngamdata, ref gdata, p) + corr + p - p * Math.Log(p + q) +
                        (q - 0.5) * r8_lnrel(ref data.lnreldata, -p / (p + q));
            }
            else
            {
                corr = r8_lgmc(ref data.lgmcdata, p) + r8_lgmc(ref data.lgmcdata, q) - r8_lgmc(ref data.lgmcdata, p + q);

                value = -0.5 * Math.Log(q) + sq2pil + corr
                        + (p - 0.5) * Math.Log(p / (p + q))
                        + q * r8_lnrel(ref data.lnreldata, -p / (p + q));
            }

            return value;
        }

        public class r8LgamsData
        {
            public r8LngamData lngamdata = new r8LngamData();
        }
        public static void r8_lgams(ref r8LgamsData data, ref r8GammaData gdata, double x, ref double algam, ref double sgngam )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LGAMS evaluates the log of |gamma(x)| and sign, for an R8 argument.
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
        //    Output, double &ALGAM, the logarithm of the absolute value of
        //    gamma ( X ).
        //
        //    Output, double &SGNGAM, the sign (+1 or -1) of gamma ( X ).
        //
        {
            int k;

            algam = r8_lngam(ref data.lngamdata, ref gdata, x);
            sgngam = 1.0;

            if (x <= 0.0)
            {
                k = (int) ((-r8_aint(x) % 2.0) + 0.1);

                if (k == 0)
                {
                    sgngam = -1.0;
                }
            }
        }

        public class r8LgicData
        {
            public double eps = 0.0;
        }
        public static double r8_lgic(ref r8LgicData data, double a, double x, double alx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LGIC evaluates the log complementary incomplete gamma function for large X.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 September 2011
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
            //    Output, double R8_LGIC, the log complementary incomplete 
            //    gamma function.
            //
        {
            double fk;
            int k;
            double p;
            double r;
            double s;
            double t;
            double value;
            double xma;
            double xpa;

            if (data.eps == 0.0)
            {
                data.eps = 0.5 * r8_mach(3);
            }

            xpa = x + 1.0 - a;
            xma = x - 1.0 - a;

            r = 0.0;
            p = 1.0;
            s = p;
            for (k = 1; k <= 300; k++)
            {
                fk = (double) (k);
                t = fk * (a - fk) * (1.0 + r);
                r = -t / ((xma + 2.0 * fk) * (xpa + 2.0 * fk) + t);
                p = r * p;
                s = s + p;
                if (Math.Abs(p) < data.eps * s)
                {
                    value = a * alx - x + Math.Log(s / xpa);
                    return value;
                }
            }

            Console.WriteLine("");
            Console.WriteLine("R8_LGIC - Fatal error!");
            Console.WriteLine("  No convergence in 300 iterations.");

            return (1);
        }

        public class r8LgitData
        {
            public double eps = 0.0;
        }
        public static double r8_lgit(ref r8LgitData data, double a, double x, double algap1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LGIT evaluates the log of Tricomi's incomplete gamma function.
            //
            //  Discussion:
            //
            //    Perron's continued fraction is used for large X and X <= A.
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
            //    Input, double ALGAP1, the logarithm of the gamma function of A+1.
            //
            //    Output, double R8_LGIT, the log of Tricomi's incomplete
            //    gamma function.
            //
        {
            double a1x;
            double ax;
            double fk;
            double hstar;
            int k;
            double p;
            double r;
            double s;
            //static double sqeps = 0.0;
            double t;
            double value;

            if (data.eps == 0.0)
            {
                data.eps = 0.5 * r8_mach(3);
                //  sqeps = sqrt ( r8_mach ( 4 ) );
            }

            if (x <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LGIT - Fatal error!");
                Console.WriteLine("  X <= 0.");
                return (1);
            }

            if (a < x)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LGIT - Fatal error!");
                Console.WriteLine("  A < X.");
                return (1);
            }

            ax = a + x;
            a1x = ax + 1.0;
            r = 0.0;
            p = 1.0;
            s = p;

            for (k = 1; k <= 200; k++)
            {
                fk = (double) (k);
                t = (a + fk) * x * (1.0 + r);
                r = t / ((ax + fk) * (a1x + fk) - t);
                p = r * p;
                s = s + p;
                if (Math.Abs(p) < data.eps * s)
                {
                    hstar = 1.0 - x * s / a1x;
                    value = -x - algap1 - Math.Log(hstar);
                    return value;
                }
            }

            Console.WriteLine("");
            Console.WriteLine("R8_LGIT - Fatal error!");
            Console.WriteLine("  No convergence after 200 iterations.");
            return (1);
        }

        public class r8LgmcData
        {
            public int nalgm = 0;
            public double xbig = 0.0;
            public double xmax = 0.0;
        }
        public static double r8_lgmc(ref r8LgmcData data, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LGMC evaluates the log gamma correction factor for an R8 argument.
            //
            //  Discussion:
            //
            //    For 10 <= X, compute the log gamma correction factor so that
            //
            //      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
            //                          + ( x - 0.5 ) * log ( x ) - x 
            //                          + r8_lgmc ( x )
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
            //    Output, double R8_LGMC, the correction factor.
            //
        {
            double[] algmcs = {
                +0.1666389480451863247205729650822,
                -0.1384948176067563840732986059135E-04,
                +0.9810825646924729426157171547487E-08,
                -0.1809129475572494194263306266719E-10,
                +0.6221098041892605227126015543416E-13,
                -0.3399615005417721944303330599666E-15,
                +0.2683181998482698748957538846666E-17,
                -0.2868042435334643284144622399999E-19,
                +0.3962837061046434803679306666666E-21,
                -0.6831888753985766870111999999999E-23,
                +0.1429227355942498147573333333333E-24,
                -0.3547598158101070547199999999999E-26,
                +0.1025680058010470912000000000000E-27,
                -0.3401102254316748799999999999999E-29,
                +0.1276642195630062933333333333333E-30
            }
            ;
            double value;

            if (data.nalgm == 0)
            {
                data.nalgm = r8_inits(algmcs, 15, r8_mach(3));
                data.xbig = 1.0 / Math.Sqrt(r8_mach(3));
                data.xmax = Math.Exp(r8_min(Math.Log(r8_mach(2) / 12.0),
                    -Math.Log(12.0 * r8_mach(1))));
            }

            if (x < 10.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LGMC - Fatal error!");
                Console.WriteLine("  X must be at least 10.");
                return (1);
            }
            else if (x < data.xbig)
            {
                value = r8_csevl(2.0 * (10.0 / x)
                                     * (10.0 / x) - 1.0, algmcs, data.nalgm) / x;
            }
            else if (x < data.xmax)
            {
                value = 1.0 / (12.0 * x);
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public class r8LiData
        {
            public double sqeps = 0.0;

        }
        
        public static double r8_li(ref r8LiData data, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LI evaluates the logarithmic integral for an R8 argument.
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
            //    Output, double R8_LI, the logarithmic integral evaluated at X.
            //
        {
            double value;

            if (data.sqeps == 0.0)
            {
                data.sqeps = Math.Sqrt(r8_mach(3));
            }

            if (x < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LI - Fatal error!");
                Console.WriteLine("  Function undefined for X <= 0.");
                return (1);
            }

            if (x == 0.0)
            {
                value = 0.0;
                return value;
            }

            if (x == 1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LI - Fatal error!");
                Console.WriteLine("  Function undefined for X = 1.");
                return (1);
            }

            if (Math.Abs(1.0 - x) < data.sqeps)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LI - Warning!");
                Console.WriteLine("  Answer less than half precision.");
                Console.WriteLine("  X is too close to 1.");
            }

            value = r8_ei(Math.Log(x));

            return value;
        }

        public class r8LngamData
        {
            public double dxrel = 0.0;
            public double xmax = 0.0;
            public r8LgmcData lgmcdata = new r8LgmcData();
        }
        public static double r8_lngam(ref r8LngamData data, ref r8GammaData gdata, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LNGAM: log of the absolute value of gamma of an R8 argument.
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
            //    Output, double R8_LNGAM, the logarithm of the absolute value of
            //    the gamma function of X.
            //
        {
            const double pi = 3.14159265358979323846264338327950;
            double sinpiy;
            const double sq2pil = 0.91893853320467274178032973640562;
            const double sqpi2l = +0.225791352644727432363097614947441;
            double value;
            double y;

            if (data.xmax == 0.0)
            {
                data.xmax = r8_mach(2) / Math.Log(r8_mach(2));
                data.dxrel = Math.Sqrt(r8_mach(4));
            }

            y = Math.Abs(x);

            if (y <= 10.0)
            {
                value = Math.Log(Math.Abs(r8_gamma(ref gdata, x)));
                return value;
            }

            if (data.xmax < y)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LNGAM - Fatal error!");
                Console.WriteLine("  Result overflows, |X| too big.");
                return (1);
            }

            if (0.0 < x)
            {
                value = sq2pil + (x - 0.5) * Math.Log(x) - x + r8_lgmc(ref data.lgmcdata, y);
                return value;
            }

            sinpiy = Math.Abs(Math.Sin(pi * y));

            if (sinpiy == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LNGAM - Fatal error!");
                Console.WriteLine("  X is a negative int.");
                return (1);
            }

            value = sqpi2l + (x - 0.5) * Math.Log(y) - x - Math.Log(sinpiy) - r8_lgmc(ref data.lgmcdata, y);

            if (Math.Abs((x - r8_aint(x - 0.5)) * value / x) < data.dxrel)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LNGAM - Warning!");
                Console.WriteLine("  Result is half precision because");
                Console.WriteLine("  X is too near a negative int.");
            }

            return value;
        }

        public class r8LnrelData
        {
            public int nlnrel = 0;
            public double xmin = 0.0;

        }
        
        public static double r8_lnrel(ref r8LnrelData data, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LNREL evaluates log ( 1 + X ) for an R8 argument.
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
            //    Output, double R8_LNREL, the value of LOG ( 1 + X ).
            //
        {
            double[] alnrcs = {
                +0.10378693562743769800686267719098E+01,
                -0.13364301504908918098766041553133,
                +0.19408249135520563357926199374750E-01,
                -0.30107551127535777690376537776592E-02,
                +0.48694614797154850090456366509137E-03,
                -0.81054881893175356066809943008622E-04,
                +0.13778847799559524782938251496059E-04,
                -0.23802210894358970251369992914935E-05,
                +0.41640416213865183476391859901989E-06,
                -0.73595828378075994984266837031998E-07,
                +0.13117611876241674949152294345011E-07,
                -0.23546709317742425136696092330175E-08,
                +0.42522773276034997775638052962567E-09,
                -0.77190894134840796826108107493300E-10,
                +0.14075746481359069909215356472191E-10,
                -0.25769072058024680627537078627584E-11,
                +0.47342406666294421849154395005938E-12,
                -0.87249012674742641745301263292675E-13,
                +0.16124614902740551465739833119115E-13,
                -0.29875652015665773006710792416815E-14,
                +0.55480701209082887983041321697279E-15,
                -0.10324619158271569595141333961932E-15,
                +0.19250239203049851177878503244868E-16,
                -0.35955073465265150011189707844266E-17,
                +0.67264542537876857892194574226773E-18,
                -0.12602624168735219252082425637546E-18,
                +0.23644884408606210044916158955519E-19,
                -0.44419377050807936898878389179733E-20,
                +0.83546594464034259016241293994666E-21,
                -0.15731559416479562574899253521066E-21,
                +0.29653128740247422686154369706666E-22,
                -0.55949583481815947292156013226666E-23,
                +0.10566354268835681048187284138666E-23,
                -0.19972483680670204548314999466666E-24,
                +0.37782977818839361421049855999999E-25,
                -0.71531586889081740345038165333333E-26,
                +0.13552488463674213646502024533333E-26,
                -0.25694673048487567430079829333333E-27,
                +0.48747756066216949076459519999999E-28,
                -0.92542112530849715321132373333333E-29,
                +0.17578597841760239233269760000000E-29,
                -0.33410026677731010351377066666666E-30,
                +0.63533936180236187354180266666666E-31
            }
            ;
            double value;

            if (data.nlnrel == 0)
            {
                data.nlnrel = r8_inits(alnrcs, 43, 0.1 * r8_mach(3));
                data.xmin = -1.0 + Math.Sqrt(r8_mach(4));
            }

            if (x <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LNREL - Fatal error!");
                Console.WriteLine("  X <= -1.");
                return (1);
            }
            else if (x < data.xmin)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LNREL - Warning!");
                Console.WriteLine("  Result is less than half precision.");
                Console.WriteLine("  X is too close to - 1.");
            }

            if (Math.Abs(x) <= 0.375)
            {
                value = x * (1.0 - x * r8_csevl(x / 0.375, alnrcs, data.nlnrel));
            }
            else
            {
                value = Math.Log(1.0 + x);
            }

            return value;
        }

        public class r8LogData
        {
            public int nterms = 0;

        }
        public static double r8_log(ref r8LogData data, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LOG evaluates the logarithm of an R8.
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
            //    Input, double X, the evaluation point.
            //
            //    Output, double R8_LOG, the logarithm of X.
            //
        {
            const double aln2 = 0.06814718055994530941723212145818;
            double[] alncen = {
                0.0,
                +0.22314355131420975576629509030983,
                +0.40546510810816438197801311546434,
                +0.55961578793542268627088850052682,
                +0.69314718055994530941723212145817
            }
            ;
            double[] alncs = {
                +0.13347199877973881561689386047187E+01,
                +0.69375628328411286281372438354225E-03,
                +0.42934039020450834506559210803662E-06,
                +0.28933847795432594580466440387587E-09,
                +0.20512517530340580901741813447726E-12,
                +0.15039717055497386574615153319999E-15,
                +0.11294540695636464284521613333333E-18,
                +0.86355788671171868881946666666666E-22,
                +0.66952990534350370613333333333333E-25,
                +0.52491557448151466666666666666666E-28,
                +0.41530540680362666666666666666666E-31
            }
            ;
            double[] center = {
                1.0,
                1.25,
                1.50,
                1.75
            }
            ;
            int n = 0;
            int ntrval;
            double t = 0;
            double t2;
            double value;
            double xn;
            double y = 0;

            if (data.nterms == 0)
            {
                data.nterms = r8_inits(alncs, 11, 28.9 * r8_mach(3));
            }

            if (x <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LOG - Fatal error!");
                Console.WriteLine("  X <= 0.");
                return (1);
            }

            r8_upak(x, ref y, ref n);

            xn = (double) (n - 1);
            y = 2.0 * y;
            ntrval = (int) (4.0 * y - 2.5);

            if (ntrval == 5)
            {
                t = ((y - 1.0) - 1.0) / (y + 2.0);
            }
            else if (ntrval < 5)
            {
                t = (y - center[ntrval - 1]) / (y + center[ntrval - 1]);
            }

            t2 = t * t;
            value = 0.625 * xn + (aln2 * xn + alncen[ntrval - 1] + 2.0 * t
                                  + t * t2 * r8_csevl(578.0 * t2 - 1.0, alncs, data.nterms));

            return value;
        }

        public static double r8_log10(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_LOG10 evaluates the logarithm, base 10, of an R8.
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
            //    Input, double X, the evaluation point.
            //
            //    Output, double R8_LOG10, the logarithm, base 10, of X.
            //
        {
            const double aloge = 0.43429448190325182765112891891661;
            double value;

            value = aloge * Math.Log(x);

            return value;
        }
    }
}