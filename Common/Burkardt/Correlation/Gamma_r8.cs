using System;
using Burkardt.Types;

namespace Burkardt.CorrelationNS
{
    public static partial class Correlation
    {
        public static double r8_lgmc(double x)

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
            int nalgm = 0;
            double value;
            double xbig = 0.0;
            double xmax = 0.0;

            if (nalgm == 0)
            {
                nalgm = r8_inits(algmcs, 15, typeMethods.r8_mach(3));
                xbig = 1.0 / Math.Sqrt(typeMethods.r8_mach(3));
                xmax = Math.Exp(Math.Min(Math.Log(typeMethods.r8_mach(2) / 12.0),
                    -Math.Log(12.0 * typeMethods.r8_mach(1))));
            }

            if (x < 10.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_LGMC - Fatal error!");
                Console.WriteLine("  X must be at least 10.");
                return (1);
            }
            else if (x < xbig)
            {
                value = r8_csevl(2.0 * (10.0 / x)
                                     * (10.0 / x) - 1.0, algmcs, nalgm) / x;
            }
            else if (x < xmax)
            {
                value = 1.0 / (12.0 * x);
            }
            else
            {
                value = 0.0;
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

            alnsml = Math.Log(typeMethods.r8_mach(1));
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

                    alnbig = Math.Log(typeMethods.r8_mach(2));
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
                            xmin = Math.Max(xmin, -xmax + 1.0);
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

        public static double r8_gamma(double x)

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
            double dxrel = 0.0;
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
            int ngcs = 0;
            double pi = 3.14159265358979323846264338327950;
            double sinpiy;
            double sq2pil = 0.91893853320467274178032973640562;
            double value;
            double xmax = 0.0;
            double xmin = 0.0;
            double xsml = 0.0;
            double y;

            if (ngcs == 0)
            {
                ngcs = r8_inits(gcs, 42, 0.1 * typeMethods.r8_mach(3));
                r8_gaml(ref xmin, ref xmax);
                xsml = Math.Exp(Math.Max(Math.Log(typeMethods.r8_mach(1)),
                    -Math.Log(typeMethods.r8_mach(2))) + 0.01);
                dxrel = Math.Sqrt(typeMethods.r8_mach(4));
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
                value = 0.9375 + r8_csevl(2.0 * y - 1.0, gcs, ngcs);

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

                    if (x < -0.5 && Math.Abs((x - Math.Truncate(x - 0.5)) / x) < dxrel)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8_GAMMA - Warning!");
                        Console.WriteLine("  X too near a negative int,");
                        Console.WriteLine("  answer is half precision.");
                    }

                    if (y < xsml)
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
                if (xmax < x)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Fatal error!");
                    Console.WriteLine("  X so big that Gamma overflows.");
                    return (1);
                }

                //
                //  Underflow.
                //
                if (x < xmin)
                {
                    value = 0.0;
                    return value;
                }

                value = Math.Exp((y - 0.5) * Math.Log(y) - y + sq2pil + r8_lgmc(y));

                if (0.0 < x)
                {
                    return value;
                }

                if (Math.Abs((x - Math.Truncate(x - 0.5)) / x) < dxrel)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Warning!");
                    Console.WriteLine("  X too near a negative int,");
                    Console.WriteLine("  answer is half precision.");
                }

                sinpiy = Math.Sin(pi * y);

                if (sinpiy == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8_GAMMA - Fatal error!");
                    Console.WriteLine("  X is a negative int.");
                    return (1);
                }

                value = -pi / (y * sinpiy * value);
            }

            return value;
        }
    }
}