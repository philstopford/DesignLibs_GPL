using System;
using Burkardt;
using Burkardt.CDFLib;
using Burkardt.TestValues;

namespace CDFLibTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdflib_test() tests cdflib().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("cdflib_test():");
            Console.WriteLine("  Test cdflib().");

            beta_inc_test();
            cdfbet_test();
            cdfbin_test();
            cdfchi_test();
            cdfchn_test();
            cdff_test();
            cdffnc_test();
            cdfgam_test();
            cdfnbn_test();
            cdfnor_test();

            cdfpoi_test();
            cdft_test();
            cumbet_test();
            cumbin_test();
            cumchi_test();
            cumchn_test();
            cumf_test();
            cumfnc_test();
            cumgam_test();
            cumnbn_test();

            cumnor_test();
            cumpoi_test();
            cumt_test();
            beta_test();
            error_f_test();
            gamma_test();
            gamma_inc_test();
            psi_test();

            Console.WriteLine("");
            Console.WriteLine("cdflib_test():");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void beta_inc_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   beta_inc_test tests BETA_INC and BETA_INC_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double ccdf_compute = 0;
            double ccdf_lookup = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int ierror = 0;
            int n_data = 0;
            double x = 0;
            double y = 0;

            Console.WriteLine("");
            Console.WriteLine("beta_inc_test");
            Console.WriteLine("  BETA_INC computes the incomplete Beta ratio.");
            Console.WriteLine("  BETA_INC_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            string cout = "    X         Y         A         B         CDF           CDF";
            cout +=
                "                                           (Lookup)      (Computed)";
            Console.WriteLine(cout);

            n_data = 0;

            for (;;)
            {
                CDF.beta_inc_values(ref n_data, ref a, ref b, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                y = 1.0 - x;

                CDF.beta_inc(a, b, x, y, ref cdf_compute, ref ccdf_compute, ref ierror);

                Console.WriteLine(x.ToString("0.######").PadLeft(10)
                    + y.ToString("0.######").PadLeft(10)
                    + a.ToString("0.######").PadLeft(10)
                    + b.ToString("0.######").PadLeft(10)
                    + cdf_lookup.ToString("0.######").PadLeft(14)
                    + cdf_compute.ToString("0.######").PadLeft(14) + "");
            }

            Console.WriteLine("");
            cout = "    X         Y         A         B         1-CDF         CCDF";
            cout +=
                "                                           (Lookup)      (Computed)";
            Console.WriteLine(cout);

            n_data = 0;

            for (;;)
            {
                CDF.beta_inc_values(ref n_data, ref a, ref b, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                y = 1.0 - x;

                CDF.beta_inc(a, b, x, y, ref cdf_compute, ref ccdf_compute, ref ierror);

                Console.WriteLine(x.ToString("0.######").PadLeft(10)
                    + y.ToString("0.######").PadLeft(10)
                    + a.ToString("0.######").PadLeft(10)
                    + b.ToString("0.######").PadLeft(10)
                    + ccdf_lookup.ToString("0.######").PadLeft(14)
                    + ccdf_compute.ToString("0.######").PadLeft(14) + "");

            }
        }

        static void cdfbet_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfbet_test tests CDFBET.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double bound = 0;
            double p = 0;
            double q = 0;
            int status = 0;
            int which = 0;
            double x = 0;
            double y = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfbet_test");
            Console.WriteLine("  CDFBET computes one missing parameter from the BETA CDF:");
            Console.WriteLine("");
            Console.WriteLine("   BETA_CDF ( (P,Q), (X,Y), A, B )");
            Console.WriteLine("");
            Console.WriteLine("      P           Q               X           Y"
                + "            A           B");
            Console.WriteLine("");

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 0.25;
                    y = 1.0 - x;
                    a = 2.0;
                    b = 3.0;
                }
                else if (which == 2)
                {
                    p = 0.261719;
                    q = 1.0 - p;
                    x = -1.0;
                    y = -1.0;
                    a = 2.0;
                    b = 3.0;
                }
                else if (which == 3)
                {
                    p = 0.261719;
                    q = 1.0 - p;
                    x = 0.25;
                    y = 1.0 - x;
                    a = -1.0;
                    b = 3.0;
                }
                else if (which == 4)
                {
                    p = 0.261719;
                    q = 1.0 - p;
                    x = 0.25;
                    y = 1.0 - x;
                    a = 2.0;
                    b = -1.0;
                }

                CDF.cdfbet(which, ref p, ref q, ref x, ref y, ref a, ref b, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFBET returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(10) + "  "
                     + q.ToString().PadLeft(10) + "  "
                     + x.ToString().PadLeft(10) + "  "
                     + y.ToString().PadLeft(10) + "  "
                     + a.ToString().PadLeft(10) + "  "
                     + b.ToString().PadLeft(10) + "");
            }
        }

        static void cdfbin_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfbin_test tests CDFBIN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double ompr = 0;
            double p = 0;
            double pr = 0;
            double q = 0;
            double s = 0;
            int status = 0;
            int which = 0;
            double xn = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfbin_test");
            Console.WriteLine("  CDFBIN computes one missing parameter from the");
            Console.WriteLine("  Binomial CDF:");
            Console.WriteLine("");
            Console.WriteLine("   BINOMIAL_CDF ( (P,Q), S, XN, (PR,OMPR) )");
            Console.WriteLine("");
            Console.WriteLine("      P           Q                S          "
                + "XN         PR         OMPR");
            Console.WriteLine("");

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    s = 5.0;
                    xn = 8.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 2)
                {
                    p = 0.067347;
                    q = 1.0 - p;
                    s = -1.0;
                    xn = 8.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 3)
                {
                    p = 0.067347;
                    q = 1.0 - p;
                    s = 5.0;
                    xn = -1.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 4)
                {
                    p = 0.067347;
                    q = 1.0 - p;
                    s = 5.0;
                    xn = 8.0;
                    pr = -1.0;
                    ompr = -1.0;
                }

                CDF.cdfbin(which, ref p, ref q, ref s, ref xn, ref pr, ref ompr, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFBIN returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(10) + "  "
                     + q.ToString().PadLeft(10) + "  "
                     + s.ToString().PadLeft(10) + "  "
                     + xn.ToString().PadLeft(10) + "  "
                     + pr.ToString().PadLeft(10) + "  "
                     + ompr.ToString().PadLeft(10) + "");
            }
        }

        static void cdfchi_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfchi_test tests CDFCHI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double df = 0;
            double p = 0;
            double q = 0;
            int status = 0;
            int which = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfchi_test");
            Console.WriteLine("  CDFCHI computes one missing parameter from the");
            Console.WriteLine("  Chi Square CDF:");
            Console.WriteLine("");
            Console.WriteLine("   CHI_CDF ( (P,Q), X, DF )");
            Console.WriteLine("");
            Console.WriteLine("      P           Q                X          DF");
            Console.WriteLine("");

            for (which = 1; which <= 3; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 5.0;
                    df = 8.0;
                }
                else if (which == 2)
                {
                    p = 0.242424;
                    q = 1.0 - p;
                    x = -1.0;
                    df = 8.0;
                }
                else if (which == 3)
                {
                    p = 0.242424;
                    q = 1.0 - p;
                    x = 5.0;
                    df = -1.0;
                }

                CDF.cdfchi(which, ref p, ref q, ref x, ref df, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFCHI returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(10) + "  "
                     + q.ToString().PadLeft(10) + "  "
                     + x.ToString().PadLeft(10) + "  "
                     + df.ToString().PadLeft(10) + "");
            }
        }

        static void cdfchn_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfchn_test tests CDFCHN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double df = 0;
            double p = 0;
            double pnonc = 0;
            double q = 0;
            int status = 0;
            int which = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfchn_test");
            Console.WriteLine("  CDFCHN computes one missing parameter from the");
            Console.WriteLine("  Chi Square CDF:");
            Console.WriteLine("");
            Console.WriteLine("   CHI_Noncentral_CDF ( (P,Q), X, DF, PNONC )");
            Console.WriteLine("");
            Console.WriteLine("     P         Q             X        DF     PNONC");
            Console.WriteLine("");

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 5.0;
                    df = 8.0;
                    pnonc = 0.5;
                }
                else if (which == 2)
                {
                    p = 0.211040;
                    q = 1.0 - p;
                    x = -1.0;
                    df = 8.0;
                    pnonc = 0.5;
                }
                else if (which == 3)
                {
                    p = 0.211040;
                    q = 1.0 - p;
                    x = 5.0;
                    df = -1.0;
                    pnonc = 0.5;
                }
                else if (which == 4)
                {
                    p = 0.211040;
                    q = 1.0 - p;
                    x = 5.0;
                    df = 8.0;
                    pnonc = -1.0;
                }

                CDF.cdfchn(which, ref p, ref q, ref x, ref df, ref pnonc, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFCHN returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(8) + "  "
                     + q.ToString().PadLeft(8) + "  "
                     + x.ToString().PadLeft(8) + "  "
                     + df.ToString().PadLeft(8) + "  "
                     + pnonc.ToString().PadLeft(8) + "");
            }
        }

        static void cdff_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdff_test tests CDFF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double dfd = 0;
            double dfn = 0;
            double f = 0;
            double p = 0;
            double q = 0;
            int status = 0;
            int which = 0;

            Console.WriteLine("");
            Console.WriteLine("cdff_test");
            Console.WriteLine("  CDFF computes one missing parameter from the");
            Console.WriteLine("  F CDF:");
            Console.WriteLine("");
            Console.WriteLine("   F_CDF ( (P,Q), F, DFN, DFD )");
            Console.WriteLine("");
            Console.WriteLine("     P         Q             F       DFN       DFD");
            Console.WriteLine("");

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = 3.0;
                }
                else if (which == 2)
                {
                    p = 0.893510;
                    q = 1.0 - p;
                    f = -1.0;
                    dfn = 8.0;
                    dfd = 3.0;
                }
                else if (which == 3)
                {
                    p = 0.893510;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = -1.0;
                    dfd = 3.0;
                }
                else if (which == 4)
                {
                    p = 0.893510;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = -1.0;
                }

                CDF.cdff(which, ref p, ref q, ref f, ref dfn, ref dfd, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFF returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(8) + "  "
                     + q.ToString().PadLeft(8) + "  "
                     + f.ToString().PadLeft(8) + "  "
                     + dfn.ToString().PadLeft(8) + "  "
                     + dfd.ToString().PadLeft(8) + "");
            }
        }

        public static void cdffnc_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdffnc_test tests CDFFNC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double dfd = 0;
            double dfn = 0;
            double f = 0;
            double p = 0;
            double pnonc = 0;
            double q = 0;
            int status = 0;
            int which = 0;

            Console.WriteLine("");
            Console.WriteLine("cdffnc_test");
            Console.WriteLine("  CDFFNC computes one missing parameter from the");
            Console.WriteLine("  noncentral F CDF:");
            Console.WriteLine("");
            Console.WriteLine("   F_noncentral_CDF ( (P,Q), F, DFN, DFD, PNONC )");
            Console.WriteLine("");
            Console.WriteLine("         P         Q         F       DFN       DFD     PNONC");
            Console.WriteLine("");

            for (which = 1; which <= 5; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = 3.0;
                    pnonc = 17.648016;
                }
                else if (which == 2)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = -1.0;
                    dfn = 8.0;
                    dfd = 3.0;
                    pnonc = 17.648016;
                }
                else if (which == 3)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = -1.0;
                    dfd = 3.0;
                    pnonc = 17.648016;
                }
                else if (which == 4)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = -1.0;
                    pnonc = 17.648016;
                }
                else if (which == 5)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = 3.0;
                    pnonc = -1.0;
                }

                CDF.cdffnc(which, ref p, ref q, ref f, ref dfn, ref dfd, ref pnonc, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFFNC returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(8) + "  "
                     + q.ToString().PadLeft(8) + "  "
                     + f.ToString().PadLeft(8) + "  "
                     + dfn.ToString().PadLeft(8) + "  "
                     + dfd.ToString().PadLeft(8) + "  "
                     + pnonc.ToString().PadLeft(8) + "");
            }
        }

        static void cdfgam_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfgam_test tests CDFGAM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double p = 0;
            double q = 0;
            double scale = 0;
            double shape = 0;
            int status = 0;
            int which = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfgam_test");
            Console.WriteLine("  CDFGAM computes one missing parameter from the");
            Console.WriteLine("  Gamma CDF:");
            Console.WriteLine("");
            Console.WriteLine("   Gamma_CDF ( (P,Q), X, SHAPE, SCALE )");
            Console.WriteLine("");
            Console.WriteLine("    P         Q              X     SHAPE     SCALE");
            Console.WriteLine("");

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 5.0;
                    shape = 8.0;
                    scale = 3.0;
                }
                else if (which == 2)
                {
                    p = 0.981998;
                    q = 1.0 - p;
                    x = -1.0;
                    shape = 8.0;
                    scale = 3.0;
                }
                else if (which == 3)
                {
                    p = 0.981998;
                    q = 1.0 - p;
                    x = 5.0;
                    shape = -1.0;
                    scale = 3.0;
                }
                else if (which == 4)
                {
                    p = 0.981998;
                    q = 1.0 - p;
                    x = 5.0;
                    shape = 8.0;
                    scale = -1.0;
                }

                CDF.cdfgam(which, ref p, ref q, ref x, ref shape, ref scale, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFGAM returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(8) + "  "
                     + q.ToString().PadLeft(9) + "  "
                     + x.ToString().PadLeft(8) + "  "
                     + shape.ToString().PadLeft(8) + "  "
                     + scale.ToString().PadLeft(8) + "");
            }
        }

        static void cdfnbn_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfnbn_test tests CDFNBN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double f = 0;
            double ompr = 0;
            double p = 0;
            double pr = 0;
            double q = 0;
            double s = 0;
            int status = 0;
            int which = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfnbn_test");
            Console.WriteLine("  CDFNBN computes one missing parameter from the");
            Console.WriteLine("  Negative_Binomial CDF:");
            Console.WriteLine("");
            Console.WriteLine("   Negative_BINOMIAL_CDF ( (P,Q), F, S, (PR,OMPR) )");
            Console.WriteLine("");
            Console.WriteLine("    P         Q               F         S       PR        OMPR");
            Console.WriteLine("");

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    f = 3.0;
                    s = 5.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 2)
                {
                    p = 0.988752;
                    q = 1.0 - p;
                    f = -1.0;
                    s = 5.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 3)
                {
                    p = 0.988752;
                    q = 1.0 - p;
                    f = 3.0;
                    s = -1.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 4)
                {
                    p = 0.988752;
                    q = 1.0 - p;
                    f = 3.0;
                    s = 5.0;
                    pr = -1.0;
                    ompr = -1.0;
                }

                CDF.cdfnbn(which, ref p, ref q, ref f, ref s, ref pr, ref ompr, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFNBN returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(8) + "  "
                     + q.ToString().PadLeft(9) + "  "
                     + f.ToString().PadLeft(8) + "  "
                     + s.ToString().PadLeft(8) + "  "
                     + pr.ToString().PadLeft(8) + "  "
                     + ompr.ToString().PadLeft(8) + "");
            }
        }

        static void cdfnor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfnor_test tests CDFNOR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double mean = 0;
            double p = 0;
            double q = 0;
            double sd = 0;
            int status = 0;
            int which = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfnor_test");
            Console.WriteLine("  CDFNOR computes one missing parameter from the");
            Console.WriteLine("  Normal CDF:");
            Console.WriteLine("");
            Console.WriteLine("   Normal_CDF ( (P,Q), X, MEAN, SD )");
            Console.WriteLine("");
            Console.WriteLine("    P         Q               X      MEAN       SD");
            Console.WriteLine("");

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 3.0;
                    mean = 5.0;
                    sd = 0.875;
                }
                else if (which == 2)
                {
                    p = 0.011135;
                    q = 1.0 - p;
                    x = -1.0;
                    mean = 5.0;
                    sd = 0.875;
                }
                else if (which == 3)
                {
                    p = 0.011135;
                    q = 1.0 - p;
                    x = 3.0;
                    mean = -1.0;
                    sd = 0.875;
                }
                else if (which == 4)
                {
                    p = 0.011135;
                    q = 1.0 - p;
                    x = 3.0;
                    mean = 5.0;
                    sd = -1.0;
                }

                CDF.cdfnor(which, ref p, ref q, ref x, ref mean, ref sd, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFNOR returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(9) + "  "
                     + q.ToString().PadLeft(8) + "  "
                     + x.ToString().PadLeft(8) + "  "
                     + mean.ToString().PadLeft(8) + "  "
                     + sd.ToString().PadLeft(8) + "");
            }

            return;
        }

        static void cdfpoi_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdfpoi_test tests CDFPOI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double p = 0;
            double q = 0;
            double s = 0;
            int status = 0;
            int which = 0;
            double xlam = 0;

            Console.WriteLine("");
            Console.WriteLine("cdfpoi_test");
            Console.WriteLine("  CDFPOI computes one missing parameter from the");
            Console.WriteLine("  Poisson CDF:");
            Console.WriteLine("");
            Console.WriteLine("   POISSON_CDF ( (P,Q), S, XLAM )");
            Console.WriteLine("");
            Console.WriteLine("     P         Q         S         XLAM");
            Console.WriteLine("");

            for (which = 1; which <= 3; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    s = 3.0;
                    xlam = 5.0;
                }
                else if (which == 2)
                {
                    p = 0.265026;
                    q = 1.0 - p;
                    s = -1.0;
                    xlam = 5.0;
                }
                else if (which == 3)
                {
                    p = 0.265026;
                    q = 1.0 - p;
                    s = 3.0;
                    xlam = -1.0;
                }

                CDF.cdfpoi(which, ref p, ref q, ref s, ref xlam, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFPOI returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(9) + "  "
                     + q.ToString().PadLeft(9) + "  "
                     + s.ToString().PadLeft(9) + "  "
                     + xlam.ToString().PadLeft(9) + "");
            }
        }

        static void cdft_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cdft_test tests CDFT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double bound = 0;
            double df = 0;
            double p = 0;
            double q = 0;
            int status = 0;
            double t = 0;
            int which = 0;

            Console.WriteLine("");
            Console.WriteLine("cdft_test");
            Console.WriteLine("  CDFT computes one missing parameter from the");
            Console.WriteLine("  T CDF:");
            Console.WriteLine("");
            Console.WriteLine("   T_CDF ( (P,Q), T, DF )");
            Console.WriteLine("");
            Console.WriteLine("    P         Q         T         DF");
            Console.WriteLine("");

            for (which = 1; which <= 3; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    t = 3.0;
                    df = 5.0;
                }
                else if (which == 2)
                {
                    p = 0.984950;
                    q = 1.0 - p;
                    t = -1.0;
                    df = 5.0;
                }
                else if (which == 3)
                {
                    p = 0.984950;
                    q = 1.0 - p;
                    t = 3.0;
                    df = -1.0;
                }

                CDF.cdft(which, ref p, ref q, ref t, ref df, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  CDFT returned STATUS = " + status + "");
                    continue;
                }

                Console.WriteLine("  "
                     + p.ToString().PadLeft(9) + "  "
                     + q.ToString().PadLeft(9) + "  "
                     + t.ToString().PadLeft(9) + "  "
                     + df.ToString().PadLeft(9) + "");
            }
        }

        static void cumbet_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumbet_test tests CUMBET, BETA_INC_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            double x = 0;
            double y = 0;

            Console.WriteLine("");
            Console.WriteLine("cumbet_test");
            Console.WriteLine("  CUMBET computes the Beta CDF");
            Console.WriteLine("  BETA_INC_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            string cout = "    X         Y         A         B         CDF           CDF";
            cout +=
                "                                           (Lookup)      (Computed)";
            Console.WriteLine(cout);

            int n_data = 0;

            for (;;)
            {
                CDF.beta_inc_values(ref n_data, ref a, ref b, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                y = 1.0 - x;

                CDF.cumbet(x, y, a, b, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine(" "
                     + x.ToString().PadLeft(9) + "  "
                     + y.ToString().PadLeft(9) + "  "
                     + a.ToString().PadLeft(9) + "  "
                     + b.ToString().PadLeft(9) + "  "
                     + cdf_lookup.ToString().PadLeft(9) + "  "
                     + cdf_compute.ToString().PadLeft(9) + "");
            }

            Console.WriteLine("");
            cout = ("    X         Y         A         B         1-CDF         CCDF");
            cout +=
                "                                           (Lookup)      (Computed)";
            Console.WriteLine(cout);

            n_data = 0;

            for (;;)
            {
                CDF.beta_inc_values(ref n_data, ref a, ref b, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                y = 1.0 - x;

                CDF.cumbet(x, y, a, b, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(9) + "  "
                     + y.ToString().PadLeft(9) + "  "
                     + a.ToString().PadLeft(9) + "  "
                     + b.ToString().PadLeft(9) + "  "
                     + ccdf_lookup.ToString().PadLeft(9) + "  "
                     + ccdf_compute.ToString().PadLeft(9) + "");
            }
        }

        static void cumbin_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumbin_test tests CUMBIN, BINOMIAL_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            double ompr = 0;
            int s = 0;
            double s_double = 0;
            double pr = 0;
            int x = 0;
            double x_double = 0;

            Console.WriteLine("");
            Console.WriteLine("cumbin_test");
            Console.WriteLine("  CUMBIN computes the Binomial CDF");
            Console.WriteLine("  BINOMIAL_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("   X   S    Pr       CDF           CDF");
            Console.WriteLine("                    (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.binomial_cdf_values(ref n_data, ref x, ref pr, ref s, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                ompr = 1.0 - pr;

                s_double = (double) s;
                x_double = (double) x;

                CDF.cumbin(s_double, x_double, pr, ompr, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + s.ToString().PadLeft(2) + "  "
                     + x.ToString().PadLeft(2) + "  "
                     + pr.ToString().PadLeft(8) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("   X   S    Pr       1-CDF         CCDF");
            Console.WriteLine("                    (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.binomial_cdf_values(ref n_data, ref x, ref pr, ref s, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                ompr = 1.0 - pr;

                s_double = (double) s;
                x_double = (double) x;

                CDF.cumbin(s_double, x_double, pr, ompr, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + s.ToString().PadLeft(2) + "  "
                     + x.ToString().PadLeft(2) + "  "
                     + pr.ToString().PadLeft(8) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumchi_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumchi_test tests CUMCHI, CHI_SQUARE_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int df = 0;
            double df_double = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumchi_test");
            Console.WriteLine("  CUMCHI computes the chi square CDF");
            Console.WriteLine("  CHI_SQUARE_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X       DF    CDF           CDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.chi_square_cdf_values(ref n_data, ref df, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                df_double = (double) df;

                CDF.cumchi(x, df_double, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + df.ToString().PadLeft(2) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    X       DF    1-CDF         CCDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.chi_square_cdf_values(ref n_data, ref df, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                df_double = (double) df;

                CDF.cumchi(x, df_double, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + df.ToString().PadLeft(2) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumchn_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumchn_test tests CUMCHN, CHI_NONCENTRAL_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int df = 0;
            double df_double = 0;
            double lambda = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumchn_test");
            Console.WriteLine("  CUMCHN computes the CDF for the noncentral chi-squared");
            Console.WriteLine("  distribution.");
            Console.WriteLine("  CHI_NONCENTRAL_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    DF    Lambda    X         CDF           CDF");
            Console.WriteLine("                             (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.chi_noncentral_cdf_values(ref n_data, ref x, ref lambda, ref df, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                df_double = (double) df;

                CDF.cumchn(x, df_double, lambda, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + df.ToString().PadLeft(6) + "  "
                     + lambda.ToString().PadLeft(8) + "  "
                     + x.ToString().PadLeft(8) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    DF    Lambda    X         1-CDF         CCDF");
            Console.WriteLine("                             (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.chi_noncentral_cdf_values(ref n_data, ref x, ref lambda, ref df, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                df_double = (double) df;

                CDF.cumchn(x, df_double, lambda, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + df.ToString().PadLeft(6) + "  "
                     + lambda.ToString().PadLeft(8) + "  "
                     + x.ToString().PadLeft(8) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
            
        }

        static void cumf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumf_test tests CUMF, F_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int dfd = 0;
            double dfd_double;
            int dfn = 0;
            double dfn_double = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumf_test");
            Console.WriteLine("  CUMF computes the F CDF");
            Console.WriteLine("  F_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X      DFN DFD    CDF           CDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.f_cdf_values(ref n_data, ref dfn, ref dfd, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                dfn_double = (double) dfn;
                dfd_double = (double) dfd;

                CDF.cumf(x, dfn_double, dfd_double, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + dfn.ToString().PadLeft(2) + "  "
                     + dfd.ToString().PadLeft(2) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");

            }

            Console.WriteLine("");
            Console.WriteLine("    X      DFN DFD    1-CDF         CCDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.f_cdf_values(ref n_data, ref dfn, ref dfd, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                dfn_double = (double) dfn;
                dfd_double = (double) dfd;

                CDF.cumf(x, dfn_double, dfd_double, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + dfn.ToString().PadLeft(2) + "  "
                     + dfd.ToString().PadLeft(2) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumfnc_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumfnc_test tests CUMFNC, F_NONCENTRAL_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int dfd = 0;
            double dfd_double = 0;
            int dfn = 0;
            double dfn_double = 0;
            double lambda = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumfnc_test");
            Console.WriteLine("  CUMFNC computes the noncentral F CDF");
            Console.WriteLine("  F_NONCENTRAL_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X      DFN DFD    LAMBDA    CDF           CDF");
            Console.WriteLine("                               (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.f_noncentral_cdf_values(ref n_data, ref dfn, ref dfd, ref lambda, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                dfn_double = (double) dfn;
                dfd_double = (double) dfd;

                CDF.cumfnc(x, dfn_double, dfd_double, lambda, ref cdf_compute,
                    ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + dfn.ToString().PadLeft(2) + "  "
                     + dfd.ToString().PadLeft(2) + "  "
                     + lambda.ToString().PadLeft(8) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    X      DFN DFD    LAMBDA    1-CDF         CCDF");
            Console.WriteLine("                               (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.f_noncentral_cdf_values(ref n_data, ref dfn, ref dfd, ref lambda, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                dfn_double = (double) dfn;
                dfd_double = (double) dfd;

                CDF.cumfnc(x, dfn_double, dfd_double, lambda, ref cdf_compute,
                    ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + dfn.ToString().PadLeft(2) + "  "
                     + dfd.ToString().PadLeft(2) + "  "
                     + lambda.ToString().PadLeft(8) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumgam_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumgam_test tests CUMGAM, GAMMA_INC_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumgam_test");
            Console.WriteLine("  CUMGAM computes the Gamma CDF");
            Console.WriteLine("  GAMMA_INC_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    A         X         CDF           CDF");
            Console.WriteLine("                        (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.gamma_inc_values(ref n_data, ref a, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                CDF.cumgam(x, a, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + a.ToString().PadLeft(8) + "  "
                     + x.ToString().PadLeft(8) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    A         X         CDF           CDF");
            Console.WriteLine("                        (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.gamma_inc_values(ref n_data, ref a, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                CDF.cumgam(x, a, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + a.ToString().PadLeft(8) + "  "
                     + x.ToString().PadLeft(8) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumnbn_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumnbn_test tests cumnbn() and negative_binomial_cdf_values()
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int f = 0;
            double f_double;
            double ompr;
            int s = 0;
            double s_double;
            double pr = 0;

            Console.WriteLine("");
            Console.WriteLine("cumnbn_test():");
            Console.WriteLine("  cumnbn() computes the Negative Binomial CDF");
            Console.WriteLine("  negative_binomial_cdf_values() looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("   F   S    Pr       CDF           CDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.negative_binomial_cdf_values(ref n_data, ref f, ref s, ref pr, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                ompr = 1.0 - pr;

                f_double = (double) f;
                s_double = (double) s;

                CDF.cumnbn(f_double, s_double, pr, ompr, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + f.ToString().PadLeft(2) + "  "
                     + s.ToString().PadLeft(2) + "  "
                     + pr.ToString().PadLeft(8) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("   F   S    Pr       1-CDF         CCDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.negative_binomial_cdf_values(ref n_data, ref f, ref s, ref pr, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                ompr = 1.0 - pr;

                f_double = (double) f;
                s_double = (double) s;

                CDF.cumnbn(f_double, s_double, pr, ompr, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + f.ToString().PadLeft(2) + "  "
                     + s.ToString().PadLeft(2) + "  "
                     + pr.ToString().PadLeft(8) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumnor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumnor_test tests CUMNOR, NORMAL_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumnor_test");
            Console.WriteLine("  CUMNOR computes the Normal CDF");
            Console.WriteLine("  NORMAL_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X         CDF           CDF");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                Normal.normal_cdf_values(ref n_data, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                CDF.cumnor(x, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    X         1-CDF         CCDF");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Normal.normal_cdf_values(ref n_data, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                CDF.cumnor(x, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumpoi_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumpoi_test tests CUMPOI, POISSON_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            double lambda = 0;
            int x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumpoi_test");
            Console.WriteLine("  CUMPOI computes the Poisson CDF");
            Console.WriteLine("  POISSON_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("     X    LAMBDA    CDF           CDF");
            Console.WriteLine("                   (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.poisson_cdf_values(ref n_data, ref lambda, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double x_double = (double) x;

                CDF.cumpoi(x_double, lambda, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(4) + "  "
                     + lambda.ToString().PadLeft(8) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("     X    LAMBDA    1-CDF         CCDF");
            Console.WriteLine("                   (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.poisson_cdf_values(ref n_data, ref lambda, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double x_double = (double) x;
                double ccdf_lookup = 1.0 - cdf_lookup;

                CDF.cumpoi(x_double, lambda, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(4) + "  "
                     + lambda.ToString().PadLeft(8) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }
        }

        static void cumt_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    cumt_test tests CUMT, STUDENT_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ccdf_compute = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int df = 0;
            double df_double;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("cumt_test");
            Console.WriteLine("  CUMT computes the Student T CDF");
            Console.WriteLine("  STUDENT_CDF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X       DF    CDF           CDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.student_cdf_values(ref n_data, ref df, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                df_double = (double) df;

                CDF.cumt(x, df_double, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + df.ToString().PadLeft(2) + "  "
                     + cdf_lookup.ToString().PadLeft(12) + "  "
                     + cdf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    X       DF    1-CDF         CCDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                CDF.student_cdf_values(ref n_data, ref df, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double ccdf_lookup = 1.0 - cdf_lookup;

                df_double = (double) df;

                CDF.cumt(x, df_double, ref cdf_compute, ref ccdf_compute);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + df.ToString().PadLeft(2) + "  "
                     + ccdf_lookup.ToString().PadLeft(12) + "  "
                     + ccdf_compute.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void beta_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    beta_test tests BETA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("beta_test");
            Console.WriteLine("  BETA evaluates the Beta function;");
            Console.WriteLine("  GAMMA_X evaluates the Gamma function.");

            double a = 2.2;
            double b = 3.7;
            double apb = a + b;

            double beta1 = CDF.beta(a, b);
            double beta2 = CDF.gamma_x(a) * CDF.gamma_x(b) / CDF.gamma_x(apb);

            Console.WriteLine("");
            Console.WriteLine("  Argument A =                   " + a + "");
            Console.WriteLine("  Argument B =                   " + b + "");
            Console.WriteLine("  Beta(A,B) =                    " + beta1 + "");
            Console.WriteLine("  (Expected value = 0.0454 )");
            Console.WriteLine("");
            Console.WriteLine("  Gamma(A)*Gamma(B)/Gamma(A+B) = " + beta2 + "");
        }

        static void error_f_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    error_f_test tests ERROR_F, ERROR_FC, ERF_VALUES..
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 November 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double erf_lookup = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("error_f_test");
            Console.WriteLine("  ERROR_F computes the error function ERF;");
            Console.WriteLine("  ERROR_FC the complementary error function ERFC.");
            Console.WriteLine("  ERF_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X         ERF           ERF");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                ErrorFunc.erf_values(ref n_data, ref x, ref erf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double erf_compute = CDF.error_f(x);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + erf_lookup.ToString().PadLeft(12) + "  "
                     + erf_compute.ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    X         ERFC          ERFC");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine("");

            int ind = 0;
            n_data = 0;

            for (;;)
            {
                ErrorFunc.erf_values(ref n_data, ref x, ref erf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double erfc_lookup = 1.0 - erf_lookup;
                double erfc_compute = CDF.error_fc(ind, x);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + erfc_lookup.ToString().PadLeft(12) + "  "
                     + erfc_compute.ToString().PadLeft(12) + "");
            }
        }

        static void gamma_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    gamma_test tests GAMMA_X.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double gamma_lookup = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("gamma_test");
            Console.WriteLine("  GAMMA_X computes the Gamma function;");
            Console.WriteLine("  tgamma() computes the Gamma function (math library)");
            Console.WriteLine("  GAMMA_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X         GAMMA         tgamma()     GAMMA_X");
            Console.WriteLine("              (Lookup)                   (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.gamma_values(ref n_data, ref x, ref gamma_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double tgamma_compute = Helpers.Gamma(x);

                double gamma_x_compute = CDF.gamma_x(x);

                Console.WriteLine("  " + x.ToString().PadLeft(8)
                    + "  " + gamma_lookup.ToString().PadLeft(12)
                    + "  " + tgamma_compute.ToString().PadLeft(12)
                    + "  " + gamma_x_compute.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void gamma_inc_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    gamma_inc_test tests GAMMA_INC, GAMMA_INC_INV.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int ierror = 0;
            double p = 0;
            double q = 0;
            int test_num = 10;
            double x2 = 0;

            double a = 3.0;
            int ind = 1;
            double x0 = 0;

            Console.WriteLine("");
            Console.WriteLine("gamma_inc_test");
            Console.WriteLine("  GAMMA_INC evaluates the incomplete Gamma ratio;");
            Console.WriteLine("  GAMMA_INC_INV inverts it.");
            Console.WriteLine("");
            Console.WriteLine("  Parameters:");
            Console.WriteLine("");
            Console.WriteLine("    A = " + a + "");
            Console.WriteLine("");
            Console.WriteLine("    X             P             Q             Inverse");
            Console.WriteLine("");

            for (int i = 0; i <= test_num; i++)
            {
                double x = (double) i / (double) test_num;

                CDF.gamma_inc(a, x, ref p, ref q, ind);

                CDF.gamma_inc_inv(a, ref x2, ref x0, ref p, ref q, ref ierror);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + p.ToString().PadLeft(12) + "  "
                     + q.ToString().PadLeft(12) + "  "
                     + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void psi_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    psi_test tests PSI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double psi_lookup = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("psi_test");
            Console.WriteLine("  PSI computes the Psi function;");
            Console.WriteLine("  PSI_VALUES looks up some values.");
            Console.WriteLine("");
            Console.WriteLine("    X         PSI           PSI");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                CDF.psi_values(ref n_data, ref x, ref psi_lookup);

                if (n_data == 0)
                {
                    break;
                }

                double psi_compute = CDF.psi(x);

                Console.WriteLine("  "
                     + x.ToString().PadLeft(8) + "  "
                     + psi_lookup.ToString().PadLeft(12) + "  "
                     + psi_compute.ToString().PadLeft(12) + "");
            }
        }
    }

    internal class TestValues
    {
    }
}