using System;
using Burkardt.AppliedStatistics;
using Burkardt.FullertonFnLib;
using Burkardt.Linpack;
using Burkardt.Values;

namespace FullertonTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FN_TEST.
        //
        //  Discussion:
        //
        //    FN_TEST tests the FN library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FN_TEST:");
        Console.WriteLine("  Test the FN library.");

        i4_mach_test();
        r8_acos_test();
        r8_acosh_test();
        r8_ai_test();
        r8_aid_test();
        r8_aint_test();
        r8_asin_test();
        r8_asinh_test();
        r8_atan_test();
        r8_atan2_test();
        r8_atanh_test();
        r8_besi0_test();
        r8_besi1_test();
        r8_besj0_test();
        r8_besj1_test();
        r8_besk0_test();
        r8_besk1_test();
        r8_besy0_test();
        r8_besy1_test();
        r8_beta_test();
        r8_betai_test();
        r8_bi_test();
        r8_bid_test();
        r8_binom_test();
        r8_cbrt_test();
        r8_chi_test();
        r8_chu_test();
        r8_ci_test();
        r8_cin_test();
        r8_cinh_test();
        r8_cos_test();
        r8_cos_deg_test();
        r8_cosh_test();
        r8_cot_test();
        r8_csevl_test();
        r8_dawson_test();
        r8_e1_test();
        r8_ei_test();
        r8_erf_test();
        r8_erfc_test();
        r8_exp_test();
        r8_fac_test();
        r8_gamic_test();
        r8_gamit_test();
        r8_gaml_test();
        r8_gamma_test();
        r8_gamr_test();
        r8_inits_test();
        r8_int_test();
        r8_lbeta_test();
        r8_lgams_test();
        r8_lgmc_test();
        r8_li_test();
        r8_lngam_test();
        r8_lnrel_test();
        r8_log_test();
        r8_log10_test();
        r8_mach_test();
        r8_pak_test();
        r8_poch_test();
        r8_psi_test();
        r8_rand_test();
        r8_randgs_test();
        r8_random_test();
        r8_ren_test();
        r8_shi_test();
        r8_si_test();
        r8_sin_test();
        r8_sin_deg_test();
        r8_sinh_test();
        r8_spence_test();
        r8_sqrt_test();
        r8_tan_test();
        r8_tanh_test();
        r8_upak_test();
        Console.WriteLine("");
        Console.WriteLine("FN_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void i4_mach_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MACH_TEST tests I4_MACH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("I4_MACH_TEST:");
        Console.WriteLine("  I4_MACH evaluates integer machine numbers.");
        Console.WriteLine("");
        Console.WriteLine("  I4_MACH(1) = the standard input unit.");
        Console.WriteLine("  I4_MACH(2) = the standard output unit.");
        Console.WriteLine("  I4_MACH(3) = the standard punch unit.");
        Console.WriteLine("  I4_MACH(4) = the standard error message unit.");
        Console.WriteLine("  I4_MACH(5) = the number of bits per integer storage unit.");
        Console.WriteLine("  I4_MACH(6) = the number of characters per integer storage unit.");
        Console.WriteLine("  I4_MACH(7) = A, the base.");
        Console.WriteLine("  I4_MACH(8) = S, the number of base A digits.");
        Console.WriteLine("  I4_MACH(9) = A^S-1, the largest integer.");
        Console.WriteLine("  I4_MACH(10) = B, the base.");
        Console.WriteLine("  I4_MACH(11) = T, the number of single precision base B digits.");
        Console.WriteLine("  I4_MACH(12) = EMIN, the smallest single precision exponent E.");
        Console.WriteLine("  I4_MACH(13) = EMAX, the largest single precision exponent E.");
        Console.WriteLine("  I4_MACH(14) = T, the number of double precision base B digits.");
        Console.WriteLine("  I4_MACH(15) = EMIN, the smallest double precision exponent E.");
        Console.WriteLine("  I4_MACH(16) = EMAX, the largest double precision exponent E.");
        Console.WriteLine("");
        Console.WriteLine("    I     I4_MACH(I)");
        Console.WriteLine("");

        for (i = 1; i <= 16; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + FullertonLib.i4_mach(i).ToString().PadLeft(12) + "");
        }

    }

    private static void r8_acos_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ACOS_TEST tests R8_ACOS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ACOS_TEST:");
        Console.WriteLine("  R8_ACOS evaluates the arccosine function.");
        Console.WriteLine("");
        Console.WriteLine("             X      ARCCOS(X)  R8_ACOS(X)         Diff");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            arc.arccos_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_acos(x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_acosh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ACOSH_TEST tests R8_ACOSH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ACOSH_TEST:");
        Console.WriteLine("  R8_ACOSH evaluates the hyperbolic arccosine function");
        Console.WriteLine("");
        Console.WriteLine("             X      ARCCOSH(X)  R8_ACOSH(X)        Diff");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            arc.arccosh_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_acosh(x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_ai_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_AI_TEST tests R8_AI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_AI_TEST:");
        Console.WriteLine("  Test R8_AI.");
        Console.WriteLine("");
        Console.WriteLine("             X   AIRY_AI(X)  R8_AI(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8AIData data = new();

        for (;;)
        {
            Airy.airy_ai_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_ai(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_aid_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_AID_TEST tests R8_AID.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_AID_TEST:");
        Console.WriteLine("  Test R8_AID.");
        Console.WriteLine("");
        Console.WriteLine("             X   AIRY_AID(X)  R8_AID(X)         Diff");

        n_data = 0;
        FullertonLib.r8AIDData data = new();

        for (;;)
        {
            Airy.airy_ai_prime_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_aid(ref data, x);

            Console.WriteLine("");
            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_aint_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_AINT_TEST tests R8_AINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_AINT_TEST:");
        Console.WriteLine("  R8_AINT rounds an R8 towards 0.");
        Console.WriteLine("");
        Console.WriteLine("             X      AINT(X)  R8_AINT(X)         Diff");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Intgr.int_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_aint(x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_asin_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ASIN_TEST tests R8_ASIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ASIN_TEST:");
        Console.WriteLine("  Test R8_ASIN.");
        Console.WriteLine("");
        Console.WriteLine("             X      ARCSIN(X)  R8_ASIN(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ASINData data = new();

        for (;;)
        {
            arc.arcsin_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_asin(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_asinh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ASINH_TEST tests R8_ASINH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ASINH_TEST:");
        Console.WriteLine("  Test R8_ASINH");
        Console.WriteLine("");
        Console.WriteLine("             X      ARCSINH(X)  R8_ASINH(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ASINHData data = new();

        for (;;)
        {
            arc.arcsinh_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_asinh(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_atan_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ATAN_TEST tests R8_ATAN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ATAN_TEST:");
        Console.WriteLine("  Test R8_ATAN.");
        Console.WriteLine("");
        Console.WriteLine("             X      ARCTAN(X)  R8_ATAN(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ATANData data = new();

        for (;;)
        {
            arc.arctan_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_atan(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_atan2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ATAN2_TEST tests R8_ATAN2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;
        double y = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ATAN2_TEST:");
        Console.WriteLine("  Test R8_ATAN2.");
        Console.WriteLine("");
        Console.WriteLine("             X             Y      ARCTAN2(Y,X)  R8_ATAN2(Y,X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ATAN2Data data = new();

        for (;;)
        {
            arc.arctan2_values(ref n_data, ref x, ref y, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_atan2(ref data, y, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + y.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_atanh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ATANH_TEST tests R8_ATANH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ATANH_TEST:");
        Console.WriteLine("  Test ARCTANH_VALUES, R4_ATANH, R8_ATANH");
        Console.WriteLine("");
        Console.WriteLine("             X      ARCTANH(X)  R8_ATANH(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ATANHData data = new();

        for (;;)
        {
            arc.arctanh_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_atanh(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besi0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESI0_TEST tests R8_BESI0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESI0_TEST:");
        Console.WriteLine("  Test R8_BESI0");
        Console.WriteLine("");
        Console.WriteLine("             X      BESI0(X)  R8_BESI0(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESI0Data data = new();

        for (;;)
        {
            Bessel.bessel_i0_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besi0(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besi1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESI1_TEST tests R8_BESI1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESI1_TEST:");
        Console.WriteLine("  Test R8_BESI1");
        Console.WriteLine("");
        Console.WriteLine("             X      BESI1(X)  R8_BESI1(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESI1Data data = new();

        for (;;)
        {
            Bessel.bessel_i1_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besi1(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besj0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESJ0_TEST tests R8_BESJ0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESJ0_TEST:");
        Console.WriteLine("  Test R8_BESJ0");
        Console.WriteLine("");
        Console.WriteLine("             X      BESJ0(X)  R8_BESJ0(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESJ0Data data = new();

        for (;;)
        {
            Bessel.bessel_j0_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besj0(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besj1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESJ1_TEST tests R8_BESJ1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESJ1_TEST:");
        Console.WriteLine("  Test R8_BESJ1");
        Console.WriteLine("");
        Console.WriteLine("             X      BESJ1(X)  R8_BESJ1(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESJ1Data data = new();

        for (;;)
        {
            Bessel.bessel_j1_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besj1(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besk0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESK0_TEST tests R8_BESK0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESK0_TEST:");
        Console.WriteLine("  Test R8_BESK0");
        Console.WriteLine("");
        Console.WriteLine("             X      BESK0(X)  R8_BESK0(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESK0Data data = new();

        for (;;)
        {
            Bessel.bessel_k0_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besk0(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besk1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESK1_TEST tests R8_BESK1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESK1_TEST:");
        Console.WriteLine("  Test R8_BESK1");
        Console.WriteLine("");
        Console.WriteLine("             X      BESK1(X)  R8_BESK1(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESK1Data data = new();

        for (;;)
        {
            Bessel.bessel_k1_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besk1(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besy0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESY0_TEST tests R8_BESY0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESY0_TEST:");
        Console.WriteLine("  Test R8_BESY0");
        Console.WriteLine("");
        Console.WriteLine("             X      BESY0(X)  R8_BESY0(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESY0Data data = new();

        for (;;)
        {
            Bessel.bessel_y0_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besy0(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_besy1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BESY1_TEST tests R8_BESY1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BESY1_TEST:");
        Console.WriteLine("  Test R8_BESY1");
        Console.WriteLine("");
        Console.WriteLine("             X      BESY1(X)  R8_BESY1(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BESY1Data data = new();

        for (;;)
        {
            Bessel.bessel_y1_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_besy1(ref globaldata, ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_beta_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BETA_TEST tests R8_BETA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx1 = 0;
        double fx2;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("R8_BETA_TEST:");
        Console.WriteLine("  Test R8_BETA.");
        Console.WriteLine("");
        Console.WriteLine("             X        BETA(A,B)  R8_BETA(A,B)       Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8GammaData gdata = new();
        FullertonLib.BesselData globaldata = new();
        FullertonLib.r8BetaData data = new();

        for (;;)
        {
            Beta.beta_values(ref n_data, ref a, ref b, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_beta(ref gdata, ref globaldata, ref data, a, b);


            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + b.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_betai_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BETAI_TEST tests R8_BETAI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BETAI_TEST:");
        Console.WriteLine("  Test R8_BETAI.");
        Console.WriteLine("");
        Console.WriteLine("             X        BETA(A,B,X)  R8_BETAI(A,B,X)       Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8Beta1Data data = new();
        FullertonLib.r8GammaData gammadata = new();

        for (;;)
        {
            Beta.beta_inc_values(ref n_data, ref a, ref b, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_betai( ref data, ref gammadata, x, a, b);

            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + b.ToString().PadLeft(14)
                                   + "  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_bi_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BI_TEST tests R8_BI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BI_TEST:");
        Console.WriteLine("  Test R8_BI.");
        Console.WriteLine("");
        Console.WriteLine("             X   AIRY_BI(X)  R8_BI(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8BiData data = new();

        for (;;)
        {
            Airy.airy_bi_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_bi(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_bid_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BID_TEST tests R8_BID.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_BID_TEST:");
        Console.WriteLine("  Test R8_BID.");
        Console.WriteLine("");
        Console.WriteLine("             X   AIRY_BID(X)  R8_BID(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8BidData data = new();

        for (;;)
        {
            Airy.airy_bi_prime_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_bid(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_binom_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_BINOM_TEST tests R8_BINOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a = 0;
        int b = 0;
        int fx1 = 0;
        double fx2;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("R8_BINOM_TEST:");
        Console.WriteLine("  Test R8_BINOM.");
        Console.WriteLine("");
        Console.WriteLine("             X        BINOM(A,B)  R8_BINOM(A,B)       Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8BinomData data = new();

        for (;;)
        {
            Binomial.binomial_values(ref n_data, ref a, ref b, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_binom(ref data, a,  b);

            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + b.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_cbrt_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CBRT_TEST tests R8_CBRT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_CBRT_TEST:");
        Console.WriteLine("  Test R8_CBRT");
        Console.WriteLine("");
        Console.WriteLine("             X      CBRT(X)  R8_CBRT(X)        Diff");
        Console.WriteLine("");

        n_data = 0;

        FullertonLib.r8CBRTData data = new();

        for (;;)
        {
            CubeRoot.cbrt_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_cbrt(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_chi_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHI_TEST tests R8_CHI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_CHI_TEST:");
        Console.WriteLine("  Test R8_CHI.");
        Console.WriteLine("");
        Console.WriteLine("             X      CHI(X)  R8_CHI(X)         Diff");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Chi.chi_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_chi(x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_chu_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHU_TEST tests R8_CHU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_CHU_TEST:");
        Console.WriteLine("  Test R8_CHU.");
        Console.WriteLine("");
        Console.WriteLine("             A               B               X     CHU(A,B,X)  R8_CHU(A,B,X)  Diff");
        Console.WriteLine("");

        n_data = 0;

        FullertonLib.r8CHUData data = new();
        FullertonLib.r8GammaData gammadata = new();

        for (;;)
        {
            Hypergeometric.hypergeometric_u_values(ref n_data, ref a, ref b, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_chu(ref data, ref gammadata, a, b, x);

            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + b.ToString().PadLeft(14)
                                   + "  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_ci_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CI_TEST tests R8_CI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_CI_TEST:");
        Console.WriteLine("  Test R8_CI.");
        Console.WriteLine("");
        Console.WriteLine("             X      CI(X)  R8_CI(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8CIData data = new();

        for (;;)
        {
            Cosine.ci_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_ci(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_cin_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CIN_TEST tests R8_CIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_CIN_TEST:");
        Console.WriteLine("  Test R8_CIN.");
        Console.WriteLine("");
        Console.WriteLine("             X      CIN(X)  R8_CIN(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8CINData data = new();

        for (;;)
        {
            Cosine.cin_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_cin(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_cinh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CINH_TEST tests R8_CINH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_CINH_TEST:");
        Console.WriteLine("  Test R8_CINH.");
        Console.WriteLine("");
        Console.WriteLine("             X      CINH(X)  R8_CINH(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8CINHData data = new();

        for (;;)
        {
            Cosine.cinh_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_cinh(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_cos_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COS_TEST tests R8_COS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_COS_TEST:");
        Console.WriteLine("  Test R8_COS.");
        Console.WriteLine("");
        Console.WriteLine("             X      COS(X)  R8_COS(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8CosData data = new();

        for (;;)
        {
            Cosine.cos_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_cos(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_cos_deg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COS_DEG_TEST tests R8_COS_DEG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_COS_DEG_TEST:");
        Console.WriteLine("  Test R8_COS_DEG.");
        Console.WriteLine("");
        Console.WriteLine("             X      COS_DEG(X)  R8_COS_DEG(X)         Diff");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Cosine.cos_degree_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_cos_deg(x);


            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_cosh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COSH_TEST tests R8_COSH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_COSH_TEST:");
        Console.WriteLine("  Test R8_COSH");
        Console.WriteLine("");
        Console.WriteLine("             X      COSH(X)  R8_COSH(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8CoshData data = new();

        for (;;)
        {
            Cosine.cosh_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_cosh(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_cot_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_COT_TEST tests R8_COT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_COT_TEST:");
        Console.WriteLine("  Test R8_COT.");
        Console.WriteLine("");
        Console.WriteLine("             X      COT(X)  R8_COT(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8CotData data = new();

        for (;;)
        {
            Cot.cot_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_cot(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_csevl_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CSEVL_TEST tests R8_CSEVL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double err;
        int i;
        int n;
        double[] expcs =  {
                2.532131755504016,
                1.130318207984970,
                0.271495339534077,
                0.044336849848664,
                0.005474240442094,
                0.000542926311914,
                0.000044977322954,
                0.000003198436462,
                0.000000199212481,
                0.000000011036772,
                0.000000000550590,
                0.000000000024980,
                0.000000000001039,
                0.000000000000040,
                0.000000000000001
            }
            ;
        double s;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_CSEVL_TEST:");
        Console.WriteLine("  R8_CSEVL evaluates a Chebyshev approximant");
        Console.WriteLine("  of N terms at a point X.");
        Console.WriteLine("");
        Console.WriteLine("  Here we use an approximant to the exponential function.");
        Console.WriteLine("  and average the absolute error at 21 points.");
        Console.WriteLine("");
        Console.WriteLine("   N    error");
        Console.WriteLine("");

        for (n = 1; n <= 12; n++)
        {
            err = 0.0;
            for (i = -10; i <= 10; i++)
            {
                x = i / 10.0;
                s = FullertonLib.r8_csevl(x, expcs, n);
                err += Math.Abs(s - Math.Exp(x));
            }

            err /= 21.0;
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + err.ToString().PadLeft(14) + "");
        }

    }

    private static void r8_dawson_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_DAWSON_TEST tests R8_DAWSON.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_DAWSON_TEST:");
        Console.WriteLine("  Test R8_DAWSON.");
        Console.WriteLine("");
        Console.WriteLine("             X      DAWSON(X)  R8_DAWSON(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8DawsonData data = new();

        for (;;)
        {
            Dawson.dawson_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_dawson(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_e1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_E1_TEST tests R8_E1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_E1_TEST:");
        Console.WriteLine("  Test R8_E1.");
        Console.WriteLine("");
        Console.WriteLine("             X      E1(X)  R8_E1(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8E1Data data = new();

        for (;;)
        {
            ExpIntegral.e1_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_e1(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_ei_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_EI_TEST tests R8_EI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_EI_TEST:");
        Console.WriteLine("  Test R8_EI.");
        Console.WriteLine("");
        Console.WriteLine("             X      EI(X)  R8_EI(X)         Diff");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            ExpIntegral.ei_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_ei(x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_erf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERF_TEST tests R8_ERF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ERF_TEST:");
        Console.WriteLine("  Test  R8_ERF.");
        Console.WriteLine("");
        Console.WriteLine("             X      ERF(X)  R8_ERF(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ErfData data = new();

        for (;;)
        {
            ErrorFunc.erf_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_erf(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_erfc_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ERFC_TEST tests R8_ERFC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_ERFC_TEST:");
        Console.WriteLine("  Test R8_ERFC.");
        Console.WriteLine("");
        Console.WriteLine("             X      ERFC(X)  R8_ERFC(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ErfCData data = new();

        for (;;)
        {
            ErrorFunc.erfc_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_erfc(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_exp_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_EXP_TEST tests R8_EXP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_EXP_TEST:");
        Console.WriteLine("  Test R8_EXP.");
        Console.WriteLine("");
        Console.WriteLine("             X      EXP(X)  R8_EXP(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ExpData data = new();

        for (;;)
        {
            Exponential.exp_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_exp(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_fac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FAC_TEST tests R8_FAC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("R8_FAC_TEST:");
        Console.WriteLine("  R8_FAC evaluates the factorial function.");
        Console.WriteLine("");
        Console.WriteLine("             N      FAC(N)  R8_FAC(N)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8FacData data = new();

        for (;;)
        {
            Factorial.r8_factorial_values(ref n_data, ref n, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_fac(ref data, n);

            Console.WriteLine("  " + n.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_gamic_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMIC_TEST tests R8_GAMIC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_GAMIC_TEST:");
        Console.WriteLine("  R8_GAMIC evaluates the incomplete gamma function.");
        Console.WriteLine("");
        Console.WriteLine("             X        GAMIC(A,X)  R8_GAMIC(A,X)       Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8GamicData data = new();
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Gamma.gamma_inc_values(ref n_data, ref a, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_gamic(ref data, ref gdata, a, x);

            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_gamit_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMIT_TEST tests R8_GAMIT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_GAMIT_TEST:");
        Console.WriteLine("  R8_GAMIT evaluates Tricomi's incomplete Gamma function.");
        Console.WriteLine("");
        Console.WriteLine("             X        GAMIT(A,X)  R8_GAMIT(A,X)       Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8GamitData data = new();
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Gamma.gamma_inc_tricomi_values(ref n_data, ref a, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_gamit(ref data, ref gdata, a, x);

            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_gaml_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAML_TEST tests R8_GAML.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double xmax = 0;
        double xmin = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_GAML_TEST:");
        Console.WriteLine("  R8_GAML returns bounds for the argument of the gamma function.");

        FullertonLib.r8_gaml(ref xmin, ref xmax);

        Console.WriteLine("");
        Console.WriteLine("  Lower limit XMIN = " + xmin + "");
        Console.WriteLine("  Upper limit XMAX = " + xmax + "");

    }

    private static void r8_gamma_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_TEST tests R8_GAMMA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_GAMMA_TEST:");
        Console.WriteLine("  R8_GAMMA evaluates the gamma function.");
        Console.WriteLine("");
        Console.WriteLine("             X      GAMMA(X)  R8_GAMMA(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Gamma.gamma_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_gamma(ref gdata, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_gamr_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMR_TEST tests R8_GAMR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        double gx = 0;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_GAMR_TEST:");
        Console.WriteLine("  R8_GAMR evaluates 1.0/Gamma(x).");
        Console.WriteLine("");
        Console.WriteLine("             X    1/GAMMA(X)  R8_GAMR(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8GamrData data = new();
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Gamma.gamma_values(ref n_data, ref x, ref gx);

            if (n_data == 0)
            {
                break;
            }

            fx1 = 1.0 / gx;
            fx2 = FullertonLib.r8_gamr(ref data, ref gdata, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_inits_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_INITS_TEST tests R8_INITS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        double[] sincs =  {
                -0.374991154955873175839919279977323464,
                -0.181603155237250201863830316158004754,
                +0.005804709274598633559427341722857921,
                -0.000086954311779340757113212316353178,
                +0.000000754370148088851481006839927030,
                -0.000000004267129665055961107126829906,
                +0.000000000016980422945488168181824792,
                -0.000000000000050120578889961870929524,
                +0.000000000000000114101026680010675628,
                -0.000000000000000000206437504424783134,
                +0.000000000000000000000303969595918706,
                -0.000000000000000000000000371357734157,
                +0.000000000000000000000000000382486123,
                -0.000000000000000000000000000000336623,
                +0.000000000000000000000000000000000256
            }
            ;
        double tol;

        Console.WriteLine("");
        Console.WriteLine("R8_INITS_TEST:");
        Console.WriteLine("  R8_INITS determines the Chebyshev interpolant degree");
        Console.WriteLine("  necessary to guarantee a desired accuracy level.");
        Console.WriteLine("");
        Console.WriteLine("  Here, we use a 15 term Chebyshev expansion for the");
        Console.WriteLine("  sine function.");
        Console.WriteLine("");
        Console.WriteLine("  Accuracy    Terms Needed");
        Console.WriteLine("");

        tol = 1.0;
        for (i = 1; i <= 18; i++)
        {
            n = FullertonLib.r8_inits(sincs, 15, tol);
            Console.WriteLine("  " + tol.ToString().PadLeft(14)
                                   + "  " + n.ToString().PadLeft(4) + "");
            tol /= 10.0;
        }

    }

    private static void r8_int_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_INT_TEST tests R8_INT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_INT_TEST:");
        Console.WriteLine("  R8_INT rounds an R8 to an integer value.");
        Console.WriteLine("");
        Console.WriteLine("             X      INT(X)  R8_INT(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8IntData data = new();

        for (;;)
        {
            Intgr.int_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_int(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_lbeta_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LBETA_TEST tests R8_LBETA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx1 = 0;
        double fx2;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("R8_LBETA_TEST:");
        Console.WriteLine("  R8_LBETA evaluates the logarithm of the Beta function.");
        Console.WriteLine("");
        Console.WriteLine("             X        LBETA(A,B)  R8_LBETA(A,B)       Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8LBetaData data = new();
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Beta.beta_log_values(ref n_data, ref a, ref b, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_lbeta(ref data, ref gdata, a, b);

            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + b.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_lgams_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LGAMS_TEST tests R8_LGAMS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2 = 0;
        int n_data;
        double slngam = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_LGAMS_TEST:");
        Console.WriteLine("  R8_LGAMS evaluates the sign of Gamma(x) and");
        Console.WriteLine("  the logarithm of the absolute value of Gamma(x).");
        Console.WriteLine("");
        Console.WriteLine("             X        LNGAM(X)  Sign(Gamma(X)) ALNGAM        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8LgamsData data = new();
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Gamma.gamma_log_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            FullertonLib.r8_lgams(ref data, ref gdata, x, ref fx2, ref slngam);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + slngam.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");

        }

    }

    private static void r8_lgmc_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LGMC_TEST tests R8_LGMC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        double gamma_log = 0;
        int n_data;
        const double r8_pi = 3.141592653589793;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_LGMC_TEST:");
        Console.WriteLine("  R8_LGMC evaluates the correction log gamma factor.");
        Console.WriteLine("  r8_lgmc(x) = log ( gamma ( x ) ) - log ( sqrt ( 2 * pi )");
        Console.WriteLine("  - ( x - 0.5 ) * log ( x ) + x");
        Console.WriteLine("");
        Console.WriteLine("             X        LGMC(X)  R8_LGMC(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8LgmcData data = new();

        for (;;)
        {
            Gamma.gamma_log_values(ref n_data, ref x, ref gamma_log);

            if (n_data == 0)
            {
                break;
            }

            switch (x)
            {
                //
                //  Function requires 10 <= x.
                //
                case >= 10.0:
                    fx1 = gamma_log - Math.Log(Math.Sqrt(2.0 * r8_pi)) - (x - 0.5) * Math.Log(x) + x;
                    fx2 = FullertonLib.r8_lgmc(ref data, x);
                    Console.WriteLine("  " + x.ToString().PadLeft(14)
                                           + "  " + fx1.ToString().PadLeft(14)
                                           + "  " + fx2.ToString().PadLeft(14)
                                           + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
                    break;
            }
        }

    }

    private static void r8_li_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LI_TEST tests R8_LI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_LI_TEST:");
        Console.WriteLine("  R8_LI evaluates the logarithmic integral");
        Console.WriteLine("");
        Console.WriteLine("             X      LI(X)  R8_LI(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8LiData data = new();

        for (;;)
        {
            Log.logarithmic_integral_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_li(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_lngam_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LNGAM_TEST tests R8_LNGAM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_LNGAM_TEST:");
        Console.WriteLine("  Test R8_LNGAM");
        Console.WriteLine("");
        Console.WriteLine("             X        LNGAM(X)  R8_LNGAM(X)        Diff");

        n_data = 0;
        FullertonLib.r8LngamData data = new();
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Gamma.gamma_log_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_lngam(ref data, ref gdata, x);

            Console.WriteLine("");
            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_lnrel_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LNREL_TEST tests R8_LNREL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL logcense.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_LNREL_TEST:");
        Console.WriteLine("  R8_LNREL evaluates ln(1+x)");
        Console.WriteLine("");
        Console.WriteLine("             X      LOG(1+X)  R8_LNREL(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8LnrelData data = new();

        for (;;)
        {
            Log.log_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x -= 1.0;

            fx2 = FullertonLib.r8_lnrel(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_log_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LOG_TEST tests R8_LOG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL logcense.
        //
        //  Modified:
        //
        //    25 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_LOG_TEST:");
        Console.WriteLine("  R8_LOG evaluates ln(x)");
        Console.WriteLine("");
        Console.WriteLine("             X      LOG(X)  R8_LOG(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8LogData data = new();

        for (;;)
        {
            Log.log_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_log(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_log10_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_LOG10_TEST tests R8_LOG10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL log10cense.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_LOG10_TEST:");
        Console.WriteLine("  Test R8_LOG10");
        Console.WriteLine("");
        Console.WriteLine("             X      LOG10(X)  R8_LOG10(X)        Diff");

        n_data = 0;

        for (;;)
        {
            Log.log10_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_log10(x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_mach_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_MACH_TEST tests R8_MACH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R8_MACH_TEST:");
        Console.WriteLine("  R8_MACH evaluates double precision machine numbers.");
        Console.WriteLine("");
        Console.WriteLine("  R8_MACH (1) = B^(EMIN-1), the smallest positive magnitude.");
        Console.WriteLine("  R8_MACH (2) = B^EMAX*(1 - B^(-T)), the largest magnitude.");
        Console.WriteLine("  R8_MACH (3) = B^(-T), the smallest relative spacing.");
        Console.WriteLine("  R8_MACH (4) = B^(1-T), the largest relative spacing.");
        Console.WriteLine("  R8_MACH (5) = LOG10(B)");
        Console.WriteLine("");
        Console.WriteLine("    I     R8_MACH(I)");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + FullertonLib.r8_mach(i).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_pak_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PAK_TEST tests R8_PAK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        int[] n_test =  {
                7, 8, 7, 7, 4, 0, -1, 0, 7, 2, 0
            }
            ;
        double x = 0;
        double y;
        double[] y_test =  {
                0.5,
                0.5,
                -0.5,
                0.75,
                0.9375,
                0.5,
                0.5,
                0.625,
                0.5048828125,
                0.7853981633974483,
                0.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("R8_PAK_TEST:");
        Console.WriteLine("  R8_PAK converts a mantissa and base 2 exponent to an R8.");
        Console.WriteLine("");
        Console.WriteLine("    Mantissa     Exponent         R8");
        Console.WriteLine("");
        FullertonLib.r8PakData data = new();

        for (i = 0; i < 11; i++)
        {
            y = y_test[i];
            n = n_test[i];

            x = FullertonLib.r8_pak(ref data, y, n);
            Console.WriteLine("  " + y.ToString().PadLeft(24)
                                   + "  " + n.ToString().PadLeft(8)
                                   + "  " + x.ToString().PadLeft(24) + "");
        }

    }

    private static void r8_poch_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_POCH_TEST tests R8_POCH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx1 = 0;
        double fx2;
        int n = 0;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_POCH_TEST:");
        Console.WriteLine("  R8_POCH evaluates the Pochhammer symbol..");
        Console.WriteLine("");
        Console.WriteLine("             X        POCH(A,X)  R8_POCH(A,X)       Diff");

        n_data = 0;
        FullertonLib.r8PochData data = new();
        FullertonLib.r8GammaData gdata = new();

        for (;;)
        {
            Factorial.r8_rise_values(ref n_data, ref a, ref n, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x = n;
            fx2 = FullertonLib.r8_poch(ref data, ref gdata, a, x);

            Console.WriteLine("  " + a.ToString().PadLeft(14)
                                   + "  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_psi_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PSI_TEST tests R8_PSI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL psicense.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_PSI_TEST:");
        Console.WriteLine("  Test R8_PSI");
        Console.WriteLine("");
        Console.WriteLine("             X      PSI(X)  R8_PSI(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8PsiData data = new();

        for (;;)
        {
            Psi.psi_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
                    
            }

            fx2 = FullertonLib.r8_psi(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_rand_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RAND_TEST tests R8_RAND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double average;
        int i;
        int[] i_value  = {
                1, 2, 3, 4, 10, 100, 1000
            }
            ;
        int ix0;
        int ix1;
        int k;
        double r;
        double[] r_value =  {
                0.0004127026,
                0.6750836372,
                0.1614754200,
                0.9086198807,
                0.5527787209,
                0.3600893021,
                0.2176990509
            }
            ;
        double s;
        double variance;

        Console.WriteLine("");
        Console.WriteLine("R8_RAND_TEST:");
        Console.WriteLine("  R8_RAND returns a random R8 value.");
        Console.WriteLine("");
        Console.WriteLine("               I       R8_RAND        Expected");
        Console.WriteLine("");
        //
        //  Start the sequence.
        //
        ix0 = 0;
        ix1 = 0;
        s = 0.0;

        k = 0;

        for (i = 1; i <= 1000; i++)
        {
            r = FullertonLib.r8_rand(s, ref ix0, ref ix1);

            if (i == i_value[k])
            {
                Console.WriteLine("  " + i.ToString().PadLeft(14)
                                       + "  " + r.ToString().PadLeft(14)
                                       + "  " + r_value[k].ToString().PadLeft(14) + "");
                k += 1;
            }
        }

        //
        //  Restart the sequence.
        //
        ix0 = 0;
        ix1 = 0;
        s = 0.0;

        average = 0.0;
        for (i = 1; i <= 1000000; i++)
        {
            r = FullertonLib.r8_rand(s, ref ix0, ref ix1);
            average += r;
        }

        average /= 1000000.0;
        Console.WriteLine("");
        Console.WriteLine("     Average =  "
                          + "  " + average.ToString().PadLeft(14)
                          + "  " + 0.5.ToString().PadLeft(14) + "");
        //
        //  Restart the sequence.
        //
        ix0 = 0;
        ix1 = 0;
        s = 0.0;

        variance = 0.0;
        for (i = 1; i <= 1000000; i++)
        {
            r = FullertonLib.r8_rand(s, ref ix0, ref ix1);
            variance += Math.Pow(r - average, 2);
        }

        variance /= 1000000.0;
        Console.WriteLine("     Variance = "
                          + "  " + variance.ToString().PadLeft(14)
                          + "  " + (1.0 / 12.0).ToString().PadLeft(14) + "");

    }

    private static void r8_randgs_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RANDGS_TEST tests R8_RANDGS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double m;
        double m2;
        double r;
        double sd;
        double sd2;
        int seed;

        m = 3.0;
        sd = 2.0;
        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("R8_RANDGS_TEST:");
        Console.WriteLine("  R8_RANDGS generates a random normal R8.");
        Console.WriteLine("  Mean =  " + m + "");
        Console.WriteLine("  Standard deviation = " + sd + "");
        Console.WriteLine("");
        Console.WriteLine("               I       R8_RANDGS");
        Console.WriteLine("");

        m2 = 0.0;
        sd2 = 0.0;

        for (i = 1; i <= 1000; i++)
        {
            r = FullertonLib.r8_randgs(m, sd, ref seed);
            m2 += r;
            sd2 += Math.Pow(m - r, 2);
            switch (i)
            {
                case <= 10:
                    Console.WriteLine("  " + i.ToString().PadLeft(14)
                                           + "  " + r.ToString().PadLeft(14) + "");
                    break;
            }
        }

        m2 /= 1000.0;
        sd2 = Math.Sqrt(sd2 / 1000.0);

        Console.WriteLine("");
        Console.WriteLine("  Sequence mean =  " + m2 + "");
        Console.WriteLine("  Sequence standard deviation = " + sd2 + "");

    }

    private static void r8_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RANDOM_TEST tests R8_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double average;
        int i;
        int[] i_value =  {
                1, 2, 3, 4, 10, 100, 1000
            }
            ;
        int ix0;
        int ix1;
        int k;
        int n = 32;
        double r;
        double[] t = new double[33];
        double variance;

        Console.WriteLine("");
        Console.WriteLine("R8_RANDOM_TEST:");
        Console.WriteLine("  R8_RANDOM returns a random R8 value.");
        Console.WriteLine("");
        Console.WriteLine("               I       R8_RANDOM");
        Console.WriteLine("");
        //
        //  Start the sequence.
        //
        ix0 = 0;
        ix1 = 0;
        FullertonLib.r8_random_init(n, ref t, ref ix0, ref ix1);

        k = 0;

        for (i = 1; i <= 1000; i++)
        {
            r = FullertonLib.r8_random(n, ref t, ref ix0, ref ix1);

            if (i == i_value[k])
            {
                Console.WriteLine("  " + i.ToString().PadLeft(14)
                                       + "  " + r.ToString().PadLeft(14) + "");
                k += 1;
            }
        }

        //
        //  Restart the sequence.
        //
        ix0 = 0;
        ix1 = 0;
        FullertonLib.r8_random_init(n, ref t, ref ix0, ref ix1);

        average = 0.0;
        for (i = 1; i <= 1000000; i++)
        {
            r = FullertonLib.r8_random(n, ref t, ref ix0, ref ix1);
            average += r;
        }

        average /= 1000000.0;
        Console.WriteLine("");
        Console.WriteLine("     Average =  "
                          + "  " + average.ToString().PadLeft(14)
                          + "  " + 0.5.ToString().PadLeft(14) + "");
        //
        //  Restart the sequence.
        //
        ix0 = 0;
        ix1 = 0;
        FullertonLib.r8_random_init(n, ref t, ref ix0, ref ix1);

        variance = 0.0;
        for (i = 1; i <= 1000000; i++)
        {
            r = FullertonLib.r8_random(n, ref t, ref ix0, ref ix1);
            variance += Math.Pow(r - average, 2);
        }

        variance /= 1000000.0;
        Console.WriteLine("     Variance = "
                          + "  " + variance.ToString().PadLeft(14)
                          + "  " + (1.0 / 12.0).ToString().PadLeft(14) + "");

    }

    private static void r8_ren_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_REN_TEST tests R8_REN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double average;
        int i;
        int[] i_value =  {
                1, 2, 3, 4, 10, 100, 1000
            }
            ;
        int k;
        double r;
        double[] r_value =  {
                0.470393,
                0.799066,
                0.883261,
                0.407667,
                0.955566,
                0.173576,
                0.0121733
            }
            ;
        int seed;
        double variance;

        Console.WriteLine("");
        Console.WriteLine("R8_REN_TEST:");
        Console.WriteLine("  R8_REN returns a random R8 value.");
        Console.WriteLine("");
        Console.WriteLine("               I       R8_REN         Expected");
        Console.WriteLine("");

        seed = 100001;
        k = 0;

        for (i = 1; i <= 1000; i++)
        {
            r = FullertonLib.r8_ren(ref seed);

            if (i == i_value[k])
            {
                Console.WriteLine("  " + i.ToString().PadLeft(14)
                                       + "  " + r.ToString().PadLeft(14)
                                       + "  " + r_value[k].ToString().PadLeft(14) + "");
                k += 1;
            }
        }

        seed = 123456789;
        average = 0.0;
        for (i = 1; i <= 1000000; i++)
        {
            r = FullertonLib.r8_ren(ref seed);
            average += r;
        }

        average /= 1000000.0;
        Console.WriteLine("");
        Console.WriteLine("     Average =  "
                          + "  " + average.ToString().PadLeft(14)
                          + "  " + 0.5.ToString().PadLeft(14) + "");

        seed = 123456789;
        variance = 0.0;
        for (i = 1; i <= 1000000; i++)
        {
            r = FullertonLib.r8_ren(ref seed);
            variance += Math.Pow(r - average, 2);
        }

        variance /= 1000000.0;
        Console.WriteLine("     Variance = "
                          + "  " + variance.ToString().PadLeft(14)
                          + "  " + (1.0 / 12.0).ToString().PadLeft(14) + "");

    }

    private static void r8_shi_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SHI_TEST tests R8_SHI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_SHI_TEST:");
        Console.WriteLine("  Test R8_SHI.");
        Console.WriteLine("");
        Console.WriteLine("             X      SHI(X)  R8_SHI(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8ShiData data = new();

        for (;;)
        {
            Shi.shi_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_shi(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_si_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SI_TEST tests R8_SI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_SI_TEST:");
        Console.WriteLine("  Test R8_SI.");
        Console.WriteLine("");
        Console.WriteLine("             X      SI(X)  R8_SI(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8SiData data = new();

        for (;;)
        {
            Sine.si_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_si(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_sin_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIN_TEST tests R8_SIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_SIN_TEST:");
        Console.WriteLine("  Test R8_SIN.");
        Console.WriteLine("");
        Console.WriteLine("             X      SIN(X)  R8_SIN(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8SinData data = new();

        for (;;)
        {
            Sine.sin_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_sin(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_sin_deg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIN_DEG_TEST tests R8_SIN_DEG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_SIN_DEG_TEST:");
        Console.WriteLine("  Test R8_SIN_DEG.");
        Console.WriteLine("");
        Console.WriteLine("             X      SIN_DEG(X)  R8_SIN_DEG(X)         Diff");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Sine.sin_degree_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_sin_deg(x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_sinh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SINH_TEST tests R8_SINH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_SINH_TEST:");
        Console.WriteLine("  Test R8_SINH");
        Console.WriteLine("");
        Console.WriteLine("             X      SINH(X)  R8_SINH(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8SinhData data = new();

        for (;;)
        {
            Sine.sinh_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_sinh(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_spence_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SPENCE_TEST tests R8_SPENCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_SPENCE_TEST:");
        Console.WriteLine("  Test R8_SPENCE");
        Console.WriteLine("");
        Console.WriteLine("             X      SPENCE(X)  R8_SPENCE(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8SpenceData data = new();

        for (;;)
        {
            Dilogarithm.dilogarithm_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_spence(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_sqrt_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SQRT_TEST tests R8_SQRT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_SQRT_TEST:");
        Console.WriteLine("  Test R8_SQRT");
        Console.WriteLine("");
        Console.WriteLine("             X      SQRT(X)  R8_SQRT(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8SqrtData data = new();

        for (;;)
        {
            Sqrt.sqrt_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_sqrt(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_tan_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TAN_TEST tests R8_TAN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_TAN_TEST:");
        Console.WriteLine("  Test R8_TAN.");
        Console.WriteLine("");
        Console.WriteLine("             X      TAN(X)  R8_TAN(X)         Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8TanData data = new();

        for (;;)
        {
            Tangent.tan_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_tan(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_tanh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_TANH_TEST tests R8_TANH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double fx2;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_TANH_TEST:");
        Console.WriteLine("  Test R8_TANH");
        Console.WriteLine("");
        Console.WriteLine("             X      TANH(X)  R8_TANH(X)        Diff");
        Console.WriteLine("");

        n_data = 0;
        FullertonLib.r8TanhData data = new();

        for (;;)
        {
            Tangent.tanh_values(ref n_data, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            fx2 = FullertonLib.r8_tanh(ref data, x);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + fx1.ToString().PadLeft(14)
                                   + "  " + fx2.ToString().PadLeft(14)
                                   + "  " + Math.Abs(fx1 - fx2).ToString().PadLeft(14) + "");
        }

    }

    private static void r8_upak_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_UPAK_TEST tests R8_UPAK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 0;
        double x = 0;
        double[] x_test =  {
                64.0,
                128.0,
                -64.0,
                96.0,
                15.0,
                0.5,
                0.25,
                0.625,
                64.625,
                3.141592653589793,
                0.0
            }
            ;
        double y = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_UPAK_TEST:");
        Console.WriteLine("  R8_UPAK converts an R8 to a mantissa and base 2 exponent.");
        Console.WriteLine("");
        Console.WriteLine("             X         Mantissa     Exponent");
        Console.WriteLine("");

        for (i = 0; i < 11; i++)
        {
            x = x_test[i];

            FullertonLib.r8_upak(x, ref y, ref n);

            Console.WriteLine("  " + x.ToString().PadLeft(14)
                                   + "  " + y.ToString().PadLeft(14)
                                   + "  " + n.ToString().PadLeft(8) + "");
        }

    }
}