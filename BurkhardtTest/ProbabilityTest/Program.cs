using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
partial class Program
{
static void Main(string[] args)
{
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PROB_TEST.
//
//  Discussion:
//
//    PROB_TEST tests the PROB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2018
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("PROB_TEST");
        Console.WriteLine("  Test the PROB library.");

        angle_cdf_test();
        angle_mean_test();
        angle_pdf_test();

        anglit_cdf_test();
        anglit_sample_test();

        arcsin_cdf_test();
        arcsin_sample_test();

        benford_cdf_test();
        benford_pdf_test();

        bernoulli_cdf_test();
        bernoulli_sample_test();

        bessel_i0_test();
        bessel_i1_test();

        beta_binomial_cdf_test();
        beta_binomial_sample_test();

        beta_cdf_test();
        beta_inc_test();
        beta_sample_test();

        binomial_cdf_test();
        binomial_sample_test();

        birthday_cdf_test();
        birthday_sample_test();

        bradford_cdf_test();
        bradford_sample_test();

        buffon_box_pdf_test();
        buffon_box_sample_test();

        buffon_pdf_test();
        buffon_sample_test();

        burr_cdf_test();
        burr_sample_test();

        cardioid_cdf_test();
        cardioid_sample_test();

        cauchy_cdf_test();
        cauchy_sample_test();

        chebyshev1_cdf_test();
        chebyshev1_sample_test();

        chi_cdf_test();
        chi_sample_test();

        chi_square_cdf_test();
        chi_square_sample_test();

        chi_square_noncentral_sample_test();

        circular_normal_01_sample_test();
        circular_normal_sample_test();

        cosine_cdf_test();
        cosine_sample_test();

        coupon_complete_pdf_test();
        coupon_sample_test();

        deranged_cdf_test();
        deranged_sample_test();

        dipole_cdf_test();
        dipole_sample_test();

        dirichlet_sample_test();
        dirichlet_pdf_test();

        dirichlet_mix_sample_test();
        dirichlet_mix_pdf_test();

        discrete_cdf_test();
        discrete_sample_test();

        disk_sample_test();

        empirical_discrete_cdf_test();
        empirical_discrete_sample_test();

        english_letter_cdf_test();

        english_sentence_length_cdf_test();
        english_sentence_length_sample_test();

        english_word_length_cdf_test();
        english_word_length_sample_test();

        erlang_cdf_test();
        erlang_sample_test();

        exponential_cdf_test();
        exponential_sample_test();

        exponential_01_cdf_test();
        exponential_01_sample_test();

        extreme_values_cdf_test();
        extreme_values_sample_test();

        f_cdf_test();
        f_sample_test();

        fermi_dirac_sample_test();

        fisher_pdf_test();

        fisk_cdf_test();
        fisk_sample_test();

        folded_normal_cdf_test();
        folded_normal_sample_test();

        frechet_cdf_test();
        frechet_sample_test();

        gamma_cdf_test();
        gamma_sample_test();

        genlogistic_cdf_test();
        genlogistic_sample_test();

        geometric_cdf_test();
        geometric_sample_test();

        gompertz_cdf_test();
        gompertz_sample_test();

        gumbel_cdf_test();
        gumbel_sample_test();

        half_normal_cdf_test();
        half_normal_sample_test();

        hypergeometric_cdf_test();
        hypergeometric_sample_test();

        i4_choose_test();
        i4_choose_log_test();
        i4_is_power_of_10_test();
        i4_uniform_ab_test();
        i4vec_uniform_ab_new_test();
        i4vec_unique_count_test();

        inverse_gaussian_cdf_test();
        inverse_gaussian_sample_test();

        laplace_cdf_test();
        laplace_sample_test();

        levy_cdf_test();

        logistic_cdf_test();
        logistic_sample_test();

        log_normal_cdf_test();
        log_normal_sample_test();

        log_series_cdf_test();
        log_series_sample_test();

        log_uniform_cdf_test();
        log_uniform_sample_test();

        lorentz_cdf_test();
        lorentz_sample_test();

        maxwell_cdf_test();
        maxwell_sample_test();

        multinomial_coef_test();
        multinomial_sample_test();
        multinomial_pdf_test();

        multinoulli_pdf_test();

        nakagami_cdf_test();
        nakagami_sample_test();

        negative_binomial_cdf_test();
        negative_binomial_sample_test();

        normal_01_cdf_test();
        normal_01_samples_test();

        normal_cdf_test();
        normal_samples_test();

        normal_truncated_ab_cdf_test();
        normal_truncated_ab_sample_test();

        normal_truncated_a_cdf_test();
        normal_truncated_a_sample_test();

        normal_truncated_b_cdf_test();
        normal_truncated_b_sample_test();

        pareto_cdf_test();
        pareto_sample_test();

        pearson_05_pdf_test();

        planck_pdf_test();
        planck_sample_test();

        poisson_cdf_test();
        poisson_sample_test();

        power_cdf_test();
        power_sample_test();

        quasigeometric_cdf_test();
        quasigeometric_sample_test();

        r8_beta_test();
        r8_ceiling_test();
        r8_error_f_test();
        r8_factorial_test();
        r8_gamma_inc_test();
        r8_gamma_log_int_test();
        r8_uniform_01_test();
        r8_zeta_test();

        rayleigh_cdf_test();
        rayleigh_sample_test();

        reciprocal_cdf_test();
        reciprocal_sample_test();

        runs_pdf_test();
        runs_sample_test();

        sech_cdf_test();
        sech_sample_test();

        semicircular_cdf_test();
        semicircular_sample_test();

        student_cdf_test();
        student_sample_test();

        student_noncentral_cdf_test();

        tfn_test();

        triangle_cdf_test();
        triangle_sample_test();

        triangular_cdf_test();
        triangular_sample_test();

        uniform_01_cdf_test();
        uniform_01_order_sample_test();
        uniform_01_sample_test();

        uniform_cdf_test();
        uniform_sample_test();

        uniform_discrete_cdf_test();
        uniform_discrete_sample_test();

        uniform_nsphere_sample_test();

        von_mises_cdf_test();
        von_mises_sample_test();

        weibull_cdf_test();
        weibull_sample_test();

        weibull_discrete_cdf_test();
        weibull_discrete_sample_test();

        zipf_cdf_test();
        zipf_sample_test();

        Console.WriteLine("");
        Console.WriteLine("PROB_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }










static void uniform_01_order_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_ORDER_SAMPLE_TEST tests UNIFORM_01_ORDER_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int n;
int seed = 123456789;
double[] x;

Console.WriteLine("");
Console.WriteLine("UNIFORM_01_ORDER_SAMPLE_TEST");
Console.WriteLine("  For the Uniform 01 Order PDF:");
Console.WriteLine("  UNIFORM_ORDER_SAMPLE samples.");

n = 10;
x = uniform_01_order_sample(n, seed);

typeMethods.r8vec_print(n, x, "  Ordered sample:");



return;
}

static void uniform_nsphere_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_NSPHERE_SAMPLE_TEST tests UNIFORM_NSPHERE_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int i;
int j;
int n;
int seed = 123456789;
double[] x;

n = 3;

Console.WriteLine("");
Console.WriteLine("UNIFORM_NSPHERE_SAMPLE_TEST");
Console.WriteLine("  UNIFORM_NSPHERE_SAMPLE samples the Uniform Nsphere distribution.");

Console.WriteLine("");
Console.WriteLine("  Dimension N of sphere = " + n + "");
Console.WriteLine("");
Console.WriteLine("  Points on the sphere:");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = uniform_nsphere_sample(n, seed);
Console.WriteLine("  " + setw(6) + i + "  ";
for (j = 0; j < n; j++)
{
Console.WriteLine(x.ToString().PadLeft(12)[j] + "  ";
}

Console.WriteLine("");

}

return;
}

static void uniform_01_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_CDF_TEST tests UNIFORM_01_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("UNIFORM_01_CDF_TEST");

Console.WriteLine("  UNIFORM_01_CDF evaluates the Uniform 01 CDF;");
Console.WriteLine("  UNIFORM_01_CDF_INV inverts the Uniform 01 CDF.");
Console.WriteLine("  UNIFORM_01_PDF evaluates the Uniform 01 PDF;");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = uniform_01_sample(seed);
pdf = uniform_01_pdf(x);
cdf = uniform_01_cdf(x);
x2 = uniform_01_cdf_inv(cdf);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void uniform_01_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_SAMPLE_TEST tests UNIFORM_01_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

int i;
double mean;
int seed = 123456789;
double variance;
double[] x = new double [SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("UNIFORM_01_SAMPLE_TEST");
Console.WriteLine("  UNIFORM_01_MEAN computes the Uniform 01 mean;");
Console.WriteLine("  UNIFORM_01_SAMPLE samples the Uniform 01 distribution;");
Console.WriteLine("  UNIFORM_01_VARIANCE computes the Uniform 01 variance.");

mean = uniform_01_mean();
variance = uniform_01_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = uniform_01_sample(seed);
}

mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void uniform_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_CDF_TEST tests UNIFORM_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("UNIFORM_CDF_TEST");
Console.WriteLine("  UNIFORM_CDF evaluates the Uniform CDF;");
Console.WriteLine("  UNIFORM_CDF_INV inverts the Uniform CDF.");
Console.WriteLine("  UNIFORM_PDF evaluates the Uniform PDF;");

a = 1.0;
b = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!uniform_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("UNIFORM_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = uniform_sample(a, b, seed);
pdf = uniform_pdf(x, a, b);
cdf = uniform_cdf(x, a, b);
x2 = uniform_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void uniform_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_SAMPLE_TEST tests UNIFORM_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double[] x = new double [SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("UNIFORM_SAMPLE_TEST");
Console.WriteLine("  UNIFORM_MEAN computes the Uniform mean;");
Console.WriteLine("  UNIFORM_SAMPLE samples the Uniform distribution;");
Console.WriteLine("  UNIFORM_VARIANCE computes the Uniform variance;");

a = 1.0;
b = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!uniform_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("UNIFORM_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = uniform_mean(a, b);
variance = uniform_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = uniform_sample(a, b, seed);
}

mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void uniform_discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_CDF_TEST tests UNIFORM_DISCRETE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
int a;
int b;
double cdf;
int i;
double pdf;
int seed = 123456789;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("UNIFORM_DISCRETE_CDF_TEST");
Console.WriteLine("  UNIFORM_DISCRETE_CDF evaluates the Uniform Discrete CDF;");
Console.WriteLine("  UNIFORM_DISCRETE_CDF_INV inverts the Uniform Discrete CDF.");
Console.WriteLine("  UNIFORM_DISCRETE_PDF evaluates the Uniform Discrete PDF;");

a = 1;
b = 6;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!uniform_discrete_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("UNIFORM_DISCRETE_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = uniform_discrete_sample(a, b, seed);
pdf = uniform_discrete_pdf(x, a, b);
cdf = uniform_discrete_cdf(x, a, b);
x2 = uniform_discrete_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void uniform_discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_SAMPLE_TEST tests UNIFORM_DISCRETE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

int a;
int b;
int i;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("UNIFORM_DISCRETE_SAMPLE_TEST");
Console.WriteLine("  UNIFORM_DISCRETE_MEAN computes the Uniform Discrete mean;");
Console.WriteLine("  UNIFORM_DISCRETE_SAMPLE samples the Uniform Discrete distribution;");
Console.WriteLine("  UNIFORM_DISCRETE_VARIANCE computes the Uniform Discrete variance;");

a = 1;
b = 6;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!uniform_discrete_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("UNIFORM_DISCRETE_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = uniform_discrete_mean(a, b);
variance = uniform_discrete_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = uniform_discrete_sample(a, b, seed);
}

mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void von_mises_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CDF_TEST tests VON_MISES_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("VON_MISES_CDF_TEST");
Console.WriteLine("  VON_MISES_CDF evaluates the Von Mises CDF;");
Console.WriteLine("  VON_MISES_CDF_INV inverts the Von Mises CDF.");
Console.WriteLine("  VON_MISES_PDF evaluates the Von Mises PDF;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!von_mises_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("VON_MISES_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = von_mises_sample(a, b, seed);
pdf = von_mises_pdf(x, a, b);
cdf = von_mises_cdf(x, a, b);
x2 = von_mises_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void von_mises_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_SAMPLE_TEST tests VON_MISES_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double[] x = new double [SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("VON_MISES_SAMPLE_TEST");
Console.WriteLine("  VON_MISES_MEAN computes the Von Mises mean;");
Console.WriteLine("  VON_MISES_SAMPLE samples the Von Mises distribution;");
Console.WriteLine("  VON_MISES_CIRCULAR_VARIANCE computes the Von Mises circular variance;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!von_mises_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("VON_MISES_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = von_mises_mean(a, b);
variance = von_mises_circular_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =              " + mean + "");
Console.WriteLine("  PDF circular variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = von_mises_sample(a, b, seed);
}

mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
variance = typeMethods.r8vec_circular_variance(SAMPLE_NUM, x);
xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =              " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =              " + mean + "");
Console.WriteLine("  Sample circular variance = " + variance + "");
Console.WriteLine("  Sample maximum =           " + xmax + "");
Console.WriteLine("  Sample minimum =           " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void weibull_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF_TEST tests WEIBULL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
double c;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("WEIBULL_CDF_TEST");
Console.WriteLine("  WEIBULL_CDF evaluates the Weibull CDF;");
Console.WriteLine("  WEIBULL_CDF_INV inverts the Weibull CDF.");
Console.WriteLine("  WEIBULL_PDF evaluates the Weibull PDF;");

a = 2.0;
b = 3.0;
c = 4.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!weibull_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("WEIBULL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = Weibull.weibull_sample(a, b, c, seed);
pdf = Weibull.weibull_pdf(x, a, b, c);
cdf = Weibull.weibull_cdf(x, a, b, c);
x2 = Weibull.weibull_cdf_inv(cdf, a, b, c);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void weibull_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_SAMPLE_TEST tests WEIBULL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double a;
double b;
double c;
int i;
double mean;
int seed = 123456789;
double variance;
double[] x = new double [SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("WEIBULL_SAMPLE_TEST");
Console.WriteLine("  WEIBULL_MEAN computes the Weibull mean;");
Console.WriteLine("  WEIBULL_SAMPLE samples the Weibull distribution;");
Console.WriteLine("  WEIBULL_VARIANCE computes the Weibull variance.");

a = 2.0;
b = 3.0;
c = 4.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!weibull_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("WEIBULL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = Weibull.weibull_mean(a, b, c);
variance = Weibull.weibull_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = Weibull.weibull_sample(a, b, c, seed);
}

mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;

# undef SAMPLE_NUM
}

static void weibull_discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_CDF_TEST tests WEIBULL_DISCRETE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
double cdf;
int i;
double pdf;
int seed = 123456789;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("WEIBULL_DISCRETE_CDF_TEST");
Console.WriteLine("  WEIBULL_DISCRETE_CDF evaluates the Weibull Discrete CDF;");
Console.WriteLine("  WEIBULL_DISCRETE_CDF_INV inverts the Weibull Discrete CDF.");
Console.WriteLine("  WEIBULL_DISCRETE_PDF evaluates the Weibull Discrete PDF;");

a = 0.5;
b = 1.5;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!weibull_discrete_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("WEIBULL_DISCRETE_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = Weibull.weibull_discrete_sample(a, b, seed);
pdf = Weibull.weibull_discrete_pdf(x, a, b);
cdf = Weibull.weibull_discrete_cdf(x, a, b);
x2 = Weibull.weibull_discrete_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void weibull_discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_SAMPLE_TEST tests WEIBULL_DISCRETE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double[] x = new double [SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("WEIBULL_DISCRETE_SAMPLE_TEST");
Console.WriteLine("  WEIBULL_DISCRETE_SAMPLE samples the Weibull Discrete distribution;");

a = 0.5;
b = 1.5;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!weibull_discrete_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("WEIBULL_DISCRETE_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = Weibull.weibull_discrete_sample(a, b, seed);
}

mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void zipf_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_CDF_TEST tests ZIPF_CDF, ZIPF_CDF_INV, ZIPF_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double cdf;
double pdf;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("ZIPF_CDF_TEST");
Console.WriteLine("  ZIPF_CDF evaluates the Zipf CDF;");
Console.WriteLine("  ZIPF_CDF_INV inverts the Zipf CDF;");
Console.WriteLine("  ZIPF_PDF evaluates the Zipf PDF;");

a = 2.0E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!zipf_check(a))
{
Console.WriteLine("");
Console.WriteLine("ZIPF_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF    CDF_INV()");
Console.WriteLine("");

for (x = 1; x <= 20; x++)
{
pdf = zipf_pdf(x, a);
cdf = zipf_cdf(x, a);
x2 = zipf_cdf_inv(a, cdf);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void zipf_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_SAMPLE_TEST tests ZIPF_MEAN, ZIPF_SAMPLE, ZIPF_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double a;
int j;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("ZIPF_SAMPLE_TEST");
Console.WriteLine("  ZIPF_MEAN computes the Zipf mean;");
Console.WriteLine("  ZIPF_SAMPLE samples the Zipf distribution;");
Console.WriteLine("  ZIPF_VARIANCE computes the Zipf variance.");

a = 4.0E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!zipf_check(a))
{
Console.WriteLine("");
Console.WriteLine("ZIPF_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = zipf_mean(a);
variance = zipf_variance(a);

Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (j = 0; j < SAMPLE_NUM; j++)
{
x[j] = zipf_sample(a, seed);
}

mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
}
}
}
}