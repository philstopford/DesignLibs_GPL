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









static void normal_01_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_TEST tests NORMAL_01_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
Console.WriteLine("NORMAL_01_CDF_TEST");
Console.WriteLine("  NORMAL_01_CDF evaluates the Normal 01 CDF;");
Console.WriteLine("  NORMAL_01_CDF_INV inverts the Normal 01 CDF.");
Console.WriteLine("  NORMAL_01_PDF evaluates the Normal 01 PDF;");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = normal_01_sample(seed);
pdf = normal_01_pdf(x);
cdf = normal_01_cdf(x);
x2 = normal_01_cdf_inv(cdf);

Console.WriteLine("  "
+ setw(24) + setprecision(16) + x + "  "
+ setw(12) + setprecision(6) + pdf + "  "
+ setw(12) + setprecision(6) + cdf + "  "
+ setw(24) + setprecision(16) + x2 + "");
}

Console.WriteLine(setprecision(6);

return;
}

static void normal_01_samples_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_SAMPLE_TEST tests NORMAL_01_SAMPLES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double mean;
int seed = 123456789;
double variance;
double[] x;
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("NORMAL_01_SAMPLES_TEST");
Console.WriteLine("  NORMAL_01_MEAN computes the Normal 01 mean;");
Console.WriteLine("  NORMAL_01_SAMPLES samples the Normal 01 distribution;");
Console.WriteLine("  NORMAL_01_VARIANCE computes the Normal 01 variance;");

mean = normal_01_mean();
variance = normal_01_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

x = normal_01_samples(SAMPLE_NUM, seed);

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

static void normal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF_TEST tests NORMAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
double cdf;
int i;
double mu;
double pdf;
int seed = 123456789;
double sigma;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("NORMAL_CDF_TEST");
Console.WriteLine("  NORMAL_CDF evaluates the Normal CDF;");
Console.WriteLine("  NORMAL_CDF_INV inverts the Normal CDF.");
Console.WriteLine("  NORMAL_PDF evaluates the Normal PDF;");

mu = 100.0;
sigma = 15.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter MU =    " + mu + "");
Console.WriteLine("  PDF parameter SIGMA = " + sigma + "");

if (!normal_check(mu, sigma))
{
Console.WriteLine("");
Console.WriteLine("NORMAL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = normal_sample(mu, sigma, seed);
pdf = normal_pdf(x, mu, sigma);
cdf = normal_cdf(x, mu, sigma);
x2 = normal_cdf_inv(cdf, mu, sigma);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void normal_samples_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_SAMPLES_TEST tests NORMAL_SAMPLES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double mean;
double mu;
int seed = 123456789;
double sigma;
double variance;
double[] x;
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("NORMAL_SAMPLES_TEST");
Console.WriteLine("  NORMAL_MEAN computes the Normal mean;");
Console.WriteLine("  NORMAL_SAMPLES samples the Normal distribution;");
Console.WriteLine("  NORMAL_VARIANCE computes the Normal variance;");

mu = 100.0;
sigma = 15.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter MU=     " + mu + "");
Console.WriteLine("  PDF parameter SIGMA = " + sigma + "");

if (!normal_check(mu, sigma))
{
Console.WriteLine("");
Console.WriteLine("NORMAL_SAMPLES_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = normal_mean(mu, sigma);
variance = normal_variance(mu, sigma);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

x = normal_samples(SAMPLE_NUM, mu, sigma, seed);

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

static void normal_truncated_ab_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_TRUNCATED_AB_CDF_TEST tests NORMAL_TRUNCATED_AB_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
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
double mu;
double pdf;
double s;
int seed;
double x;
double x2;

a = 50.0;
b = 150.0;
mu = 100.0;
s = 25.0;
seed = 123456789;

Console.WriteLine("");
Console.WriteLine("NORMAL_TRUNCATED_AB_CDF_TEST");
Console.WriteLine("  NORMAL_TRUNCATED_AB_CDF evaluates the Normal Truncated AB CDF.");
Console.WriteLine("  NORMAL_TRUNCATED_AB_CDF_INV inverts the Normal Truncated AB CDF.");
Console.WriteLine("  NORMAL_TRUNCATED_AB_PDF evaluates the Normal Truncated AB PDF.");
Console.WriteLine("");
Console.WriteLine("  The parent normal distribution has");
Console.WriteLine("    mean =               " + mu + "");
Console.WriteLine("    standard deviation = " + s + "");
Console.WriteLine("  The parent distribution is truncated to");
Console.WriteLine("  the interval [" + a + "," + b + "]");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = normal_truncated_ab_sample(mu, s, a, b, seed);

pdf = normal_truncated_ab_pdf(x, mu, s, a, b);

cdf = normal_truncated_ab_cdf(x, mu, s, a, b);

x2 = normal_truncated_ab_cdf_inv(cdf, mu, s, a, b);

Console.WriteLine("  " + setw(14) + x
+ "  " + pdf.ToString().PadLeft(14)
+ "  " + cdf.ToString().PadLeft(14)
+ "  " + setw(14) + x2 + "");
}

return;
}

static void normal_truncated_ab_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_TRUNCATED_AB_SAMPLE_TEST tests NORMAL_TRUNCATED_AB_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
int i;
double mean;
double mu;
double s;
int sample_num = 1000;
int seed;
double variance;
double[] x;
double xmax;
double xmin;

a = 50.0;
b = 150.0;
mu = 100.0;
s = 25.0;
seed = 123456789;

Console.WriteLine("");
Console.WriteLine("NORMAL_TRUNCATED_AB_SAMPLE_TEST");
Console.WriteLine("  NORMAL_TRUNCATED_AB_MEAN computes the Normal Truncated AB mean;");
Console.WriteLine("  NORMAL_TRUNCATED_AB_SAMPLE samples the Normal Truncated AB distribution;");
Console.WriteLine("  NORMAL_TRUNCATED_AB_VARIANCE computes the Normal Truncated AB variance.");
Console.WriteLine("");
Console.WriteLine("  The parent normal distribution has");
Console.WriteLine("    mean =               " + mu + "");
Console.WriteLine("    standard deviation = " + s + "");
Console.WriteLine("  The parent distribution is truncated to");
Console.WriteLine("  the interval [" + a + "," + b + "]");

mean = normal_truncated_ab_mean(mu, s, a, b);

variance = normal_truncated_ab_variance(mu, s, a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean      =               " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

x = (double*) malloc(sample_num * sizeof(double));

for (i = 0; i < sample_num; i++)
{
x[i] = normal_truncated_ab_sample(mu, s, a, b, seed);
}

mean = typeMethods.r8vec_mean(sample_num, x);
variance = typeMethods.r8vec_variance(sample_num, x);
xmax = typeMethods.r8vec_max(sample_num, x);
xmin = typeMethods.r8vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");



return;
}

static void normal_truncated_a_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_TRUNCATED_A_CDF_TEST tests NORMAL_TRUNCATED_A_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double cdf;
int i;
double mu;
double pdf;
double s;
int seed;
double x;
double x2;

a = 50.0;
mu = 100.0;
s = 25.0;
seed = 123456789;

Console.WriteLine("");
Console.WriteLine("NORMAL_TRUNCATED_A_CDF_TEST");
Console.WriteLine("  NORMAL_TRUNCATED_A_CDF evaluates the Normal Truncated A CDF.");
Console.WriteLine("  NORMAL_TRUNCATED_A_CDF_INV inverts the Normal Truncated A CDF.");
Console.WriteLine("  NORMAL_TRUNCATED_A_PDF evaluates the Normal Truncated A PDF.");
Console.WriteLine("");
Console.WriteLine("  The parent normal distribution has");
Console.WriteLine("    mean =               " + mu + "");
Console.WriteLine("    standard deviation = " + s + "");
Console.WriteLine("  The parent distribution is truncated to");
Console.WriteLine("  the interval [" + a + ",+oo]");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = normal_truncated_a_sample(mu, s, a, seed);

pdf = normal_truncated_a_pdf(x, mu, s, a);

cdf = normal_truncated_a_cdf(x, mu, s, a);

x2 = normal_truncated_a_cdf_inv(cdf, mu, s, a);

Console.WriteLine("  " + setw(14) + x
+ "  " + pdf.ToString().PadLeft(14)
+ "  " + cdf.ToString().PadLeft(14)
+ "  " + setw(14) + x2 + "");
}

return;
}

static void normal_truncated_a_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_TRUNCATED_A_SAMPLE_TEST tests NORMAL_TRUNCATED_A_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
int i;
double mean;
double mu;
double s;
int sample_num = 1000;
int seed;
double variance;
double[] x;
double xmax;
double xmin;

a = 50.0;
mu = 100.0;
s = 25.0;
seed = 123456789;

Console.WriteLine("");
Console.WriteLine("NORMAL_TRUNCATED_A_SAMPLE_TEST");
Console.WriteLine("  NORMAL_TRUNCATED_A_MEAN computes the Normal Truncated A mean;");
Console.WriteLine("  NORMAL_TRUNCATED_A_SAMPLE samples the Normal Truncated A distribution;");
Console.WriteLine("  NORMAL_TRUNCATED_A_VARIANCE computes the Normal Truncated A variance.");
Console.WriteLine("");
Console.WriteLine("  The parent normal distribution has");
Console.WriteLine("    mean =               " + mu + "");
Console.WriteLine("    standard deviation = " + s + "");
Console.WriteLine("  The parent distribution is truncated to");
Console.WriteLine("  the interval [" + a + ",+oo]");

mean = normal_truncated_a_mean(mu, s, a);

variance = normal_truncated_a_variance(mu, s, a);

Console.WriteLine("");
Console.WriteLine("  PDF mean      =               " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

x = (double*) malloc(sample_num * sizeof(double));

for (i = 0; i < sample_num; i++)
{
x[i] = normal_truncated_a_sample(mu, s, a, seed);
}

mean = typeMethods.r8vec_mean(sample_num, x);
variance = typeMethods.r8vec_variance(sample_num, x);
xmax = typeMethods.r8vec_max(sample_num, x);
xmin = typeMethods.r8vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");



return;
}

static void normal_truncated_b_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_TRUNCATED_B_CDF_TEST tests NORMAL_TRUNCATED_B_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double b;
double cdf;
int i;
double mu;
double pdf;
double s;
int seed;
double x;
double x2;

b = 150.0;
mu = 100.0;
s = 25.0;
seed = 123456789;

Console.WriteLine("");
Console.WriteLine("NORMAL_TRUNCATED_B_CDF_TEST");
Console.WriteLine("  NORMAL_TRUNCATED_B_CDF evaluates the Normal Truncated B CDF.");
Console.WriteLine("  NORMAL_TRUNCATED_B_CDF_INV inverts the Normal Truncated B CDF.");
Console.WriteLine("  NORMAL_TRUNCATED_B_PDF evaluates the Normal Truncated B PDF.");
Console.WriteLine("");
Console.WriteLine("  The parent normal distribution has");
Console.WriteLine("    mean =               " + mu + "");
Console.WriteLine("    standard deviation = " + s + "");
Console.WriteLine("  The parent distribution is truncated to");
Console.WriteLine("  the interval [-oo," + b + "]");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = normal_truncated_b_sample(mu, s, b, seed);

pdf = normal_truncated_b_pdf(x, mu, s, b);

cdf = normal_truncated_b_cdf(x, mu, s, b);

x2 = normal_truncated_b_cdf_inv(cdf, mu, s, b);

Console.WriteLine("  " + setw(14) + x
+ "  " + pdf.ToString().PadLeft(14)
+ "  " + cdf.ToString().PadLeft(14)
+ "  " + setw(14) + x2 + "");
}

return;
}

static void normal_truncated_b_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_TRUNCATED_B_SAMPLE_TEST tests NORMAL_TRUNCATED_B_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double b;
int i;
double mean;
double mu;
double s;
int sample_num = 1000;
int seed;
double variance;
double[] x;
double xmax;
double xmin;

b = 150.0;
mu = 100.0;
s = 25.0;
seed = 123456789;

Console.WriteLine("");
Console.WriteLine("NORMAL_TRUNCATED_B_SAMPLE_TEST");
Console.WriteLine("  NORMAL_TRUNCATED_B_MEAN computes the Normal Truncated B mean;");
Console.WriteLine("  NORMAL_TRUNCATED_B_SAMPLE samples the Normal Truncated B distribution;");
Console.WriteLine("  NORMAL_TRUNCATED_B_VARIANCE computes the Normal Truncated B variance.");
Console.WriteLine("");
Console.WriteLine("  The parent normal distribution has");
Console.WriteLine("    mean =               " + mu + "");
Console.WriteLine("    standard deviation = " + s + "");
Console.WriteLine("  The parent distribution is truncated to");
Console.WriteLine("  the interval [-oo," + b + "]");

mean = normal_truncated_b_mean(mu, s, b);

variance = normal_truncated_b_variance(mu, s, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean      =               " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

x = (double*) malloc(sample_num * sizeof(double));

for (i = 0; i < sample_num; i++)
{
x[i] = normal_truncated_b_sample(mu, s, b, seed);
}

mean = typeMethods.r8vec_mean(sample_num, x);
variance = typeMethods.r8vec_variance(sample_num, x);
xmax = typeMethods.r8vec_max(sample_num, x);
xmin = typeMethods.r8vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");



return;
}

static void pareto_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_CDF_TEST tests PARETO_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
Console.WriteLine("PARETO_CDF_TEST");
Console.WriteLine("  PARETO_CDF evaluates the Pareto CDF;");
Console.WriteLine("  PARETO_CDF_INV inverts the Pareto CDF.");
Console.WriteLine("  PARETO_PDF evaluates the Pareto PDF;");

a = 0.5;
b = 5.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!pareto_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("PARETO_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = pareto_sample(a, b, seed);
pdf = pareto_pdf(x, a, b);
cdf = pareto_cdf(x, a, b);
x2 = pareto_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void pareto_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_SAMPLE_TEST tests PARETO_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
Console.WriteLine("PARETO_SAMPLE_TEST");
Console.WriteLine("  PARETO_MEAN computes the Pareto mean;");
Console.WriteLine("  PARETO_SAMPLE samples the Pareto distribution;");
Console.WriteLine("  PARETO_VARIANCE computes the Pareto variance;");

a = 0.5;
b = 5.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!pareto_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("PARETO_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = pareto_mean(a, b);
variance = pareto_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = pareto_sample(a, b, seed);
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

static void pearson_05_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_PDF_TEST tests PEARSON_05_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
double c;
double pdf;
double x;

Console.WriteLine("");
Console.WriteLine("PEARSON_05_PDF");
Console.WriteLine("  PEARSON_05_PDF evaluates the Pearson 05 PDF.");

x = 5.0;

a = 1.0;
b = 2.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + a + "");
Console.WriteLine("  PDF parameter B = " + b + "");
Console.WriteLine("  PDF parameter C = " + c + "");

if (!pearson_05_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("PEARSON_05_PDF - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

pdf = pearson_05_pdf(x, a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF argument X =  " + x + "");
Console.WriteLine("  PDF value =       " + pdf + "");

return;
}

static void planck_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_PDF_TEST tests PLANCK_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
int i;
double pdf;
int seed = 123456789;
double x;

Console.WriteLine("");
Console.WriteLine("PLANCK_PDF_TEST");
Console.WriteLine("  PLANCK_PDF evaluates the Planck PDF.");
Console.WriteLine("  PLANCK_SAMPLE samples the Planck PDF.");

a = 2.0E+00;
b = 3.0E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + a + "");
Console.WriteLine("  PDF parameter B = " + b + "");

if (!planck_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("PLANCK_PDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = planck_sample(a, b, seed);

pdf = planck_pdf(x, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "");
}

return;
}

static void planck_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_SAMPLE_TEST tests PLANCK_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
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
Console.WriteLine("PLANCK_SAMPLE_TEST");
Console.WriteLine("  PLANCK_MEAN computes the Planck mean;");
Console.WriteLine("  PLANCK_SAMPLE samples the Planck distribution;");
Console.WriteLine("  PLANCK_VARIANCE computes the Planck variance;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!planck_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("PLANCK_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = planck_mean(a, b);
variance = planck_variance(a, b);

Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = planck_sample(a, b, seed);
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

static void poisson_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF_TEST tests POISSON_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double cdf;
int i;
double pdf;
int seed = 123456789;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("POISSON_CDF_TEST");
Console.WriteLine("  POISSON_CDF evaluates the Poisson CDF;");
Console.WriteLine("  POISSON_CDF_INV inverts the Poisson CDF.");
Console.WriteLine("  POISSON_PDF evaluates the Poisson PDF;");

a = 10.0E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!poisson_check(a))
{
Console.WriteLine("");
Console.WriteLine("POISSON_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = poisson_sample(a, seed);
pdf = poisson_pdf(x, a);
cdf = poisson_cdf(x, a);
x2 = poisson_cdf_inv(cdf, a);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void poisson_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_SAMPLE_TEST tests POISSON_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

double a;
int i;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("POISSON_SAMPLE_TEST");
Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");
Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");
Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");

a = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");

if (!poisson_check(a))
{
Console.WriteLine("");
Console.WriteLine("POISSON_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = poisson_mean(a);
variance = poisson_variance(a);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = poisson_sample(a, seed);
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

static void power_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    POWER_CDF_TEST tests POWER_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2016
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
Console.WriteLine("POWER_CDF_TEST");
Console.WriteLine("  POWER_CDF evaluates the Power CDF;");
Console.WriteLine("  POWER_CDF_INV inverts the Power CDF.");
Console.WriteLine("  POWER_PDF evaluates the Power PDF;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!power_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("POWER_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = power_sample(a, b, seed);
pdf = power_pdf(x, a, b);
cdf = power_cdf(x, a, b);
x2 = power_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void power_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    POWER_SAMPLE_TEST tests POWER_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2016
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
Console.WriteLine("POWER_SAMPLE_TEST");
Console.WriteLine("  POWER_MEAN computes the Power mean;");
Console.WriteLine("  POWER_SAMPLE samples the Power distribution;");
Console.WriteLine("  POWER_VARIANCE computes the Power variance;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!power_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("POWER_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = power_mean(a, b);
variance = power_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = power_sample(a, b, seed);
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

static void quasigeometric_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_CDF_TEST tests QUASIGEOMETRIC_CDF.
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
double cdf;
int i;
double pdf;
int seed = 123456789;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("QUASIGEOMETRIC_CDF_TEST");
Console.WriteLine("  QUASIGEOMETRIC_CDF evaluates the Quasigeometric CDF;");
Console.WriteLine("  QUASIGEOMETRIC_CDF_INV inverts the Quasigeometric CDF.");
Console.WriteLine("  QUASIGEOMETRIC_PDF evaluates the Quasigeometric PDF;");

a = 0.4825;
b = 0.5893;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + a + "");
Console.WriteLine("  PDF parameter B = " + b + "");

if (!quasigeometric_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("QUASIGEOMETRIC_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = quasigeometric_sample(a, b, seed);
pdf = quasigeometric_pdf(x, a, b);
cdf = quasigeometric_cdf(x, a, b);
x2 = quasigeometric_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void quasigeometric_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_SAMPLE_TEST tests QUASIGEOMETRIC_SAMPLE.
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
int j;
double mean;
int sample_num = 1000;
int seed = 123456789;
double variance;
int* x;
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("QUASIGEOMETRIC_SAMPLE_TEST");
Console.WriteLine("  QUASIGEOMETRIC_MEAN computes the Quasigeometric mean;");
Console.WriteLine("  QUASIGEOMETRIC_SAMPLE samples the Quasigeometric distribution;");
Console.WriteLine("  QUASIGEOMETRIC_VARIANCE computes the Quasigeometric variance.");

a = 0.4825;
b = 0.5893;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + a + "");
Console.WriteLine("  PDF parameter B = " + b + "");

if (!quasigeometric_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("QUASIGEOMETRIC_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = quasigeometric_mean(a, b);
variance = quasigeometric_variance(a, b);

Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

x = new int[sample_num];

for (j = 0; j < sample_num; j++)
{
x[j] = quasigeometric_sample(a, b, seed);
}

mean = typeMethods.i4vec_mean(sample_num, x);
variance = typeMethods.i4vec_variance(sample_num, x);
xmax = typeMethods.i4vec_max(sample_num, x);
xmin = typeMethods.i4vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");



return;
}

static void r8_beta_test()

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
//    04 March 2016
//
//  Author:
//
//    John Burkardt
//
{
double fxy1;
double fxy2;
int n_data;
streamsize ss;
double x;
double y;
//
//  Save the current precision.
//
ss = cout.precision();

Console.WriteLine("");
Console.WriteLine("R8_BETA_TEST:");
Console.WriteLine("  R8_BETA evaluates the Beta function.");
Console.WriteLine("");
Console.WriteLine("      X              Y         BETA(X,Y)         R8_BETA(X,Y)");
Console.WriteLine("                               tabulated         computed.");
Console.WriteLine("");

n_data = 0;

for (;;)
{
beta_values(n_data, x, y, fxy1);

if (n_data == 0)
{
break;
}

fxy2 = r8_beta(x, y);

Console.WriteLine("  " + x.ToString().PadLeft(12)
+ "  " + y.ToString().PadLeft(12)
+ "  " + setw(24) + setprecision(16) + fxy1
+ "  " + setw(24) + setprecision(16) + fxy2 + "");
}

//
//  Restore the default precision.
//
cout.precision(ss);

return;
}

static void r8_ceiling_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_CEILING_TEST tests R8_CEILING.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
int i;
int ival;
double rval;

Console.WriteLine("");
Console.WriteLine("R8_CEILING_TEST");
Console.WriteLine("  R8_CEILING rounds an R8 up.");
Console.WriteLine("");
Console.WriteLine("       X           R8_CEILING(X)");
Console.WriteLine("");

for (i = -6; i <= 6; i++)
{
rval = (double) (i) / 5.0;
ival = r8_ceiling(rval);
Console.WriteLine("  "
+ setw(14) + rval + "  "
+ setw(6) + ival + "");
}

return;
}

static void r8_error_f_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERROR_F_TEST tests R8_ERROR_F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int i;
int seed;
double x;
double y;
double z;

Console.WriteLine("");
Console.WriteLine("R8_ERROR_F_TEST");
Console.WriteLine("  R8_ERROR_F evaluates ERF(X).");
Console.WriteLine("");
Console.WriteLine("X   -> Y = R8_ERROR_F(X) -> Z = R8_ERROR_F_INVERSE(Y)");
Console.WriteLine("");

seed = 123456789;

x = 1.0;

for (i = 1; i <= 20; i++)
{
x = normal_01_sample(seed);
y = r8_error_f(x);
z = r8_error_f_inverse(y);
Console.WriteLine("  " + setw(14) + x
+ "  " + setw(14) + y
+ "  " + setw(14) + z + "");
}

return;
}

static void r8_factorial_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_TEST tests R8_FACTORIAL.
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
double f;
int i;

Console.WriteLine("");
Console.WriteLine("R8_FACTORIAL_TEST");
Console.WriteLine("  R8_FACTORIAL evaluates the factorial function.");
Console.WriteLine("");
Console.WriteLine("    I                R8_FACTORIAL(I)");
Console.WriteLine("");

for (i = 0; i <= 20; i++)
{
f = r8_factorial(i);

Console.WriteLine("  "
+ setw(4) + i + "  "
+ setw(24) + f + "");
}

return;
}

static void r8_gamma_inc_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_INC_TEST tests R8_GAMMA_INC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
double a;
double fx;
double fx2;
int n_data;
double x;

Console.WriteLine("");
Console.WriteLine("R8_GAMMA_INC_TEST:");
Console.WriteLine("  R8_GAMMA_INC evaluates the normalized incomplete Gamma");
Console.WriteLine("  function.");
Console.WriteLine("");
Console.WriteLine("   A      X       Exact F       R8_GAMMA_INC(A,X)");
Console.WriteLine("");

n_data = 0;

for (;;)
{
gamma_inc_values(n_data, a, x, fx);

if (n_data == 0)
{
break;
}

fx2 = r8_gamma_inc(a, x);

Console.WriteLine("  "
+ setw(8) + a + "  "
+ setw(8) + x + "  "
+ setw(16) + fx + "  "
+ setw(16) + fx2 + "");
}

return;
}

static void r8_gamma_log_int_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG_INT_TEST tests R8_GAMMA_LOG_INT;
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
double g;
int i;

Console.WriteLine("");
Console.WriteLine("R8_GAMMA_LOG_INT_TEST");
Console.WriteLine("  R8_GAMMA_LOG_INT evaluates the logarithm of the gamma function");
Console.WriteLine("  for integer argument.");

Console.WriteLine("");
Console.WriteLine("       I     R8_GAMMA_LOG_INT(I)");
Console.WriteLine("");

for (i = 1; i <= 20; i++)
{
g = r8_gamma_log_int(i);

Console.WriteLine("  "
+ setw(6) + i + "  "
+ setw(12) + g + "");
}

return;
}

static void r8_uniform_01_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 1000

int i;
double max;
double mean;
double min;
int seed = 123456789;
double x[N];
double variance;

Console.WriteLine("");
Console.WriteLine("R8_UNIFORM_01_TEST");
Console.WriteLine("  R8_UNIFORM_01 samples a uniform random distribution in [0,1].");
Console.WriteLine("  distributed random numbers.");
Console.WriteLine("  Using initial random number seed = " + seed + "");

for (i = 0; i < N; i++)
{
x[i] = r8_uniform_01(seed);
}

Console.WriteLine("");
Console.WriteLine("  First few values:");
Console.WriteLine("");
for (i = 0; i < 10; i++)
{
Console.WriteLine("  " + setw(6) + i
+ "  " + setw(14) + x[i] + "");
}

min = typeMethods.r8vec_min(N, x);
max = typeMethods.r8vec_max(N, x);
mean = typeMethods.r8vec_mean(N, x);
variance = typeMethods.r8vec_variance(N, x);

Console.WriteLine("");
Console.WriteLine("  Number of samples was " + N + "");
Console.WriteLine("  Minimum value was " + min + "");
Console.WriteLine("  Maximum value was " + max + "");
Console.WriteLine("  Average value was " + mean + "");
Console.WriteLine("  Variance was      " + variance + "");

return;
# undef N
}

static void r8_zeta_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_ZETA_TEST tests R8_ZETA.
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
int i;
double p;
double v;

Console.WriteLine("");
Console.WriteLine("R8_ZETA_TEST");
Console.WriteLine("  R8_ZETA estimates the Zeta function.");

Console.WriteLine("");
Console.WriteLine("       P     R8_Zeta(P)");
Console.WriteLine("");
for (i = 1; i <= 25; i++)
{
p = (double) (i);
v = r8_zeta(p);
Console.WriteLine("  " + setw(6) + p
+ "  " + setw(14) + v + "");
}

Console.WriteLine("");
for (i = 0; i <= 8; i++)
{
p = 3.0 + (double) (i) / 8.0;
v = r8_zeta(p);
Console.WriteLine("  " + setw(6) + p
+ "  " + setw(14) + v + "");
}

return;
}

static void rayleigh_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CDF_TEST tests RAYLEIGH_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
double a;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("RAYLEIGH_CDF_TEST");
Console.WriteLine("  RAYLEIGH_CDF evaluates the Rayleigh CDF;");
Console.WriteLine("  RAYLEIGH_CDF_INV inverts the Rayleigh CDF.");
Console.WriteLine("  RAYLEIGH_PDF evaluates the Rayleigh PDF;");

a = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!rayleigh_check(a))
{
Console.WriteLine("");
Console.WriteLine("RAYLEIGH_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = rayleigh_sample(a, seed);
pdf = rayleigh_pdf(x, a);
cdf = rayleigh_cdf(x, a);
x2 = rayleigh_cdf_inv(cdf, a);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void rayleigh_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_SAMPLE_TEST tests RAYLEIGH_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
double[] x = new double [SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("RAYLEIGH_SAMPLE_TEST");
Console.WriteLine("  RAYLEIGH_MEAN computes the Rayleigh mean;");
Console.WriteLine("  RAYLEIGH_SAMPLE samples the Rayleigh distribution;");
Console.WriteLine("  RAYLEIGH_VARIANCE computes the Rayleigh variance.");

a = 2.0E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!rayleigh_check(a))
{
Console.WriteLine("");
Console.WriteLine("RAYLEIGH_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = rayleigh_mean(a);
variance = rayleigh_variance(a);

Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (j = 0; j < SAMPLE_NUM; j++)
{
x[j] = rayleigh_sample(a, seed);
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

static void reciprocal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_CDF_TEST tests RECIPROCAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2016
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
Console.WriteLine("RECIPROCAL_CDF_TEST");
Console.WriteLine("  RECIPROCAL_CDF evaluates the Reciprocal CDF;");
Console.WriteLine("  RECIPROCAL_CDF_INV inverts the Reciprocal CDF.");
Console.WriteLine("  RECIPROCAL_PDF evaluates the Reciprocal PDF;");

a = 1.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!reciprocal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("RECIPROCAL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = reciprocal_sample(a, b, seed);
pdf = reciprocal_pdf(x, a, b);
cdf = reciprocal_cdf(x, a, b);
x2 = reciprocal_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void reciprocal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_SAMPLE_TEST tests RECIPROCAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
Console.WriteLine("RECIPROCAL_SAMPLE_TEST");
Console.WriteLine("  RECIPROCAL_MEAN computes the Reciprocal mean;");
Console.WriteLine("  RECIPROCAL_SAMPLE samples the Reciprocal distribution;");
Console.WriteLine("  RECIPROCAL_VARIANCE computes the Reciprocal variance;");

a = 1.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!reciprocal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("RECIPROCAL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = reciprocal_mean(a, b);
variance = reciprocal_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = reciprocal_sample(a, b, seed);
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

static void runs_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_PDF_TEST tests RUNS_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int m;
int n;
double pdf;
double pdf_total;
int r;

Console.WriteLine("");
Console.WriteLine("RUNS_PDF_TEST");
Console.WriteLine("  RUNS_PDF evaluates the Runs PDF;");
Console.WriteLine("");
Console.WriteLine("  M is the number of symbols of one kind,");
Console.WriteLine("  N is the number of symbols of the other kind,");
Console.WriteLine("  R is the number of runs (sequences of one symbol)");
Console.WriteLine("");
Console.WriteLine("         M         N         R      PDF");
Console.WriteLine("");

m = 6;

for (n = 0; n <= 9; n++)
{
Console.WriteLine("");
pdf_total = 0.0;

for (r = 1; r <= 2 * i4_min(m, n) + 2; r++)
{
pdf = runs_pdf(m, n, r);

Console.WriteLine("  " + setw(8) + m
  + "  " + setw(8) + n
  + "  " + setw(8) + r
  + "  " + pdf.ToString().PadLeft(14) + "");

pdf_total = pdf_total + pdf;
}

Console.WriteLine("  " + setw(8) + m
+ "  " + "        "
+ "  " + "        "
+ "  " + pdf.ToString().PadLeft(14)_total + "");

}

return;
}

static void runs_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_SAMPLE_TEST tests RUNS_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int SAMPLE_NUM = 1000;

int i;
int m;
double mean;
int n;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("RUNS_SAMPLE_TEST");
Console.WriteLine("  RUNS_MEAN computes the Runs mean;");
Console.WriteLine("  RUNS_SAMPLE samples the Runs distribution;");
Console.WriteLine("  RUNS_VARIANCE computes the Runs variance");

m = 10;
n = 5;

Console.WriteLine("");
Console.WriteLine("  PDF parameter M = " + m + "");
Console.WriteLine("  PDF parameter N = " + n + "");

mean = runs_mean(m, n);
variance = runs_variance(m, n);

Console.WriteLine("  PDF mean =        " + mean + "");
Console.WriteLine("  PDF variance =    " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = runs_sample(m, n, seed);
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

static void sech_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    SECH_CDF_TEST tests SECH_CDF.
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
Console.WriteLine("SECH_CDF_TEST");
Console.WriteLine("  SECH_CDF evaluates the Sech CDF;");
Console.WriteLine("  SECH_CDF_INV inverts the Sech CDF.");
Console.WriteLine("  SECH_PDF evaluates the Sech PDF;");

a = 3.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!sech_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("SECH_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = sech_sample(a, b, seed);
pdf = sech_pdf(x, a, b);
cdf = sech_cdf(x, a, b);
x2 = sech_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void sech_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    SECH_SAMPLE_TEST tests SECH_SAMPLE.
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
Console.WriteLine("SECH_SAMPLE_TEST");
Console.WriteLine("  SECH_MEAN computes the Sech mean;");
Console.WriteLine("  SECH_SAMPLE samples the Sech distribution;");
Console.WriteLine("  SECH_VARIANCE computes the Sech variance;");

a = 3.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!sech_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("SECH_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = sech_mean(a, b);
variance = sech_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = sech_sample(a, b, seed);
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

static void semicircular_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_CDF_TEST tests SEMICIRCULAR_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2016
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
Console.WriteLine("SEMICIRCULAR_CDF_TEST");
Console.WriteLine("  SEMICIRCULAR_CDF evaluates the Semicircular CDF;");
Console.WriteLine("  SEMICIRCULAR_CDF_INV inverts the Semicircular CDF.");
Console.WriteLine("  SEMICIRCULAR_PDF evaluates the Semicircular PDF;");

a = 3.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!semicircular_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("SEMICIRCULAR_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = semicircular_sample(a, b, seed);
pdf = semicircular_pdf(x, a, b);
cdf = semicircular_cdf(x, a, b);
x2 = semicircular_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void semicircular_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_SAMPLE_TEST tests SEMICIRCULAR_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2016
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
Console.WriteLine("SEMICIRCULAR_SAMPLE_TEST");
Console.WriteLine("  SEMICIRCULAR_MEAN computes the Semicircular mean;");
Console.WriteLine("  SEMICIRCULAR_SAMPLE samples the Semicircular distribution;");
Console.WriteLine("  SEMICIRCULAR_VARIANCE computes the Semicircular variance;");

a = 3.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!semicircular_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("SEMICIRCULAR_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = semicircular_mean(a, b);
variance = semicircular_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = semicircular_sample(a, b, seed);
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

static void student_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CDF_TEST tests STUDENT_CDF.
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
double c;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;

Console.WriteLine("");
Console.WriteLine("STUDENT_CDF_TEST");
Console.WriteLine("  STUDENT_CDF evaluates the Student CDF;");
Console.WriteLine("  STUDENT_PDF evaluates the Student PDF;");
Console.WriteLine("  STUDENT_SAMPLE samples the Student PDF;");

a = 0.5;
b = 2.0;
c = 6.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + a + "");
Console.WriteLine("  PDF parameter B = " + b + "");
Console.WriteLine("  PDF parameter C = " + c + "");

if (!student_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("STUDENT_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameter values are illegal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = student_sample(a, b, c, seed);
pdf = student_pdf(x, a, b, c);
cdf = student_cdf(x, a, b, c);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "");
}

return;
}

static void student_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_SAMPLE_TEST tests STUDENT_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2016
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
Console.WriteLine("STUDENT_SAMPLE_TEST");
Console.WriteLine("  STUDENT_MEAN evaluates the Student mean;");
Console.WriteLine("  STUDENT_SAMPLE samples the Student PDF;");
Console.WriteLine("  STUDENT_VARIANCE computes the Student variance;");

a = 0.5;
b = 2.0;
c = 6.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + a + "");
Console.WriteLine("  PDF parameter B = " + b + "");
Console.WriteLine("  PDF parameter C = " + c + "");

if (!student_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("STUDENT_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameter values are illegal.");
return;
}

mean = student_mean(a, b, c);
variance = student_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = student_sample(a, b, c, seed);
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

static void student_noncentral_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_NONCENTRAL_CDF_TEST tests STUDENT_NONCENTRAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double b;
double cdf;
int idf;
double x;

Console.WriteLine("");
Console.WriteLine("STUDENT_NONCENTRAL_CDF_TEST");
Console.WriteLine("  STUDENT_NONCENTRAL_CDF evaluates the Student Noncentral CDF;");

x = 0.50;
idf = 10;
b = 1.0;

cdf = student_noncentral_cdf(x, idf, b);

Console.WriteLine("");
Console.WriteLine("  PDF argument X =              " + x + "");
Console.WriteLine("  PDF parameter IDF =           " + idf + "");
Console.WriteLine("  PDF parameter B =             " + b + "");
Console.WriteLine("  CDF value =                   " + cdf + "");

return;
}

static void tfn_test()

//****************************************************************************80
//
//  Purpose:
//
//    TFN_TEST tests TFN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double h;
int n_data;
double t;
double t2;

Console.WriteLine("");
Console.WriteLine("TFN_TEST");
Console.WriteLine("  TFN evaluates Owen's T function;");
Console.WriteLine("");
Console.WriteLine("      H             A           T(H,A)          Exact");
Console.WriteLine("");

n_data = 0;

for (;;)
{
owen_values(n_data, h, a, t);

if (n_data <= 0)
{
break;
}

t2 = tfn(h, a);

Console.WriteLine("  "
+ setw(14) + h + "  "
+ setw(14) + a + "  "
+ setw(14) + t2 + "  "
+ setw(14) + t + "");
}

return;
}

static void triangle_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CDF_TEST tests TRIANGLE_CDF.
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
double c;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("TRIANGLE_CDF_TEST");
Console.WriteLine("  TRIANGLE_CDF evaluates the Triangle CDF;");
Console.WriteLine("  TRIANGLE_CDF_INV inverts the Triangle CDF.");
Console.WriteLine("  TRIANGLE_PDF evaluates the Triangle PDF;");

a = 1.0;
b = 3.0;
c = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!triangle_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("TRIANGLE_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = triangle_sample(a, b, c, seed);
pdf = triangle_pdf(x, a, b, c);
cdf = triangle_cdf(x, a, b, c);
x2 = triangle_cdf_inv(cdf, a, b, c);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void triangle_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE_TEST tests TRIANGLE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2016
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
Console.WriteLine("TRIANGLE_SAMPLE_TEST");
Console.WriteLine("  TRIANGLE_MEAN computes the Triangle mean;");
Console.WriteLine("  TRIANGLE_SAMPLE samples the Triangle distribution;");
Console.WriteLine("  TRIANGLE_VARIANCE computes the Triangle variance;");

a = 1.0;
b = 3.0;
c = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =        " + a + "");
Console.WriteLine("  PDF parameter B =        " + b + "");
Console.WriteLine("  PDF parameter C =        " + c + "");

if (!triangle_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("TRIANGLE_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = triangle_mean(a, b, c);
variance = triangle_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = triangle_sample(a, b, c, seed);
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

static void triangular_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_CDF_TEST tests TRIANGULAR_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2016
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
Console.WriteLine("TRIANGULAR_CDF_TEST");
Console.WriteLine("  TRIANGULAR_CDF evaluates the Triangular CDF;");
Console.WriteLine("  TRIANGULAR_CDF_INV inverts the Triangular CDF.");
Console.WriteLine("  TRIANGULAR_PDF evaluates the Triangular PDF;");

a = 1.0;
b = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!triangular_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("TRIANGULAR_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = triangular_sample(a, b, seed);
pdf = triangular_pdf(x, a, b);
cdf = triangular_cdf(x, a, b);
x2 = triangular_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ x.ToString().PadLeft(12) + "  "
+ pdf.ToString().PadLeft(12) + "  "
+ cdf.ToString().PadLeft(12) + "  "
+ x2.ToString().PadLeft(12) + "");
}

return;
}

static void triangular_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_SAMPLE_TEST tests TRIANGULAR_SAMPLE.
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
Console.WriteLine("TRIANGULAR_SAMPLE_TEST");
Console.WriteLine("  TRIANGULAR_MEAN computes the Triangular mean;");
Console.WriteLine("  TRIANGULAR_SAMPLE samples the Triangular distribution;");
Console.WriteLine("  TRIANGULAR_VARIANCE computes the Triangular variance;");

a = 1.0;
b = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!triangular_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("TRIANGULAR_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = triangular_mean(a, b);
variance = triangular_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = triangular_sample(a, b, seed);
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