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













static void disk_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DISK_SAMPLE_TEST tests DISK_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2016
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

double a;
double b;
double c;
int j;
double* mean;
int seed = 123456789;
double variance;
double x[2 * SAMPLE_NUM];
double* xmax;
double* xmin;
double* y;

Console.WriteLine("");
Console.WriteLine("DISK_SAMPLE");
Console.WriteLine("  DISK_MEAN returns the Disk mean;");
Console.WriteLine("  DISK_SAMPLE samples the Disk distribution;");
Console.WriteLine("  DISK_VARIANCE returns the Disk variance;");

a = 10.0;
b = 4.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

mean = disk_mean(a, b, c);
variance = disk_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  Disk mean ="
+ "  " + setw(12) + mean[0]
+ "  " + setw(12) + mean[1] + "");
Console.WriteLine("  Disk variance = " + setw(12) + variance + "");

delete[] mean;

for (j = 0; j < SAMPLE_NUM; j++)
{
y = disk_sample(a, b, c, seed);
x[0 + j * 2] = y[0];
x[1 + j * 2] = y[1];
delete[] y;
}

variance = 0.0;
for (j = 0; j < SAMPLE_NUM; j++)
{
variance = variance + pow(x[0 + j * 2] - a, 2) + pow(x[1 + j * 2] - b, 2);
}

variance = variance / (double) (SAMPLE_NUM);


mean = r8row_mean(2, SAMPLE_NUM, x);
xmax = r8row_max(2, SAMPLE_NUM, x);
xmin = r8row_min(2, SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     "
+ setw(12) + mean[0] + "  "
+ setw(12) + mean[1] + "");
Console.WriteLine("  Sample variance = "
+ setw(12) + variance + "");
Console.WriteLine("  Sample maximum =  "
+ setw(12) + xmax[0] + "  "
+ setw(12) + xmax[1] + "");
Console.WriteLine("  Sample minimum =  "
+ setw(12) + xmin[0] + "  "
+ setw(12) + xmin[1] + "");

delete[] mean;
delete[] xmax;
delete[] xmin;

return;
# undef SAMPLE_NUM
}

static void empirical_discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_CDF_TEST tests EMPIRICAL_DISCRETE_CDF.
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
# define A 6

double b[A] =  {
1.0, 1.0, 3.0, 2.0, 1.0, 2.0
}
;
double c[A] =  {
0.0, 1.0, 2.0, 4.5, 6.0, 10.0
}
;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("EMPIRICAL_DISCRETE_CDF_TEST");
Console.WriteLine("  EMPIRICAL_DISCRETE_CDF evaluates the Empirical Discrete CDF;");
Console.WriteLine("  EMPIRICAL_DISCRETE_CDF_INV inverts the Empirical Discrete CDF.");
Console.WriteLine("  EMPIRICAL_DISCRETE_PDF evaluates the Empirical Discrete PDF;");

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + A + "");
r8vec_print(A, b, "  PDF parameter B = ");
r8vec_print(A, c, "  PDF parameter C = ");

if (!empirical_discrete_check(A, b, c))
{
Console.WriteLine("");
Console.WriteLine("EMPIRICAL_DISCRETE_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = empirical_discrete_sample(A, b, c, seed);
pdf = empirical_discrete_pdf(x, A, b, c);
cdf = empirical_discrete_cdf(x, A, b, c);
x2 = empirical_discrete_cdf_inv(cdf, A, b, c);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
# undef A
}

static void empirical_discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_SAMPLE_TEST tests EMPIRICAL_DISCRETE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define A 6
# define SAMPLE_NUM 1000

double b[A] =  {
1.0, 1.0, 3.0, 2.0, 1.0, 2.0
}
;
double c[A] =  {
0.0, 1.0, 2.0, 4.5, 6.0, 10.0
}
;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("EMPIRICAL_DISCRETE_SAMPLE_TEST");
Console.WriteLine("  EMPIRICAL_DISCRETE_MEAN computes the Empirical Discrete mean;");
Console.WriteLine("  EMPIRICAL_DISCRETE_SAMPLE samples the Empirical Discrete distribution;");
Console.WriteLine("  EMPIRICAL_DISCRETE_VARIANCE computes the Empirical Discrete variance.");

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + A + "");
r8vec_print(A, b, "  PDF parameter B = ");
r8vec_print(A, c, "  PDF parameter C = ");

if (!empirical_discrete_check(A, b, c))
{
Console.WriteLine("");
Console.WriteLine("EMPIRICAL_DISCRETE_SAMPLE_TEST- Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = empirical_discrete_mean(A, b, c);
variance = empirical_discrete_variance(A, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = empirical_discrete_sample(A, b, c, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef A
# undef SAMPLE_NUM
}

static void english_letter_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_LETTER_CDF_TEST tests ENGLISH_LETTER_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2016
//
//  Author:
//
//    John Burkardt
//
{
char c;
char c2;
double cdf;
int i;
double pdf;
int seed;

Console.WriteLine("");
Console.WriteLine("ENGLISH_LETTER_CDF_TEST");
Console.WriteLine("  ENGLISH_LETTER_CDF evaluates the English Letter CDF;");
Console.WriteLine("  ENGLISH_LETTER_CDF_INV inverts the English Letter CDF.");
Console.WriteLine("  ENGLISH_LETTER_PDF evaluates the English Letter PDF;");

seed = 123456789;

Console.WriteLine("");
Console.WriteLine("   C              PDF             CDF    CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
c = english_letter_sample(seed);
pdf = english_letter_pdf(c);
cdf = english_letter_cdf(c);
c2 = english_letter_cdf_inv(cdf);

Console.WriteLine("  '" + c + "'"
+ "  " + setw(14) + pdf
+ "  " + setw(14) + cdf
+ "  '" + c2 + "'");
}

return;
}

static void english_sentence_length_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_CDF_TEST tests ENGLISH_SENTENCE_LENGTH_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
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
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("ENGLISH_SENTENCE_LENGTH_CDF_TEST");
Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_CDF evaluates the English Sentence Length CDF;");
Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the English Sentence Length CDF.");
Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_PDF evaluates the English Sentence Length PDF;");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = english_sentence_length_sample(seed);

pdf = english_sentence_length_pdf(x);

cdf = english_sentence_length_cdf(x);

x2 = english_sentence_length_cdf_inv(cdf);

Console.WriteLine("  " + setw(12) + x
+ "  " + setw(12) + pdf
+ "  " + setw(12) + cdf
+ "  " + setw(12) + x2 + "");
}

return;
}

static void english_sentence_length_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_SAMPLE_TEST tests ENGLISH_SENTENCE_LENGTH_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

int i;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("ENGLISH_SENTENCE_LENGTH_SAMPLE_TEST");
Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_MEAN computes the English Sentence Length mean;");
Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_SAMPLE samples the English Sentence Length distribution;");
Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_VARIANCE computes the English Sentence Length variance.");

mean = english_sentence_length_mean();
variance = english_sentence_length_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = english_sentence_length_sample(seed);
}

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void english_word_length_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_CDF_TEST tests ENGLISH_WORD_LENGTH_CDF.
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
double cdf;
int i;
double pdf;
int seed = 123456789;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("ENGLISH_WORD_LENGTH_CDF_TEST");
Console.WriteLine("  ENGLISH_WORD_LENGTH_CDF evaluates the English Word LengthCDF;");
Console.WriteLine("  ENGLISH_WORD_LENGTH_CDF_INV inverts the English Word LengthCDF.");
Console.WriteLine("  ENGLISH_WORD_LENGTH_PDF evaluates the English Word LengthPDF;");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = english_word_length_sample(seed);

pdf = english_word_length_pdf(x);

cdf = english_word_length_cdf(x);

x2 = english_word_length_cdf_inv(cdf);

Console.WriteLine("  " + setw(12) + x
+ "  " + setw(12) + pdf
+ "  " + setw(12) + cdf
+ "  " + setw(12) + x2 + "");
}

return;
}

static void english_word_length_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_SAMPLE_TEST tests ENGLISH_WORD_LENGTH_SAMPLE.
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
# define SAMPLE_NUM 1000

int i;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("ENGLISH_WORD_LENGTH_SAMPLE_TEST");
Console.WriteLine("  ENGLISH_WORD_LENGTH_MEAN computes the English Word Lengthmean;");
Console.WriteLine("  ENGLISH_WORD_LENGTH_SAMPLE samples the English Word Lengthdistribution;");
Console.WriteLine("  ENGLISH_WORD_LENGTH_VARIANCE computes the English Word Lengthvariance.");

mean = english_word_length_mean();
variance = english_word_length_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = english_word_length_sample(seed);
}

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void erlang_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_CDF_TEST tests ERLANG_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
int c;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("ERLANG_CDF_TEST");
Console.WriteLine("  ERLANG_CDF evaluates the Erlang CDF;");
Console.WriteLine("  ERLANG_CDF_INV inverts the Erlang CDF.");
Console.WriteLine("  ERLANG_PDF evaluates the Erlang PDF;");

a = 1.0;
b = 2.0;
c = 3;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!erlang_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("ERLANG_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = erlang_sample(a, b, c, seed);
pdf = erlang_pdf(x, a, b, c);
cdf = erlang_cdf(x, a, b, c);
x2 = erlang_cdf_inv(cdf, a, b, c);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void erlang_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_SAMPLE_TEST tests ERLANG_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int c;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("ERLANG_SAMPLE_TEST");
Console.WriteLine("  ERLANG_MEAN computes the Erlang mean;");
Console.WriteLine("  ERLANG_SAMPLE samples the Erlang distribution;");
Console.WriteLine("  ERLANG_VARIANCE computes the Erlang variance;");

a = 1.0;
b = 2.0;
c = 3;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!erlang_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("ERLANG_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = erlang_mean(a, b, c);
variance = erlang_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = erlang_sample(a, b, c, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void exponential_01_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_CDF_TEST tests EXPONENTIAL_01_CDF.
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
Console.WriteLine("EXPONENTIAL_01_CDF_TEST");
Console.WriteLine("  EXPONENTIAL_01_CDF evaluates the Exponential 01 CDF;");
Console.WriteLine("  EXPONENTIAL_01_CDF_INV inverts the Exponential 01 CDF.");
Console.WriteLine("  EXPONENTIAL_01_PDF evaluates the Exponential 01 PDF;");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = exponential_01_sample(seed);
pdf = exponential_01_pdf(x);
cdf = exponential_01_cdf(x);
x2 = exponential_01_cdf_inv(cdf);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void exponential_01_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_SAMPLE_TEST tests EXPONENTIAL_01_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2016
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("EXPONENTIAL_01_SAMPLE_TEST");
Console.WriteLine("  EXPONENTIAL_01_MEAN computes the Exponential 01 mean;");
Console.WriteLine("  EXPONENTIAL_01_SAMPLE samples the Exponential 01 distribution;");
Console.WriteLine("  EXPONENTIAL_01_VARIANCE computes the Exponential 01 variance.");

mean = exponential_01_mean();
variance = exponential_01_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = exponential_01_sample(seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void exponential_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF_TEST tests EXPONENTIAL_CDF.
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
Console.WriteLine("EXPONENTIAL_CDF_TEST");
Console.WriteLine("  EXPONENTIAL_CDF evaluates the Exponential CDF;");
Console.WriteLine("  EXPONENTIAL_CDF_INV inverts the Exponential CDF.");
Console.WriteLine("  EXPONENTIAL_PDF evaluates the Exponential PDF;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!exponential_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("EXPONENTIAL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = exponential_sample(a, b, seed);
pdf = exponential_pdf(x, a, b);
cdf = exponential_cdf(x, a, b);
x2 = exponential_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void exponential_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_SAMPLE_TEST tests EXPONENTIAL_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("EXPONENTIAL_SAMPLE_TEST");
Console.WriteLine("  EXPONENTIAL_MEAN computes the Exponential mean;");
Console.WriteLine("  EXPONENTIAL_SAMPLE samples the Exponential distribution;");
Console.WriteLine("  EXPONENTIAL_VARIANCE computes the Exponential variance;");

a = 1.0;
b = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!exponential_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("EXPONENTIAL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = exponential_mean(a, b);
variance = exponential_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = exponential_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void extreme_values_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CDF_TEST tests EXTREME_VALUES_CDF.
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
Console.WriteLine("EXTREME_VALUES_CDF_TEST");
Console.WriteLine("  EXTREME_VALUES_CDF evaluates the Extreme Values CDF;");
Console.WriteLine("  EXTREME_VALUES_CDF_INV inverts the Extreme Values CDF.");
Console.WriteLine("  EXTREME_VALUES_PDF evaluates the Extreme Values PDF;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!extreme_values_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("EXTREME_VALUES_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = extreme_values_sample(a, b, seed);
pdf = extreme_values_pdf(x, a, b);
cdf = extreme_values_cdf(x, a, b);
x2 = extreme_values_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void extreme_values_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_SAMPLE_TEST tests EXTREME_VALUES_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("EXTREME_VALUES_SAMPLE_TEST");
Console.WriteLine("  EXTREME_VALUES_MEAN computes the Extreme Values mean;");
Console.WriteLine("  EXTREME_VALUES_SAMPLE samples the Extreme Values distribution;");
Console.WriteLine("  EXTREME_VALUES_VARIANCE computes the Extreme Values variance;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!extreme_values_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("EXTREME_VALUES_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = extreme_values_mean(a, b);
variance = extreme_values_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = extreme_values_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void f_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    F_CDF_TEST tests F_CDF.
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
int m;
int n;
double pdf;
int seed = 123456789;
double x;

Console.WriteLine("");
Console.WriteLine("F_CDF_TEST");
Console.WriteLine("  F_CDF evaluates the F CDF;");
Console.WriteLine("  F_PDF evaluates the F PDF;");
Console.WriteLine("  F_SAMPLE samples the F PDF;");

m = 1;
n = 1;

Console.WriteLine("");
Console.WriteLine("  PDF parameter M = " + m + "");
Console.WriteLine("  PDF parameter N = " + n + "");

if (!f_check(m, n))
{
Console.WriteLine("");
Console.WriteLine("F_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameter values are illegal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = f_sample(m, n, seed);
pdf = f_pdf(x, m, n);
cdf = f_cdf(x, m, n);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "");
}

return;
}

static void f_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    F_SAMPLE_TEST tests F_SAMPLE.
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
# define SAMPLE_NUM 1000

int i;
int m;
double mean;
int n;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("F_SAMPLE_TEST");
Console.WriteLine("  F_MEAN computes the F mean;");
Console.WriteLine("  F_SAMPLE samples the F distribution;");
Console.WriteLine("  F_VARIANCE computes the F variance;");

m = 8;
n = 6;

Console.WriteLine("");
Console.WriteLine("  PDF parameter M = " + m + "");
Console.WriteLine("  PDF parameter N = " + n + "");

if (!f_check(m, n))
{
Console.WriteLine("");
Console.WriteLine("F_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = f_mean(m, n);
variance = f_variance(m, n);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = f_sample(m, n, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void fermi_dirac_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    FERMI_DIRAC_SAMPLE_TEST tests FERMI_DIRAC_SAMPLE.
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
double mean;
int sample_num = 10000;
int seed;
int test;
double u;
double u_test[7] =  {
1.0, 2.0, 4.0, 8.0, 16.0,
32.0, 1.0
}
;
double v;
double v_test[7] =  {
1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 0.25
}
;
double variance;
double x[10000];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("FERMI_DIRAC_SAMPLE_TEST");
Console.WriteLine("  FERMI_DIRAC_SAMPLE samples the Fermi Dirac distribution.");

for (test = 0; test < 7; test++)
{
u = u_test[test];
v = v_test[test];
seed = 123456789;

Console.WriteLine("");
Console.WriteLine("  U =          " + u + "");
Console.WriteLine("  V =          " + v + "");

for (i = 0; i < sample_num; i++)
{
x[i] = fermi_dirac_sample(u, v, seed);
}

mean = r8vec_mean(sample_num, x);
variance = r8vec_variance(sample_num, x);
xmax = r8vec_max(sample_num, x);
xmin = r8vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + mean + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");
}

return;
}

static void fisher_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    FISHER_PDF_TEST tests FISHER_PDF.
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
int j;
double kappa;
double mu[3];
int n = 10;
double pdf;
int seed;
int test;
int test_num = 3;
double* x;

Console.WriteLine("");
Console.WriteLine("FISHER_SAMPLE_TEST");
Console.WriteLine("  FISHER_PDF evaluates the Fisher PDF.");

for (test = 1; test <= test_num; test++)
{
if (test == 1)
{
kappa = 0.0;
mu[0] = 1.0;
mu[1] = 0.0;
mu[2] = 0.0;
}
else if (test == 2)
{
kappa = 0.5;
mu[0] = 1.0;
mu[1] = 0.0;
mu[2] = 0.0;
}
else if (test == 3)
{
kappa = 10.0;
mu[0] = 1.0;
mu[1] = 0.0;
mu[2] = 0.0;
}

Console.WriteLine("");
Console.WriteLine("  PDF parameters:");
Console.WriteLine("    Concentration parameter KAPPA = " + kappa + "");
Console.WriteLine("    Direction MU(1:3) = "
+ "  " + mu[0]
+ "  " + mu[1]
+ "  " + mu[2] + "");

Console.WriteLine("");
Console.WriteLine("      X                         PDF");
Console.WriteLine("");

seed = 123456789;

for (j = 0; j < n; j++)
{
x = fisher_sample(kappa, mu, 1, seed);

pdf = fisher_pdf(x, kappa, mu);

Console.WriteLine("  " + setw(10) + x[0]
  + "  " + setw(10) + x[1]
  + "  " + setw(10) + x[2]
  + "  " + setw(14) + pdf + "");

delete[] x;
}

}

return;
}

static void fisk_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    FISK_CDF_TEST tests FISK_CDF.
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
double b;
double c;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("FISK_CDF_TEST");
Console.WriteLine("  FISK_CDF evaluates the Fisk CDF;");
Console.WriteLine("  FISK_CDF_INV inverts the Fisk CDF.");
Console.WriteLine("  FISK_PDF evaluates the Fisk PDF;");

a = 1.0;
b = 2.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!fisk_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("FISK_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = fisk_sample(a, b, c, seed);
pdf = fisk_pdf(x, a, b, c);
cdf = fisk_cdf(x, a, b, c);
x2 = fisk_cdf_inv(cdf, a, b, c);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void fisk_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    FISK_SAMPLE_TEST tests FISK_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
double c;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("FISK_SAMPLE_TEST");
Console.WriteLine("  FISK_MEAN computes the Fisk mean;");
Console.WriteLine("  FISK_SAMPLE samples the Fisk distribution;");
Console.WriteLine("  FISK_VARIANCE computes the Fisk variance;");

a = 1.0;
b = 2.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!fisk_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("FISK_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = fisk_mean(a, b, c);
variance = fisk_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = fisk_sample(a, b, c, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void folded_normal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_CDF_TEST tests FOLDED_NORMAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
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
Console.WriteLine("FOLDED_NORMAL_CDF_TEST");
Console.WriteLine("  FOLDED_NORMAL_CDF evaluates the Folded Normal CDF;");
Console.WriteLine("  FOLDED_NORMAL_CDF_INV inverts the Folded Normal CDF.");
Console.WriteLine("  FOLDED_NORMAL_PDF evaluates the Folded Normal PDF;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!folded_normal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("FOLDED_NORMAL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = folded_normal_sample(a, b, seed);
pdf = folded_normal_pdf(x, a, b);
cdf = folded_normal_cdf(x, a, b);
x2 = folded_normal_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void folded_normal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_SAMPLE_TEST tests FOLDED_NORMAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("FOLDED_NORMAL_SAMPLE_TEST");
Console.WriteLine("  FOLDED_NORMAL_MEAN computes the Folded Normal mean;");
Console.WriteLine("  FOLDED_NORMAL_SAMPLE samples the Folded Normal distribution;");
Console.WriteLine("  FOLDED_NORMAL_VARIANCE computes the Folded Normal variance;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!folded_normal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("FOLDED_NORMAL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = folded_normal_mean(a, b);
variance = folded_normal_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = folded_normal_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void frechet_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_CDF_TEST tests FRECHET_CDF.
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
double alpha;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("FRECHET_CDF_TEST");
Console.WriteLine("  FRECHET_CDF evaluates the Frechet CDF;");
Console.WriteLine("  FRECHET_CDF_INV inverts the Frechet CDF.");
Console.WriteLine("  FRECHET_PDF evaluates the Frechet PDF;");

alpha = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter ALPHA =  " + alpha + "");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = frechet_sample(alpha, seed);
pdf = frechet_pdf(x, alpha);
cdf = frechet_cdf(x, alpha);
x2 = frechet_cdf_inv(cdf, alpha);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void frechet_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_SAMPLE_TEST tests FRECHET_SAMPLE.
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
# define SAMPLE_NUM 1000

double alpha;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("FRECHET_SAMPLE_TEST");
Console.WriteLine("  FRECHET_MEAN computes the Frechet mean;");
Console.WriteLine("  FRECHET_SAMPLE samples the Frechet distribution;");
Console.WriteLine("  FRECHET_VARIANCE computes the Frechet variance;");

alpha = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter ALPHA =  " + alpha + "");

mean = frechet_mean(alpha);
variance = frechet_variance(alpha);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = frechet_sample(alpha, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void gamma_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_CDF_TEST tests GAMMA_CDF.
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
double b;
double c;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;

Console.WriteLine("");
Console.WriteLine("GAMMA_CDF_TEST");
Console.WriteLine("  GAMMA_CDF evaluates the Gamma CDF;");
Console.WriteLine("  GAMMA_PDF evaluates the Gamma PDF;");
Console.WriteLine("  GAMMA_SAMPLE samples the Gamma PDF;");

a = 1.0;
b = 1.5;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A = " + a + "");
Console.WriteLine("  PDF parameter B = " + b + "");
Console.WriteLine("  PDF parameter B = " + c + "");

if (!gamma_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("GAMMA_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameter values are illegal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = gamma_sample(a, b, c, seed);
pdf = gamma_pdf(x, a, b, c);
cdf = gamma_cdf(x, a, b, c);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "");
}

return;
}

static void gamma_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_SAMPLE_TEST tests GAMMA_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
double c;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("GAMMA_SAMPLE_TEST");
Console.WriteLine("  GAMMA_MEAN computes the Gamma mean;");
Console.WriteLine("  GAMMA_SAMPLE samples the Gamma distribution;");
Console.WriteLine("  GAMMA_VARIANCE computes the Gamma variance;");

a = 1.0;
b = 3.0;
c = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!gamma_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("GAMMA_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = gamma_mean(a, b, c);
variance = gamma_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = gamma_sample(a, b, c, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void genlogistic_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_CDF_TEST tests GENLOGISTIC_CDF.
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
Console.WriteLine("GENLOGISTIC_CDF_TEST");
Console.WriteLine("  GENLOGISTIC_CDF evaluates the Genlogistic CDF;");
Console.WriteLine("  GENLOGISTIC_CDF_INV inverts the Genlogistic CDF.");
Console.WriteLine("  GENLOGISTIC_PDF evaluates the Genlogistic PDF;");

a = 1.0;
b = 2.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!genlogistic_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("GENLOGISTIC_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = genlogistic_sample(a, b, c, seed);
pdf = genlogistic_pdf(x, a, b, c);
cdf = genlogistic_cdf(x, a, b, c);
x2 = genlogistic_cdf_inv(cdf, a, b, c);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void genlogistic_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_SAMPLE_TEST tests GENLOGISTIC_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
double c;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("GENLOGISTIC_SAMPLE_TEST");
Console.WriteLine("  GENLOGISTIC_MEAN computes the Genlogistic mean;");
Console.WriteLine("  GENLOGISTIC_SAMPLE samples the Genlogistic distribution;");
Console.WriteLine("  GENLOGISTIC_VARIANCE computes the Genlogistic variance;");

a = 1.0;
b = 2.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!genlogistic_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("GENLOGISTIC_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = genlogistic_mean(a, b, c);
variance = genlogistic_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = genlogistic_sample(a, b, c, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void geometric_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CDF_TEST tests GEOMETRIC_CDF.
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
double cdf;
int i;
double pdf;
int seed = 123456789;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("GEOMETRIC_CDF_TEST");
Console.WriteLine("  GEOMETRIC_CDF evaluates the Geometric CDF;");
Console.WriteLine("  GEOMETRIC_CDF_INV inverts the Geometric CDF.");
Console.WriteLine("  GEOMETRIC_PDF evaluates the Geometric PDF;");

a = 0.25E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!geometric_check(a))
{
Console.WriteLine("");
Console.WriteLine("GEOMETRIC_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = geometric_sample(a, seed);
pdf = geometric_pdf(x, a);
cdf = geometric_cdf(x, a);
x2 = geometric_cdf_inv(cdf, a);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void geometric_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_SAMPLE_TEST tests GEOMETRIC_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
int j;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("GEOMETRIC_SAMPLE_TEST");
Console.WriteLine("  GEOMETRIC_MEAN computes the Geometric mean;");
Console.WriteLine("  GEOMETRIC_SAMPLE samples the Geometric distribution;");
Console.WriteLine("  GEOMETRIC_VARIANCE computes the Geometric variance.");

a = 0.25E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!geometric_check(a))
{
Console.WriteLine("");
Console.WriteLine("GEOMETRIC_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = geometric_mean(a);
variance = geometric_variance(a);

Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (j = 0; j < SAMPLE_NUM; j++)
{
x[j] = geometric_sample(a, seed);
}

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void gompertz_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_CDF_TEST tests GOMPERTZ_CDF
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
double a;
double b;
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("GOMPERTZ_CDF_TEST");
Console.WriteLine("  GOMPERTZ_CDF evaluates the Gompertz CDF;");
Console.WriteLine("  GOMPERTZ_CDF_INV inverts the Gompertz CDF.");
Console.WriteLine("  GOMPERTZ_PDF evaluates the Gompertz PDF;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!gompertz_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("GOMPERTZ_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = gompertz_sample(a, b, seed);
pdf = gompertz_pdf(x, a, b);
cdf = gompertz_cdf(x, a, b);
x2 = gompertz_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void gompertz_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_SAMPLE_TEST tests GOMPERTZ_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("GOMPERTZ_SAMPLE_TEST");
Console.WriteLine("  GOMPERTZ_MEAN computes the Gompertz mean;");
Console.WriteLine("  GOMPERTZ_SAMPLE samples the Gompertz distribution;");
Console.WriteLine("  GOMPERTZ_VARIANCE computes the Gompertz variance;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!gompertz_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("GOMPERTZ_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = gompertz_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void gumbel_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_CDF_TEST tests GUMBEL_CDF.
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
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("GUMBEL_CDF_TEST");
Console.WriteLine("  GUMBEL_CDF evaluates the Gumbel CDF;");
Console.WriteLine("  GUMBEL_CDF_INV inverts the Gumbel CDF.");
Console.WriteLine("  GUMBEL_PDF evaluates the Gumbel PDF;");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = gumbel_sample(seed);
pdf = gumbel_pdf(x);
cdf = gumbel_cdf(x);
x2 = gumbel_cdf_inv(cdf);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void gumbel_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_SAMPLE_TEST tests GUMBEL_SAMPLE.
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
# define SAMPLE_NUM 1000

int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("GUMBEL_SAMPLE_TEST");
Console.WriteLine("  GUMBEL_MEAN computes the Gumbel mean;");
Console.WriteLine("  GUMBEL_SAMPLE samples the Gumbel distribution;");
Console.WriteLine("  GUMBEL_VARIANCE computes the Gumbel variance.");

mean = gumbel_mean();
variance = gumbel_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = gumbel_sample(seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void half_normal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_CDF_TEST tests HALF_NORMAL_CDF.
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
Console.WriteLine("HALF_NORMAL_CDF_TEST");
Console.WriteLine("  HALF_NORMAL_CDF evaluates the Half Normal CDF;");
Console.WriteLine("  HALF_NORMAL_CDF_INV inverts the Half Normal CDF.");
Console.WriteLine("  HALF_NORMAL_PDF evaluates the Half Normal PDF;");

a = 0.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!half_normal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("HALF_NORMAL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = half_normal_sample(a, b, seed);
pdf = half_normal_pdf(x, a, b);
cdf = half_normal_cdf(x, a, b);
x2 = half_normal_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void half_normal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_SAMPLE_TEST tests HALF_NORMAL_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("HALF_NORMAL_SAMPLE_TEST");
Console.WriteLine("  HALF_NORMAL_MEAN computes the Half Normal mean;");
Console.WriteLine("  HALF_NORMAL_SAMPLE samples the Half Normal distribution;");
Console.WriteLine("  HALF_NORMAL_VARIANCE computes the Half Normal variance;");

a = 0.0;
b = 10.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!half_normal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("HALF_NORMAL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = half_normal_mean(a, b);
variance = half_normal_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = half_normal_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void hypergeometric_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_TEST tests HYPERGEOMETRIC_CDF.
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
double cdf;
int l;
int m;
int n;
double pdf;
int x;

Console.WriteLine("");
Console.WriteLine("HYPERGEOMETRIC_CDF_TEST");
Console.WriteLine("  HYPERGEOMETRIC_CDF evaluates the Hypergeometric CDF.");
Console.WriteLine("  HYPERGEOMETRIC_PDF evaluates the Hypergeometric PDF.");

n = 10;
m = 7;
l = 100;

Console.WriteLine("");
Console.WriteLine("  Total number of balls L =         " + l + "");
Console.WriteLine("  Number of white balls M =         " + m + "");
Console.WriteLine("  Number of balls taken N =         " + n + "");

if (!hypergeometric_check(n, m, l))
{
Console.WriteLine("");
Console.WriteLine("HYPERGEOMETRIC_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

x = 7;

pdf = hypergeometric_pdf(x, n, m, l);

cdf = hypergeometric_cdf(x, n, m, l);

Console.WriteLine("  PDF argument X =                " + x + "");
Console.WriteLine("  PDF value =                   = " + pdf + "");
Console.WriteLine("  CDF value =                   = " + cdf + "");

return;
}

static void hypergeometric_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_SAMPLE_TEST tests HYPERGEOMETRIC_SAMPLE.
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
# define SAMPLE_NUM 1000

int j;
int l;
int m;
double mean;
int n;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("HYPERGEOMETRIC_SAMPLE_TEST");
Console.WriteLine("  HYPERGEOMETRIC_MEAN computes the Hypergeometric mean;");
Console.WriteLine("  HYPERGEOMETRIC_SAMPLE samples the Hypergeometric distribution;");
Console.WriteLine("  HYPERGEOMETRIC_VARIANCE computes the Hypergeometric variance.");

n = 10;
m = 7;
l = 100;

Console.WriteLine("");
Console.WriteLine("  Total number of balls L =         " + l + "");
Console.WriteLine("  Number of white balls M =         " + m + "");
Console.WriteLine("  Number of balls taken N =         " + n + "");

if (!hypergeometric_check(n, m, l))
{
Console.WriteLine("");
Console.WriteLine("HYPERGEOMETRIC_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = hypergeometric_mean(n, m, l);
variance = hypergeometric_variance(n, m, l);

Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (j = 0; j < SAMPLE_NUM; j++)
{
x[j] = hypergeometric_sample(n, m, l, seed);
}

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void i4_choose_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE_TEST tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
int cnk;
int k;
int n;

Console.WriteLine("");
Console.WriteLine("I4_CHOOSE_TEST");
Console.WriteLine("  I4_CHOOSE evaluates C(N,K).");
Console.WriteLine("");
Console.WriteLine("       N       K     CNK");

for (n = 0; n <= 4; n++)
{
Console.WriteLine("");
for (k = 0; k <= n; k++)
{
cnk = i4_choose(n, k);

cout + "  "
+ setw(6) + n + "  "
+ setw(6) + k + "  "
+ setw(6) + cnk + "");
}
}

return;
}

static void i4_choose_log_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE_LOG_TEST tests I4_CHOOSE_LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2016
//
//  Author:
//
//    John Burkardt
//
{
int cnk;
double elcnk;
int k;
double lcnk;
int n;

Console.WriteLine("");
Console.WriteLine("I4_CHOOSE_LOG_TEST");
Console.WriteLine("  I4_CHOOSE_LOG evaluates log(C(N,K)).");
Console.WriteLine("");
Console.WriteLine("       N       K            lCNK           elCNK     CNK");

for (n = 0; n <= 4; n++)
{
Console.WriteLine("");
for (k = 0; k <= n; k++)
{
lcnk = i4_choose_log(n, k);
elcnk = exp(lcnk);
cnk = i4_choose(n, k);

Console.WriteLine("  " + setw(6) + n
  + "  " + setw(6) + k
  + "  " + setw(14) + lcnk
  + "  " + setw(14) + elcnk
  + "  " + setw(6) + cnk + "");
}
}

return;
}

static void i4_is_power_of_10_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_IS_POWER_OF_10_TEST tests I4_IS_POWER_OF_10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2016
//
//  Author:
//
//    John Burkardt
//
{
int i;

Console.WriteLine("");
Console.WriteLine("I4_IS_POWER_OF_10_TEST");
Console.WriteLine("  I4_IS_POWER_OF_10 reports whether an I4 is a power of 10.");
Console.WriteLine("");
Console.WriteLine("  I     I4_IS_POWER_OF_10(I)");
Console.WriteLine("");

for (i = 97; i <= 103; i++)
{
Console.WriteLine("  " + setw(6) + i
+ "  " + i4_is_power_of_10(i) + "");
}

return;
}

static void i4_uniform_ab_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
int a = -100;
int b = 200;
int i;
int j;
int seed = 123456789;

Console.WriteLine("");
Console.WriteLine("I4_UNIFORM_AB_TEST");
Console.WriteLine("  I4_UNIFORM_AB computes pseudorandom values");
Console.WriteLine("  in an interval [A,B].");

Console.WriteLine("");
Console.WriteLine("  The lower endpoint A = " + a + "");
Console.WriteLine("  The upper endpoint B = " + b + "");
Console.WriteLine("  The initial seed is " + seed + "");
Console.WriteLine("");

for (i = 1; i <= 20; i++)
{
j = i4_uniform_ab(a, b, seed);

Console.WriteLine("  " + setw(8) + i
+ "  " + setw(8) + j + "");
}

return;
}

static void i4vec_uniform_ab_new_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM_AB_NEW_TEST tests I4VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
int a = -100;
int b = 200;
int n = 20;
int seed = 123456789;
int* v;

Console.WriteLine("");
Console.WriteLine("I4VEC_UNIFORM_AB_NEW_TEST");
Console.WriteLine("  I4VEC_UNIFORM_AB_NEW computes pseudorandom values");
Console.WriteLine("  in an interval [A,B].");

Console.WriteLine("");
Console.WriteLine("  The lower endpoint A = " + a + "");
Console.WriteLine("  The upper endpoint B = " + b + "");
Console.WriteLine("  The initial seed is " + seed + "");
Console.WriteLine("");

v = i4vec_uniform_ab_new(n, a, b, seed);

i4vec_print(n, v, "  The random vector:");

delete[] v;

return;
}

static void i4vec_unique_count_test()

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIQUE_COUNT_TEST tests I4VEC_UNIQUE_COUNT.
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
int* a;
int a_hi;
int a_lo;
int a_unique;
int n;
int seed;

Console.WriteLine("");
Console.WriteLine("I4VEC_UNIQUE_COUNT_TEST");
Console.WriteLine("  I4VEC_UNIQUE_COUNT counts unique entries in an I4VEC.");

n = 20;
a_lo = 0;
a_hi = n;
seed = 123456789;

a = i4vec_uniform_ab_new(n, a_lo, a_hi, seed);

i4vec_print(n, a, "  Array:");

a_unique = i4vec_unique_count(n, a);

Console.WriteLine("");
Console.WriteLine("  Number of unique entries is " + a_unique + "");

delete[] a;

return;
}

static void inverse_gaussian_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_CDF_TEST tests INVERSE_GAUSSIAN_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2016
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

Console.WriteLine("");
Console.WriteLine("INVERSE_GAUSSIAN_CDF_TEST");
Console.WriteLine("  INVERSE_GAUSSIAN_CDF evaluates the Inverse Gaussian CDF;");
Console.WriteLine("  INVERSE_GAUSSIAN_PDF evaluates the Inverse Gaussian PDF;");

a = 5.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!inverse_gaussian_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("INVERSE_GAUSSIAN_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = inverse_gaussian_sample(a, b, seed);
pdf = inverse_gaussian_pdf(x, a, b);
cdf = inverse_gaussian_cdf(x, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "");
}

return;
}

static void inverse_gaussian_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_SAMPLE_TEST tests INVERSE_GAUSSIAN_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("INVERSE_GAUSSIAN_SAMPLE_TEST");
Console.WriteLine("  INVERSE_GAUSSIAN_MEAN computes the Inverse Gaussian mean;");
Console.WriteLine("  INVERSE_GAUSSIAN_SAMPLE samples the Inverse Gaussian distribution;");
Console.WriteLine("  INVERSE_GAUSSIAN_VARIANCE computes the Inverse Gaussian variance;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!inverse_gaussian_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("INVERSE_GAUSSIAN_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = inverse_gaussian_mean(a, b);
variance = inverse_gaussian_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = inverse_gaussian_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void laplace_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CDF_TEST tests LAPLACE_CDF.
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
Console.WriteLine("LAPLACE_CDF_TEST");
Console.WriteLine("  LAPLACE_CDF evaluates the Laplace CDF;");
Console.WriteLine("  LAPLACE_CDF_INV inverts the Laplace CDF.");
Console.WriteLine("  LAPLACE_PDF evaluates the Laplace PDF;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!laplace_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LAPLACE_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = laplace_sample(a, b, seed);
pdf = laplace_pdf(x, a, b);
cdf = laplace_cdf(x, a, b);
x2 = laplace_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void laplace_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_SAMPLE_TEST tests LAPLACE_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("LAPLACE_SAMPLE_TEST");
Console.WriteLine("  LAPLACE_MEAN computes the Laplace mean;");
Console.WriteLine("  LAPLACE_SAMPLE samples the Laplace distribution;");
Console.WriteLine("  LAPLACE_VARIANCE computes the Laplace variance;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!laplace_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LAPLACE_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = laplace_mean(a, b);
variance = laplace_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = laplace_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void levy_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LEVY_CDF_TEST tests LEVY_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2016
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
Console.WriteLine("LEVY_CDF_TEST");
Console.WriteLine("  LEVY_CDF evaluates the Levy CDF;");
Console.WriteLine("  LEVY_CDF_INV inverts the Levy CDF.");
Console.WriteLine("  LEVY_PDF evaluates the Levy PDF;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = levy_sample(a, b, seed);
pdf = levy_pdf(x, a, b);
cdf = levy_cdf(x, a, b);
x2 = levy_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void logistic_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CDF_TEST tests LOGISTIC_CDF.
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
Console.WriteLine("LOGISTIC_CDF_TEST");
Console.WriteLine("  LOGISTIC_CDF evaluates the Logistic CDF;");
Console.WriteLine("  LOGISTIC_CDF_INV inverts the Logistic CDF.");
Console.WriteLine("  LOGISTIC_PDF evaluates the Logistic PDF;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!logistic_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LOGISTIC_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = logistic_sample(a, b, seed);
pdf = logistic_pdf(x, a, b);
cdf = logistic_cdf(x, a, b);
x2 = logistic_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void logistic_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_SAMPLE_TEST tests LOGISTIC_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("LOGISTIC_SAMPLE_TEST");
Console.WriteLine("  LOGISTIC_MEAN computes the Logistic mean;");
Console.WriteLine("  LOGISTIC_SAMPLE samples the Logistic distribution;");
Console.WriteLine("  LOGISTIC_VARIANCE computes the Logistic variance;");

a = 2.0;
b = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!logistic_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LOGISTIC_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = logistic_mean(a, b);
variance = logistic_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = logistic_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void log_normal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF_TEST tests LOG_NORMAL_CDF.
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
Console.WriteLine("LOG_NORMAL_CDF_TEST");
Console.WriteLine("  LOG_NORMAL_CDF evaluates the Log Normal CDF;");
Console.WriteLine("  LOG_NORMAL_CDF_INV inverts the Log Normal CDF.");
Console.WriteLine("  LOG_NORMAL_PDF evaluates the Log Normal PDF;");

a = 10.0;
b = 2.25;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!log_normal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LOG_NORMAL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = log_normal_sample(a, b, seed);
pdf = log_normal_pdf(x, a, b);
cdf = log_normal_cdf(x, a, b);
x2 = log_normal_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void log_normal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_SAMPLE_TEST tests LOG_NORMAL_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("LOG_NORMAL_SAMPLE_TEST");
Console.WriteLine("  LOG_NORMAL_MEAN computes the Log Normal mean;");
Console.WriteLine("  LOG_NORMAL_SAMPLE samples the Log Normal distribution;");
Console.WriteLine("  LOG_NORMAL_VARIANCE computes the Log Normal variance;");

a = 1.0;
b = 2.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!normal_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LOG_NORMAL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = log_normal_mean(a, b);
variance = log_normal_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = log_normal_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void log_series_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CDF_TEST tests LOG_SERIES_CDF.
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
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("LOG_SERIES_CDF_TEST");
Console.WriteLine("  LOG_SERIES_CDF evaluates the Log Series CDF;");
Console.WriteLine("  LOG_SERIES_CDF_INV inverts the Log Series CDF.");
Console.WriteLine("  LOG_SERIES_PDF evaluates the Log Series PDF;");

a = 0.25E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!log_series_check(a))
{
Console.WriteLine("");
Console.WriteLine("LOG_SERIES_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = log_series_sample(a, seed);
pdf = log_series_pdf(x, a);
cdf = log_series_cdf(x, a);
x2 = log_series_cdf_inv(cdf, a);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void log_series_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_SAMPLE_TEST tests LOG_SERIES_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
int j;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("LOG_SERIES_SAMPLE_TEST");
Console.WriteLine("  LOG_SERIES_MEAN computes the Log Series mean;");
Console.WriteLine("  LOG_SERIES_SAMPLE samples the Log Series distribution;");
Console.WriteLine("  LOG_SERIES_VARIANCE computes the Log Series variance.");

a = 0.25E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!log_series_check(a))
{
Console.WriteLine("");
Console.WriteLine("LOG_SERIES_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = log_series_mean(a);
variance = log_series_variance(a);

Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (j = 0; j < SAMPLE_NUM; j++)
{
x[j] = log_series_sample(a, seed);
}

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void log_uniform_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_CDF_TEST tests LOG_UNIFORM_CDF.
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
Console.WriteLine("LOG_UNIFORM_CDF_TEST");
Console.WriteLine("  LOG_UNIFORM_CDF evaluates the Log Uniform CDF;");
Console.WriteLine("  LOG_UNIFORM_CDF_INV inverts the Log Uniform CDF.");
Console.WriteLine("  LOG_UNIFORM_PDF evaluates the Log Uniform PDF;");

a = 2.0;
b = 20.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!log_uniform_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LOG_UNIFORM_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = log_uniform_sample(a, b, seed);
pdf = log_uniform_pdf(x, a, b);
cdf = log_uniform_cdf(x, a, b);
x2 = log_uniform_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void log_uniform_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_SAMPLE_TEST tests LOG_UNIFORM_SAMPLE;
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("LOG_UNIFORM_SAMPLE_TEST");
Console.WriteLine("  LOG_UNIFORM_MEAN computes the Log Uniform mean;");
Console.WriteLine("  LOG_UNIFORM_SAMPLE samples the Log Uniform distribution;");

a = 2.0;
b = 20.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!log_uniform_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("LOG_UNIFORM_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = log_uniform_mean(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = log_uniform_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void lorentz_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_CDF_TEST tests LORENTZ_CDF.
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
double cdf;
int i;
double pdf;
int seed = 123456789;
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("LORENTZ_CDF_TEST");
Console.WriteLine("  LORENTZ_CDF evaluates the Lorentz CDF;");
Console.WriteLine("  LORENTZ_CDF_INV inverts the Lorentz CDF.");
Console.WriteLine("  LORENTZ_PDF evaluates the Lorentz PDF;");

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = lorentz_sample(seed);
pdf = lorentz_pdf(x);
cdf = lorentz_cdf(x);
x2 = lorentz_cdf_inv(cdf);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void lorentz_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_SAMPLE_TEST tests LORENTZ_SAMPLE.
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
# define SAMPLE_NUM 1000

int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("LORENTZ_SAMPLE_TEST");
Console.WriteLine("  LORENTZ_MEAN computes the Lorentz mean;");
Console.WriteLine("  LORENTZ_SAMPLE samples the Lorentz distribution;");
Console.WriteLine("  LORENTZ_VARIANCE computes the Lorentz variance.");

mean = lorentz_mean();
variance = lorentz_variance();

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = lorentz_sample(seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void maxwell_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_CDF_TEST tests MAXWELL_CDF.
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
Console.WriteLine("MAXWELL_CDF_TEST8");
Console.WriteLine("  MAXWELL_CDF evaluates the Maxwell CDF;");
Console.WriteLine("  MAXWELL_CDF_INV inverts the Maxwell CDF.");
Console.WriteLine("  MAXWELL_PDF evaluates the Maxwell PDF;");

a = 2.0E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!maxwell_check(a))
{
Console.WriteLine("");
Console.WriteLine("MAXWELL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = maxwell_sample(a, seed);
pdf = maxwell_pdf(x, a);
cdf = maxwell_cdf(x, a);
x2 = maxwell_cdf_inv(cdf, a);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void maxwell_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_SAMPLE_TEST tests MAXWELL_SAMPLE.
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
# define SAMPLE_NUM 1000

double a;
int j;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
double xmax;
double xmin;

Console.WriteLine("");
Console.WriteLine("MAXWELL_SAMPLE_TEST");
Console.WriteLine("  MAXWELL_MEAN computes the Maxwell mean;");
Console.WriteLine("  MAXWELL_SAMPLE samples the Maxwell distribution;");
Console.WriteLine("  MAXWELL_VARIANCE computes the Maxwell variance.");

a = 2.0E+00;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =             " + a + "");

if (!maxwell_check(a))
{
Console.WriteLine("");
Console.WriteLine("MAXWELL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = maxwell_mean(a);
variance = maxwell_variance(a);

Console.WriteLine("  PDF mean =                    " + mean + "");
Console.WriteLine("  PDF variance =                " + variance + "");

for (j = 0; j < SAMPLE_NUM; j++)
{
x[j] = maxwell_sample(a, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
}

static void multinomial_coef_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_TEST tests MULTINOMIAL_COEF1, MULTINOMIAL_COEF2.
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
# define MAXFACTOR 5

int factor[MAXFACTOR];
int i;
int j;
int n;
int ncomb1;
int ncomb2;
int nfactor;

Console.WriteLine("");
Console.WriteLine("MULTINOMIAL_TEST");
Console.WriteLine("  MULTINOMIAL_COEF1 computes multinomial");
Console.WriteLine("  coefficients using the Gamma function;");
Console.WriteLine("  MULTINOMIAL_COEF2 computes multinomial");
Console.WriteLine("  coefficients directly.");

Console.WriteLine("");
Console.WriteLine("  Line 10 of the BINOMIAL table:");
Console.WriteLine("");

n = 10;
nfactor = 2;

for (i = 0; i <= n; i++)
{
factor[0] = i;
factor[1] = n - i;

ncomb1 = multinomial_coef1(nfactor, factor);

ncomb2 = multinomial_coef2(nfactor, factor);

Console.WriteLine("  "
+ setw(2) + factor[0] + "  "
+ setw(2) + factor[1] + "  "
+ setw(5) + ncomb1 + "  "
+ setw(5) + ncomb2 + "");
}

Console.WriteLine("");
Console.WriteLine("  Level 5 of the TRINOMIAL coefficients:");

n = 5;
nfactor = 3;

for (i = 0; i <= n; i++)
{
factor[0] = i;

Console.WriteLine("");

for (j = 0; j <= n - factor[0]; j++)
{
factor[1] = j;
factor[2] = n - factor[0] - factor[1];

ncomb1 = multinomial_coef1(nfactor, factor);

ncomb2 = multinomial_coef2(nfactor, factor);

Console.WriteLine("  "
  + setw(2) + factor[0] + "  "
  + setw(2) + factor[1] + "  "
  + setw(2) + factor[2] + "  "
  + setw(5) + ncomb1 + "  "
  + setw(5) + ncomb2 + "");
}
}

return;
# undef MAXFACTOR
}

static void multinomial_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_SAMPLE_TEST tests MULTINOMIAL_SAMPLE.
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
# define B 3
# define SAMPLE_NUM 1000

int a;
double c[B] =  {
0.125, 0.500, 0.375
}
;
int i;
int j;
double* mean;
int seed = 123456789;
double* variance;
int x[B * SAMPLE_NUM];
int* xmax;
int* xmin;
int* y;

Console.WriteLine("");
Console.WriteLine("MULTINOMIAL_SAMPLE_TEST");
Console.WriteLine("  MULTINOMIAL_MEAN computes the Multinomial mean;");
Console.WriteLine("  MULTINOMIAL_SAMPLE samples the Multinomial distribution;");
Console.WriteLine("  MULTINOMIAL_VARIANCE computes the Multinomial variance;");

a = 5;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + B + "");
r8vec_print(B, c, "  PDF parameter C:");

if (!multinomial_check(a, B, c))
{
Console.WriteLine("");
Console.WriteLine("MULTINOMIAL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = multinomial_mean(a, B, c);
variance = multinomial_variance(a, B, c);
r8vec_print(B, mean, "  PDF mean:");
r8vec_print(B, variance, "  PDF variance:");

delete[] mean;
delete[] variance;

for (j = 0; j < SAMPLE_NUM; j++)
{
y = multinomial_sample(a, B, c, seed);
for (i = 0; i < B; i++)
{
x[i + j * B] = y[i];
}

delete[] y;
}

mean = i4row_mean(B, SAMPLE_NUM, x);
variance = i4row_variance(B, SAMPLE_NUM, x);
xmax = i4row_max(B, SAMPLE_NUM, x);
xmin = i4row_min(B, SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("");
Console.WriteLine("  Component Mean, Variance, Min, Max:");
Console.WriteLine("");

for (i = 0; i < B; i++)
{
Console.WriteLine("  "
+ setw(6) + i + 1 + "  "
+ setw(12) + mean[i] + "  "
+ setw(12) + variance[i] + "  "
+ setw(12) + xmin[i] + "  "
+ setw(12) + xmax[i] + "");
}

delete[] mean;
delete[] variance;
delete[] xmax;
delete[] xmin;

return;
# undef B
# undef SAMPLE_NUM
}

static void multinomial_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_PDF_TEST tests MULTINOMIAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define B 3

int a;
double c[B] =  {
0.1, 0.5, 0.4
}
;
double pdf;
int x[B] =  {
0, 2, 3
}
;

Console.WriteLine("");
Console.WriteLine("MULTINOMIAL_PDF_TEST");
Console.WriteLine("  MULTINOMIAL_PDF evaluates the Multinomial PDF;");

a = 5;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + B + "");
r8vec_print(B, c, "  PDF parameter C:");

if (!multinomial_check(a, B, c))
{
Console.WriteLine("");
Console.WriteLine("MULTINOMIAL_PDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

i4vec_print(B, x, "  PDF argument X:");

pdf = multinomial_pdf(x, a, B, c);

Console.WriteLine("");
Console.WriteLine("  PDF value = " + pdf + "");

return;
# undef B
# undef SAMPLE_NUM
}

static void multinoulli_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOULLLI_PDF_TEST tests MULTINOULLI_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2018
//
//  Author:
//
//    John Burkardt
//
{
int n = 5;
double pdf;
int seed;
double* theta;
double theta_sum;
int x;

Console.WriteLine("");
Console.WriteLine("MULTINOULLI_PDF_TEST");
Console.WriteLine("  MULTINOULLI_PDF evaluates the Multinoulli PDF.");

seed = 123456789;
theta = r8vec_uniform_01_new(n, seed);
theta_sum = 0.0;
for (x = 0; x < n; x++)
{
theta_sum = theta_sum + theta[x];
}

for (x = 0; x < n; x++)
{
theta[x] = theta[x] / theta_sum;
}

Console.WriteLine("");
Console.WriteLine("   X     pdf(X)");
Console.WriteLine("");
for (x = -1; x <= n; x++)
{
pdf = multinoulli_pdf(x, n, theta);
Console.WriteLine("  " + setw(2) + x
+ "  " + setw(14) + pdf + "");
}

delete[] theta;

return;
}

static void nakagami_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_CDF_TEST tests NAKAGAMI_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2016
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
double x;
double x2;

Console.WriteLine("");
Console.WriteLine("NAKAGAMI_CDF_TEST");
Console.WriteLine("  NAKAGAMI_CDF evaluates the Nakagami CDF;");
Console.WriteLine("  NAKAGAMI_CDF_INV inverts the Nakagami CDF;");
Console.WriteLine("  NAKAGAMI_PDF evaluates the Nakagami PDF;");

a = 1.0;
b = 2.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!nakagami_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("NAKAGAMI_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF         CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = a + b * sqrt((double) (i) / c / 10.0);
pdf = nakagami_pdf(x, a, b, c);
cdf = nakagami_cdf(x, a, b, c);
x2 = nakagami_cdf_inv(cdf, a, b, c);

Console.WriteLine("  " + setw(12) + x
+ "  " + setw(12) + pdf
+ "  " + setw(12) + cdf
+ "  " + setw(12) + x2 + "");
}

return;
}

static void nakagami_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_SAMPLE_TEST tests NAKAGAMI_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
{
double a;
double b;
double c;
double mean;
double variance;

Console.WriteLine("");
Console.WriteLine("NAKAGAMI_SAMPLE_TEST");
Console.WriteLine("  NAKAGAMI_MEAN evaluates the Nakagami mean;");
Console.WriteLine("  NAKAGAMI_VARIANCE evaluates the Nakagami variance;");

a = 1.0;
b = 2.0;
c = 3.0;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");
Console.WriteLine("  PDF parameter C =      " + c + "");

if (!nakagami_check(a, b, c))
{
Console.WriteLine("");
Console.WriteLine("NAKAGAMI_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = nakagami_mean(a, b, c);
variance = nakagami_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =      " + mean + "");
Console.WriteLine("  PDF variance =  " + variance + "");

return;
}

static void negative_binomial_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_TEST tests NEGATIVE_BINOMIAL_CDF.
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
int a;
double b;
double cdf;
int i;
double pdf;
int seed = 123456789;
int x;
int x2;

Console.WriteLine("");
Console.WriteLine("NEGATIVE_BINOMIAL_CDF_TEST");
Console.WriteLine("  NEGATIVE_BINOMIAL_CDF evaluates the Negative Binomial CDF;");
Console.WriteLine("  NEGATIVE_BINOMIAL_CDF_INV inverts the Negative Binomial CDF.");
Console.WriteLine("  NEGATIVE_BINOMIAL_PDF evaluates the Negative Binomial PDF;");

a = 2;
b = 0.25;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!negative_binomial_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("NEGATIVE_BINOMIAL_CDF_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

Console.WriteLine("");
Console.WriteLine("       X            PDF           CDF            CDF_INV");
Console.WriteLine("");

for (i = 1; i <= 10; i++)
{
x = negative_binomial_sample(a, b, seed);
pdf = negative_binomial_pdf(x, a, b);
cdf = negative_binomial_cdf(x, a, b);
x2 = negative_binomial_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
}

return;
}

static void negative_binomial_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_SAMPLE_TEST tests NEGATIVE_BINOMIAL_SAMPLE.
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
# define SAMPLE_NUM 1000

int a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
int x[SAMPLE_NUM];
int xmax;
int xmin;

Console.WriteLine("");
Console.WriteLine("NEGATIVE_BINOMIAL_SAMPLE_TEST");
Console.WriteLine("  NEGATIVE_BINOMIAL_MEAN computes the Negative Binomial mean;");
Console.WriteLine("  NEGATIVE_BINOMIAL_SAMPLE samples the Negative Binomial distribution;");
Console.WriteLine("  NEGATIVE_BINOMIAL_VARIANCE computes the Negative Binomial variance;");

a = 2;
b = 0.75;

Console.WriteLine("");
Console.WriteLine("  PDF parameter A =      " + a + "");
Console.WriteLine("  PDF parameter B =      " + b + "");

if (!negative_binomial_check(a, b))
{
Console.WriteLine("");
Console.WriteLine("NEGATIVE_BINOMIAL_SAMPLE_TEST - Fatal error!");
Console.WriteLine("  The parameters are not legal.");
return;
}

mean = negative_binomial_mean(a, b);
variance = negative_binomial_variance(a, b);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = negative_binomial_sample(a, b, seed);
}

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

return;
# undef SAMPLE_NUM
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
# define SAMPLE_NUM 1000

double mean;
int seed = 123456789;
double variance;
double* x;
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

delete[] x;

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double mean;
double mu;
int seed = 123456789;
double sigma;
double variance;
double* x;
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

delete[] x;

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
+ "  " + setw(14) + pdf
+ "  " + setw(14) + cdf
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
double* x;
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

mean = r8vec_mean(sample_num, x);
variance = r8vec_variance(sample_num, x);
xmax = r8vec_max(sample_num, x);
xmin = r8vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

delete[] x;

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
+ "  " + setw(14) + pdf
+ "  " + setw(14) + cdf
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
double* x;
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

mean = r8vec_mean(sample_num, x);
variance = r8vec_variance(sample_num, x);
xmax = r8vec_max(sample_num, x);
xmin = r8vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

delete[] x;

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
+ "  " + setw(14) + pdf
+ "  " + setw(14) + cdf
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
double* x;
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

mean = r8vec_mean(sample_num, x);
variance = r8vec_variance(sample_num, x);
xmax = r8vec_max(sample_num, x);
xmin = r8vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

delete[] x;

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

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

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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

mean = i4vec_mean(sample_num, x);
variance = i4vec_variance(sample_num, x);
xmax = i4vec_max(sample_num, x);
xmin = i4vec_min(sample_num, x);

Console.WriteLine("");
Console.WriteLine("  Sample size =     " + sample_num + "");
Console.WriteLine("  Sample mean =     " + mean + "");
Console.WriteLine("  Sample variance = " + variance + "");
Console.WriteLine("  Sample maximum =  " + xmax + "");
Console.WriteLine("  Sample minimum =  " + xmin + "");

delete[] x;

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

Console.WriteLine("  " + setw(12) + x
+ "  " + setw(12) + y
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

min = r8vec_min(N, x);
max = r8vec_max(N, x);
mean = r8vec_mean(N, x);
variance = r8vec_variance(N, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
int j;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
  + "  " + setw(14) + pdf + "");

pdf_total = pdf_total + pdf;
}

Console.WriteLine("  " + setw(8) + m
+ "  " + "        "
+ "  " + "        "
+ "  " + setw(14) + pdf_total + "");

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
# define SAMPLE_NUM 1000

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

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
double c;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
double c;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
double* x;

Console.WriteLine("");
Console.WriteLine("UNIFORM_01_ORDER_SAMPLE_TEST");
Console.WriteLine("  For the Uniform 01 Order PDF:");
Console.WriteLine("  UNIFORM_ORDER_SAMPLE samples.");

n = 10;
x = uniform_01_order_sample(n, seed);

r8vec_print(n, x, "  Ordered sample:");

delete[] x;

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
double* x;

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
Console.WriteLine(setw(12) + x[j] + "  ";
}

Console.WriteLine("");
delete[] x;
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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

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

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_circular_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
x = weibull_sample(a, b, c, seed);
pdf = weibull_pdf(x, a, b, c);
cdf = weibull_cdf(x, a, b, c);
x2 = weibull_cdf_inv(cdf, a, b, c);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
double c;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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

mean = weibull_mean(a, b, c);
variance = weibull_variance(a, b, c);

Console.WriteLine("");
Console.WriteLine("  PDF mean =     " + mean + "");
Console.WriteLine("  PDF variance = " + variance + "");

for (i = 0; i < SAMPLE_NUM; i++)
{
x[i] = weibull_sample(a, b, c, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
x = weibull_discrete_sample(a, b, seed);
pdf = weibull_discrete_pdf(x, a, b);
cdf = weibull_discrete_cdf(x, a, b);
x2 = weibull_discrete_cdf_inv(cdf, a, b);

Console.WriteLine("  "
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

double a;
double b;
int i;
double mean;
int seed = 123456789;
double variance;
double x[SAMPLE_NUM];
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
x[i] = weibull_discrete_sample(a, b, seed);
}

mean = r8vec_mean(SAMPLE_NUM, x);
variance = r8vec_variance(SAMPLE_NUM, x);
xmax = r8vec_max(SAMPLE_NUM, x);
xmin = r8vec_min(SAMPLE_NUM, x);

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
+ setw(12) + x + "  "
+ setw(12) + pdf + "  "
+ setw(12) + cdf + "  "
+ setw(12) + x2 + "");
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
# define SAMPLE_NUM 1000

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

mean = i4vec_mean(SAMPLE_NUM, x);
variance = i4vec_variance(SAMPLE_NUM, x);
xmax = i4vec_max(SAMPLE_NUM, x);
xmin = i4vec_min(SAMPLE_NUM, x);

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