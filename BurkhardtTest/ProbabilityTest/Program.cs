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
        }
    }
}