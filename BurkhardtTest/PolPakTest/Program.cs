using System;

namespace PolPakTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    polpak_test() tests polpak().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("polpak_test():");
            Console.WriteLine("  Test polpak().");

            gudermannianTest.agud_test();
            alignTest.align_enum_test();
            bellTest.bell_test();
            bellTest.bell_poly_coef_test();
            benfordTest.benford_test();
            bernoulliTest.bernoulli_number_test();
            bernoulliTest.bernoulli_number2_test();
            bernoulliTest.bernoulli_number3_test();
            bernoulliTest.bernoulli_poly_test();
            bernoulliTest.bernoulli_poly2_test();
            bernsteinTest.bernstein_poly_test();
            bernsteinTest.bpab_test();
            cardanTest.cardan_poly_test();
            cardanTest.cardan_poly_coef_test();
            cardinalTest.cardinal_cos_test();
            cardinalTest.cardinal_sin_test();
            catalanTest.catalan_test();
            catalanTest.catalan_row_next_test();
            charlierTest.charlier_test();
            chebyshevTest.cheby_t_poly_test();
            chebyshevTest.cheby_t_poly_coef_test();
            chebyshevTest.cheby_t_poly_zero_test();
            chebyshevTest.cheby_u_poly_test();
            chebyshevTest.cheby_u_poly_coef_test();
            chebyshevTest.cheby_u_poly_zero_test();
            chebyshevTest.chebyshev_discrete_test();
            collatzTest.collatz_count_test();
            collatzTest.collatz_count_max_test();
            combTest.comb_row_next_test();
            commulTest.commul_test();
            compSymmPolyTest.complete_symmetric_poly_test();
            cospowerTest.cos_power_int_test();
            delannoyTest.delannoy_test();
            eulerTest.euler_number_test();
            eulerTest.euler_number2_test();
            eulerTest.euler_poly_test();
            eulerTest.eulerian_test();
            hofStadterTest.f_hofstadter_test();
            fibonacciTest.fibonacci_direct_test();
            fibonacciTest.fibonacci_floor_test();
            fibonacciTest.fibonacci_recursive_test();
            hofStadterTest.g_hofstadter_test();
            gegenbauerTest.gegenbauer_poly_test();
            hermiteTest.gen_hermite_poly_test();
            laguerreTest.gen_laguerre_poly_test();
            gudermannianTest.gud_test();
            hailTest.hail_test();
            hofStadterTest.h_hofstadter_test();
            hermiteTest.hermite_poly_phys_test();
            hermiteTest.hermite_poly_phys_coef_test();
            i4Test.i4_choose_test();
            i4Test.i4_factor_test();
            i4Test.i4_factorial_test();
            i4Test.i4_factorial2_test();
            i4Test.i4_is_fibonacci_test();
            i4Test.i4_is_triangular_test();
            i4Test.i4_partition_distinct_count_test();
            i4Test.i4_to_triangle_lower_test();
            jacobiTest.jacobi_poly_test();
            jacobiTest.jacobi_symbol_test();
            krawtchoukTest.krawtchouk_test();
            laguerreTest.laguerre_associated_test();
            laguerreTest.laguerre_poly_test();
            laguerreTest.laguerre_poly_coef_test();
            legendreTest.legendre_poly_test();
            legendreTest.legendre_poly_coef_test();
            legendreTest.legendre_associated_test();
            legendreTest.legendre_associated_normalized_test();
            legendreTest.legendre_function_q_test();
            legendreTest.legendre_symbol_test();
            lerchTest.lerch_test();
            lgammaTest.lgamma_test();
            meixnerTest.meixner_test();
            mertensTest.mertens_test();
            moebiusTest.moebius_test();
            motzkinTest.motzkin_test();
            normalTest.normal_01_cdf_inverse_test();
            omegaTest.omega_test();
            pentagonTest.pentagon_num_test();
            primeTest.phi_test();
            planepartitionTest.plane_partition_num_test();
            bernoulliTest.poly_bernoulli_test();
            polynomialTest.poly_coef_count_test();
            primeTest.prime_test();
            pyramidTest.pyramid_num_test();
            pyramidTest.pyramid_square_num_test();
            r8Test.r8_agm_test();
            r8Test.r8_beta_test();
            r8Test.r8_choose_test();
            r8Test.r8_cube_root_test();
            r8Test.r8_erf_test();
            r8Test.r8_erf_inverse_test();
            r8Test.r8_euler_constant_test();
            r8Test.r8_factorial_test();
            r8Test.r8_factorial_log_test();
            r8Test.r8_gamma_test();
            r8Test.r8_hyper_2f1_test();
            r8Test.r8_psi_test();
            r8Test.r8poly_degree_test();
            r8Test.r8poly_print_test();
            r8Test.r8poly_value_horner_test();
            sigmaTest.sigma_test();
            simplexTest.simplex_num_test();
            sinpowerTest.sin_power_int_test();
            slicesTest.slices_test();
            sphericalharmonicTest.spherical_harmonic_test();
            stirlingTest.stirling1_test();
            stirlingTest.stirling2_test();
            tauTest.tau_test();
            terahedronTest.tetrahedron_num_test();
            triangleTest.triangle_num_test();
            triangleTest.triangle_lower_to_i4_test();
            tribonacciTest.tribonacci_direct_test();
            tribonacciTest.tribonacci_recursive_test();
            tribonacciTest.tribonacci_roots_test();
            polynomialTest.trinomial_test();
            hofStadterTest.v_hofstadter_test();
            vibonacciTest.vibonacci_test();
            zeckendorfTest.zeckendorf_test();
            zernikeTest.zernike_poly_test();
            zernikeTest.zernike_poly_coef_test();
            zetaTest.zeta_m1_test();
            zetaTest.zeta_naive_test();

            Console.WriteLine("");
            Console.WriteLine("polpak_test():");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}