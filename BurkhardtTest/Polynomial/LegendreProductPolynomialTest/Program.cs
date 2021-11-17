using System;

namespace LegendreProductPolynomialTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LPP_TEST.
        //
        //  Discussion:
        //
        //    LPP_TEST tests the LEGENDRE_PRODUCT_POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LPP_TEST:");
        Console.WriteLine("  Test the LEGENDRE_PRODUCT_POLYNOMIAL library.");

        i4Test.i4_choose_test ( );
        i4Test.i4_uniform_ab_test ( );

        i4Test.i4vec_permute_test ( );
        i4Test.i4vec_print_test ( );
        i4Test.i4vec_sort_heap_index_a_test ( );
        i4Test.i4vec_sum_test ( );
        i4Test.i4vec_uniform_ab_new_test ( );

        r8Test.r8vec_permute_test ( );
        r8Test.r8vec_print_test ( );
        r8Test.r8vec_uniform_ab_new_test ( );

        r8Test.r8mat_print_test ( );
        r8Test.r8mat_print_some_test ( );
        r8Test.r8mat_uniform_ab_new_test ( );

        permTest.perm_uniform_test ( );

        compTest.comp_enum_test ( );
        compTest.comp_next_grlex_test ( );
        compTest.comp_random_grlex_test ( );
        compTest.comp_rank_grlex_test ( );
        compTest.comp_unrank_grlex_test ( );

        monoTest.mono_next_grlex_test ( );
        monoTest.mono_print_test ( );
        monoTest.mono_rank_grlex_test ( );
        monoTest.mono_unrank_grlex_test ( );
        monoTest.mono_upto_enum_test ( );
        monoTest.mono_upto_next_grlex_test ( ); 
        monoTest.mono_upto_random_test ( );

        polyTest.polynomial_compress_test ( );
        polyTest.polynomial_print_test ( );
        polyTest.polynomial_sort_test ( );
        polyTest.polynomial_value_test ( );

        lpTest.lp_coefficients_test ( );
        lpTest.lp_value_test ( );
        lpTest.lp_values_test ( );

        lppTest.lpp_to_polynomial_test ( );
        lppTest.lpp_value_test ( );

        Console.WriteLine("");
        Console.WriteLine("LPP_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}