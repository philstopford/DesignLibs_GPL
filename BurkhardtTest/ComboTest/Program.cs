using System;

namespace ComboTest;

internal partial class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for COMBO_TEST.
        //
        //  Discussion:
        //
        //    COMBO_TEST tests the COMBO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("COMBO_TEST");

        Console.WriteLine("  Test the COMBO library.");

        backtrack_test();

        bal_seq_check_test();
        bal_seq_enum_test();
        bal_seq_rank_test();
        bal_seq_successor_test();
        bal_seq_to_tableau_test();
        bal_seq_unrank_test();

        bell_numbers_test();

        cycle_check_test();
        cycle_to_perm_test();

        dist_enum_test();
        dist_next_test();

        edge_check_test();
        edge_degree_test();
        edge_enum_test();

        gray_code_check_test();
        gray_code_enum_test();
        gray_code_rank_test();
        gray_code_successor_test();
        gray_code_unrank_test();

        i4_choose_test();
        i4_factorial_test();
        i4_fall_test();
        i4_uniform_ab_test();

        i4vec_backtrack_test();
        i4vec_dot_product_test();
        i4vec_part1_new_test();
        i4vec_part2_test();
        i4vec_part2_new_test();
        i4vec_search_binary_a_test();
        i4vec_search_binary_d_test();
        i4vec_sort_insert_a_test();
        i4vec_sort_insert_d_test();
        i4vec_uniform_ab_new_test();

        knapsack_01_test();
        knapsack_reorder_test();
        knapsack_rational_test();

        ksubset_colex_check_test();
        ksubset_colex_rank_test();
        ksubset_colex_successor_test();
        ksubset_colex_unrank_test();
        ksubset_enum_test();
        ksubset_lex_check_test();
        ksubset_lex_rank_test();
        ksubset_lex_successor_test();
        ksubset_lex_unrank_test();
        ksubset_revdoor_rank_test();
        ksubset_revdoor_successor_test();
        ksubset_revdoor_unrank_test();

        marriage_test();

        mountain_test();

        npart_enum_test();
        npart_rsf_lex_random_test();
        npart_rsf_lex_rank_test();
        npart_rsf_lex_successor_test();
        npart_rsf_lex_unrank_test();
        npart_sf_lex_successor_test();
        npart_table_test();

        part_enum_test();
        part_rsf_check_test();
        part_sf_check_test();
        part_sf_conjugate_test();
        part_sf_majorize_test();
        part_successor_test();
        part_table_test();

        partition_greedy_test();

        partn_enum_test();
        partn_sf_check_test();
        partn_successor_test();

        perm_check_test();
        perm_enum_test();
        perm_inv_test();
        perm_lex_rank_test();
        perm_lex_successor_test();
        perm_lex_unrank_test();
        perm_mul_test();
        perm_parity_test();
        perm_print_test();
        perm_random_test();
        perm_tj_rank_test();
        perm_tj_successor_test();
        perm_tj_unrank_test();
        perm_to_cycle_test();

        pruefer_check_test();
        pruefer_enum_test();
        pruefer_rank_test();
        pruefer_successor_test();
        pruefer_to_tree_test();
        pruefer_unrank_test();

        queens_test();

        r8_choose_test();
        r8_gamma_log_test();

        r8vec_backtrack_test();

        rgf_check_test();
        rgf_enum_test();
        rgf_g_table_test();
        rgf_rank_test();
        rgf_successor_test();
        rgf_to_setpart_test();
        rgf_unrank_test();

        setpart_check_test();
        setpart_enum_test();
        setpart_to_rgf_test();

        stirling_numbers1_test();
        stirling_numbers2_test();

        subset_check_test();
        subset_colex_rank_test();
        subset_colex_successor_test();
        subset_colex_unrank_test();
        subset_complement_test();
        subset_distance_test();
        subset_enum_test();
        subset_intersect_test();
        subset_lex_rank_test();
        subset_lex_successor_test();
        subset_lex_unrank_test();
        subset_union_test();
        subset_weight_test();
        subset_xor_test();

        subsetsumswap_test();

        tableau_check_test();
        tableau_enum_test();
        tableau_to_bal_seq_test();

        tree_check_test();
        tree_enum_test();
        tree_rank_test();
        tree_successor_test();
        tree_to_pruefer_test();
        tree_unrank_test();

        Console.WriteLine("");
        Console.WriteLine("COMBO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}