namespace TestValues
{
    public static class Polyomino
    {

        public static void polyomino_chiral_count_values(ref int n_data, ref int order,
                ref long number)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_CHIRAL_COUNT_VALUES counts chiral polyominoes (allowing holes).
            //
            //  Discussion:
            //
            //    Polyominoes are connected planar shapes formed by adjoining unit squares.
            //
            //    The number of unit squares in a polyomino is its order.
            //
            //    If we do not ignore reflections, but ignore rotations when comparing 
            //    then we are considering the class of "chiral" polyominoes.  In that case,
            //    for instance, there are 18 fixed polyominoes of order 5.
            //
            //    As the order increases, the number of polyominoes grows very rapidly.
            //    The list offered here goes no further than order 28, but the later
            //    numbers in the list are too large to represent as 32 byte integers. 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Solomon Golomb,
            //    Polyominoes: Puzzles, Patterns, Problems, and Packings,
            //    Princeton University Press, 1996,
            //    ISBN: 9780691024448
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 
            //    before the first call.  On each call, the routine increments N_DATA by 1, 
            //    and returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int ORDER, the order of a polyomino.
            //
            //    Output, long long ref int NUMBER, the number of chiral polyominos 
            //    of this order.
            //
        {
            int N_MAX = 31;

            long[] number_vec =
            {
                1L,
                1L,
                1L,
                2L,
                7L,
                18L,
                60L,
                196L,
                704L,
                2500L,
                9189L,
                33896L,
                126759L,
                476270L,
                1802312L,
                6849777L,
                26152418L,
                100203194L,
                385221143L,
                1485200848L,
                5741256764L,
                22245940545L,
                86383382827L,
                336093325058L,
                1309998125640L,
                5114451441106L,
                19998172734786L,
                78306011677182L,
                307022182222506L,
                1205243866707468L,
                4736694001644862L
            };
            int[] order_vec =
            {
                0,
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15,
                16, 17, 18, 19, 20,
                21, 22, 23, 24, 25,
                26, 27, 28, 29, 30
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                order = 0;
                number = 0;
            }
            else
            {
                order = order_vec[n_data - 1];
                number = number_vec[n_data - 1];
            }
        }

        public static void polyomino_fixed_count_values(ref int n_data, ref int order,
                ref long number)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_FIXED_COUNT_VALUES counts fixed polyominoes (allowing holes).
            //
            //  Discussion:
            //
            //    Polyominoes are connected planar shapes formed by adjoining unit squares.
            //
            //    The number of unit squares in a polyomino is its order.
            //
            //    If we do not ignore reflections and rotations when comparing polyominoes,
            //    then we are considering the class of "fixed" polyominoes.  In that case,
            //    for instance, there are 65 fixed polyominoes of order 5.
            //
            //    As the order increases, the number of polyominoes grows very rapidly.
            //    The list offered here goes no further than order 28, but the later
            //    numbers in the list are too large to represent as 32 byte integers. 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Solomon Golomb,
            //    Polyominoes: Puzzles, Patterns, Problems, and Packings,
            //    Princeton University Press, 1996,
            //    ISBN: 9780691024448
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 
            //    before the first call.  On each call, the routine increments N_DATA by 1, 
            //    and returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int ORDER, the order of a polyomino.
            //
            //    Output, long long ref int NUMBER, the number of fixed polyominos 
            //    of this order.
            //
        {
            int N_MAX = 29;

            long[] number_vec =
            {
                1L,
                1L,
                2L,
                6L,
                19L,
                63L,
                216L,
                760L,
                2725L,
                9910L,
                36446L,
                135268L,
                505861L,
                1903890L,
                7204874L,
                27394666L,
                104592937L,
                400795844L,
                1540820542L,
                5940738676L,
                22964779660L,
                88983512783L,
                345532572678L,
                1344372335524L,
                5239988770268L,
                20457802016011L,
                79992676367108L,
                313224032098244L,
                1228088671826973L
            };
            int[] order_vec =
            {
                0,
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15,
                16, 17, 18, 19, 20,
                21, 22, 23, 24, 25,
                26, 27, 28
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                order = 0;
                number = 0;
            }
            else
            {
                order = order_vec[n_data - 1];
                number = number_vec[n_data - 1];
            }
        }

        public static void polyomino_free_count_values(ref int n_data, ref int order,
                ref long number)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYOMINO_FREE_COUNT_VALUES counts free polyominoes (allowing holes).
            //
            //  Discussion:
            //
            //    Polyominoes are connected planar shapes formed by adjoining unit squares.
            //
            //    The number of unit squares in a polyomino is its order.
            //
            //    If we ignore reflections and rotations when comparing polyominoes,
            //    then we are considering the class of "free" polyominoes.  In that case,
            //    for instance, there are just 12 free polyominoes of order 5, the
            //    so called "pentominoes".
            //
            //    As the order increases, the number of polyominoes grows very rapidly.
            //    The list offered here goes no further than order 28, but the later
            //    numbers in the list are too large to represent as 32 byte integers. 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Solomon Golomb,
            //    Polyominoes: Puzzles, Patterns, Problems, and Packings,
            //    Princeton University Press, 1996,
            //    ISBN: 9780691024448
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 
            //    before the first call.  On each call, the routine increments N_DATA by 1, 
            //    and returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int ORDER, the order of a polyomino.
            //
            //    Output, long long ref int NUMBER, the number of free polyominos of 
            //    this order.
            //
        {
            int N_MAX = 29;

            long[] number_vec =
            {
                1L,
                1L,
                1L,
                2L,
                5L,
                12L,
                35L,
                108L,
                369L,
                1285L,
                4655L,
                17073L,
                63600L,
                238591L,
                901971L,
                3426576L,
                13079255L,
                50107909L,
                192622052L,
                742624232L,
                2870671950L,
                11123060678L,
                43191857688L,
                168047007728L,
                654999700403L,
                2557227044764L,
                9999088822075L,
                39153010938487L,
                153511100594603L
            };
            int[] order_vec =
            {
                0,
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15,
                16, 17, 18, 19, 20,
                21, 22, 23, 24, 25,
                26, 27, 28
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                order = 0;
                number = 0;
            }
            else
            {
                order = order_vec[n_data - 1];
                number = number_vec[n_data - 1];
            }
        }

    }
}