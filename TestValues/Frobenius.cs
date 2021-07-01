namespace TestValues
{
    public static class Frobenius
    {

        public static void frobenius_number_data_values(ref int n_data, int order, ref int[] c, ref int f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FROBENIUS_NUMBER_DATA_VALUES returns data for the Frobenius problem.
            //
            //  Discussion:
            //
            //    The user should first call FROBENIUS_NUMBER_ORDER_VALUES to get the
            //    order or size of the "C" vector that will be returned by this routine.
            //
            //    The Frobenius number of order N and data C is the solution of the
            //    Frobenius coin sum problem for N coin denominations C(1) through C(N).
            //
            //    The Frobenius coin sum problem assumes the existence of
            //    N coin denominations, and asks for the largest value that cannot
            //    be formed by any combination of coins of these denominations.
            //
            //    The coin denominations are assumed to be distinct positive integers.
            //
            //    For general N, this problem is fairly difficult to handle.
            //
            //    For N = 2, it is known that:
            //
            //    * if C1 and C2 are not relatively prime, then
            //      there are infinitely large values that cannot be formed.
            //
            //    * otherwise, the largest value that cannot be formed is
            //      C1 * C2 - C1 - C2, and that exactly half the values between
            //      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
            //
            //    As a simple example, if C1 = 2 and C2 = 7, then the largest
            //    unrepresentable value is 5, and there are (5+1)/2 = 3
            //    unrepresentable values, namely 1, 3, and 5.
            //
            //    For a general N, and a set of coin denominations C1, C2,, CN,
            //    the Frobenius number F(N, C(1:N) ) is defined as the largest value
            //    B for which the equation
            //
            //      C1*X1 + C2*X2 + ... + CN*XN = B
            //
            //    has no nonnegative integer solution X(1:N).
            //
            //    In Mathematica, the Frobenius number can be determined by
            //
            //      FrobeniusNumber[ {C1,...,CN} ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Gerard Cornuejols, Regina Urbaniak, Robert Weismantel, Laurence Wolsey,
            //    Decomposition of Integer Programs and of Generating Sets,
            //    in Algorithms - ESA 97,
            //    Lecture Notes in Computer Science 1284,
            //    edited by R Burkard, G Woeginger,
            //    Springer, 1997, pages 92-103.
            //
            //    Robert Owens,
            //    An Algorithm to Solve the Frobenius Problem,
            //    Mathematics Magazine,
            //    Volume 76, Number 4, October 2003, 264-275.
            //
            //    James Sylvester,
            //    Question 7382,
            //    Mathematical Questions with their Solutions,
            //    Educational Times,
            //    Volume 41, page 21, 1884.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input, ref int N_DATA.  Unlike most other routines in this
            //    library, this routine assumes that N_DATA has already been processed by a call
            //    to FROBENIUS_NUMBER_ORDER_VALUE.  Therefore, this routine will return the
            //    next set of data as long as N_DATA is in the expected range.
            //
            //    Input, int ORDER, the order of the problem.
            //
            //    Output, int C[ORDER], the denominations of the problem.
            //
            //    Output, ref int F, the value of the function.
            //
        {
            int CVEC_MAX = 77;
            int N_MAX = 19;

            int[] c_vec =
            {
                2, 5,
                3, 17,
                4, 19,
                5, 13,
                12, 11,
                99, 100,
                6, 9, 20,
                5, 17, 23,
                137, 251, 256,
                31, 41, 47, 61,
                271, 277, 281, 283,
                10, 18, 26, 33, 35,
                34, 37, 38, 40, 43,
                12223, 12224, 36674, 61119, 85569,
                1000, 1476, 3764, 4864, 4871, 7773,
                12228, 36679, 36682, 46908, 61139, 73365,
                12137, 36405, 24269, 36407, 84545, 60683,
                13211, 13212, 39638, 66060, 52864, 79268, 92482,
                13429, 26850, 26855, 40280, 40281, 53711, 53714, 67141
            };
            int[] f_vec =
            {
                3,
                31,
                53,
                47,
                109,
                9701,
                43,
                41,
                4948,
                240,
                13022,
                67,
                175,
                89643481,
                47350,
                89716838,
                58925134,
                104723595,
                45094583
            };
            int i;
            int v_data = 0;

            if (n_data < 1 || N_MAX < n_data)
            {
                n_data = 0;
                v_data = 0;
                for (i = 0; i < order; i++)
                {
                    c[i] = 0;
                }

                f = 0;
            }
            else
            {
                for (i = 0; i < order; i++)
                {
                    c[i] = c_vec[v_data + i];
                }

                v_data = v_data + order;
                if (n_data == N_MAX)
                {
                    v_data = 0;
                }

                f = f_vec[n_data - 1];
            }
        }

        public static void frobenius_number_order_values(ref int n_data, ref int order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FROBENIUS_NUMBER_ORDER_VALUES returns orders of the Frobenius problem.
            //
            //  Discussion:
            //
            //    This routine returns the order or size of a Frobenius problem.
            //    To get the corresponding data, call FROBENIUS_NUMBER_DATA_VALUES.
            //
            //    The Frobenius number of order N and data C is the solution of a Frobenius
            //    coin sum problem for N coin denominations C(1) through C(N).
            //
            //    The Frobenius coin sum problem assumes the existence of
            //    N coin denominations, and asks for the largest value that cannot
            //    be formed by any combination of coins of these denominations.
            //
            //    The coin denominations are assumed to be distinct positive integers.
            //
            //    For general order N, this problem is fairly difficult to handle.
            //
            //    For order N = 2, it is known that:
            //
            //    * if C1 and C2 are not relatively prime, then
            //      there are infinitely large values that cannot be formed.
            //
            //    * otherwise, the largest value that cannot be formed is
            //      C1 * C2 - C1 - C2, and that exactly half the values between
            //      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
            //
            //    As a simple example, if C1 = 2 and C2 = 7, then the largest
            //    unrepresentable value is 5, and there are (5+1)/2 = 3
            //    unrepresentable values, namely 1, 3, and 5.
            //
            //    For a general N, and a set of coin denominations C1, C2,, CN,
            //    the Frobenius number F(N, C(1:N) ) is defined as the largest value
            //    B for which the equation
            //
            //      C1*X1 + C2*X2 + ... + CN*XN = B
            //
            //    has no nonnegative integer solution X(1:N).
            //
            //    In Mathematica, the Frobenius number can be determined by
            //
            //      FrobeniusNumber[ {C1,...,CN} ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Gerard Cornuejols, Regina Urbaniak, Robert Weismantel, Laurence Wolsey,
            //    Decomposition of Integer Programs and of Generating Sets,
            //    in Algorithms - ESA 97,
            //    Lecture Notes in Computer Science 1284,
            //    edited by R Burkard, G Woeginger,
            //    Springer, 1997, pages 92-103,
            //    LC: QA76.9.A43.E83.
            //
            //    Robert Owens,
            //    An Algorithm to Solve the Frobenius Problem,
            //    Mathematics Magazine,
            //    Volume 76, Number 4, October 2003, 264-275.
            //
            //    James Sylvester,
            //    Question 7382,
            //    Mathematical Questions with their Solutions,
            //    Educational Times,
            //    Volume 41, page 21, 1884.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
            //    before the first call.  On each call, the routine increments N_DATA by 1,
            //    and returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int ORDER, the order of a Frobenius problem.
            //
        {
            int N_MAX = 19;

            int[] order_vec =
            {
                2,
                2,
                2,
                2,
                2,
                2,
                3,
                3,
                3,
                4,
                4,
                5,
                5,
                5,
                6,
                6,
                6,
                7,
                8
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
            }
            else
            {
                order = order_vec[n_data - 1];
            }
        }

        public static void frobenius_number_order2_values(ref int n_data, ref int c1, ref int c2, ref int f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FROBENIUS_NUMBER_ORDER2_VALUES returns values of the order 2 Frobenius number.
            //
            //  Discussion:
            //
            //    The Frobenius number of order N is the solution of the Frobenius
            //    coin sum problem for N coin denominations.
            //
            //    The Frobenius coin sum problem assumes the existence of
            //    N coin denominations, and asks for the largest value that cannot
            //    be formed by any combination of coins of these denominations.
            //
            //    The coin denominations are assumed to be distinct positive integers.
            //
            //    For general N, this problem is fairly difficult to handle.
            //
            //    For N = 2, it is known that:
            //
            //    * if C1 and C2 are not relatively prime, then
            //      there are infinitely large values that cannot be formed.
            //
            //    * otherwise, the largest value that cannot be formed is
            //      C1 * C2 - C1 - C2, and that exactly half the values between
            //      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
            //
            //    As a simple example, if C1 = 2 and C2 = 7, then the largest
            //    unrepresentable value is 5, and there are (5+1)/2 = 3
            //    unrepresentable values, namely 1, 3, and 5.
            //
            //    For a general N, and a set of coin denominations C1, C2,, CN,
            //    the Frobenius number F(N, C(1:N) ) is defined as the largest value
            //    B for which the equation
            //
            //      C1*X1 + C2*X2 + ... + CN*XN = B
            //
            //    has no nonnegative integer solution X(1:N).
            //
            //    In Mathematica, the Frobenius number can be determined by
            //
            //      FrobeniusNumber[ {C1,...,CN} ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 November 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    James Sylvester,
            //    Question 7382,
            //    Mathematical Questions with their Solutions,
            //    Educational Times,
            //    Volume 41, page 21, 1884.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int C1, &C2, the parameters of the function.
            //
            //    Output, ref int F, the value of the function.
            //
        {
            int N_MAX = 6;

            int[] c1_vec =
            {
                2,
                3,
                4,
                5,
                12,
                99
            };
            int[] c2_vec =
            {
                5,
                17,
                19,
                13,
                11,
                100
            };
            int[] f_vec =
            {
                3,
                31,
                53,
                47,
                109,
                9701
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                c1 = 0;
                c2 = 0;
                f = 0;
            }
            else
            {
                c1 = c1_vec[n_data - 1];
                c2 = c2_vec[n_data - 1];
                f = f_vec[n_data - 1];
            }
        }

    }
}