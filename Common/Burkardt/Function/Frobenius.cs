using Burkardt.Types;

namespace Burkardt.Function;

public static class Frobenius
{
    public static int frobenius_number_order2(int c1, int c2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FROBENIUS_NUMBER_ORDER2 returns the Frobenius number for order 2.
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
        //    For a general N, and a set of coin denominations C1, C2, ..., CN,
        //    the Frobenius number F(N, C(1:N) ) is defined as the largest value
        //    B for which the equation
        //
        //      C1*X1 + C2*X2 + ... + CN*XN = B
        //
        //    has no nonnegative integer solution X(1:N).
        //
        //    In the Mathematica Package "NumberTheory", the Frobenius number
        //    can be determined by
        //
        //    <<NumberTheory`Frobenius`
        //    FrobeniusF[ {C1,...,CN} ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2007
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
        //    Input, int C1, C2, the coin denominations. C1 and C2
        //    should be positive and relatively prime.
        //
        //    Output, int FROBENIUS_NUMBER_ORDER2, the Frobenius number of (C1,C2).
        //
    {
        const int i4_huge = 2147483647;
        int value;

        switch (c1)
        {
            case <= 0:
                value = i4_huge;
                break;
            default:
            {
                switch (c2)
                {
                    case <= 0:
                        value = i4_huge;
                        break;
                    default:
                    {
                        if (typeMethods.i4_gcd(c1, c2) != 1)
                        {
                            value = i4_huge;
                        }
                        else
                        {
                            value = c1 * c2 - c1 - c2;
                        }

                        break;
                    }
                }

                break;
            }
        }

        return value;
    }

}