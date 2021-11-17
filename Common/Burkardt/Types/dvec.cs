using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int[] get_seq(int d, int norm, int seq_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GET_SEQ generates all positive integer D-vectors that sum to NORM.
        //
        //  Discussion:
        //
        //    This function computes a list, in reverse dictionary order, of
        //    all D-vectors of positive values that sum to NORM.
        //
        //    For example, call get_seq ( 3, 5, 6, fs ) returns
        //
        //      3  1  1
        //      2  2  1
        //      2  1  2
        //      1  3  1
        //      1  2  2
        //      1  1  3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 December 2012
        //
        //  Author:
        //
        //    Original MATLAB version by Florian Heiss, Viktor Winschel.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //  Parameters:
        //
        //    Input, int D, the dimension.
        //    1 <= D.
        //
        //    Input, int NORM, the value that each row must sum to.
        //    D <= NORM.
        //
        //    Input, int SEQ_NUM, the number of rows of FS.
        //
        //    Output, int GET_SEQ[SEQ_NUM*D].  Each row of FS represents 
        //    one vector with all elements positive and summing to NORM.
        //
    {
        int a;
        int c;
        int[] fs;
        int i;
        int row;
        int[] seq;

        seq = new int[d];
        fs = new int[seq_num * d];

        if (norm < d)
        {
            Console.WriteLine("");
            Console.WriteLine("GET_SEQ - Fatal error!");
            Console.WriteLine("  NORM = " + norm + " < D = " + "");
            return null;
        }

        for (i = 0; i < d; i++)
        {
            seq[i] = 0;
        }

        //
        //  The algorithm is written to work with vectors whose minimum value is
        //  allowed to be zero.  So we subtract D from NORM at the beginning and
        //  then increment the result vectors by 1 at the end!
        //
        a = norm - d;
        seq[0] = a;

        row = 0;
        for (i = 0; i < d; i++)
        {
            fs[row + i * seq_num] = seq[i] + 1;
        }

        c = 0;

        while (seq[d - 1] < a)
        {
            if (c == d - 1)
            {
                for (i = c - 1; 0 <= i; i--)
                {
                    c = i;
                    if (seq[i] != 0)
                    {
                        break;
                    }
                }
            }

            seq[c] -= 1;
            c += 1;
            seq[c] = a;
            for (i = 0; i < c; i++)
            {
                seq[c] -= seq[i];
            }

            if (c < d - 1)
            {
                for (i = c + 1; i < d; i++)
                {
                    seq[i] = 0;
                }
            }

            row += 1;
            for (i = 0; i < d; i++)
            {
                fs[row + i * seq_num] = seq[i] + 1;
            }
        }

        return fs;
    }

    public static void dvec_add(int n, int[] dvec1, int[] dvec2, ref int[] dvec3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVEC_ADD adds two (signed) decimal vectors.
        //
        //  Discussion:
        //
        //    A DVEC is an integer vector of decimal digits, intended to
        //    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
        //    is the coefficient of 10^(N-2), and DVEC(N) contains sign
        //    information.  It is 0 if the number is positive, and 9 if
        //    the number is negative.
        //
        //  Example:
        //
        //    N = 4
        //
        //      DVEC1     +   DVEC2     =   DVEC3
        //
        //    ( 0 0 1 7 ) + ( 0 1 0 4 ) = ( 0 0 1 2 1 )
        //
        //          17    +       104   =         121
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the length of the vectors.
        //
        //    Input, int DVEC1[N], DVEC2[N], the vectors to be added.
        //
        //    Output, int DVEC3[N], the sum of the two input vectors.
        //
    {
        int base_ = 10;
        int i;

        for (i = 0; i < n; i++)
        {
            dvec3[i] = dvec1[i] + dvec2[i];
        }

        for (i = 0; i < n; i++)
        {
            while (base_ <= dvec3[i])
            {
                dvec3[i] -= base_;
                if (i < n - 1)
                {
                    dvec3[i + 1] += 1;
                }
                else
                {
                    Console.WriteLine("dec_add: Overflow error!");
                    return;
                }
            }
        }
    }

    public static void dvec_complementx(int n, int[] dvec1, ref int[] dvec2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVEC_COMPLEMENTX computes the ten's complement of a decimal vector.
        //
        //  Discussion:
        //
        //    A DVEC is an integer vector of decimal digits, intended to
        //    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
        //    is the coefficient of 10^(N-2), and DVEC(N) contains sign
        //    information.  It is 0 if the number is positive, and 9 if
        //    the number is negative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the length of the vectors.
        //
        //    Input, int DVEC1[N], the vector to be complemented.
        //
        //    Output, int DVEC2[N], the complemented vector.
        //
    {
        int base_ = 10;
        int[] dvec3;
        int[] dvec4;
        int i;

        dvec3 = new int[n];
        dvec4 = new int[n];

        for (i = 0; i < n; i++)
        {
            dvec3[i] = base_ - 1 - dvec1[i];
        }

        dvec4[0] = 1;
        for (i = 1; i < n; i++)
        {
            dvec4[i] = 0;
        }

        dvec_add(n, dvec3, dvec4, ref dvec2);
    }

    public static void dvec_mul(int n, int[] dvec1, int[] dvec2, ref int[] dvec3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVEC_MUL computes the product of two decimal vectors.
        //
        //  Discussion:
        //
        //    A DVEC is an integer vector of decimal digits, intended to
        //    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
        //    is the coefficient of 10^(N-2), and DVEC(N) contains sign
        //    information.  It is 0 if the number is positive, and 9 if
        //    the number is negative.
        //
        //    Since the user may want to make calls like
        //
        //      dvec_mul ( n, dvec1, dvec1, dvec3 )
        //    or even
        //      dvec_mul ( n, dvec1, dvec1, dvec1 )
        //
        //    we need to copy the arguments, work on them, and then copy out the result.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the length of the vectors.
        //
        //    Input, int DVEC1[N], DVEC2[N], the vectors to be multiplied.
        //
        //    Output, int DVEC3[N], the product of the two input vectors.
        //
    {
        int base_ = 10;
        int carry;
        int[] dveca;
        int[] dvecb;
        int[] dvecc;
        int i;
        int j;
        int product_sign;

        dveca = new int[n];
        dvecb = new int[n];
        dvecc = new int[n];
        //
        //  Copy the input.
        //
        for (i = 0; i < n; i++)
        {
            dveca[i] = dvec1[i];
        }

        for (i = 0; i < n; i++)
        {
            dvecb[i] = dvec2[i];
        }

        //
        //  Record the sign of the product.
        //  Make the factors positive.
        //
        product_sign = 1;

        if (dveca[n - 1] != 0)
        {
            product_sign = -product_sign;
            dvec_complementx(n, dveca, ref dveca);
        }

        if (dvecb[n - 1] != 0)
        {
            product_sign = -product_sign;
            dvec_complementx(n, dvecb, ref dvecb);
        }

        for (i = 0; i < n; i++)
        {
            dvecc[i] = 0;
        }

        //
        //  Multiply.
        //
        for (i = 1; i <= n - 1; i++)
        {
            for (j = i; j <= n - 1; j++)
            {
                dvecc[j - 1] += dveca[i - 1] * dvecb[j - i];
            }
        }

        //
        //  Take care of carries.
        //  Unlike the DVEC_ADD routine, we do NOT allow carries into the
        //  N-th position.
        //
        for (i = 0; i < n - 1; i++)
        {
            carry = dvecc[i] / base_;
            dvecc[i] -= carry * base_;

            if (i < n - 2)
            {
                dvecc[i + 1] += carry;
            }
        }

        switch (product_sign)
        {
            //
            //  Take care of the sign of the product.
            //
            case < 0:
                dvec_complementx(n, dvecc, ref dvecc);
                break;
        }

        //
        //  Copy the output.
        //
        for (i = 0; i < n; i++)
        {
            dvec3[i] = dvecc[i];
        }
    }

    public static void dvec_print(int n, int[] dvec, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVEC_PRINT prints a decimal integer vector, with an optional title.
        //
        //  Discussion:
        //
        //    A DVEC is an integer vector of decimal digits, intended to
        //    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
        //    is the coefficient of 10^(N-2), and DVEC(N) contains sign
        //    information.  It is 0 if the number is positive, and 9 if
        //    the number is negative.
        //
        //    The vector is printed "backwards", that is, the first entry
        //    printed is DVEC(N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, int DVEC[N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;
        string cout = "";

        switch (title.Length)
        {
            case > 0:
                Console.WriteLine("");
                Console.WriteLine(title + "");
                Console.WriteLine("");
                break;
        }

        cout += dvec[n - 1] switch
        {
            9 => "-",
            _ => "+"
        };

        for (i = n - 2; 0 <= i; i--)
        {
            cout += dvec[i];
        }

        Console.WriteLine(cout);

    }

    public static void dvec_sub(int n, int[] dvec1, int[] dvec2, ref int[] dvec3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVEC_SUB subtracts two decimal vectors.
        //
        //  Discussion:
        //
        //    A DVEC is an integer vector of decimal digits, intended to
        //    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
        //    is the coefficient of 10^(N-2), and DVEC(N) contains sign
        //    information.  It is 0 if the number is positive, and 9 if
        //    the number is negative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the length of the vectors.
        //
        //    Input, int DVEC1[N], DVEC2[N]), the vectors to be subtracted.
        //
        //    Output, int DVEC3[N], the value of DVEC1 - DVEC2.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            dvec3[i] = dvec2[i];
        }

        dvec_complementx(n, dvec3, ref dvec3);

        dvec_add(n, dvec1, dvec3, ref dvec3);

    }

    public static int dvec_to_i4(int n, ref int[] dvec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVEC_TO_I4 makes an integer from a (signed) decimal vector.
        //
        //  Discussion:
        //
        //    A DVEC is an integer vector of decimal digits, intended to
        //    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
        //    is the coefficient of 10^(N-2), and DVEC(N) contains sign
        //    information.  It is 0 if the number is positive, and 9 if
        //    the number is negative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //
        //    Input, int DVEC[N], the decimal vector.
        //
        //    Output, int DVEC_TO_I4, the integer.
        //
    {
        int base_ = 10;
        int[] dvec2;
        int i;
        int i_sign;
        int i4;

        dvec2 = new int[n];

        for (i = 0; i < n; i++)
        {
            dvec2[i] = dvec[i];
        }

        i_sign = 1;

        if (dvec2[n - 1] == base_ - 1)
        {
            i_sign = -1;
            dvec_complementx(n - 1, dvec2, ref dvec2);
        }

        i4 = 0;
        for (i = n - 2; 0 <= i; i--)
        {
            i4 = base_ * i4 + dvec2[i];
        }

        i4 = i_sign * i4;

        return i4;
    }
}