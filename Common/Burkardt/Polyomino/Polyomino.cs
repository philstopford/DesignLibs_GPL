using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.PolyominoNS;

public static class Polyomino
{
    public static void polyomino_condense(int mp, int np, int[] p, ref int mq, ref int nq, ref int[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_CONDENSE condenses a polyomino.
        //
        //  Discussion:
        //
        //    A polyomino is a shape formed by connecting unit squares edgewise.
        //
        //    A polyomino can be represented by an MxN matrix, whose entries are
        //    1 for squares that are part of the polyomino, and 0 otherwise.
        //
        //    This program is given an MxN matrix that is meant to represent a 
        //    polyomino.  It first replaces all nonzero entries by the value 1.
        //    It then "condenses" the matrix, if possible, by removing initial and
        //    final rows and columns that are entirely zero.
        //
        //    While this procedure might save a slight amount of space, its purpose
        //    is to simplify the task of manipulating polyominos, embedding them in
        //    larger shapes, and detecting whether two polyominos describe the same
        //    shape.
        //
        //    It is entirely possible, and usual, that the output quantities are
        //    simply copies of the input quantities.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MP, NP, the number of rows and columns in the representation
        //    of the polyomino P.
        //
        //    Input, int P[MP*NP], a matrix of 0's and 1's representing the 
        //    polyomino.  
        //
        //    Output, int &MQ, &NQ, the number of rows and columns of the
        //    condensed polyomino.
        //
        //    Output, int *Q[MQ*NQ], the representation of the condensed
        //    polyomino.
        //
    {
        int i;
        int j;
        //
        //  Discard nonsense.
        //
        if (mp <= 0 || np <= 0)
        {
            mq = 0;
            nq = 0;
            q = null;
            return;
        }

        //
        //  Seek first and last nonzero rows, columns.
        //
        int i_min = -1;
        for (i = 0; i < mp; i++)
        {
            for (j = 0; j < np; j++)
            {
                if (p[i + j * mp] == 0)
                {
                    continue;
                }

                i_min = i;
                break;
            }

            if (i_min != -1)
            {
                break;
            }
        }

        switch (i_min)
        {
            //
            //  If I_MIN = -1, then we have a null matrix.
            //
            case -1:
                mq = 0;
                nq = 0;
                q = null;
                return;
        }

        int i_max = mp;
        for (i = mp - 1; 0 <= i; i--)
        {
            for (j = 0; j < np; j++)
            {
                if (p[i + j * mp] == 0)
                {
                    continue;
                }

                i_max = i;
                break;
            }

            if (i_max != mp)
            {
                break;
            }
        }

        int j_min = -1;
        for (j = 0; j < np; j++)
        {
            for (i = 0; i < mp; i++)
            {
                if (p[i + j * mp] == 0)
                {
                    continue;
                }

                j_min = j;
                break;
            }

            if (j_min != -1)
            {
                break;
            }
        }

        int j_max = np;
        for (j = np - 1; 0 <= j; j--)
        {
            for (i = 0; i < mp; i++)
            {
                if (p[i + j * mp] == 0)
                {
                    continue;
                }

                j_max = j;
                break;
            }

            if (j_max != np)
            {
                break;
            }
        }

        //
        //  Measure the nonzero block.
        //
        mq = i_max + 1 - i_min;
        nq = j_max + 1 - j_min;
        q = new int [mq * nq];
        //
        //  Copy the nonzero block.
        //
        for (j = 0; j < nq; j++)
        {
            for (i = 0; i < mq; i++)
            {
                if (p[i + i_min + (j + j_min) * mp] != 0)
                {
                    q[i + j * mq] = 1;
                }
                else
                {
                    q[i + j * mq] = 0;
                }
            }
        }
    }

    public static int[] polyomino_embed_list(int mr, int nr, int[] r, int mp, int np, int[] p,
            int number)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_EMBED_LIST lists the polyomino embeddings in a region.
        //
        //  Discusion:
        //
        //    A region R is a subset of an MRxNR grid of squares.
        //
        //    A polyomino P is a subset of an MPxNP grid of squares.
        //
        //    Both objects are represented by binary matrices, with the property that
        //    there are no initial or final zero rows or columns.
        //
        //    For this computation, we regard P as a "fixed" polyomino; in other words,
        //    no reflections or rotations will be allowed.
        //
        //    An "embedding" of P into R is an offset (MI,NJ) such that 
        //      P(I,J) = R(I+MI,J+NJ) 
        //      for 1 <= I <= MP, 1 <= J <= NP, and 
        //      for 0 <= MI <= MR-MP, 0 <= MJ <= NR-NP.
        //    We can detect an embedding simply by taking what amounts to a kind of
        //    dot product of P with a corresponding subregion of R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 May 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MR, NR, the number of rows and columns in the representation
        //    of the region R.
        //
        //    Input, int R(MR,NR), a matrix of 0's and 1's representing the 
        //    region.
        //
        //    Input, int MP, NP, the number of rows and columns in the representation
        //    of the polyomino P.
        //
        //    Input, int P[MP*NP], a matrix of 0's and 1's representing the 
        //    polyomino.
        //
        //    Input, int NUMBER, the number of embeddings.
        //
        //    Output, int POLYOMINO_EMBED_LIST[NUMBER*2], for each embedding, the I 
        //    and J offsets applied to the polyomino P.
        //
    {
        int mi;

        int[] list = new int [number * 2];
        //
        //  Count the 1's in P.
        //
        int pr = typeMethods.i4mat_sum(mp, np, ref p);
        //
        //  For each possible (I,J) coordinate of the upper left corner of a subset of R,
        //  see if it matches P.
        //
        int k = 0;
        for (mi = 0; mi <= mr - mp; mi++)
        {
            int nj;
            for (nj = 0; nj <= nr - np; nj++)
            {
                int srp = 0;
                int j;
                for (j = 0; j < np; j++)
                {
                    int i;
                    for (i = 0; i < mp; i++)
                    {
                        srp += p[i + j * mp] * r[i + mi + (j + nj) * mr];
                    }
                }

                if (srp != pr)
                {
                    continue;
                }

                list[k + 0 * number] = mi;
                list[k + 1 * number] = nj;
                k += 1;
            }
        }

        return list;
    }

    public static int polyomino_embed_number(int mr, int nr, int[] r, int mp, int np, int[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_EMBED_NUMBER counts the number of polyomino embeddings in a region.
        //
        //  Discusion:
        //
        //    A region R is a subset of an MRxNR grid of squares.
        //
        //    A polyomino P is a subset of an MPxNP grid of squares.
        //
        //    Both objects are represented by binary matrices, with the property that
        //    there are no initial or final zero rows or columns.
        //
        //    For this computation, we regard P as a "fixed" polyomino; in other words,
        //    no reflections or rotations will be allowed.
        //
        //    An "embedding" of P into R is an offset (MI,NJ) such that 
        //      P(I,J) = R(I+MI,J+NJ) 
        //      for 1 <= I <= MP, 1 <= J <= NP, and 
        //      for 0 <= MI <= MR-MP, 0 <= MJ <= NR-NP.
        //    We can detect an embedding simply by taking what amounts to a kind of
        //    dot product of P with a corresponding subregion of R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 May 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MR, NR, the number of rows and columns in the representation
        //    of the region R.
        //
        //    Input, int R[MR*NR], a matrix of 0's and 1's representing the 
        //    region.
        //
        //    Input, int MP, NP, the number of rows and columns in the representation
        //    of the polyomino P.
        //
        //    Input, int P[MP*NP], a matrix of 0's and 1's representing the 
        //    polyomino.
        //
        //    Output, int POLYOMINO_EMBED_NUMBER, the number of distinct embeddings of 
        //    P into R.
        //
    {
        int mi;

        int number = 0;
        //
        //  Count the 1's in P.
        //
        int pr = typeMethods.i4mat_sum(mp, np, ref p);
        //
        //  For each possible (I,J) coordinate of the upper left corner of a subset of R,
        //  see if it matches P.
        //
        for (mi = 0; mi <= mr - mp; mi++)
        {
            int nj;
            for (nj = 0; nj <= nr - np; nj++)
            {
                int srp = 0;
                int j;
                for (j = 0; j < np; j++)
                {
                    int i;
                    for (i = 0; i < mp; i++)
                    {
                        srp += p[i + j * mp] * r[i + mi + (j + nj) * mr];
                    }
                }

                if (srp == pr)
                {
                    number += 1;
                }
            }
        }

        return number;
    }

    public static void polyomino_enumerate_chiral(ref int n_data, ref int order,
            ref long number)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_ENUMERATE_CHIRAL counts chiral polyominoes (allowing holes).
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 
        //    before the first call.  On each call, the routine increments N_DATA by 1, 
        //    and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &ORDER, the order of a polyomino.
        //
        //    Output, long long int &NUMBER, the number of chiral polyominos 
        //    of this order.
        //
    {
        const int N_MAX = 31;

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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static void polyomino_enumerate_fixed(ref int n_data, ref int order,
            ref long number)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_ENUMERATE_FIXED counts fixed polyominoes (allowing holes).
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 
        //    before the first call.  On each call, the routine increments N_DATA by 1, 
        //    and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &ORDER, the order of a polyomino.
        //
        //    Output, long long int &NUMBER, the number of fixed polyominos 
        //    of this order.
        //
    {
        const int N_MAX = 29;

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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static void polyomino_enumerate_free(ref int n_data, ref int order,
            ref long number)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_ENUMERATE_FREE counts free polyominoes (allowing holes).
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 
        //    before the first call.  On each call, the routine increments N_DATA by 1, 
        //    and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &ORDER, the order of a polyomino.
        //
        //    Output, long long int &NUMBER, the number of free polyominos of 
        //    this order.
        //
    {
        const int N_MAX = 29;

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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static int[] polyomino_index(int m, int n, int[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_INDEX assigns an index to each nonzero entry of a polyomino.
        //
        //  Example:
        //
        //    P = 
        //      1 0 1 1
        //      1 1 1 0
        //      0 1 1 0
        //
        //    PIN =
        //      1 0 2 3
        //      4 5 6 0
        //      0 7 8 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array that
        //    represents the polyomino.
        //
        //    Input, int P[M*N], the polyomino.  It is assumed that every entry
        //    is a 0 or a 1.
        //
        //    Output, int *PIN[M*N], the index of each nonzero entry.
        //
    {
        int i;

        int[] pin = new int[m * n];

        int k = 0;
        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                if (p[i + j * m] != 0)
                {
                    k += 1;
                    pin[i + j * m] = k;
                }
                else
                {
                    pin[i + j * m] = 0;
                }
            }
        }

        return pin;
    }

    public static void polyomino_lp_write(string filename, string label, int m, int n, int[] a,
            int[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_LP_WRITE writes an LP file for the polyomino problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILENAME, the output filename.
        //
        //    Input, string LABEL, the problem title.
        //
        //    Input, int M, the number of equations
        //
        //    Input, int N, the number of variables.
        //
        //    Input, int A[M*N], the coefficients.
        //
        //    Input, int B[M], the right hand sides.
        //
    {
        int i;
        int j;
        List<string> output = new();
        string cout;

        output.Add(label + "");
        output.Add("");

        output.Add("Maximize");
        output.Add("  Obj: 0");

        output.Add("Subject to");

        for (i = 0; i < m; i++)
        {
            cout = "";
            bool first = true;
            for (j = 0; j < n; j++)
            {
                if (a[i + j * m] == 0)
                {
                    continue;
                }

                switch (a[i + j * m])
                {
                    case < 0:
                        cout += " -";
                        break;
                    default:
                    {
                        switch (first)
                        {
                            case false:
                                cout += " +";
                                break;
                        }

                        break;
                    }
                }

                cout += Math.Abs(a[i + j * m]) switch
                {
                    1 => " x" + j + 1,
                    _ => " " + Math.Abs(a[i + j * m]) + " x" + j + 1
                };

                first = false;
            }

            output.Add(cout + " = " + b[i] + "");
        }

        output.Add("Binary");
        cout = " ";
        for (j = 0; j < n; j++)
        {
            cout += " x" + j + 1;
        }

        output.Add(cout);

        output.Add("End");

        try
        {
            File.WriteAllLines(filename, output);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("POLYOMINO_LP_WRITE - Error!");
            Console.WriteLine("  Could not open the output file.");
        }
    }

    public static void polyomino_print(int m, int n, int[] p, string label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_PRINT prints a polyomino.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the representation
        //    of the polyomino P.
        //
        //    Input, int P[M*N], a matrix of 0's and 1's representing the 
        //    polyomino.  The matrix should be "tight", that is, there should be a
        //    1 in row 1, and in column 1, and in row M, and in column N.
        //
        //    Input, string LABEL, a title for the polyomino.
        //
    {
        Console.WriteLine("");
        Console.WriteLine(label + "");
        Console.WriteLine("");
        if (m <= 0 || n <= 0)
        {
            Console.WriteLine("  [ Null matrix ]");
        }
        else
        {
            string cout = "";
            int i;
            for (i = 0; i < m; i++)
            {
                int j;
                for (j = 0; j < n; j++)
                {
                    cout += " " + p[i + j * m];
                }

                Console.WriteLine(cout);
                cout = "";
            }
        }
    }

    public static void polyomino_transform(int m, int n, int[] p, int rotate, int reflect,
            ref int mq, ref int nq, ref int[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_TRANSFORM transforms a polyomino.
        //
        //  Discussion:
        //
        //    A polyomino can be rotated or reflected.
        //
        //    This program is given a polyomino and returns the resulting polyomino
        //    after the specified reflection and rotation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the representation
        //    of the polyomino P.
        //
        //    Input, int P[M*N], a matrix of 0's and 1's representing the 
        //    polyomino.  The matrix should be "tight", that is, there should be a
        //    1 in row 1, and in column 1, and in row M, and in column N.
        //
        //    Input, int ROTATE, is 0, 1, 2, or 3, the number of 90 degree
        //    counterclockwise rotations to be applied.
        //
        //    Input, int REFLECT, is 0 or 1.  If it is 1, then each row of the
        //    polyomino matrix is to be reflected before any rotations are performed.
        //
        //    Output, int &MQ, &NQ, the number of rows and columns of the
        //    representation of the transformed polyomino
        //
        //    Output, int Q[MQ*NQ], the representation of the transformed
        //    polyomino.
        //
    {
        int k;

        mq = m;
        nq = n;

        reflect %= 2;

        q = typeMethods.i4mat_copy_new(m, n, p);

        switch (reflect)
        {
            case 1:
                typeMethods.i4mat_flip_cols(mq, nq, ref q);
                break;
        }

        rotate %= 4;

        for (k = 1; k <= rotate; k++)
        {
            typeMethods.i4mat_transpose(mq, nq, ref q);
            int r = mq;
            int s = nq;
            mq = s;
            nq = r;
            typeMethods.i4mat_flip_rows(mq, nq, ref q);
        }
    }
}