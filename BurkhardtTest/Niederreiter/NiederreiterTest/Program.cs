using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.NiederreiterNS;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace NiederreiterTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for NIEDERREITER_TEST.
            //
            //  Discussion:
            //
            //    NIEDERREITER_TEST tests the NIEDERREITER library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            gfarit();
            gfplys();
            niederreiter();
        }


        static void gfplys()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for GFPLYS.
            //
            //  Discussion:
            //
            //    GFPLYS writes out data about irreducible polynomials.
            //
            //    The program calculates irreducible polynomials for various
            //    finite fields, and writes them out to the file "gfplys.txt".
            //
            //    Finite field arithmetic is carried out with the help of
            //    precalculated addition and multiplication tables found on
            //    the file "gfarit.txt".  This file should have been computed
            //    and written by the program GFARIT.
            //
            //    The format of the irreducible polynomials on the output file is
            //
            //      Q
            //      d1   a(1)  a(2) ... a(d1)
            //      d2   b(1)  b(2) ... b(d2)
            //      ...
            //
            //    where 
            //
            //      Q is the order of the field, 
            //      d1 is the degree of the first irreducible polynomial, 
            //      a(1), a(2), ..., a(d1) are its coefficients.
            //
            //    Polynomials stored as arrays have the coefficient of degree N in 
            //    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
            //    DEG is just to remind us of this last fact.  A polynomial which is
            //    identically 0 is given degree -1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 September 2007
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
            //    Harald Niederreiter.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Paul Bratley, Bennett Fox, Harald Niederreiter,
            //    Algorithm 738: 
            //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 20, Number 4, 1994, pages 494-495.
            //
        {
            string output_filename = "gfplys.txt";
            List<string> output = new List<string>();
            int q_init;
            Polynomial.PLY ply = new Polynomial.PLY();

            Console.WriteLine("");
            Console.WriteLine("GFPLYS:");
            Console.WriteLine("  C++ version");
            Console.WriteLine("");
            Console.WriteLine("  A program to compute a set of irreducible");
            Console.WriteLine("  polynomials over fields of certain orders Q.");
            Console.WriteLine("");

            for (q_init = 2; q_init <= Polynomial.PLY.Q_MAX; q_init++)
            {
                irred(ref ply, ref output, q_init);
            }

            try
            {
                File.WriteAllLines(output_filename, output);
            }
            catch (Exception e)
            {
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("GFPLYS:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static int find(int n, int[] tab, int i, int tab_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIND seeks the value N in the range TAB(I) to TAB(TAB_MAX).
            //
            //  Discussion:
            //
            //    The vector TAB does not have to be sorted or have any other
            //    special properties.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 September 2007
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
            //    Harald Niederreiter.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Paul Bratley, Bennett Fox, Harald Niederreiter,
            //    Algorithm 738: 
            //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 20, Number 4, 1994, pages 494-495.
            //
            //  Parameters:
            //
            //    Input, int N, the value being sought.
            //
            //    Input, int TAB[], the table to be searched.
            //
            //    Input, int I, TAB_MAX, the first and last entries of
            //    TAB to be examined.
            //
            //    Output, int FIND, is the index ( between I and TAB_MAX) of the 
            //    entry in TAB that is equal to N, or else -1 if no such value
            //    was found.
            //
        {
            int j;
            int value;

            value = -1;

            if (tab[tab_max - 1] < n)
            {
                return value;
            }

            for (j = i; j <= tab_max; j++)
            {
                if (tab[j - 1] == n)
                {
                    value = j;
                    return value;
                }
            }

            return value;
        }

        static void irred(ref Polynomial.PLY ply, ref List<string> output, int q_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    IRRED computes and writes out a set of irreducible polynomials.
            //
            //  Discussion:
            //
            //    We find the irreducible polynomials using a sieve.  
            //
            //    Polynomials stored as arrays have the coefficient of degree n in 
            //    POLY(N), and the degree of the polynomial in POLY(-1).  The parameter
            //    DEG is just to remind us of this last fact.  A polynomial which is
            //    identically 0 is given degree -1.
            //
            //    Note that the value of NPOL controls the number of polynomials
            //    computed, and hence the maximum spatial dimension for the
            //    subsequence Niederreiter sequences, JVB, 07 June 2010.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2010
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
            //    Harald Niederreiter.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Paul Bratley, Bennett Fox, Harald Niederreiter,
            //    Algorithm 738: 
            //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 20, Number 4, 1994, pages 494-495.
            //
            //  Parameters:
            //
            //    Input, ofstream &OUTPUT, a reference to the output stream.
            //
            //    Input, int Q_INIT, the order of the field.
            //
            //  Local Parameters:
            //
            //    Local, int SIEVE_MAX, the size of the sieve.  
            //
            //    Array MONPOL holds monic polynomials.
            //
            //    Array SIEVE says whether the polynomial is still OK.
            //
            //    Local, int NPOLS, the number of irreducible polynomials to
            //    be calculated for a given field.
            //
        {
            int SIEVE_MAX = 400;

            int i;
            int j;
            int k;
            int l;
            int[] monpol = new int[SIEVE_MAX];
            int n;
            int npols = 50;
            int[] pi;
            int[] pj;
            int[] pk;
            bool[] sieve = new bool[SIEVE_MAX];

            if (q_init <= 1 || Polynomial.PLY.Q_MAX < q_init)
            {
                Console.WriteLine("");
                Console.WriteLine("IRRED - Fatal error!");
                Console.WriteLine("  Bad value of Q = " + q_init + "");
                return;
            }

            ply.P = typeMethods.i4_characteristic(q_init);
            //
            //  If no field of order Q_INIT exists, there is nothing to do.
            //
            if (ply.P <= 0)
            {
                return;
            }

            Console.WriteLine("  IRRED setting up case for Q = " + q_init + "");
            //
            //  Set up the field arithmetic tables.
            //  (Note that SETFLD sets Q = q_init!)
            //
            setfld(ref ply, q_init);
            //
            //  Set up the sieve containing only monic polynomials.
            //
            i = 0;
            j = 1;
            k = ply.Q;

            for (n = 1; n <= SIEVE_MAX; n++)
            {
                i = i + 1;

                if (i == j)
                {
                    i = k;
                    j = 2 * k;
                    k = ply.Q * k;
                }

                monpol[n - 1] = i;
                sieve[n - 1] = true;
            }

            //
            //  Write out the irreducible polynomials as they are found.
            //
            n = 0;
            output.Add(ply.Q.ToString().PadLeft(3) + "");

            for (i = 1; i <= SIEVE_MAX; i++)
            {
                if (sieve[i - 1])
                {
                    pi = Polynomial.itop(ref ply, monpol[i - 1]);
                    k = pi[0];
                    string cout = k.ToString().PadLeft(3);
                    for (l = 0; l <= k; l++)
                    {
                        cout += pi[l + 1].ToString().PadLeft(3);
                    }

                    output.Add(cout);
                    n = n + 1;

                    if (n == npols)
                    {
                        return;
                    }

                    for (j = i; j <= SIEVE_MAX; j++)
                    {
                        pj = Polynomial.itop(ref ply, monpol[j - 1]);
                        pk = Polynomial.plymul(ref ply, pi, pj);

                        k = find(Polynomial.ptoi(ref ply, pk), monpol, j, SIEVE_MAX);

                        if (k != -1)
                        {
                            sieve[k - 1] = false;
                        }
                    }
                }
            }

            Console.WriteLine("");
            Console.WriteLine("IRRED - Warning!");
            Console.WriteLine("  The sieve size SIEVE_MAX is too small.");
            Console.WriteLine("  Number of irreducible polynomials found: " + n + "");
            Console.WriteLine("  Number needed: " + npols + "");
        }


        static void gfarit()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for GFARIT.
            //
            //  Discussion:
            //
            //    GFARIT writes the arithmetic tables called "gfarit.txt".
            //
            //    The program calculates addition and multiplication tables
            //    for arithmetic in finite fields, and writes them out to
            //    the file "gfarit.txt".  Tables are only calculated for fields
            //    of prime-power order Q, the other cases being trivial.
            //
            //    For each value of Q, the file contains first Q, then the
            //    addition table, and lastly the multiplication table.
            //
            //    After "gfarit.txt" has been set up, run GFPLYS to set up 
            //    the file "gfplys.txt".  That operation requires reading 
            //    "gfarit.txt".  
            //
            //    The files "gfarit.txt" and "gfplys.txt" should be saved 
            //    for future use.  
            //
            //    Thus, a user needs to run GFARIT and GFPLYS just once,
            //    before running the set of programs associated with GENIN.  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 September 2007
            //
            //  Author:
            //
            //    Paul Bratley, Bennet Fox, Harald Niederreiter.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Paul Bratley, Bennett Fox, Harald Niederreiter,
            //    Algorithm 738: 
            //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 20, Number 4, 1994, pages 494-495.
            //
        {
            string output_filename = "gfarit.txt";
            List<string> output = new List<string>();
            int q_init;
            Polynomial.PLY ply = new Polynomial.PLY();


            Console.WriteLine("");
            Console.WriteLine("GFARIT:");
            Console.WriteLine("");
            Console.WriteLine("  A program which computes a set of arithmetic");
            Console.WriteLine("  tables, and writes them to a file.");
            Console.WriteLine("");
            Console.WriteLine("  Tables will be created for fields of prime or prime power order");
            Console.WriteLine("  Q between 2 and " + Polynomial.PLY.Q_MAX + ".");
            Console.WriteLine("");


            for (q_init = 2; q_init <= Polynomial.PLY.Q_MAX; q_init++)
            {
                gftab(ref ply, ref output, q_init);
            }

            try
            {
                File.WriteAllLines(output_filename, output);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("GFARIT - Fatal error!");
                Console.WriteLine("  Could not open the output file: \"" + output_filename + "\"");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("GFARIT:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static void niederreiter()
        {
            int base_;
            int dim_num;

            Console.WriteLine("");
            Console.WriteLine("NIEDERREITER_TEST");
            Console.WriteLine("  Test the NIEDERREITER routines.");

            base_ = 2;
            test01(base_);

            base_ = 3;
            test01(base_);

            base_ = 13;
            test01(base_);

            base_ = 2;
            test02(base_);

            base_ = 3;
            test02(base_);

            base_ = 2;
            dim_num = 20;
            test03(base_, dim_num);

            base_ = 2;
            dim_num = 29;
            test03(base_, dim_num);

            Console.WriteLine("");
            Console.WriteLine("NIEDERREITER_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void gftab(ref Polynomial.PLY ply, ref List<string> output, int q_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GFTAB computes and writes data for a particular field size Q_INIT.
            //
            //  Discussion:
            //
            //    A polynomial with coefficients A(*) in the field of order Q
            //    can also be stored in an integer I, with
            //
            //      I = AN*Q**N + ... + A0.
            //
            //    Polynomials stored as arrays have the
            //    coefficient of degree n in POLY(N), and the degree of the
            //    polynomial in POLY(-1).  The parameter DEG is just to remind
            //    us of this last fact.  A polynomial which is identically 0
            //    is given degree -1.
            //
            //    IRRPLY holds irreducible polynomials for constructing
            //    prime-power fields.  IRRPLY(-2,I) says which field this
            //    row is used for, and then the rest of the row is a
            //    polynomial (with the degree in IRRPLY(-1,I) as usual).
            //    The chosen irreducible poly is copied into MODPLY for use.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 September 2007
            //
            //  Author:
            //
            //    Paul Bratley, Bennet Fox, Harald Niederreiter.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Paul Bratley, Bennett Fox, Harald Niederreiter,
            //    Algorithm 738: 
            //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 20, Number 4, 1994, pages 494-495.
            //
            //  Parameters:
            //
            //    Input, ofstream &OUTPUT, a reference to the output stream.
            //
            //    Input, int Q_INIT, the order of the field for which the
            //    addition and multiplication tables are needed.
            //
        {
            int[,] gfadd = new int [Polynomial.PLY.Q_MAX, Polynomial.PLY.Q_MAX];
            int[,] gfmul = new int [Polynomial.PLY.Q_MAX, Polynomial.PLY.Q_MAX];
            int i;
            int[,] irrply =
            {
                {4, 2, 1, 1, 1, 0, 0, 0},
                {8, 3, 1, 1, 0, 1, 0, 0},
                {9, 2, 1, 0, 1, 0, 0, 0},
                {16, 4, 1, 1, 0, 0, 1, 0},
                {25, 2, 2, 0, 1, 0, 0, 0},
                {27, 3, 1, 2, 0, 1, 0, 0},
                {32, 5, 1, 0, 1, 0, 0, 1},
                {49, 2, 1, 0, 1, 0, 0, 0}
            };
            int j;
            int[] modply = new int[Polynomial.PLY.DEG_MAX + 2];
            int[] pi;
            int[] pj;
            int[] pk;
            int[] pl;

            if (q_init <= 1 || Polynomial.PLY.Q_MAX < q_init)
            {
                Console.WriteLine("");
                Console.WriteLine("GFTAB - Fatal error!");
                Console.WriteLine("  Bad value of Q_INIT.");
                return;
            }

            ply.P = typeMethods.i4_characteristic(q_init);
            //
            //  If QIN is not a prime power, we are not interested.
            //
            if (ply.P == 0 || ply.P == q_init)
            {
                return;
            }

            Console.WriteLine("  GFTAB computing table for Q = " + q_init
                                                                 + "  with characteristic P = " + ply.P + ".");
            //
            //  Otherwise, we set up the elements of the common /FIELD/
            //  ready to do arithmetic mod P, the characteristic of Q_INIT.
            //
            setfld(ref ply, q_init);
            //
            //  Next find a suitable irreducible polynomial and copy it to array MODPLY.
            //
            i = 1;

            while (irrply[i - 1, -2 + 2] != q_init)
            {
                i = i + 1;
            }

            for (j = -1; j <= irrply[i - 1, -1 + 2]; j++)
            {
                modply[j + 1] = irrply[i - 1, j + 2];
            }

            for (j = irrply[i - 1, -1 + 2] + 1; j <= Polynomial.PLY.DEG_MAX; j++)
            {
                modply[j + 1] = 0;
            }

            //
            //  Deal with the trivial cases.
            //
            for (i = 0; i < q_init; i++)
            {
                gfadd[i, 0] = i;
                gfadd[0, i] = i;
                gfmul[i, 0] = 0;
                gfmul[0, i] = 0;
            }

            for (i = 1; i < q_init; i++)
            {
                gfmul[i, 1] = i;
                gfmul[1, i] = i;
            }

            //
            //  Now deal with the rest.  Each integer from 1 to Q-1
            //  is treated as a polynomial with coefficients handled mod P.
            //  Multiplication of polynomials is mod MODPLY.
            //
            pl = new int[Polynomial.PLY.DEG_MAX + 2];

            for (i = 1; i < q_init; i++)
            {
                pi = Polynomial.itop(i, ply.P);

                for (j = 1; j <= i; j++)
                {
                    pj = Polynomial.itop(j, ply.P);
                    pk = Polynomial.plyadd(ref ply, pi, pj);
                    gfadd[i, j] = Polynomial.ptoi(pk, ply.P);
                    gfadd[j, i] = gfadd[i, j];

                    if (1 < i && 1 < j)
                    {
                        pk = Polynomial.plymul(ref ply, pi, pj);
                        Polynomial.plydiv(ref ply, pk, modply, ref pj, ref pl);
                        gfmul[i, j] = Polynomial.ptoi(pl, ply.P);
                        gfmul[j, i] = gfmul[i, j];
                    }
                }
            }

            //
            //  Write out the tables.
            //
            output.Add(" " + q_init + "");

            for (i = 0; i < q_init; i++)
            {
                string cout = "";
                for (j = 0; j < q_init; j++)
                {
                    cout += " " + gfadd[i, j];
                }

                output.Add(cout);
            }

            for (i = 0; i < q_init; i++)
            {
                string cout = "";
                for (j = 0; j < q_init; j++)
                {
                    cout += " " + gfmul[i, j];
                }

                output.Add(cout);
            }
        }

        static void setfld(ref Polynomial.PLY ply, int q_init)

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    SETFLD sets up the arithmetic tables for a finite field.
            //
            //  Discussion:
            //
            //    This subroutine sets up addition, multiplication, and
            //    subtraction tables for the finite field of order QIN.
            //
            //    A polynomial with coefficients A(*) in the field of order Q
            //    can also be stored in an integer I, with
            //
            //      I = AN*Q**N + ... + A0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 September 2007
            //
            //  Author:
            //
            //    Paul Bratley, Bennet Fox, Harald Niederreiter.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Paul Bratley, Bennett Fox, Harald Niederreiter,
            //    Algorithm 738: 
            //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
            //    ACM Transactions on Mathematical Software,
            //    Volume 20, Number 4, 1994, pages 494-495.
            //
            //  Parameters:
            //
            //    Input, int Q_INIT, the order of the field.
            //
        {
            int i;
            string[] input;
            string input_filename = "gfarit.txt";
            int j;
            int n;

            if (q_init <= 1 || Polynomial.PLY.Q_MAX < q_init)
            {
                Console.WriteLine("");
                Console.WriteLine("SETFLD - Fatal error!");
                Console.WriteLine("  Bad value of Q = " + q_init + "");
                return;
            }

            ply.Q = q_init;
            ply.P = typeMethods.i4_characteristic(ply.Q);

            if (ply.P == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("SETFLD - Fatal error!");
                Console.WriteLine("  There is no field of order Q = " + ply.Q + "");
                return;
            }

            //
            //  Set up to handle a field of prime or prime-power order.
            //  Calculate the addition and multiplication tables.
            //
            for (i = 0; i < ply.P; i++)
            {
                for (j = 0; j < ply.P; j++)
                {
                    ply.add[i, j] = (i + j) % ply.P;
                }
            }

            for (i = 0; i < ply.P; i++)
            {
                for (j = 0; j < ply.P; j++)
                {
                    ply.mul[i, j] = (i * j) % ply.P;
                }
            }

            //
            //  Use the addition table to set the subtraction table.
            //
            for (i = 0; i < ply.P; i++)
            {
                for (j = 0; j < ply.P; j++)
                {
                    ply.sub[ply.add[i, j], i] = j;
                }
            }
        }


        static void test01(int base_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests NIEDERREITER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int BASE, the base_ to use in the computation.
            //    BASE should be a prime, or a power of a prime.
            //
        {
            const int dim_max = 4;

            int dim;
            int dim_num;
            int i;
            double[] r = new double[dim_max];
            int seed;
            int seed_in;
            int seed_out;
            Niederreiter.NiederReiterData data = new Niederreiter.NiederReiterData();

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  NIEDERREITER computes the next element of ");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we call NIEDERREITER repeatedly.");
            Console.WriteLine("");
            Console.WriteLine("  Using base_ BASE =      " + base_ + "");

            for (dim_num = 2; dim_num <= dim_max; dim_num++)
            {
                seed = 0;

                Console.WriteLine("");
                Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");
                Console.WriteLine("");
                Console.WriteLine("    Seed    Seed     Niederreiter");
                Console.WriteLine("      In     Out");
                Console.WriteLine("");
                for (i = 0; i <= 110; i++)
                {
                    seed_in = seed;
                    Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                    seed_out = seed;
                    if (i <= 11 || 95 <= i)
                    {
                        string cout = "  " + seed_in.ToString().PadLeft(8)
                                           + "  " + seed_out.ToString().PadLeft(8);
                        for (dim = 0; dim < dim_num; dim++)
                        {
                            cout += "  " + r[dim].ToString().PadLeft(10);
                        }

                        Console.WriteLine(cout);
                    }
                    else if (i == 12)
                    {
                        Console.WriteLine("......................");
                    }
                }
            }

        }

        static void test02(int base_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests NIEDERREITER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int BASE, the base_ to use in the computation.
            //    BASE should be a prime, or a power of a prime.
            //
        {
            const int dim_num = 3;

            int dim;
            int i;
            double[] r = new double[dim_num];
            int seed;
            int seed_in;
            int seed_out;
            string cout = "";
            Niederreiter.NiederReiterData data = new Niederreiter.NiederReiterData();

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  NIEDERREITER computes the next element of");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we demonstrate how the SEED can be");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");

            Console.WriteLine("");
            Console.WriteLine("  Using base_ BASE =           " + base_ + "");
            Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");

            seed = 0;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                            + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + r[dim].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump ahead by increasing SEED:");
            Console.WriteLine("");

            seed = 100;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            for (i = 1; i <= 5; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                            + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + r[dim].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump back by decreasing SEED:");
            Console.WriteLine("");

            seed = 3;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                            + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + r[dim].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump ahead by increasing SEED:");
            Console.WriteLine("");

            seed = 98;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            cout = "";
            for (i = 1; i <= 5; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                            + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + r[dim].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

        }

        static void test03(int base_, int dim_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests NIEDERREITER.
            //
            //  Discussion:
            //
            //    Simply verify that a few terms of a sequence of given dimension
            //    can be computed.  Most recently, the NIEDERREITER code was set
            //    up to handle up to dimension 50...we think.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int BASE, the base_ to use in the computation.
            //    BASE should be a prime, or a power of a prime.
            //
            //    Input, int DIM, the spatial dimension.
            //
        {
            int dim;
            int i;
            double[] r;
            int seed;
            int seed_in;
            int seed_out;
            Niederreiter.NiederReiterData data = new Niederreiter.NiederReiterData();

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  NIEDERREITER computes the next element of");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply generate ten elements in a given base_");
            Console.WriteLine("  and dimension.");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");

            Console.WriteLine("");
            Console.WriteLine("  Using base_ BASE =           " + base_ + "");
            Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");

            seed = 0;
            r = new double[dim_num];

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                            + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + r[dim].ToString().PadLeft(10);
                    if (((dim + 1) % 5 == 0) && (dim != dim_num))
                    {
                        Console.WriteLine(cout);
                        cout = "                    ";
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }
}