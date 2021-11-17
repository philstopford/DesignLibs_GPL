using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.IO;

public static class HB
{
    public static void hb_exact_read(string[] input, ref int inputIndex, int nrow, int nrhs, int rhscrd,
            char[] rhsfmt, char[] rhstyp, ref double[] exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_EXACT_READ reads the exact solution vectors in an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, the unit from which data is read.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int RHSCRD, the number of lines in the file for
        //    right hand sides.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    Output, double EXACT[NROW*NRHS], the exact solution vectors.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int khi;
        int klo;
        char[] line = new char[255];
        int line_num;
        int m = 0;
        int r = 0;
        char[] s;
        int w = 0;

        switch (rhscrd)
        {
            case > 0:
            {
                switch (rhstyp[2])
                {
                    case 'X':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);

                        line_num = 1 + (nrow * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            line = input[inputIndex].ToCharArray();
                            inputIndex++;
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrow * nrhs);

                            khi = 0;
                            for (j = jlo; j <= jhi; j++)
                            {
                                klo = khi + 1;
                                khi = Math.Min(klo + w - 1, line.Length);
                                s = typeMethods.s_substring(line, klo, khi);
                                exact[j - 1] = Convert.ToDouble(s);
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static void hb_exact_write(ref List<string> output, int nrow, int nrhs, int rhscrd,
            char[] rhsfmt, char[] rhstyp, double[] exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_EXACT_WRITE writes the exact solution vectors to an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, List<string> &OUTPUT, the unit to which data is written.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int RHSCRD, the number of lines in the file for
        //    right hand sides.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    Input, double EXACT[NROW*NRHS], the exact solution vectors.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int line_num;
        int m = 0;
        int r = 0;
        int w = 0;

        switch (rhscrd)
        {
            case > 0:
            {
                switch (rhstyp[2])
                {
                    case 'X':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);
                        line_num = 1 + (nrow * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            string cout = "";
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrow * nrhs);
                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += exact[j - 1].ToString().PadLeft(w);
                            }

                            output.Add(cout);
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static void hb_file_read(string inputFile, ref char[] title, ref char[] key, ref int totcrd,
            ref int ptrcrd, ref int indcrd, ref int valcrd, ref int rhscrd, ref char[] mxtype, ref int nrow,
            ref int ncol, ref int nnzero, ref int neltvl, ref char[] ptrfmt, ref char[] indfmt, ref char[] valfmt,
            ref char[] rhsfmt, ref char[] rhstyp, ref int nrhs, ref int nrhsix, ref int[] colptr,
            ref int[] rowind, ref double[] values, ref double[] rhsval, ref int[] rhsptr, ref int[] rhsind,
            ref double[] rhsvec, ref double[] guess, ref double[] exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_FILE_READ reads an HB file.
        //
        //  Discussion:
        //
        //    This routine reads all the information from an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, the unit from which the data is read.
        //
        //    Output, char *TITLE, a 72 character title for the matrix.
        //
        //    Output, char *KEY, an 8 character identifier for the matrix.
        //
        //    Output, int *TOTCRD, the total number of lines of data.
        //
        //    Output, int *PTRCRD, the number of input lines for pointers.
        //
        //    Output, int *INDCRD, the number of input lines for row indices.
        //
        //    Output, int *VALCRD, the number of input lines for numerical values.
        //
        //    Output, int *RHSCRD, the number of input lines for right hand sides.
        //
        //    Output, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Output, int *NROW, the number of rows or variables.
        //
        //    Output, int *NCOL, the number of columns or elements.
        //
        //    Output, int *NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Output, int *NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Output, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Output, char *INDFMT, the 16 character format for reading indices.
        //
        //    Output, char *VALFMT, the 20 character format for reading values.
        //
        //    Output, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Output, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //
        //    Output, int *NRHS, the number of right hand sides.
        //
        //    Output, int *NRHSIX, the number of entries of storage for right
        //    hand side values, in the case where RHSTYP[0] = 'M' and
        //    MXTYPE[2] = 'A'.
        //
        //    Output, int COLPTR[NCOL+1], COLPTR[I-1] points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Output, int ROWIND[NNZERO], the row index of each item.
        //
        //    If MXTYPE[2] == 'F':
        //
        //      Output, int ROWIND[NELTVL], the row index of each item.
        //
        //    If RHSTYP[0] == 'F':
        //
        //      Output, double RHSVAL[NROW*NRHS], contains NRHS dense right hand
        //      side vectors.
        //
        //      Output, int RHSPTR[], is not used.
        //
        //      Output, int RHSIND[], is not used.
        //
        //      Output, int RHSVEC[], is not used.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'A':
        //
        //      Output, double RHSVAL[], is not used.
        //
        //      Output, int RHSPTR[NRHS+1], RHSPTR[I-1] points to the location of
        //      the first entry of right hand side I in the sparse right hand
        //      side vector.
        //
        //      Output, int RHSIND[NRHSIX], indicates, for each entry of
        //      RHSVEC, the corresponding row index.
        //
        //      Output, double RHSVEC[NRHSIX], contains the value of the right hand
        //      side entries.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'E':
        //
        //      Output, double RHSVAL[NNZERO*NRHS], contains NRHS unassembled
        //      finite element vector right hand sides.
        //
        //      Output, int RHSPTR[], is not used.
        //
        //      Output, int RHSIND[], is not used.
        //
        //      Output, double RHSVEC[], is not used.
        //
        //    Output, double GUESS[NROW*NRHS], the starting guess vectors.
        //
        //    Output, double EXACT[NROW*NRHS], the exact solution vectors.
        //
    {
        //
        //  Read the header block.
        //
        int inputIndex = 0;
        string[] input;
        try
        {
            input = File.ReadAllLines(inputFile);
        }
        catch (Exception e)
        {
            Console.WriteLine(e);
            return;
        }

        hb_header_read(input, ref inputIndex, ref title, ref key, ref totcrd, ref ptrcrd, ref indcrd,
            ref valcrd, ref rhscrd, ref mxtype, ref nrow, ref ncol, ref nnzero, ref neltvl, ref ptrfmt, ref indfmt,
            ref valfmt, ref rhsfmt, ref rhstyp, ref nrhs, ref nrhsix);
        colptr = ptrcrd switch
        {
            //
            //  Read the matrix structure.
            //
            > 0 => new int[ncol + 1],
            _ => colptr
        };

        switch (indcrd)
        {
            case > 0:
                switch (mxtype[2])
                {
                    case 'A':
                        rowind = new int[nnzero];
                        break;
                    case 'E':
                        rowind = new int[neltvl];
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("HB_FILE_READ - Fatal error!");
                        Console.WriteLine("  Illegal value of MXTYPE character 3!");
                        return;
                }

                break;
        }

        hb_structure_read(input, ref inputIndex, ncol, mxtype, nnzero, neltvl,
            ptrcrd, ptrfmt, indcrd, indfmt, ref colptr, ref rowind);
        switch (valcrd)
        {
            //
            //  Read the matrix values.
            //
            case > 0:
                switch (mxtype[2])
                {
                    case 'A':
                        values = new double[nnzero];
                        break;
                    case 'E':
                        values = new double[neltvl];
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("HB_FILE_READ - Fatal error!");
                        Console.WriteLine("  Illegal value of MXTYPE character 3!");
                        return;
                }

                hb_values_read(input, ref inputIndex, valcrd, mxtype, nnzero, neltvl,
                    valfmt, ref values);
                break;
        }

        switch (rhscrd)
        {
            //
            //  Read the right hand sides.
            //
            case > 0:
            {
                switch (rhstyp[0])
                {
                    case 'F':
                        rhsval = new double[nrow * nrhs];
                        break;
                    case 'M' when mxtype[2] == 'A':
                        rhsptr = new int[nrhs + 1];


                        rhsind = new int[nrhsix];


                        rhsvec = new double[nrhsix];
                        break;
                    case 'M' when mxtype[2] == 'E':
                        rhsval = new double[nnzero * nrhs];
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("HB_FILE_READ - Fatal error!");
                        Console.WriteLine("  Illegal combination of RHSTYP character 1");
                        Console.WriteLine("  and MXTYPE character 3!");
                        return;
                }

                hb_rhs_read(input, ref inputIndex, nrow, nnzero, nrhs, nrhsix,
                    rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, ref rhsval,
                    ref rhsind, ref rhsptr, ref rhsvec);
                switch (rhstyp[1])
                {
                    //
                    //  Read the starting guesses.
                    //
                    case 'G':
                        guess = new double[nrow * nrhs];

                        hb_guess_read(input, ref inputIndex, nrow, nrhs, rhscrd, rhsfmt, rhstyp, ref guess);
                        break;
                }

                switch (rhstyp[2])
                {
                    //
                    //  Read the exact solutions.
                    //
                    case 'X':
                        exact = new double[nrow * nrhs];

                        hb_exact_read(input, ref inputIndex, nrow, nrhs, rhscrd, rhsfmt, rhstyp, ref exact);
                        break;
                }

                break;
            }
        }
    }

    public static void hb_file_write(string outputFile, char[] title, char[] key, int totcrd,
            int ptrcrd, int indcrd, int valcrd, int rhscrd, char[] mxtype, int nrow,
            int ncol, int nnzero, int neltvl, char[] ptrfmt, char[] indfmt, char[] valfmt,
            char[] rhsfmt, char[] rhstyp, int nrhs, int nrhsix, int[] colptr,
            int[] rowind, double[] values, double[] rhsval, int[] rhsptr, int[] rhsind,
            double[] rhsvec, double[] guess, double[] exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_FILE_WRITE writes an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ofsteam &OUTPUT, the unit to which data is written.
        //
        //    Input, char *TITLE, a 72 character title for the matrix.
        //
        //    Input, char *KEY, an 8 character identifier for the matrix.
        //
        //    Input, int TOTCRD, the total number of lines of data.
        //
        //    Input, int PTRCRD, the number of input lines for pointers.
        //
        //    Input, int INDCRD, the number of input lines for row indices.
        //
        //    Input, int VALCRD, the number of input lines for numerical values.
        //
        //    Input, int RHSCRD, the number of input lines for right hand sides.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NCOL, the number of columns or elements.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Input, char *INDFMT, the 16 character format for reading indices.
        //
        //    Input, char *VALFMT, the 20 character format for reading values.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int NRHSIX, the number of row indices (set to 0
        //    in the case of unassembled matrices.)  Ignored if NRHS = 0.
        //
        //    Input, int COLPTR[NCOL+1], COLPTR(I) points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Input, int ROWIND[NNZERO], the row index of each item.
        //
        //      Input, double VALUES[NNZERO], the nonzero values of the matrix.
        //
        //    If MXTYPE[2] == 'E':
        //
        //      Input, int ROWIND[NELTVL], the row index of each item.
        //
        //      Input, double VALUES[NELTVL], the nonzero values of the matrix.
        //
        //    If RHSTYP[0] == 'F':
        //
        //      Input, double RHSVAL[NROW*NRHS], contains NRHS dense right hand
        //      side vectors.
        //
        //      Input, int RHSPTR[], is not used.
        //
        //      Input, int RHSIND[], is not used.
        //
        //      Input, double RHSVEC[], is not used.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'A':
        //
        //      Input, double RHSVAL[], is not used.
        //
        //      Input, int RHSPTR[NRHS+1], RHSPTR(I) points to the location of
        //      the first entry of right hand side I in the sparse right hand
        //      side vector.
        //
        //      Input, int RHSIND[NRHSIX], indicates, for each entry of
        //      RHSVEC, the corresponding row index.
        //
        //      Input, double RHSVEC[NRHSIX], contains the value of the right hand
        //      side entries.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'E':
        //
        //      Input, double RHSVAL[NNZERO*NRHS], contains NRHS unassembled
        //      finite element vector right hand sides.
        //
        //      Input, int RHSPTR[], is not used.
        //
        //      Input, int RHSIND[], is not used.
        //
        //      Input, double RHSVEC[], is not used.
        //
        //    Input, double GUESS[NROW*NRHS], the starting guess vectors.
        //
        //    Input, double EXACT[NROW*NRHS], the exact solution vectors.
        //
    {
        List<string> output = new();
        //
        //  Write the header block.
        //
        hb_header_write(ref output, title, key, totcrd, ptrcrd, indcrd,
            valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt,
            valfmt, rhsfmt, rhstyp, nrhs, nrhsix);
        //
        //  Write the matrix structure.
        //
        hb_structure_write(ref output, ncol, mxtype, nnzero, neltvl,
            ptrfmt, indfmt, colptr, rowind);
        //
        //  Write the matrix values.
        //
        hb_values_write(ref output, valcrd, mxtype, nnzero, neltvl,
            valfmt, values);
        //
        //  Write the right hand sides.
        //
        hb_rhs_write(ref output, nrow, nnzero, nrhs, nrhsix,
            rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval,
            rhsind, rhsptr, rhsvec);
        //
        //  Write the starting guesses.
        //
        hb_guess_write(ref output, nrow, nrhs, rhscrd, rhsfmt, rhstyp, guess);
        //
        //  Write the exact solutions.
        //
        hb_exact_write(ref output, nrow, nrhs, rhscrd, rhsfmt, rhstyp, exact);

        try
        {
            File.WriteAllLines(outputFile, output);
        }
        catch (Exception e)
        {
            Console.WriteLine(e);
            throw;
        }

    }

    public static void hb_guess_read(string[] input, ref int inputIndex, int nrow, int nrhs, int rhscrd,
            char[] rhsfmt, char[] rhstyp, ref double[] guess)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_GUESS_READ reads the starting guess vectors in an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, the unit from which data is read.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int RHSCRD, the number of lines in the file for
        //    right hand sides.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    Output, double GUESS[NROW*NRHS], the starting guess vectors.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int khi;
        int klo;
        char[] line = new char[255];
        int line_num;
        int m = 0;
        int r = 0;
        char[] s;
        int w = 0;

        switch (rhscrd)
        {
            case > 0:
            {
                switch (rhstyp[1])
                {
                    case 'G':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);

                        line_num = 1 + (nrow * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            line = input[inputIndex].ToCharArray();
                            inputIndex++;
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrow * nrhs);

                            khi = 0;
                            for (j = jlo; j <= jhi; j++)
                            {
                                klo = khi + 1;
                                khi = Math.Min(klo + w - 1, line.Length);
                                s = typeMethods.s_substring(line, klo, khi);
                                guess[j - 1] = Convert.ToDouble(s);
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }

    }

    public static void hb_guess_write(ref List<string> output, int nrow, int nrhs, int rhscrd,
            char[] rhsfmt, char[] rhstyp, double[] guess)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_GUESS_WRITE writes the starting guess vectors to an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, List<string> &OUTPUT, the unit to which data is written.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int RHSCRD, the number of lines in the file for
        //    right hand sides.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    Input, double GUESS[NROW*NRHS], the starting guess vectors.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int line_num;
        int m = 0;
        int r = 0;
        int w = 0;

        switch (rhscrd)
        {
            case > 0:
            {
                switch (rhstyp[1])
                {
                    case 'G':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);
                        line_num = 1 + (nrow * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            string cout = "";
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrow * nrhs);

                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += guess[j - 1].ToString().PadLeft(w);
                            }

                            output.Add(cout);
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static void hb_header_print(char[] title, char[] key, int totcrd, int ptrcrd,
            int indcrd, int valcrd, int rhscrd, char[] mxtype, int nrow, int ncol,
            int nnzero, int neltvl, char[] ptrfmt, char[] indfmt, char[] valfmt,
            char[] rhsfmt, char[] rhstyp, int nrhs, int nrhsix)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_HEADER_PRINT prints the header of an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, char *TITLE, a 72 character title for the matrix.
        //
        //    Input, char *KEY, an 8 character identifier for the matrix.
        //
        //    Input, int TOTCRD, the total number of lines of data.
        //
        //    Input, int PTRCRD, the number of input lines for pointers.
        //
        //    Input, int INDCRD, the number of input lines for row indices.
        //
        //    Input, int VALCRD, the number of input lines for numerical values.
        //
        //    Input, int RHSCRD, the number of input lines for right hand sides.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NCOL, the number of columns or elements.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Input, char *INDFMT, the 16 character format for reading indices.
        //
        //    Input, char *VALFMT, the 20 character format for reading values.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int NRHSIX, the number of entries of storage for right
        //    hand side values, in the case where RHSTYP[0] = 'M' and
        //    MXTYPE[2] = 'A'.
        //
    {
        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        Console.WriteLine("  TOTCRD = " + totcrd + "");
        Console.WriteLine("  PTRCRD = " + ptrcrd + "");
        Console.WriteLine("  INDCRD = " + indcrd + "");
        Console.WriteLine("  VALCRD = " + valcrd + "");
        Console.WriteLine("  RHSCRD = " + rhscrd + "");
        Console.WriteLine("");
        Console.WriteLine("  KEY =    '" + key + "'.");
        Console.WriteLine("  MXTYPE = '" + mxtype + "'.");
        Console.WriteLine("  RHSTYP = '" + rhstyp + "'.");
        Console.WriteLine("");
        Console.WriteLine("  NROW =   " + nrow + "");
        Console.WriteLine("  NCOL =   " + ncol + "");
        Console.WriteLine("  NNZERO = " + nnzero + "");
        Console.WriteLine("  NELTVL = " + neltvl + "");
        Console.WriteLine("  NRHS =   " + nrhs + "");
        Console.WriteLine("  NRHSIX = " + nrhsix + "");
        Console.WriteLine("");
        Console.WriteLine("  PTRFMT = '" + ptrfmt + "'.");
        Console.WriteLine("  INDFMT = '" + indfmt + "'.");
        Console.WriteLine("  VALFMT = '" + valfmt + "'.");
        Console.WriteLine("  RHSFMT = '" + rhsfmt + "'.");
    }

    public static void hb_header_read(string[] input, ref int inputIndex, ref char[] title, ref char[] key,
            ref int totcrd,
            ref int ptrcrd, ref int indcrd, ref int valcrd, ref int rhscrd, ref char[] mxtype, ref int nrow,
            ref int ncol, ref int nnzero, ref int neltvl, ref char[] ptrfmt, ref char[] indfmt, ref char[] valfmt,
            ref char[] rhsfmt, ref char[] rhstyp, ref int nrhs, ref int nrhsix)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_HEADER_READ reads the header of an HB file.
        //
        //  Discussion:
        //
        //    The user should already have opened the file, and positioned it
        //    to the first record.
        //
        //    Thanks to Alexander Zilyakov for pointing out that the line 
        //      if ( 0 < rhscrd )
        //    should be 
        //      if ( 0 < *rhscrd )
        //    15 August 2016.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, the unit from which data is read.
        //
        //    Output, char *TITLE, a 72 character title for the matrix.
        //
        //    Output, char *KEY, an 8 character identifier for the matrix.
        //
        //    Output, int *TOTCRD, the total number of lines of data.
        //
        //    Output, int *PTRCRD, the number of input lines for pointers.
        //
        //    Output, int *INDCRD, the number of input lines for row indices.
        //
        //    Output, int *VALCRD, the number of input lines for numerical values.
        //
        //    Output, int *RHSCRD, the number of input lines for right hand sides.
        //
        //    Output, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Output, int *NROW, the number of rows or variables.
        //
        //    Output, int *NCOL, the number of columns or elements.
        //
        //    Output, int *NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Output, int *NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Output, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Output, char *INDFMT, the 16 character format for reading indices.
        //
        //    Output, char *VALFMT, the 20 character format for reading values.
        //
        //    Output, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Output, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //
        //    Output, int *NRHS, the number of right hand sides.
        //
        //    Output, int *NRHSIX, the number of entries of storage for right
        //    hand side values, in the case where RHSTYP[0] = 'M' and
        //    MXTYPE[2] = 'A'.
        //
    {
        char[] field;
        char[] line = new char[255];
        //
        //  Read line 1.
        //
        try
        {
            line = input[inputIndex].ToArray();
            inputIndex++;
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("HB_HEADER_READ - Fatal error!");
            Console.WriteLine("  I/O error reading header line " + inputIndex + ".");
            return;
        }

        title = typeMethods.s_substring(line, 1, 72);
        typeMethods.s_trim(ref title);

        key = typeMethods.s_substring(line, 73, 80);
        typeMethods.s_trim(ref key);
        //
        //  Read line 2.
        //
        try
        {
            line = input[inputIndex].ToArray();
            inputIndex++;
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("HB_HEADER_READ - Fatal error!");
            Console.WriteLine("  I/O error reading header line " + inputIndex + ".");
            return;
        }

        field = typeMethods.s_substring(line, 1, 14);
        totcrd = Convert.ToInt32(field);

        field = typeMethods.s_substring(line, 15, 28);
        ptrcrd = Convert.ToInt32(field);

        field = typeMethods.s_substring(line, 29, 42);
        indcrd = Convert.ToInt32(field);

        field = typeMethods.s_substring(line, 43, 56);
        valcrd = Convert.ToInt32(field);

        field = typeMethods.s_substring(line, 57, 70);
        rhscrd = Convert.ToInt32(field);
        //
        //  Read line 3.
        //
        try
        {
            line = input[inputIndex].ToArray();
            inputIndex++;
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("HB_HEADER_READ - Fatal error!");
            Console.WriteLine("  I/O error reading header line " + inputIndex + ".");
            return;
        }

        mxtype = typeMethods.s_substring(line, 1, 3);
        typeMethods.s_trim(ref mxtype);

        field = typeMethods.s_substring(line, 15, 28);
        nrow = Convert.ToInt32(field);

        field = typeMethods.s_substring(line, 29, 42);
        ncol = Convert.ToInt32(field);

        field = typeMethods.s_substring(line, 43, 56);
        nnzero = Convert.ToInt32(field);

        field = typeMethods.s_substring(line, 57, 70);
        neltvl = Convert.ToInt32(field);
        //
        //  Read line 4.
        //
        try
        {
            line = input[inputIndex].ToArray();
            inputIndex++;
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("HB_HEADER_READ - Fatal error!");
            Console.WriteLine("  I/O error reading header line " + inputIndex + ".");
            return;
        }

        ptrfmt = typeMethods.s_substring(line, 1, 16);
        typeMethods.s_trim(ref ptrfmt);

        indfmt = typeMethods.s_substring(line, 17, 32);
        typeMethods.s_trim(ref indfmt);

        valfmt = typeMethods.s_substring(line, 33, 52);
        typeMethods.s_trim(ref valfmt);

        rhsfmt = typeMethods.s_substring(line, 53, 72);
        typeMethods.s_trim(ref rhsfmt);
        switch (rhscrd)
        {
            //
            //  Read line 5.
            //
            case > 0:
                try
                {
                    line = input[inputIndex].ToArray();
                    inputIndex++;
                }
                catch
                {
                    Console.WriteLine("");
                    Console.WriteLine("HB_HEADER_READ - Fatal error!");
                    Console.WriteLine("  I/O error reading header line " + inputIndex + ".");
                    return;
                }

                rhstyp = typeMethods.s_substring(line, 1, 3);
                typeMethods.s_trim(ref rhstyp);

                field = typeMethods.s_substring(line, 15, 28);
                nrhs = Convert.ToInt32(field);

                field = typeMethods.s_substring(line, 29, 42);
                nrhsix = Convert.ToInt32(field);
                break;
            default:
                rhstyp = null;
                nrhs = 0;
                nrhsix = 0;
                break;
        }

    }

    public static void hb_header_write(ref List<string> output, char[] title, char[] key, int totcrd,
            int ptrcrd, int indcrd, int valcrd, int rhscrd, char[] mxtype, int nrow,
            int ncol, int nnzero, int neltvl, char[] ptrfmt, char[] indfmt, char[] valfmt,
            char[] rhsfmt, char[] rhstyp, int nrhs, int nrhsix)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_HEADER_WRITE writes the header of an HB file.
        //
        //  Discussion:
        //
        //    The user should already have opened the file, and positioned it
        //    to the first record.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, List<string> &OUTPUT, the unit to which data is written.
        //
        //    Input, char *TITLE, a 72 character title for the matrix.
        //
        //    Input, char *KEY, an 8 character identifier for the matrix.
        //
        //    Input, int TOTCRD, the total number of lines of data.
        //
        //    Input, int PTRCRD, the number of input lines for pointers.
        //
        //    Input, int INDCRD, the number of input lines for row indices.
        //
        //    Input, int VALCRD, the number of input lines for numerical values.
        //
        //    Input, int RHSCRD, the number of input lines for right hand sides.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NCOL, the number of columns or elements.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Input, char *INDFMT, the 16 character format for reading indices.
        //
        //    Input, char *VALFMT, the 20 character format for reading values.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int NRHSIX, the number of row indices (set to 0
        //    in the case of unassembled matrices.)  Ignored if NRHS = 0.
        //
    {
        output.Add(string.Join("", title).PadRight(72)
                   + string.Join("", key).PadRight(8)
                   + "");

        output.Add(totcrd.ToString().PadRight(14)
                   + string.Join("", ptrcrd).PadRight(14)
                   + string.Join("", indcrd).PadRight(14)
                   + string.Join("", valcrd).PadRight(14)
                   + string.Join("", rhscrd).PadRight(14) + "");

        output.Add(
            string.Join("", mxtype).PadRight(3)
            + "           "
            + string.Join("", nrow).PadRight(14)
            + string.Join("", ncol).PadRight(14)
            + string.Join("", nnzero).PadRight(14)
            + string.Join("", neltvl).PadRight(14) + "");

        output.Add(
            string.Join("", ptrfmt).PadRight(16)
            + string.Join("", indfmt).PadRight(16)
            + string.Join("", valfmt).PadRight(20)
            + string.Join("", rhsfmt).PadRight(20)
            + "");

        switch (rhscrd)
        {
            case > 0:
                output.Add(string.Join("", rhstyp).PadRight(3)
                           + "           "
                           + string.Join("", nrhs).PadRight(14)
                           + string.Join("", nrhsix).PadRight(14) + "");
                break;
        }
    }

    public static double[] hb_matvec_a_mem(int nrow, int ncol, int nnzero, int nrhs,
            int[] colptr, int[] rowind, double[] values, double[] exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_MATVEC_A_MEM multiplies an assembled Harwell Boeing matrix times a vector.
        //
        //  Discussion:
        //
        //    In this "A_MEM" version of the routine, the matrix is assumed to be in
        //    "assembled" form, and all the data is assumed to be small enough
        //    to reside completely in memory; the matrix and multiplicand vectors
        //    are assumed to have been read into memory before this routine is called.
        //
        //    It is assumed that MXTYPE(3:3) = 'A', that is, that the matrix is
        //    stored in the "assembled" format.
        //
        //    Also, the storage used for the vectors X and the products A*X
        //    corresponds to RHSTYP(1:1) = 'F', that is, the "full" storage mode
        //    for vectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NCOL, the number of columns or elements.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int COLPTR[NCOL+1], COLPTR(I) points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    Input, int ROWIND[NNZERO], the row index of each item.
        //
        //    Input, double VALUES[NNZERO], the nonzero values of the matrix.
        //
        //    Input, double EXACT[NCOL*NRHS], contains NRHS dense vectors.
        //
        //    Output, double HB_MATVEC_A_MEM[NROW*NRHS], the product vectors A*X.
        //
    {
        int column;
        int k;
        double[] rhsval;
        int rhs;
        int row;

        rhsval = new double[nrow * nrhs];
        //
        //  Zero out the result vectors.
        //
        for (rhs = 1; rhs <= nrhs; rhs++)
        {
            for (row = 1; row <= nrow; row++)
            {
                rhsval[row - 1 + (rhs - 1) * nrow] = 0.0E+00;
            }
        }

        //
        //  For each column J of the matrix,
        //
        for (column = 1; column <= ncol; column++)
        {
            //
            //  For nonzero entry K
            //
            for (k = colptr[column - 1]; k <= colptr[column] - 1; k++)
            {
                row = rowind[k - 1];
                //
                //  For each right hand side vector:
                //
                //    B(I,1:NRHS) = B(I,1:NRHS) + A(I,J) * X(J,1:NRHS)
                //
                for (rhs = 1; rhs <= nrhs; rhs++)
                {
                    rhsval[row - 1 + (rhs - 1) * nrow] += values[k - 1] * exact[column - 1 + (rhs - 1) * ncol];
                }
            }
        }

        return rhsval;
    }

    public static void hb_rhs_read(string[] input, ref int inputIndex, int nrow, int nnzero, int nrhs, int nrhsix,
            int rhscrd, char[] ptrfmt, char[] indfmt, char[] rhsfmt, char[] mxtype,
            char[] rhstyp, ref double[] rhsval, ref int[] rhsind, ref int[] rhsptr, ref double[] rhsvec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_RHS_READ reads the right hand side information in an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, the unit from which data is read.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int NRHSIX, the number of entries of storage for right
        //    hand side values, in the case where RHSTYP[0] = 'M' and
        //    MXTYPE[2] = 'A'.
        //
        //    Input, int RHSCRD, the number of lines in the file for
        //    right hand sides.
        //
        //    Input, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Input, char *INDFMT, the 16 character format for reading indices.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    If RHSTYP[0] == 'F':
        //
        //      Output, double RHSVAL[NROW*NRHS], contains NRHS dense right hand
        //      side vectors.
        //
        //      Output, int RHSPTR[], is not used.
        //
        //      Output, int RHSIND[], is not used.
        //
        //      Output, int RHSVEC[], is not used.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'A':
        //
        //      Output, double RHSVAL[], is not used.
        //
        //      Output, int RHSPTR[NRHS+1], RHSPTR[I-1] points to the location of
        //      the first entry of right hand side I in the sparse right hand
        //      side vector.
        //
        //      Output, int RHSIND[NRHSIX], indicates, for each entry of
        //      RHSVEC, the corresponding row index.
        //
        //      Output, double RHSVEC[NRHSIX], contains the value of the right hand
        //      side entries.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'E':
        //
        //      Output, double RHSVAL[NNZERO*NRHS], contains NRHS unassembled
        //      finite element vector right hand sides.
        //
        //      Output, int RHSPTR[], is not used.
        //
        //      Output, int RHSIND[], is not used.
        //
        //      Output, double RHSVEC[], is not used.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int khi;
        int klo;
        char[] line = new char[255];
        int line_num;
        int m = 0;
        int r = 0;
        char[] s;
        int w = 0;
        switch (rhscrd)
        {
            //
            //  Read the right hand sides.
            //    case F                             = "full" or "dense";
            //    case not F + matrix storage is "A" = sparse pointer RHS
            //    case not F + matrix storage is "E" = finite element RHS
            //
            case > 0:
                switch (rhstyp[0])
                {
                    //
                    //  Dense right hand sides:
                    //
                    case 'F':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);

                        line_num = 1 + (nrow * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            line = input[inputIndex].ToCharArray();
                            inputIndex++;
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrow * nrhs);

                            khi = 0;
                            for (j = jlo; j <= jhi; j++)
                            {
                                klo = khi + 1;
                                khi = Math.Min(klo + w - 1, line.Length);
                                s = typeMethods.s_substring(line, klo, khi);
                                rhsval[j - 1] = Convert.ToDouble(s);
                            }
                        }

                        break;
                    }
                    //
                    //  Sparse right-hand sides stored like the matrix.
                    //  Read pointer array, indices, and values.
                    //
                    case 'M' when mxtype[2] == 'A':
                    {
                        typeMethods.s_to_format(ptrfmt, ref r, ref code, ref w, ref m);

                        line_num = 1 + (nrhs + 1 - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            line = input[inputIndex].ToCharArray();
                            inputIndex++;
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrhs + 1);

                            khi = 0;
                            for (j = jlo; j <= jhi; j++)
                            {
                                klo = khi + 1;
                                khi = Math.Min(klo + w - 1, line.Length);
                                s = typeMethods.s_substring(line, klo, khi);
                                rhsptr[j - 1] = Convert.ToInt32(s);
                            }
                        }

                        typeMethods.s_to_format(indfmt, ref r, ref code, ref w, ref m);

                        line_num = 1 + (nrhsix - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            line = input[inputIndex].ToCharArray();
                            inputIndex++;
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nnzero);

                            khi = 0;
                            for (j = jlo; j <= jhi; j++)
                            {
                                klo = khi + 1;
                                khi = Math.Min(klo + w - 1, line.Length);
                                s = typeMethods.s_substring(line, klo, khi);
                                rhsind[j - 1] = Convert.ToInt32(s);
                            }
                        }

                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);

                        line_num = 1 + (nrhsix - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            line = input[inputIndex].ToCharArray();
                            inputIndex++;
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrhsix);

                            khi = 0;
                            for (j = jlo; j <= jhi; j++)
                            {
                                klo = khi + 1;
                                khi = Math.Min(klo + w - 1, line.Length);
                                s = typeMethods.s_substring(line, klo, khi);
                                rhsvec[j - 1] = Convert.ToDouble(s);
                            }
                        }

                        break;
                    }
                    //
                    //  Sparse right hand sides in finite element format.
                    //
                    case 'M' when mxtype[2] == 'E':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);

                        line_num = 1 + (nnzero * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            line = input[inputIndex].ToCharArray();
                            inputIndex++;
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nnzero * nrhs);

                            khi = 0;
                            for (j = jlo; j <= jhi; j++)
                            {
                                klo = khi + 1;
                                khi = Math.Min(klo + w - 1, line.Length);
                                s = typeMethods.s_substring(line, klo, khi);
                                rhsval[j - 1] = Convert.ToDouble(s);
                            }
                        }

                        break;
                    }
                    case 'M':
                        Console.WriteLine("");
                        Console.WriteLine("HB_RHS_READ - Fatal error!");
                        Console.WriteLine("  Illegal value of MXTYPE character 3!");
                        break;
                    //
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("HB_RHS_READ - Fatal error!");
                        Console.WriteLine("  Illegal value of RHSTYP character 1!");
                        break;
                }

                break;
        }
    }

    public static void hb_rhs_write(ref List<string> output, int nrow, int nnzero, int nrhs, int nrhsix,
            int rhscrd, char[] ptrfmt, char[] indfmt, char[] rhsfmt, char[] mxtype,
            char[] rhstyp, double[] rhsval, int[] rhsind, int[] rhsptr, double[] rhsvec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_RHS_WRITE writes the right hand side information to an HB file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, List<string> &OUTPUT, the unit to which data is written.
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int NRHSIX, the number of entries of storage for right
        //    hand side values, in the case where RHSTYP[0] = 'M' and
        //    MXTYPE[2] = 'A'.
        //
        //    Input, int RHSCRD, the number of lines in the file for
        //    right hand sides.
        //
        //    Input, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Input, char *INDFMT, the 16 character format for reading indices.
        //
        //    Input, char *RHSFMT, the 20 character format for reading values
        //    of the right hand side.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, char *RHSTYP, the 3 character right hand side type.
        //    First character is F for full storage or M for same as matrix.
        //    Second character is G if starting "guess" vectors are supplied.
        //    Third character is X if exact solution vectors are supplied.
        //    Ignored if NRHS = 0.
        //
        //    If RHSTYP[0] == 'F':
        //
        //      Input, double RHSVAL[NROW*NRHS], contains NRHS dense right hand
        //      side vectors.
        //
        //      Input, int RHSPTR[], is not used.
        //
        //      Input, int RHSIND[], is not used.
        //
        //      Input, double RHSVEC[], is not used.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'A':
        //
        //      Input, double RHSVAL[], is not used.
        //
        //      Input, int RHSPTR[NRHS+1], RHSPTR(I) points to the location of
        //      the first entry of right hand side I in the sparse right hand
        //      side vector.
        //
        //      Input, int RHSIND[NRHSIX], indicates, for each entry of
        //      RHSVEC, the corresponding row index.
        //
        //      Input, double RHSVEC[NRHSIX], contains the value of the right hand
        //      side entries.
        //
        //    If RHSTYP[0] = 'M' and MXTYPE[2] = 'E':
        //
        //      Input, double RHSVAL[NNZERO*NRHS], contains NRHS unassembled
        //      finite element vector right hand sides.
        //
        //      Input, int RHSPTR[], is not used.
        //
        //      Input, int RHSIND[], is not used.
        //
        //      Input, double RHSVEC[], is not used.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int line_num;
        int m = 0;
        int r = 0;
        int w = 0;
        switch (rhscrd)
        {
            //
            //  Read the right hand sides.
            //    case F                             = "full" or "dense";
            //    case not F + matrix storage is "A" = sparse pointer RHS
            //    case not F + matrix storage is "E" = finite element RHS
            //
            case > 0:
                switch (rhstyp[0])
                {
                    //
                    //  Dense right hand sides:
                    //
                    case 'F':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);
                        line_num = 1 + (nrow * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            string cout = "";
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrow * nrhs);

                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += rhsval[j - 1].ToString().PadLeft(w);
                            }

                            output.Add(cout);
                        }

                        break;
                    }
                    //
                    //  Sparse right-hand sides stored like the matrix.
                    //  Read pointer array, indices, and values.
                    //
                    case 'M' when mxtype[2] == 'A':
                    {
                        typeMethods.s_to_format(ptrfmt, ref r, ref code, ref w, ref m);
                        line_num = 1 + (nrhs + 1 - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            string cout = "";
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrhs + 1);

                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += rhsptr[j - 1].ToString().PadLeft(w);
                            }

                            output.Add(cout);
                        }

                        typeMethods.s_to_format(indfmt, ref r, ref code, ref w, ref m);
                        line_num = 1 + (nrhsix - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            string cout = "";
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrhsix);

                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += rhsind[j - 1].ToString().PadLeft(w);
                            }

                            output.Add(cout);
                        }

                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);
                        line_num = 1 + (nrhsix - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            string cout = "";
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nrhsix);

                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += rhsvec[j - 1].ToString().PadLeft(w);
                            }

                            output.Add(cout);
                        }

                        break;
                    }
                    //
                    //  Sparse right hand sides in finite element format.
                    //
                    case 'M' when mxtype[2] == 'E':
                    {
                        typeMethods.s_to_format(rhsfmt, ref r, ref code, ref w, ref m);
                        line_num = 1 + (nnzero * nrhs - 1) / r;

                        jhi = 0;
                        for (i = 1; i <= line_num; i++)
                        {
                            string cout = "";
                            jlo = jhi + 1;
                            jhi = Math.Min(jlo + r - 1, nnzero * nrhs);

                            for (j = jlo; j <= jhi; j++)
                            {
                                cout += rhsval[j - 1].ToString().PadLeft(w);
                            }

                            output.Add(cout);
                        }

                        break;
                    }
                    case 'M':
                        Console.WriteLine("");
                        Console.WriteLine("HB_RHS_WRITE - Fatal error!");
                        Console.WriteLine("  Illegal value of MXTYPE character 3!");
                        break;
                }

                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HB_RHS_WRITE - Fatal error!");
                Console.WriteLine("  Illegal value of RHSTYP character 1!");
                break;
        }
    }

    public static void hb_structure_print(int ncol, char[] mxtype, int nnzero, int neltvl,
            int[] colptr, int[] rowind)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_STRUCTURE_PRINT prints the structure of an HB matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, int NCOL, the number of columns.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, int COLPTR[NCOL+1], COLPTR[I-1] points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Input, int ROWIND[NNZERO], the row index of each item.
        //
        //    If MXTYPE[2] == 'F':
        //
        //      Input, int ROWIND[NELTVL], the row index of each item.
        //
    {
        int j;
        int k;
        int khi;
        int klo;

        switch (mxtype[2])
        {
            case 'A':
            {
                Console.WriteLine("");
                Console.WriteLine("Column Begin   End   ----------------------------------------");
                Console.WriteLine("");
                for (j = 1; j <= Math.Min(ncol, 10); j++)
                {
                    if (colptr[j] - 1 < colptr[j - 1])
                    {
                        Console.WriteLine(j.ToString().PadLeft(6) + "   EMPTY");
                    }
                    else
                    {
                        for (klo = colptr[j - 1]; klo <= colptr[j] - 1; klo += 10)
                        {
                            string cout = "";
                            khi = Math.Min(klo + 9, colptr[j] - 1);
                            if (klo == colptr[j - 1])
                            {
                                cout += j.ToString().PadLeft(6)
                                        + colptr[j - 1].ToString().PadLeft(6)
                                        + (colptr[j] - 1).ToString().PadLeft(6) + "   ";
                            }

                            for (k = klo; k <= khi; k++)
                            {
                                cout += rowind[k - 1].ToString().PadLeft(4);
                            }

                            Console.WriteLine(cout);
                        }
                    }
                }

                Console.WriteLine("                     ----------------------------------------");
                break;
            }
            case 'E':
                Console.WriteLine("");
                Console.WriteLine("Column Begin   End   ----------------------------------------");
                Console.WriteLine("                     ----------------------------------------");

                Console.WriteLine("");
                Console.WriteLine("  I haven't thought about how to print an");
                Console.WriteLine("  unassembled matrix yet!");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HB_STRUCTURE_PRINT - Fatal error!");
                Console.WriteLine("  Illegal value of MXTYPE character #3 = " + mxtype[2] + "");
                break;
        }

    }

    public static void hb_structure_read(string[] input, ref int inputIndex, int ncol, char[] mxtype, int nnzero,
            int neltvl, int ptrcrd, char[] ptrfmt, int indcrd, char[] indfmt,
            ref int[] colptr, ref int[] rowind)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_STRUCTURE_READ reads the structure of an HB matrix.
        //
        //  Discussion:
        //
        //    The user should already have opened the file, and positioned it
        //    to just after the header records.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, the unit from which data is read.
        //
        //    Input, int NCOL, the number of columns.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, int PTRCRD, the number of input lines for pointers.
        //
        //    Input, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Input, int INDCRD, the number of input lines for indices.
        //
        //    Input, char *INDFMT, the 16 character format for reading indices.
        //
        //    Output, int COLPTR[NCOL+1], COLPTR[I-1] points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Output, int ROWIND[NNZERO], the row index of each item.
        //
        //    If MXTYPE[2] == 'F':
        //
        //      Output, int ROWIND[NELTVL], the row index of each item.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int khi;
        int klo;
        char[] line = new char[255];
        int line_num;
        int m = 0;
        int number;
        int r = 0;
        char[] s;
        int w = 0;

        typeMethods.s_to_format(ptrfmt, ref r, ref code, ref w, ref m);

        line_num = mxtype[2] switch
        {
            'A' => 1 + (ncol + 1 - 1) / r,
            _ => 1 + (ncol - 1) / r
        };

        jhi = 0;
        for (i = 1; i <= line_num; i++)
        {
            line = input[inputIndex].ToCharArray();
            inputIndex++;
            jlo = jhi + 1;
            jhi = mxtype[2] switch
            {
                'A' => Math.Min(jlo + r - 1, ncol + 1),
                _ => Math.Min(jlo + r - 1, ncol)
            };

            khi = 0;
            for (j = jlo; j <= jhi; j++)
            {
                klo = khi + 1;
                khi = Math.Min(klo + w - 1, line.Length);
                s = typeMethods.s_substring(line, klo, khi);
                colptr[j - 1] = Convert.ToInt32(s);
            }
        }

        switch (mxtype[2])
        {
            case 'A':
            {
                typeMethods.s_to_format(indfmt, ref r, ref code, ref w, ref m);

                line_num = 1 + (nnzero - 1) / r;

                jhi = 0;
                for (i = 1; i <= line_num; i++)
                {
                    line = input[inputIndex].ToCharArray();
                    inputIndex++;
                    jlo = jhi + 1;
                    jhi = Math.Min(jlo + r - 1, nnzero);

                    khi = 0;
                    for (j = jlo; j <= jhi; j++)
                    {
                        klo = khi + 1;
                        khi = Math.Min(klo + w - 1, line.Length);
                        s = typeMethods.s_substring(line, klo, khi);
                        rowind[j - 1] = Convert.ToInt32(s);
                    }
                }

                break;
            }
            case 'E':
            {
                typeMethods.s_to_format(indfmt, ref r, ref code, ref w, ref m);

                number = colptr[ncol - 1] - colptr[0];
                line_num = 1 + (number - 1) / r;

                jhi = 0;
                for (i = 1; i <= line_num; i++)
                {
                    line = input[inputIndex].ToCharArray();
                    inputIndex++;
                    jlo = jhi + 1;
                    jhi = Math.Min(jlo + r - 1, number);

                    khi = 0;
                    for (j = jlo; j <= jhi; j++)
                    {
                        klo = khi + 1;
                        khi = Math.Min(klo + w - 1, line.Length);
                        s = typeMethods.s_substring(line, klo, khi);
                        rowind[j - 1] = Convert.ToInt32(s);
                    }
                }

                break;
            }
            default:
                Console.WriteLine("");
                Console.WriteLine("HB_STRUCTURE_READ - Fatal error!");
                Console.WriteLine("  Illegal value of MXTYPE character 3.");
                break;
        }

    }

    public static void hb_structure_write(ref List<string> output, int ncol, char[] mxtype,
            int nnzero, int neltvl, char[] ptrfmt, char[] indfmt, int[] colptr,
            int[] rowind)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_STRUCTURE_WRITE writes the structure of an HB matrix.
        //
        //  Discussion:
        //
        //    If the user is creating an HB file, then the user should
        //    already have opened the file, and written the header records.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, List<string> &OUTPUT, the unit to which data is written.
        //
        //    Input, int NCOL, the number of columns.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, char *PTRFMT, the 16 character format for reading pointers.
        //
        //    Input, char *INDFMT, the 16 character format for reading indices.
        //
        //    Input, int COLPTR[NCOL+1], COLPTR[I-1] points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Input, int ROWIND[NNZERO], the row index of each item.
        //
        //    If MXTYPE[2] == 'F':
        //
        //      Input, int ROWIND[NELTVL], the row index of each item.
        //
    {
        char code = ' ';
        int j;
        int m = 0;
        int number;
        int r = 0;
        int w = 0;

        typeMethods.s_to_format(ptrfmt, ref r, ref code, ref w, ref m);

        string cout = "";
        for (j = 1; j <= ncol + 1; j++)
        {
            cout += colptr[j - 1].ToString().PadLeft(w);
            switch (j % r)
            {
                case 0:
                    output.Add(cout);
                    cout = "";
                    break;
            }
        }

        if ((ncol + 1) % r != 0)
        {
            output.Add(cout);
        }

        typeMethods.s_to_format(indfmt, ref r, ref code, ref w, ref m);

        cout = "";
        switch (mxtype[2])
        {
            case 'A':
            {
                for (j = 1; j <= nnzero; j++)
                {
                    cout += rowind[j - 1].ToString().PadLeft(w);
                    switch (j % r)
                    {
                        case 0:
                            output.Add(cout);
                            cout = "";
                            break;
                    }
                }

                if (nnzero % r != 0)
                {
                    output.Add(cout);
                }

                break;
            }
            case 'E':
            {
                number = colptr[ncol - 1] - colptr[0];

                cout = "";
                for (j = 1; j <= number; j++)
                {
                    cout += rowind[j - 1].ToString().PadLeft(w);
                    switch (j % r)
                    {
                        case 0:
                            output.Add(cout);
                            cout = "";
                            break;
                    }
                }

                if (number % r != 0)
                {
                    output.Add(cout);
                }

                break;
            }
            default:
                Console.WriteLine("");
                Console.WriteLine("HB_STRUCTURE_WRITE - Fatal error!");
                Console.WriteLine("  Illegal value of MXTYPE character 3.");
                break;
        }

    }

    public static int[] hb_ua_colind(int ncol, int[] colptr, int nnzero)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_UA_COLUMN_INDEX creates a column index for an unsymmetric assembled matrix.
        //
        //  Discussion:
        //
        //    It is assumed that the input data corresponds to a Harwell-Boeing
        //    matrix which is unsymmetric, and assembled.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, integer NCOL, the number of columns.
        //
        //    Input, integer COLPTR[NCOL+1], COLPTR(I) points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    Input, int NNZERO, the number of nonzeros.
        //
        //    Output, int HB_UA_COLIND[NNZERO], the column index of each matrix value.
        //
    {
        int[] colind;
        int i;
        int j;

        colind = new int[nnzero];

        for (i = 1; i <= ncol; i++)
        {
            for (j = colptr[i - 1]; j <= colptr[i] - 1; j++)
            {
                colind[j - 1] = i;
            }
        }

        return colind;
    }

    public static void hb_values_print(int ncol, int[] colptr, char[] mxtype, int nnzero,
            int neltvl, double[] values)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_VALUES_PRINT prints the values of an HB matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, int NCOL, the number of columns.
        //
        //    Input, int COLPTR[NCOL+1], COLPTR[I-1] points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Input, double VALUES[NNZERO], the nonzero values of the matrix.
        //
        //    If MXTYPE[2] == 'E':
        //
        //      Input, double VALUES[NELTVL], the nonzero values of the matrix.
        //
    {
        int j;
        int k;
        int khi;
        int klo;

        switch (mxtype[2])
        {
            case 'A':
            {
                Console.WriteLine("");
                Console.WriteLine("Column Begin   End   ----------------------------------------");
                Console.WriteLine("");

                for (j = 1; j <= ncol; j++)
                {
                    switch (j)
                    {
                        case > 5 when j < ncol:
                            continue;
                    }

                    if (j == ncol && 6 < ncol)
                    {
                        Console.WriteLine("Skipping intermediate columns...)");
                    }

                    if (colptr[j] - 1 < colptr[j - 1])
                    {
                        Console.WriteLine(j.ToString().PadLeft(6) + "   EMPTY");
                    }
                    else
                    {
                        for (klo = colptr[j - 1]; klo <= colptr[j] - 1; klo += 5)
                        {
                            string cout = "";
                            khi = Math.Min(klo + 4, colptr[j] - 1);
                            if (klo == colptr[j - 1])
                            {
                                cout += j.ToString().PadLeft(5)
                                        + colptr[j - 1].ToString().PadLeft(5)
                                        + (colptr[j] - 1).ToString().PadLeft(5) + "   ";
                            }

                            for (k = klo; k <= khi; k++)
                            {
                                cout += values[k - 1].ToString().PadLeft(12);
                            }

                            Console.WriteLine(cout);
                        }
                    }
                }

                Console.WriteLine("                     ----------------------------------------");
                break;
            }
            case 'E':
                Console.WriteLine("");
                Console.WriteLine("Column Begin   End   ----------------------------------------");
                Console.WriteLine("                     ----------------------------------------");

                Console.WriteLine("");
                Console.WriteLine("I haven't thought about how to print an");
                Console.WriteLine("unassembled matrix yet!");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HB_VALUES_PRINT - Fatal error!");
                Console.WriteLine("  Illegal value of MXTYPE character 3 = " + mxtype[2] + "");
                break;
        }
    }

    public static void hb_values_read(string[] input, ref int inputIndex, int valcrd, char[] mxtype, int nnzero,
            int neltvl, char[] valfmt, ref double[] values)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_VALUES_READ reads the values of an HB matrix.
        //
        //  Discussion:
        //
        //    The user should already have opened the file, and positioned it
        //    to just after the header and structure records.
        //
        //    Values are contained in an HB file if the VALCRD parameter
        //    is nonzero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, ifstream &INPUT, the unit from which data is read.
        //
        //    Input, int VALCRD, the number of input lines for numerical values.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, char *VALFMT, the 20 character format for reading values.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Output, double VALUES[NNZERO], the nonzero values of the matrix.
        //
        //    If MXTYPE[2] == 'E':
        //
        //      Output, double VALUES[NELTVL], the nonzero values of the matrix.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int khi;
        int klo;
        char[] line = new char[255];
        int line_num;
        int m = 0;
        int r = 0;
        char[] s;
        int w = 0;

        typeMethods.s_to_format(valfmt, ref r, ref code, ref w, ref m);

        switch (valcrd)
        {
            //
            //  Read the matrix values.
            //    case "A" = assembled;
            //    case "E" = unassembled finite element matrices.
            //
            case > 0:
            {
                switch (mxtype[2])
                {
                    case 'A':
                        line_num = 1 + (nnzero - 1) / r;
                        break;
                    case 'E':
                        line_num = 1 + (neltvl - 1) / r;
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("HB_VALUES_READ - Fatal error!");
                        Console.WriteLine("  Illegal value of MXTYPE character 3.");
                        return;
                }

                jhi = 0;
                for (i = 1; i <= line_num; i++)
                {
                    line = input[inputIndex].ToCharArray();
                    inputIndex++;
                    jlo = jhi + 1;
                    jhi = mxtype[2] switch
                    {
                        'A' => Math.Min(jlo + r - 1, nnzero),
                        _ => Math.Min(jlo + r - 1, neltvl)
                    };

                    khi = 0;
                    for (j = jlo; j <= jhi; j++)
                    {
                        klo = khi + 1;
                        khi = Math.Min(klo + w - 1, line.Length);
                        s = typeMethods.s_substring(line, klo, khi);
                        values[j - 1] = Convert.ToDouble(s);
                    }
                }

                break;
            }
        }
    }

    public static void hb_values_write(ref List<string> output, int valcrd, char[] mxtype,
            int nnzero, int neltvl, char[] valfmt, double[] values)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_VALUES_WRITE writes the values of an HB matrix.
        //
        //  Discussion:
        //
        //    If the user is creating an HB file, then the user should already
        //    have opened the file, and written the header and structure records.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, List<string> &OUTPUT, the unit to which data is written.
        //
        //    Input, int VALCRD, the number of input lines for numerical values.
        //
        //    Input, char *MXTYPE, the 3 character matrix type.
        //    First character is R for Real, C for complex, P for pattern only.
        //    Second character is S for symmetric, U for unsymmetric, H for
        //      Hermitian, Z for skew symmetric, R for rectangular.
        //    Third character is A for assembled and E for unassembled
        //      finite element matrices.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  In the case of unassembled finite
        //    element matrices, in which the right hand side vectors are also
        //    stored as unassembled finite element vectors, this is the total
        //    number of entries in a single unassembled right hand side vector.
        //
        //    Input, int NELTVL, the number of finite element matrix entries,
        //    set to 0 in the case of assembled matrices.
        //
        //    Input, char *VALFMT, the 20 character format for reading values.
        //
        //    If MXTYPE[2] == 'A':
        //
        //      Input, double VALUES[NNZERO], the nonzero values of the matrix.
        //
        //    If MXTYPE[2] == 'E':
        //
        //      Input, double VALUES[NELTVL], the nonzero values of the matrix.
        //
    {
        char code = ' ';
        int i;
        int j;
        int jhi;
        int jlo;
        int line_num;
        int m = 0;
        int r = 0;
        int w = 0;

        switch (valcrd)
        {
            case > 0:
            {
                typeMethods.s_to_format(valfmt, ref r, ref code, ref w, ref m);

                switch (mxtype[2])
                {
                    case 'A':
                        line_num = 1 + (nnzero - 1) / r;
                        break;
                    case 'E':
                        line_num = 1 + (neltvl - 1) / r;
                        break;
                    default:
                        Console.WriteLine("");
                        Console.WriteLine("HB_VALUES_WRITE - Fatal error!");
                        Console.WriteLine("  Illegal value of MXTYPE character 3.");
                        return;
                }

                jhi = 0;
                for (i = 1; i <= line_num; i++)
                {
                    string cout = "";
                    jlo = jhi + 1;
                    jhi = mxtype[2] switch
                    {
                        'A' => Math.Min(jlo + r - 1, nnzero),
                        _ => Math.Min(jlo + r - 1, neltvl)
                    };

                    for (j = jlo; j <= jhi; j++)
                    {
                        cout += values[j - 1].ToString().PadLeft(w);
                    }

                    output.Add(cout);
                }

                break;
            }
        }
    }

    public static double[] hb_vecmat_a_mem(int nrow, int ncol, int nnzero, int nrhs,
            int[] colptr, int[] rowind, double[] values, double[] exact)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HB_VECMAT_A_MEM multiplies a vector times an assembled Harwell Boeing matrix.
        //
        //  Discussion:
        //
        //    In this "A_MEM" version of the routine, the matrix is assumed to be in
        //    "assembled" form, and all the data is assumed to be small enough
        //    to reside completely in memory; the matrix and multiplicand vectors
        //    are assumed to have been read into memory before this routine is called.
        //
        //    It is assumed that MXTYPE(3:3) = 'A', that is, that the matrix is
        //    stored in the "assembled" format.
        //
        //    Also, the storage used for the vectors X and the products A*X
        //    corresponds to RHSTYP(1:1) = 'F', that is, the "full" storage mode
        //    for vectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Iain Duff, Roger Grimes, John Lewis,
        //    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
        //    October 1992.
        //
        //  Parameters:
        //
        //    Input, int NROW, the number of rows or variables.
        //
        //    Input, int NCOL, the number of columns or elements.
        //
        //    Input, int NNZERO.  In the case of assembled sparse matrices,
        //    this is the number of nonzeroes.  
        //
        //    Input, int NRHS, the number of right hand sides.
        //
        //    Input, int COLPTR[NCOL+1], COLPTR(I) points to the location of
        //    the first entry of column I in the sparse matrix structure.
        //
        //    Input, int ROWIND[NNZERO], the row index of each item.
        //
        //    Input, double VALUES[NNZERO], the nonzero values of the matrix.
        //
        //    Input, double EXACT[NROW*NRHS], contains NRHS dense vectors.
        //
        //    Output, double HB_VECMAT_A_MEM[NCOL*NRHS], the product vectors A'*X.
        //
    {
        int column;
        int k;
        double[] rhsval;
        int rhs;
        int row;

        rhsval = new double[ncol * nrhs];
        //
        //  Zero out the result vectors.
        //
        for (rhs = 1; rhs <= nrhs; rhs++)
        {
            for (column = 1; column <= ncol; column++)
            {
                rhsval[column - 1 + (rhs - 1) * ncol] = 0.0E+00;
            }
        }

        //
        //  For each column J of the matrix,
        //
        for (column = 1; column <= ncol; column++)
        {
            //
            //  For nonzero entry K
            //
            for (k = colptr[column - 1]; k <= colptr[column] - 1; k++)
            {
                row = rowind[k - 1];
                //
                //  For each right hand side vector:
                //
                //    B(J,1:NRHS) = B(J,1:NRHS) + X(I,1:NRHS) * A(I,J)
                //
                for (rhs = 1; rhs <= nrhs; rhs++)
                {
                    rhsval[column - 1 + (rhs - 1) * ncol] += values[k - 1] * exact[row - 1 + (rhs - 1) * nrow];
                }
            }
        }

        return rhsval;
    }
}