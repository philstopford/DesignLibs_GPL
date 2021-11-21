namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void invmod(ref int[] mat, ref int[] imat, int[] rmod, int[] cmod, int nrow,
            ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVMOD inverts a matrix using modulo arithmetic.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 December 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Roger Payne.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Payne,
        //    Inversion of matrices with contents subject to modulo arithmetic,
        //    Applied Statistics,
        //    Volume 46, Number 2, 1997, pages 295-298.
        //
        //  Parameters:
        //
        //    Input/output, int MAT[NROW*NROW].
        //    On input, the matrix to be inverted.
        //    On output, the product of the input matrix and IMAT.
        //
        //    Output, int IMAT[NROW*NROW], the inverse matrix.  
        //    If IFAULT = -1 on output, then IMAT is only a left inverse.
        //
        //    Input, int RMOD[NROW], the modulus for values in each row.
        //
        //    Input, int CMOD[NROW], the modulus for values 
        //    in each column.
        //
        //    Input, int NROW, the order of the matrix.
        //
        //    Output, int &IFAULT, an error flag.
        //    0, no error was detected.
        //    -1, only a left inverse could be formed.
        //    1, the matrix contains elements that are negative, or too large.
        //    2, the matrix contains nonzero elements in mixed modulus positions.
        //    3, the matrix cannot be inverted.
        //
    {
        int ir;
        //
        //  Check that elements in 'mixed-moduli' positions are all zero.
        //
        int n = 0;
        for (int i = 1; i <= nrow; i++)
        {
            for (int j = 1; j <= nrow; j++)
            {
                n += 1;

                if (rmod[i - 1] != cmod[j - 1] && 0 < mat[n - 1])
                {
                    ifault = 2;
                    return;
                }

                if (rmod[i - 1] >= mat[n - 1] && mat[n - 1] >= 0)
                {
                    continue;
                }

                ifault = 1;
                return;
            }
        }

        n = 0;
        for (int i = 1; i <= nrow; i++)
        {
            for (int j = 1; j <= nrow; j++)
            {
                n += 1;
                imat[n - 1] = 0;
            }
        }

//
//  Sort rows and columns into ascending order of moduli
//
        int[] rsort = new int[nrow];
        int[] csort = new int[nrow];

        msort(ref mat, ref imat, ref rmod, ref cmod, ref rsort, ref csort, nrow);
//
//  Complete initialization of inverse matrix 
//
        for (n = 1; n <= nrow * nrow; n = n + nrow + 1)
        {
            imat[n - 1] = 1;
        }

//
//  Invert the matrix.
//
        for (ir = 1; ir <= nrow; ir++)
        {
            int kir = (ir - 1) * nrow;

            int k;
            int kjr;
            switch (mat[kir + ir - 1])
            {
                case 0:
                {
                    //
//  Find a row JR below IR such that K(JR,IR)>0
//
                    bool all_zero = true;

                    for (kjr = kir + nrow + ir; kjr <= nrow * nrow; kjr += nrow)
                    {
                        if (0 >= mat[kjr - 1])
                        {
                            continue;
                        }

                        all_zero = false;
                        break;
                    }

                    switch (all_zero)
                    {
    
//
                        //  Column IR contains all zeros in rows IR or below:
                        //  look for a row above with zeros to left of column IR 
                        //  and K(JR,IR)>0
                        //
                        case true:
                        {
                            for (kjr = ir; kjr <= kir; kjr += nrow)
                            {
                                if (0 >= mat[kjr - 1])
                                {
                                    continue;
                                }

                                for (int i = kjr - ir + 1; i < kjr; i++)
                                {
                                    switch (mat[i - 1])
                                    {
                                        case > 0:
                                            ifault = 3;
                                            return;
                                    }
                                }

                                all_zero = false;
                                break;
                            }

                            break;
                        }
                    }

                    switch (all_zero)
                    {
    
//
                        //  Column IR contains all zeros
                        //
                        case true:
                            continue;
                    }

                    //
                    //  Switch row JR with row IR
                    //
                    kjr -= ir;

                    for (int i = 1; i <= nrow; i++)
                    {
                        k = mat[kir + i - 1];
                        mat[kir + i - 1] = mat[kjr + i - 1];
                        mat[kjr + i - 1] = k;

                        k = imat[kir + i - 1];
                        imat[kir + i - 1] = imat[kjr + i - 1];
                        imat[kjr + i - 1] = k;
                    }

                    break;
                }
            }

            //
            //  Find a multiplier N such that N*MAT(IR,IR)=1 mod(P{IR})
            //
            k = mat[kir + ir - 1];
            for (n = 1; n < rmod[ir - 1]; n++)
            {
                if (n * k % rmod[ir - 1] == 1)
                {
                    break;
                }
            }

            switch (n)
            {
                //
                //  Multiply row IR by N.
                //
                case > 1:
                {
                    for (int i = kir + 1; i <= ir * nrow; i++)
                    {
                        mat[i - 1] *= n;
                        imat[i - 1] *= n;
                    }

                    break;
                }
            }

            //
            //  Subtract MAT(JR,IR) * row IR from each row JR
            //
            for (kjr = 0; kjr < nrow * nrow; kjr += nrow)
            {
                n = rmod[ir - 1] - mat[kjr + ir - 1];
                if (kjr == kir || n == 0)
                {
                    continue;
                }

                for (int i = 1; i <= nrow; i++)
                {
                    mat[kjr + i - 1] = (mat[kjr + i - 1] + n * mat[kir + i - 1]) % cmod[i - 1];
                    imat[kjr + i - 1] = (imat[kjr + i - 1] + n * imat[kir + i - 1]) % cmod[i - 1];
                }
            }

        }

        //
        //  Check inversion was possible - that result has
        //  non-zero elements only on diagonal.
        //
        ifault = 0;
        //
        //  If we encounter a zero diagonal element, then only a left inverse
        //  will be formed.
        //
        for (n = 1; n <= nrow * nrow; n = n + nrow + 1)
        {
            ifault = mat[n - 1] switch
            {
                0 => -1,
                _ => ifault
            };

            mat[n - 1] = -mat[n - 1];
        }

        for (n = 1; n <= nrow * nrow; n++)
        {
            switch (mat[n - 1])
            {
                case > 0:
                    ifault = 3;
                    return;
            }
        }

        for (n = 1; n <= nrow * nrow; n = n + nrow + 1)
        {
            mat[n - 1] = -mat[n - 1];
        }

        //
        //  Unsort the rows and columns back into their original order.
        //
        musort(ref mat, ref imat, ref rmod, ref cmod, ref rsort, ref csort, nrow);

    }

    public static void msort(ref int[] mat, ref int[] imat, ref int[] rmod, ref int[] cmod, ref int[] rsort,
        ref int[] csort, int nrow)

//****************************************************************************80
//
//  Purpose:
//
//    MSORT sorts matrix rows and columns in ascending order of moduli.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    Original FORTRAN77 version by Roger Payne.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Payne,
//    Inversion of matrices with contents subject to modulo arithmetic,
//    Applied Statistics,
//    Volume 46, Number 2, 1997, pages 295-298.
//
//  Parameters:
//
//    Input/output, int MAT[NROW*NROW].
//    On output, the matrix has been sorted.
//
//    Ignoreput, int IMAT[NROW*NROW].  
//    This quantity is ignored.
//
//    Input/output, int RMOD[NROW], the modulus for values in 
//    each row.  On output, these have been rearranged according to the sorting.
//
//    Input/output, int CMOD[NROW], the modulus for values in 
//    each column.  On output, these have been rearranged according to the 
//    sorting.
//
//    Output, int RSORT[NROW], the sorted row indices.
//
//    Output, int CSORT[NROW], the sorted column indices.
//
//    Input, int NROW, the order of the matrix.
//
    {
        int i;
        int irc;
        int j;
        int jrc;
        int p;
//
//  Initialize row and column addresses.
//
        for (i = 1; i <= nrow; i++)
        {
            rsort[i - 1] = i;
            csort[i - 1] = i;
        }

//
//  Sort the rows.
//
        for (irc = 1; irc <= nrow; irc++)
        {
//
//  Find the next row.
//
            jrc = irc;
            p = rmod[irc - 1];

            for (i = irc + 1; i <= nrow; i++)
            {
                if (rmod[i - 1] >= p)
                {
                    continue;
                }

                p = rmod[i - 1];
                jrc = i;
            }

            if (irc == jrc)
            {
                continue;
            }

            i = rmod[irc - 1];
            rmod[irc - 1] = rmod[jrc - 1];
            rmod[jrc - 1] = i;

            i = rsort[irc - 1];
            rsort[irc - 1] = rsort[jrc - 1];
            rsort[jrc - 1] = i;
//
//  Switch the rows.
//
            int kirc = (irc - 1) * nrow;
            int kjrc = (jrc - 1) * nrow;

            for (j = 1; j <= nrow; j++)
            {
                i = mat[kirc + j - 1];
                mat[kirc + j - 1] = mat[kjrc + j - 1];
                mat[kjrc + j - 1] = i;
            }
        }

//
//  Sort the columns.
//
        for (irc = 1; irc <= nrow; irc++)
        {
//
//  Find the next column.
//
            jrc = irc;
            p = cmod[irc - 1];

            for (i = irc + 1; i <= nrow; i++)
            {
                if (cmod[i - 1] >= p)
                {
                    continue;
                }

                p = cmod[i - 1];
                jrc = i;
            }

            if (irc == jrc)
            {
                continue;
            }

            i = cmod[irc - 1];
            cmod[irc - 1] = cmod[jrc - 1];
            cmod[jrc - 1] = i;

            i = csort[irc - 1];
            csort[irc - 1] = csort[jrc - 1];
            csort[jrc - 1] = i;
//
//  Switch the columns.
//
            for (j = 0; j < nrow * nrow; j += nrow)
            {
                i = mat[irc + j - 1];
                mat[irc + j - 1] = mat[jrc + j - 1];
                mat[jrc + j - 1] = i;
            }
        }
    }

    public static void musort(ref int[] mat, ref int[] imat, ref int[] rmod, ref int[] cmod, ref int[] rsort,
        ref int[] csort, int nrow)

//****************************************************************************80
//
//  Purpose:
//
//    MUSORT unsorts the inverse matrix rows and columns into the original order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    Original FORTRAN77 version by Roger Payne.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Payne,
//    Inversion of matrices with contents subject to modulo arithmetic,
//    Applied Statistics,
//    Volume 46, Number 2, 1997, pages 295-298.
//
//  Parameters:
//
//    Input/output, int MAT[NROW*NROW].
//    On output, the matrix has been "unsorted".
//
//    Input/output, int IMAT[NROW*NROW].
//    On output, the matrix has been "unsorted".
//
//    Input/output, int RMOD[NROW], the modulus for values in 
//    each row.  On output, these have been restored to their original ordering.
//
//    Input/output, int CMOD[NROW], the modulus for values in
//    each column.  On output, these have been restored to their original 
//    ordering.
//
//    Input/output, int RSORT[NROW], the sorted row indices.
//
//    Input/output, int CSORT[NROW], the sorted column indices.
//
//    Input, int NROW, the order of the matrix.
//
    {
        int i;
        int irc;
        int j;
        int jrc;
        int kirc;
        int kjrc;
//
//  Sort rows of inverse (= columns of original).
//
        for (irc = 1; irc <= nrow; irc++)
        {
//
//  Find next row.
//
            if (csort[irc - 1] == irc)
            {
                continue;
            }

            for (jrc = irc + 1; jrc <= nrow; jrc++)
            {
                if (csort[jrc - 1] == irc)
                {
                    break;
                }
            }

            i = cmod[irc - 1];
            cmod[irc - 1] = cmod[jrc - 1];
            cmod[jrc - 1] = i;

            i = csort[irc - 1];
            csort[irc - 1] = csort[jrc - 1];
            csort[jrc - 1] = i;
//
//  Switch rows.
//
            kirc = (irc - 1) * nrow;
            kjrc = (jrc - 1) * nrow;

            for (j = 1; j <= nrow; j++)
            {
                i = imat[kirc + j - 1];
                imat[kirc + j - 1] = imat[kjrc + j - 1];
                imat[kjrc + j - 1] = i;
            }
        }

//
//  Sort the columns of the inverse (= rows of original).
//
        for (irc = 1; irc <= nrow; irc++)
        {
//
//  Find the next column.
//
            if (rsort[irc - 1] != irc)
            {
                for (jrc = irc + 1; jrc <= nrow; jrc++)
                {
                    if (rsort[jrc - 1] == irc)
                    {
                        break;
                    }
                }

                i = rmod[irc - 1];
                rmod[irc - 1] = rmod[jrc - 1];
                rmod[jrc - 1] = i;

                i = rsort[irc - 1];
                rsort[irc - 1] = rsort[jrc - 1];
                rsort[jrc - 1] = i;
//
//  Switch the columns of IMAT.
//
                for (j = 0; j < nrow * nrow; j += nrow)
                {
                    i = imat[irc + j - 1];
                    imat[irc + j - 1] = imat[jrc + j - 1];
                    imat[jrc + j - 1] = i;
                }

//
//  Switch the diagonal elements of MAT (others are zero).
//
                kirc = (irc - 1) * nrow + irc;
                kjrc = (jrc - 1) * nrow + jrc;

                i = mat[kirc - 1];
                mat[kirc - 1] = mat[kjrc - 1];
                mat[kjrc - 1] = i;
            }
        }
    }
}