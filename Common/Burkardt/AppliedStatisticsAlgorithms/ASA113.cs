using System;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static void swap(double[] varval, ref int[] klass, ref int[] clsize, int in_, int ik, int iv,
                                ref double critvl, ref int ntrans, ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SWAP interchanges objects between different classes to improve a criterion.
        //
        //  Discussion:
        //
        //    This routine is given a classification of objects, including the
        //    number of objects in each class, and the current value of some criterion
        //    which is desired to be minimized.
        //
        //    The routine calculates the change in criterion for all possible swaps,
        //    that is, operations in which two objects in different classes exchange 
        //    places. Each swap that would result in a lowering of the criterion is 
        //    executed, and the related quantities are updated.
        //
        //    When no more advantageous swaps can be found, the routine returns.
        //
        //    The routine relies on a user-supplied routine, CRSWAP, to report the
        //    expected change in the criterion for a given swap, and to carry
        //    out that transfer if requested.
        //
        //    The variables CLASS and CRITVL have been added to the argument list
        //    of CRSWAP.
        //
        //    Also, the order of the two classes "L" and "M" was interchanged in
        //    the call to CRSWAP.  The original order was counterintuitive.
        //
        //    Sinced CLASS is a reserved keyword in C++, the variable originally
        //    named "CLASS" has been unoriginally renamed "KLASS".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Colin Banfield, LC Bassill.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Colin Banfield, LC Bassill,
        //    Algorithm AS 113:
        //    A transfer for non-hierarchichal classification,
        //    Applied Statistics,
        //    Volume 26, Number 2, 1977, pages 206-210.
        //
        //  Parameters:
        //
        //    Input, double VARVAL(IN,IV), the data values.  There are 
        //    IN objects, each having spatial dimension IV.
        //
        //    Input/output, int KLASS(IN), the classification of 
        //    each object.
        //
        //    Input/output, int CLSIZE(IK), the number of objects 
        //    in each class.
        //
        //    Input, int IN, the number of objects.
        //
        //    Input, int IK, the number of classes.
        //
        //    Input, int IV, the number of spatial dimensions, 
        //    or variates, of the objects.
        //
        //    Input/output, double *CRITVL, the current value of the criterion.
        //
        //    Output, int *NTRANS, the number of transfers executed.
        //
        //    Output, int *IFAULT, error indicator.
        //    0, no error detected.
        //    1, the number of classes was less than 2.
        //    2, the number of objects was less than the number of classes.
        //
        {
            double eps = 1.0E-38;

            if (ik <= 1)
            {
                ifault = 1;
                return;
            }

            if ( in_ <= ik )
            {
                ifault = 2;
                return;
            }

            ifault = 0;
            int icount = 0;
            ntrans = 0;
            int itop = ( in_ *( in_ -1 ) ) / 2;

            int i = 1;

            for (;;)
            {
                i = i + 1;

                if (itop <= icount)
                {
                    break;
                }

                if ( in_ < i )
                {
                    i = 1;
                    continue;
                }

                int l = klass[i - 1];
                int it = i - 1;
                //
                //  Test the swap of object I from class M to L, 
                //  and object J from class L to M.
                //
                for (int j = 1; j <= it; j++)
                {
                    icount = icount + 1;
                    int m = klass[j - 1];

                    if (l != j)
                    {
                        if (clsize[l - 1] != 1 || clsize[m - 1] != 1)
                        {
                            int iswitch = 1;
                            double inc = crswap(varval, klass, clsize, in_, ik, iv, critvl,
                                i, j, l, m, iswitch);

                            if (inc < -eps)
                            {
                                critvl = critvl + inc;
                                icount = 0;

                                iswitch = 2;
                                crswap(varval, klass, clsize, in_, ik, iv, critvl,
                                    i, j, l, m, iswitch);

                                ntrans = ntrans + 1;
                                klass[i - 1] = m;
                                klass[j - 1] = l;
                                l = m;
                            }
                        }
                    }
                }
            }
        }
        
        public static void trnsfr ( double[] varval, ref int[] klass, ref int[] clsize, int in_, int ik, 
                                    int iv, ref double critvl, ref int ntrans, ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRNSFR transfers objects between classes to improve a criterion.
        //
        //  Discussion:
        //
        //    This routine is given a classification of objects, including the
        //    number of objects in each class, and the current value of some criterion
        //    which is desired to be minimized.
        //
        //    The routine calculates the change in criterion for all possible transfers
        //    of any object from its current class to a different class.  Each transfer
        //    that would result in a lowering of the criterion is executed, and the
        //    related quantities are updated.
        //
        //    When no more advantageous transfers can be found, the routine returns.
        //
        //    The routine relies on a user-supplied routine, CRTRAN, to report the
        //    expected change in the criterion for a given transfer, and to carry
        //    out that transfer if requested.
        //
        //    The variables CLASS and CRITVL have been added to the argument list
        //    of CRTRAN.
        //
        //    Also, the order of the two classes "L" and "M" was interchanged in
        //    the call to CRTRAN.  The original order was counterintuitive.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Colin Banfield, LC Bassill.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Colin Banfield, LC Bassill,
        //    Algorithm AS 113:
        //    A transfer for non-hierarchichal classification,
        //    Applied Statistics,
        //    Volume 26, Number 2, 1977, pages 206-210.
        //
        //  Parameters:
        //
        //    Input, double VARVAL(IN,IV), the data values.  There are IN 
        //    objects, each having spatial dimension IV.
        //
        //    Input/output, int KLASS(IN), the classification of 
        //    each object.
        //
        //    Input/output, int CLSIZE(IK), the number of objects in 
        //    each class.
        //
        //    Input, int IN, the number of objects.
        //
        //    Input, int IK, the number of classes.
        //
        //    Input, int IV, the number of spatial dimensions, or 
        //    variates, of the objects.
        //
        //    Input/output, double *CRITVL, the current value of the criterion.
        //
        //    Output, int *NTRANS, the number of transfers executed.
        //
        //    Output, int *IFAULT, error indicator.
        //    0, no error detected.
        //    1, the number of classes was less than 2.
        //    2, the number of objects was less than the number of classes.
        //
        {
            double eps = 1.0E-38;

            if (ik <= 1)
            {
                ifault = 1;
                return;
            }

            if ( in_ <= ik )
            {
                ifault = 2;
                return;
            }

            ifault = 0;
            ntrans = 0;
            int i = 0;
            int icount = 0;

            for (;;)
            {
                i = i + 1;

                if ( in_ <= icount )
                {
                    break;
                }

                if ( in_ < i )
                {
                    i = 0;
                    icount = 0;
                    continue;
                }

                int m = klass[i - 1];
                if (clsize[m - 1] <= 1)
                {
                    icount = icount + 1;
                    continue;
                }

                double inco = -eps;
                int lo = m;
                //
                //  Test the transfer of object I from class M to class L.
                //
                int iswitch;
                int l;
                for (l = 1; l <= ik; l++)
                {
                    if (l != m)
                    {
                        iswitch = 1;
                        double inc = crtran(varval, klass, clsize, in_, ik, iv, critvl,
                            i, m, l, iswitch);
                        //
                        //  Remember the values of L and INC.
                        //
                        if (inc < inco)
                        {
                            lo = l;
                            inco = inc;
                        }
                    }
                }

                icount = icount + 1;
                //
                //  Execute the transfer of object I from class M to class LO.
                //
                if (lo != m)
                {
                    l = lo;
                    critvl = critvl + inco;
                    icount = 0;

                    iswitch = 2;
                    crtran(varval, klass, clsize, in_, ik, iv, critvl,
                        i, m, l, iswitch);

                    ntrans = ntrans + 1;
                    klass[i - 1] = l;
                    clsize[l - 1] = clsize[l - 1] + 1;
                    clsize[m - 1] = clsize[m - 1] - 1;
                }
            }
        }

        public static double crswap ( double[] a, int[] c, int[] c_size, int m, int k,
                                int n, double critvl, int i1, int i2, int c1, int c2, int iswitch )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CRSWAP determines the effect of swapping two objects.
        //
        //  Discussion:
        //
        //    This computation is very inefficient.  It is only set up so that we
        //    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
        //    ASA 136.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Colin Banfield, LC Bassill,
        //    Algorithm AS 113:
        //    A transfer for non-hierarchichal classification,
        //    Applied Statistics,
        //    Volume 26, Number 2, 1977, pages 206-210.
        //
        //  Parameters:
        //
        //    Input, double A(M,N), the data values.  There are M objects,
        //    each having spatial dimension N.
        //
        //    Input, int C(M), the classification of each object.
        //
        //    Input, int C_SIZE(K), the number of objects in each class.
        //
        //    Input, int M, the number of objects.
        //
        //    Input, int K, the number of classes.
        //
        //    Input, int N, the number of spatial dimensions, or variates,
        //    of the objects.
        //
        //    Input, double *CRITVL, the current value of the criterion.
        //
        //    Input, int I1, I2, the objects to be swapped.
        //
        //    Input, int C1, C2, the current classes of objects I1 and I2.
        //
        //    Input, int ISWITCH:
        //    1, indicates that I1 and I2 should be temporarily swapped, the
        //       change in CRITVL should be computed, and then I1 and I2 restored.
        //    2, indicates that I1 and I2 will be swapped.
        //
        //    Output, double CRSWAP, the change to CRITVL that would occur if I1 and
        //    I2 were swapped.  This is only computed for ISWITCH = 1.
        //
        {
            int ci;
            double inc;

            if (iswitch == 2)
            {
                inc = 0.0;
                return inc;
            }

            double[] c_center = new double[k * n];
            //
            //  Move object I1 from class C1 to class C2.
            //  Move object I2 from class C2 to class C1.
            //
            c[i1 - 1] = c2;
            c[i2 - 1] = c1;
            //
            //  Define the critical value as the sum of the squares of the distances
            //  of the points to their cluster center.
            //
            for (int i = 1; i <= k; i++)
            {
                c_size[i - 1] = 0;
                for (int j = 1; j <= n; j++)
                {
                    c_center[i - 1 + (j - 1) * k] = 0.0;
                }
            }

            for (int i = 1; i <= m; i++)
            {
                ci = c[i - 1];
                c_size[ci - 1] = c_size[ci - 1] + 1;
                for (int j = 1; j <= n; j++)
                {
                    c_center[ci - 1 + (j - 1) * k] = c_center[ci - 1 + (j - 1) * k] + a[i - 1 + (j - 1) * m];
                }
            }

            for (int i = 1; i <= k; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    c_center[i - 1 + (j - 1) * k] = c_center[i - 1 + (j - 1) * k] / (double) (c_size[i - 1]);
                }
            }

            double critvl_new = 0.0;

            for (int i = 1; i <= m; i++)
            {
                ci = c[i - 1];
                for (int j = 1; j <= n; j++)
                {
                    critvl_new = critvl_new
                                 + Math.Pow(a[i - 1 + (j - 1) * m] - c_center[ci - 1 + (j - 1) * k], 2);
                }
            }

            inc = critvl_new - critvl;
            //
            //  Move object I1 from class C2 to class C1.
            //  Move object I2 from class C1 to class C2.
            //
            c[i1 - 1] = c1;
            c[i2 - 1] = c2;

            return inc;
        }

        public static double crtran ( double[] a, int[] c, int[] c_size, int m, int k, int n,
                                double critvl, int i1, int c1, int c2, int iswitch )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CRTRAN determines the effect of moving an object to another class.
        //
        //  Discussion:
        //
        //    This computation is very inefficient.  It is only set up so that we
        //    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
        //    ASA 136.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Colin Banfield, LC Bassill,
        //    Algorithm AS 113:
        //    A transfer for non-hierarchichal classification,
        //    Applied Statistics,
        //    Volume 26, Number 2, 1977, pages 206-210.
        //
        //  Parameters:
        //
        //    Input, double AL(M,N), the data values.  There are M objects,
        //    each having spatial dimension N.
        //
        //    Input, int C(M), the classification of each object.
        //
        //    Input, int C_SIZE(K), the number of objects in each class.
        //
        //    Input, int M, the number of objects.
        //
        //    Input, int K, the number of classes.
        //
        //    Input, int N, the number of spatial dimensions, or variates,
        //    of the objects.
        //
        //    Input, double *CRITVL, the current value of the criterion.
        //
        //    Input, int I1, the object to be transferred.
        //
        //    Input, int C1, C2, the current class of object I1, and the
        //    class to which it may be transferred.
        //
        //    Input, int ISWITCH:
        //    1, indicates that I1 should be temporarily transferred, the
        //       change in CRITVL should be computed, and then I1 restored.
        //    2, indicates that I1 will be permanently transferred.
        //
        //    Output, double CRTRAN, the change to CRITVL that would occur if I1 were
        //    transferred from class C1 to C2.  This is only computed for ISWITCH = 1.
        //
        {
            int ci;
            double inc;

            if (iswitch == 2)
            {
                inc = 0.0;
                return inc;
            }

            double[] c_center = new double[k * n];
            //
            //  Move object I from class C1 to class C2.
            //
            c[i1 - 1] = c2;
            c_size[c1 - 1] = c_size[c1 - 1] - 1;
            c_size[c2 - 1] = c_size[c2 - 1] + 1;
            //
            //  Define the critical value as the sum of the squares of the distances
            //  of the points to their cluster center.
            //
            for (int i = 1; i <= k; i++)
            {
                c_size[i - 1] = 0;
                for (int j = 1; j <= n; j++)
                {
                    c_center[i - 1 + (j - 1) * k] = 0.0;
                }
            }

            for (int i = 1; i <= m; i++)
            {
                ci = c[i - 1];
                c_size[ci - 1] = c_size[ci - 1] + 1;
                for (int j = 1; j <= n; j++)
                {
                    c_center[ci - 1 + (j - 1) * k] = c_center[ci - 1 + (j - 1) * k] + a[i - 1 + (j - 1) * m];
                }
            }

            for (int i = 1; i <= k; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    c_center[i - 1 + (j - 1) * k] = c_center[i - 1 + (j - 1) * k] / (double) (c_size[i - 1]);
                }
            }

            double critvl_new = 0.0;

            for (int i = 1; i <= m; i++)
            {
                ci = c[i - 1];
                for (int j = 1; j <= n; j++)
                {
                    critvl_new = critvl_new
                                 + Math.Pow(a[i - 1 + (j - 1) * m] - c_center[ci - 1 + (j - 1) * k], 2);
                }
            }

            inc = critvl_new - critvl;
            //
            //  Move object I1 from class C2 to class C1.
            //
            c[i1 - 1] = c1;
            c_size[c1 - 1] = c_size[c1 - 1] + 1;
            c_size[c2 - 1] = c_size[c2 - 1] - 1;

            return inc;
        }


    }
}