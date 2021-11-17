using System;
using System.Collections.Generic;
using System.IO;

namespace Tiler3D;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TILER_3D illustrates the use of 3D blending.
        //
        //  Discussion:
        //
        //    This main program works in a 3D rectangular space we can think
        //    of as "UVW" space.  We are interested in the space of data values
        //    bounded by (umin,vmin,wmin) and (umax,vmax,wmax), and we plan to
        //    divide this up, using indices (I,J,K), into NI by NJ by NK sub-boxes.
        //
        //    The code below considers each sub-box indexed by (I,J,K) and determines
        //    the values (u0,v0,w0) and (u1,v1,w1) that characterize its corners.
        //    The coordinates of this sub-box and the coordinates of the big box
        //    are then passed to SUB_BOX_TILER_3D.
        //
        //    The picture would be a LOT more interesting if the boundary were
        //    a bit more wiggly, there were more sub-boxes, and the object in
        //    each sub-box had more parts.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Gordon and Charles Hall,
        //    Construction of Curvilinear Coordinate Systems and Application to
        //    Mesh Generation,
        //    International Journal of Numerical Methods in Engineering,
        //    Volume 7, pages 461-477, 1973.
        //
        //    Joe Thompson, Bharat Soni, Nigel Weatherill,
        //    Handbook of Grid Generation,
        //    CRC Press, 1999.
        //
    {
        string file_out_name = "tiler_3d.tri";
        List<string> file_out = new();
        int i;
        int j;
        int k;
        int ni = 3;
        int nj = 3;
        int nk = 3;
        double u0;
        double u1;
        double umax = 30.0;
        double umin = 150.0;
        double v0;
        double v1;
        double vmax = 5.0;
        double vmin = 1.0;
        double w0;
        double w1;
        double wmax = 30.0;
        double wmin = -30.0;


        Console.WriteLine("");
        Console.WriteLine("TILER_3D");
        Console.WriteLine("  A simple example of transfinite interpolation in 3D.");
        //
        //  Write the first line of the output file, which is the number
        //  of triangles in the output 3D shape.  (In our simple example,
        //  each sub-box will contain a tetrahedron.)
        //


        file_out.Add(4 * ni * nj * nk + "");
        //
        //  Consider items with index (I,*,*):
        //
        for (i = 1; i <= ni; i++)
        {

            u0 = ((ni + 1 - i) * umin + (i - 1) * umax)
                 / ni;

            u1 = ((ni - i) * umin + i * umax)
                 / ni;
            //
            //  Consider items with index (I,J,*):
            //
            for (j = 1; j <= nj; j++)
            {

                v0 = ((nj - j + 1) * vmin
                      + (j - 1) * vmax)
                     / nj;

                v1 = ((nj - j) * vmin
                      + j * vmax)
                     / nj;
                //
                //  Consider items with index (I,J,K):
                //
                for (k = 1; k <= nk; k++)
                {

                    w0 = ((nk - k + 1) * wmin
                          + (k - 1) * wmax)
                         / nk;

                    w1 = ((nk - k) * wmin
                          + k * wmax)
                         / nk;
                    //
                    //  Fill sub-box (I,J,K) with the 3-D "tile".
                    //   
                    sub_box_tiler_3d(ref file_out, umin, vmin, wmin, umax, vmax, wmax,
                        u0, v0, w0, u1, v1, w1);
                }
            }
        }

        try
        {
            File.WriteAllLines(file_out_name, file_out);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TILER_3D:");
            Console.WriteLine("  Could not open the output file.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("TILER_3D:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void boundary_3d(double umin, double vmin, double wmin, double umax,
            double vmax, double wmax, double u, double v, double w, ref double x,
            ref double y, ref double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOUNDARY_3D returns the (X,Y,Z) coordinates of a point (U,V,W).
        //
        //  Discussion:
        //
        //    In this example, a single formula describes the mapping
        //
        //      (U,V,W) => (X,Y,Z).
        //
        //    It is more common (and more interesting) for the formula
        //    to depend on which face of the boundary is being considered.
        //
        //    This routine is only called for points (U,V,W) where one of the 
        //    values U, V and W is "extreme", so there are generally six cases
        //    to consider, for general boundaries.  The coding has been set
        //    up with this in mind.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Gordon and Charles Hall,
        //    Construction of Curvilinear Coordinate Systems and Application to
        //    Mesh Generation,
        //    International Journal of Numerical Methods in Engineering,
        //    Volume 7, pages 461-477, 1973.
        //
        //    Joe Thompson, Bharat Soni, Nigel Weatherill,
        //    Handbook of Grid Generation,
        //    CRC Press, 1999.
        //
        //  Parameters:
        //
        //    Input, double UMIN, VMIN, WMIN, UMAX, VMAX, WMAX, the (U,V,W) coordinates
        //    of the two opposite corners of the big box.
        //
        //    Input, double U, V, W, the (U,V,W) coordinates of a point in the big box.
        //
        //    Output, double *X, *Y, *Z, the (X,Y,Z) coordinates of the point.
        //
    {
        double DEG2RAD = Math.PI / 180.0;

        double psi;
        double theta;

        theta = u * DEG2RAD;
        psi = w * DEG2RAD;

        if (u == umin)
        {
            x = v * Math.Cos(theta) * Math.Cos(psi);
            y = v * Math.Sin(theta) * Math.Cos(psi);
            z = v * Math.Sin(psi);
        }
        else if (u == umax)
        {
            x = v * Math.Cos(theta) * Math.Cos(psi);
            y = v * Math.Sin(theta) * Math.Cos(psi);
            z = v * Math.Sin(psi);
        }
        else if (v == vmin)
        {
            x = v * Math.Cos(theta) * Math.Cos(psi);
            y = v * Math.Sin(theta) * Math.Cos(psi);
            z = v * Math.Sin(psi);
        }
        else if (v == vmax)
        {
            x = v * Math.Cos(theta) * Math.Cos(psi);
            y = v * Math.Sin(theta) * Math.Cos(psi);
            z = v * Math.Sin(psi);
        }
        else if (w == wmin)
        {
            x = v * Math.Cos(theta) * Math.Cos(psi);
            y = v * Math.Sin(theta) * Math.Cos(psi);
            z = v * Math.Sin(psi);
        }
        else if (w == wmax)
        {
            x = v * Math.Cos(theta) * Math.Cos(psi);
            y = v * Math.Sin(theta) * Math.Cos(psi);
            z = v * Math.Sin(psi);
        }
        else
        {
            Console.WriteLine("");
            Console.WriteLine("BOUNDARY_3D - Fatal error!");
            Console.WriteLine("  Illegal input value of (U,V,W).");
        }
    }

    private static void sub_box_tiler_3d(ref List<string> file_out, double umin,
            double vmin, double wmin, double umax, double vmax, double wmax,
            double u0, double v0, double w0, double u1, double v1, double w1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUB_BOX_TILER_3D "tiles" a 3D sub-box with a given pattern.
        //
        //  Discussion:
        //
        //    This routine knows the (U,V,W) coordinates of the big box and the
        //    sub box, and knows the shape of the object to be place in the sub-box.
        //    It uses transfinite interpolation to put the shape in the box.
        //    This requires that, for each point of the shape to be mapped, the 
        //    (X,Y,Z) coordinates be evaluated for points on the surface of
        //    the big box, namely, at 8 corners, 12 edges, and 6 faces.  These
        //    values are then blended to give a sensible (X,Y,Z) coordinate for
        //    the point belonging to the shape.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Gordon and Charles Hall,
        //    Construction of Curvilinear Coordinate Systems and Application to
        //    Mesh Generation,
        //    International Journal of Numerical Methods in Engineering,
        //    Volume 7, pages 461-477, 1973.
        //
        //    Joe Thompson, Bharat Soni, Nigel Weatherill,
        //    Handbook of Grid Generation,
        //    CRC Press, 1999.
        //
        //  Parameters:
        //
        //    Input, ofstream reffile_out, a call by reference for the output file.
        //
        //    Input, double UMIN, VMIN, WMIN, UMAX, VMAX, WMAX, the (U,V,W) coordinates
        //    of the two opposite corners of the big box.
        //
        //    Input, double U0, V0, W0, U1, V1, W1, the (U,V,W) coordinates of the
        //    two oppositie corners of the sub-box.
        //
    {
        int NPOINT = 4;
        double R0 = 0.0;
        double R1 = 1.0;
        double S0 = 0.0;
        double S1 = 1.0;
        double T0 = 0.0;
        double T1 = 1.0;

        int i;
        double r = 0;
        double[] r_tab = new double[NPOINT];
        double s = 0;
        double[] s_tab = new double[NPOINT];
        double t = 0;
        double[] t_tab = new double[NPOINT];
        double u = 0;
        double v = 0;
        double w = 0;
        double[] x = new double[NPOINT];
        double x000 = 0;
        double x00t = 0;
        double x001 = 0;
        double x010 = 0;
        double x011 = 0;
        double x01t = 0;
        double x0s0 = 0;
        double x0s1 = 0;
        double x0st = 0;
        double x100 = 0;
        double x101 = 0;
        double x10t = 0;
        double x110 = 0;
        double x111 = 0;
        double x11t = 0;
        double x1s0 = 0;
        double x1s1 = 0;
        double x1st = 0;
        double xr00 = 0;
        double xr01 = 0;
        double xr0t = 0;
        double xr10 = 0;
        double xr11 = 0;
        double xr1t = 0;
        double xrs0 = 0;
        double xrs1 = 0;
        double[] y = new double[NPOINT];
        double y000 = 0;
        double y00t = 0;
        double y001 = 0;
        double y010 = 0;
        double y011 = 0;
        double y01t = 0;
        double y0s0 = 0;
        double y0s1 = 0;
        double y0st = 0;
        double y100 = 0;
        double y101 = 0;
        double y10t = 0;
        double y110 = 0;
        double y111 = 0;
        double y11t = 0;
        double y1s0 = 0;
        double y1s1 = 0;
        double y1st = 0;
        double yr00 = 0;
        double yr01 = 0;
        double yr0t = 0;
        double yr10 = 0;
        double yr11 = 0;
        double yr1t = 0;
        double yrs0 = 0;
        double yrs1 = 0;
        double[] z = new double[NPOINT];
        double z000 = 0;
        double z00t = 0;
        double z001 = 0;
        double z010 = 0;
        double z011 = 0;
        double z01t = 0;
        double z0s0 = 0;
        double z0s1 = 0;
        double z0st = 0;
        double z100 = 0;
        double z101 = 0;
        double z10t = 0;
        double z110 = 0;
        double z111 = 0;
        double z11t = 0;
        double z1s0 = 0;
        double z1s1 = 0;
        double z1st = 0;
        double zr00 = 0;
        double zr01 = 0;
        double zr0t = 0;
        double zr10 = 0;
        double zr11 = 0;
        double zr1t = 0;
        double zrs0 = 0;
        double zrs1 = 0;
        //
        //  Here are the (R,S,T) coordinates of the tetrahedron that we can think
        //  of as the "tile" we need to place in each sub-box.
        //
        //  The (R,S,T) coordinate system is assumed to range from 0 to 1.
        //
        r_tab[0] = 0.2;
        s_tab[0] = 0.2;
        t_tab[0] = 0.2;

        r_tab[1] = 0.8;
        s_tab[1] = 0.2;
        t_tab[1] = 0.2;

        r_tab[2] = 0.5;
        s_tab[2] = 0.8;
        t_tab[2] = 0.2;

        r_tab[3] = 0.5;
        s_tab[3] = 0.5;
        t_tab[3] = 0.8;
        //
        //  Compute the (X,Y,Z) coordinates of the corners of the (U,V,W) box.
        //  These really only need to be computed once ever.
        //
        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umin, vmin, wmin, ref x000, ref y000, ref z000);

        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umax, vmin, wmin, ref x100, ref y100, ref z100);

        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umin, vmax, wmin, ref x010, ref y010, ref z010);

        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umax, vmax, wmin, ref x110, ref y110, ref z110);

        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umin, vmin, wmax, ref x001, ref y001, ref z001);

        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umax, vmin, wmax, ref x101, ref y101, ref z101);

        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umin, vmax, wmax, ref x011, ref y011, ref z011);

        boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
            umax, vmax, wmax, ref x111, ref y111, ref z111);
        //
        //  Now figure out the (X,Y,Z) coordinates of the tile point with
        //  given (R,S,T) coordinates.  This depends on the positions of 
        //  all sorts of points on the corners, edges, and faces of the big box.
        //
        for (i = 0; i < NPOINT; i++)
        {
            //
            //  Get the (R,S,T) coordinates of point I.
            //
            r = r_tab[i];
            s = s_tab[i];
            t = t_tab[i];
            //
            //  Get the corresponding point (U,V,W) in the rectangular space.
            //
            u = ((R1 - r) * u0
                 + (r - R0) * u1)
                / (R1 - R0);

            v = ((S1 - s) * v0
                 + (s - S0) * v1)
                / (S1 - S0);

            w = ((T1 - t) * w0
                 + (t - T0) * w1)
                / (T1 - T0);
            //
            //  Evaluate (X,Y,Z) on the 12 edges "near" (U,V,W).
            //
            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umin, vmin, w, ref x00t, ref y00t, ref z00t);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umin, vmax, w, ref x01t, ref y01t, ref z01t);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umax, vmin, w, ref x10t, ref y10t, ref z10t);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umax, vmax, w, ref x11t, ref y11t, ref z11t);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umin, v, wmin, ref x0s0, ref y0s0, ref z0s0);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umin, v, wmax, ref x0s1, ref y0s1, ref z0s1);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umax, v, wmin, ref x1s0, ref y1s0, ref z1s0);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umax, v, wmax, ref x1s1, ref y1s1, ref z1s1);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, vmin, wmin, ref xr00, ref yr00, ref zr00);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, vmin, wmax, ref xr01, ref yr01, ref zr01);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, vmax, wmin, ref xr10, ref yr10, ref zr10);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, vmax, wmax, ref xr11, ref yr11, ref zr11);
            //
            //  Evaluate (X,Y,Z) on the six faces near (U,V,W).
            //
            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umin, v, w, ref x0st, ref y0st, ref z0st);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                umax, v, w, ref x1st, ref y1st, ref z1st);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, vmin, w, ref xr0t, ref yr0t, ref zr0t);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, vmax, w, ref xr1t, ref yr1t, ref zr1t);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, v, wmin, ref xrs0, ref yrs0, ref zrs0);

            boundary_3d(umin, vmin, wmin, umax, vmax, wmax,
                u, v, wmax, ref xrs1, ref yrs1, ref zrs1);
            //
            //  Now figure out the location of 
            //    point I => (R,S,T) => (U,V,W) => (X(I),Y(I),Z(I)).
            //
            x[i] = (
                (umax - u) * (vmax - v) * (wmax - w) * x000
                - (umax - u) * (vmax - v) * (wmax - wmin) * x00t
                + (umax - u) * (vmax - v) * (w - wmin) * x001
                - (umax - u) * (vmax - vmin) * (wmax - w) * x0s0
                + (umax - u) * (vmax - vmin) * (wmax - wmin) * x0st
                - (umax - u) * (vmax - vmin) * (w - wmin) * x0s1
                + (umax - u) * (v - vmin) * (wmax - w) * x010
                - (umax - u) * (v - vmin) * (wmax - wmin) * x01t
                + (umax - u) * (v - vmin) * (w - wmin) * x011
                - (umax - umin) * (vmax - v) * (wmax - w) * xr00
                + (umax - umin) * (vmax - v) * (wmax - wmin) * xr0t
                - (umax - umin) * (vmax - v) * (w - wmin) * xr01
                + (umax - umin) * (vmax - vmin) * (wmax - w) * xrs0
                + (umax - umin) * (vmax - vmin) * (w - wmin) * xrs1
                - (umax - umin) * (v - vmin) * (wmax - w) * xr10
                + (umax - umin) * (v - vmin) * (wmax - wmin) * xr1t
                - (umax - umin) * (v - vmin) * (w - wmin) * xr11
                + (u - umin) * (vmax - v) * (wmax - w) * x100
                - (u - umin) * (vmax - v) * (wmax - wmin) * x10t
                + (u - umin) * (vmax - v) * (w - wmin) * x101
                - (u - umin) * (vmax - vmin) * (wmax - w) * x1s0
                + (u - umin) * (vmax - vmin) * (wmax - wmin) * x1st
                - (u - umin) * (vmax - vmin) * (w - wmin) * x1s1
                + (u - umin) * (v - vmin) * (wmax - w) * x110
                - (u - umin) * (v - vmin) * (wmax - wmin) * x11t
                + (u - umin) * (v - vmin) * (w - wmin) * x111
            ) / ((umax - umin) * (vmax - vmin) * (wmax - wmin));

            y[i] = (
                (umax - u) * (vmax - v) * (wmax - w) * y000
                - (umax - u) * (vmax - v) * (wmax - wmin) * y00t
                + (umax - u) * (vmax - v) * (w - wmin) * y001
                - (umax - u) * (vmax - vmin) * (wmax - w) * y0s0
                + (umax - u) * (vmax - vmin) * (wmax - wmin) * y0st
                - (umax - u) * (vmax - vmin) * (w - wmin) * y0s1
                + (umax - u) * (v - vmin) * (wmax - w) * y010
                - (umax - u) * (v - vmin) * (wmax - wmin) * y01t
                + (umax - u) * (v - vmin) * (w - wmin) * y011
                - (umax - umin) * (vmax - v) * (wmax - w) * yr00
                + (umax - umin) * (vmax - v) * (wmax - wmin) * yr0t
                - (umax - umin) * (vmax - v) * (w - wmin) * yr01
                + (umax - umin) * (vmax - vmin) * (wmax - w) * yrs0
                + (umax - umin) * (vmax - vmin) * (w - wmin) * yrs1
                - (umax - umin) * (v - vmin) * (wmax - w) * yr10
                + (umax - umin) * (v - vmin) * (wmax - wmin) * yr1t
                - (umax - umin) * (v - vmin) * (w - wmin) * yr11
                + (u - umin) * (vmax - v) * (wmax - w) * y100
                - (u - umin) * (vmax - v) * (wmax - wmin) * y10t
                + (u - umin) * (vmax - v) * (w - wmin) * y101
                - (u - umin) * (vmax - vmin) * (wmax - w) * y1s0
                + (u - umin) * (vmax - vmin) * (wmax - wmin) * y1st
                - (u - umin) * (vmax - vmin) * (w - wmin) * y1s1
                + (u - umin) * (v - vmin) * (wmax - w) * y110
                - (u - umin) * (v - vmin) * (wmax - wmin) * y11t
                + (u - umin) * (v - vmin) * (w - wmin) * y111
            ) / ((umax - umin) * (vmax - vmin) * (wmax - wmin));

            z[i] = (
                (umax - u) * (vmax - v) * (wmax - w) * z000
                - (umax - u) * (vmax - v) * (wmax - wmin) * z00t
                + (umax - u) * (vmax - v) * (w - wmin) * z001
                - (umax - u) * (vmax - vmin) * (wmax - w) * z0s0
                + (umax - u) * (vmax - vmin) * (wmax - wmin) * z0st
                - (umax - u) * (vmax - vmin) * (w - wmin) * z0s1
                + (umax - u) * (v - vmin) * (wmax - w) * z010
                - (umax - u) * (v - vmin) * (wmax - wmin) * z01t
                + (umax - u) * (v - vmin) * (w - wmin) * z011
                - (umax - umin) * (vmax - v) * (wmax - w) * zr00
                + (umax - umin) * (vmax - v) * (wmax - wmin) * zr0t
                - (umax - umin) * (vmax - v) * (w - wmin) * zr01
                + (umax - umin) * (vmax - vmin) * (wmax - w) * zrs0
                + (umax - umin) * (vmax - vmin) * (w - wmin) * zrs1
                - (umax - umin) * (v - vmin) * (wmax - w) * zr10
                + (umax - umin) * (v - vmin) * (wmax - wmin) * zr1t
                - (umax - umin) * (v - vmin) * (w - wmin) * zr11
                + (u - umin) * (vmax - v) * (wmax - w) * z100
                - (u - umin) * (vmax - v) * (wmax - wmin) * z10t
                + (u - umin) * (vmax - v) * (w - wmin) * z101
                - (u - umin) * (vmax - vmin) * (wmax - w) * z1s0
                + (u - umin) * (vmax - vmin) * (wmax - wmin) * z1st
                - (u - umin) * (vmax - vmin) * (w - wmin) * z1s1
                + (u - umin) * (v - vmin) * (wmax - w) * z110
                - (u - umin) * (v - vmin) * (wmax - wmin) * z11t
                + (u - umin) * (v - vmin) * (w - wmin) * z111
            ) / ((umax - umin) * (vmax - vmin) * (wmax - wmin));

        }

        file_out.Add(x[0] + "  " + y[0] + "  " + z[0] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[1] + "  " + y[1] + "  " + z[1] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[2] + "  " + y[2] + "  " + z[2] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");

        file_out.Add(x[1] + "  " + y[1] + "  " + z[1] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[0] + "  " + y[0] + "  " + z[0] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[3] + "  " + y[3] + "  " + z[3] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");

        file_out.Add(x[2] + "  " + y[2] + "  " + z[2] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[3] + "  " + y[3] + "  " + z[3] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[0] + "  " + y[0] + "  " + z[0] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");

        file_out.Add(x[3] + "  " + y[3] + "  " + z[3] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[2] + "  " + y[2] + "  " + z[2] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");
        file_out.Add(x[1] + "  " + y[1] + "  " + z[1] + "  "
                     + 0.0 + "  " + 0.0 + "  " + 0.0 + "");

    }
}