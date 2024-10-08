﻿/*
** SGI FREE SOFTWARE LICENSE B (Version 2.0, Sept. 18, 2008) 
** Copyright (C) 2011 Silicon Graphics, Inc.
** All Rights Reserved.
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
** of the Software, and to permit persons to whom the Software is furnished to do so,
** subject to the following conditions:
** 
** The above copyright notice including the dates of first publication and either this
** permission notice or a reference to http://oss.sgi.com/projects/FreeB/ shall be
** included in all copies or substantial portions of the Software. 
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
** INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
** PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL SILICON GRAPHICS, INC.
** BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
** TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
** OR OTHER DEALINGS IN THE SOFTWARE.
** 
** Except as contained in this notice, the name of Silicon Graphics, Inc. shall not
** be used in advertising or otherwise to promote the sale, use or other dealings in
** this Software without prior written authorization from Silicon Graphics, Inc.
*/
/*
** Original Author: Eric Veach, July 1994.
** libtess2: Mikko Mononen, http://code.google.com/p/libtess2/.
** LibTessDotNet: Remi Gillig, https://github.com/speps/LibTessDotNet
*/

using System;
using System.Diagnostics;

namespace LibTessDotNet.Double;

internal static class Geom
{
    public static bool IsWindingInside(WindingRule rule, int n)
    {
        return rule switch
        {
            WindingRule.EvenOdd => (n & 1) == 1,
            WindingRule.NonZero => n != 0,
            WindingRule.Positive => n > 0,
            WindingRule.Negative => n < 0,
            WindingRule.AbsGeqTwo => n is >= 2 or <= -2,
            _ => throw new Exception("Wrong winding rule")
        };
    }

    public static bool VertCCW(MeshUtils.Vertex u, MeshUtils.Vertex v, MeshUtils.Vertex w)
    {
        return u._s * (v._t - w._t) + v._s * (w._t - u._t) + w._s * (u._t - v._t) >= 0.0f;
    }
    public static bool VertEq(MeshUtils.Vertex lhs, MeshUtils.Vertex rhs)
    {
        return Math.Abs(lhs._s - rhs._s) < float.Epsilon && Math.Abs(lhs._t - rhs._t) < float.Epsilon;
    }
    public static bool VertLeq(MeshUtils.Vertex lhs, MeshUtils.Vertex rhs)
    {
        return lhs._s < rhs._s || Math.Abs(lhs._s - rhs._s) < float.Epsilon && lhs._t <= rhs._t;
    }

    /// <summary>
    /// Given three vertices u,v,w such that VertLeq(u,v) && VertLeq(v,w),
    /// evaluates the t-coord of the edge uw at the s-coord of the vertex v.
    /// Returns v->t - (uw)(v->s), ie. the signed distance from uw to v.
    /// If uw is vertical (and thus passes through v), the result is zero.
    /// 
    /// The calculation is extremely accurate and stable, even when v
    /// is very close to u or w.  In particular if we set v->t = 0 and
    /// let r be the negated result (this evaluates (uw)(v->s)), then
    /// r is guaranteed to satisfy MIN(u->t,w->t) &lt;= r &lt;= MAX(u->t,w->t).
    /// </summary>
    public static double EdgeEval(MeshUtils.Vertex u, MeshUtils.Vertex v, MeshUtils.Vertex w)
    {
        Debug.Assert(VertLeq(u, v) && VertLeq(v, w));

        double gapL = v._s - u._s;
        double gapR = w._s - v._s;

        return (gapL + gapR) switch
        {
            > 0.0f when gapL < gapR => v._t - u._t + (u._t - w._t) * (gapL / (gapL + gapR)),
            > 0.0f => v._t - w._t + (w._t - u._t) * (gapR / (gapL + gapR)),
            _ => 0
        };
    }

    /// <summary>
    /// Returns a number whose sign matches EdgeEval(u,v,w) but which
    /// is cheaper to evaluate. Returns &gt; 0, == 0 , or &lt; 0
    /// as v is above, on, or below the edge uw.
    /// </summary>
    public static double EdgeSign(MeshUtils.Vertex u, MeshUtils.Vertex v, MeshUtils.Vertex w)
    {
        Debug.Assert(VertLeq(u, v) && VertLeq(v, w));

        double gapL = v._s - u._s;
        double gapR = w._s - v._s;

        return (gapL + gapR) switch
        {
            > 0.0f => (v._t - w._t) * gapL + (v._t - u._t) * gapR,
            _ => 0
        };
    }

    public static bool TransLeq(MeshUtils.Vertex lhs, MeshUtils.Vertex rhs)
    {
        return lhs._t < rhs._t || Math.Abs(lhs._t - rhs._t) < float.Epsilon && lhs._s <= rhs._s;
    }

    public static double TransEval(MeshUtils.Vertex u, MeshUtils.Vertex v, MeshUtils.Vertex w)
    {
        Debug.Assert(TransLeq(u, v) && TransLeq(v, w));

        double gapL = v._t - u._t;
        double gapR = w._t - v._t;

        return (gapL + gapR) switch
        {
            > 0.0f when gapL < gapR => v._s - u._s + (u._s - w._s) * (gapL / (gapL + gapR)),
            > 0.0f => v._s - w._s + (w._s - u._s) * (gapR / (gapL + gapR)),
            _ => 0
        };
    }

    public static double TransSign(MeshUtils.Vertex u, MeshUtils.Vertex v, MeshUtils.Vertex w)
    {
        Debug.Assert(TransLeq(u, v) && TransLeq(v, w));

        double gapL = v._t - u._t;
        double gapR = w._t - v._t;

        return (gapL + gapR) switch
        {
            > 0.0f => (v._s - w._s) * gapL + (v._s - u._s) * gapR,
            _ => 0
        };
    }

    public static bool EdgeGoesLeft(MeshUtils.Edge e)
    {
        return VertLeq(e._Dst, e._Org);
    }

    public static bool EdgeGoesRight(MeshUtils.Edge e)
    {
        return VertLeq(e._Org, e._Dst);
    }

    public static double VertL1dist(MeshUtils.Vertex u, MeshUtils.Vertex v)
    {
        return Math.Abs(u._s - v._s) + Math.Abs(u._t - v._t);
    }

    public static void AddWinding(MeshUtils.Edge eDst, MeshUtils.Edge eSrc)
    {
        eDst._winding += eSrc._winding;
        eDst._Sym._winding += eSrc._Sym._winding;
    }

    public static double Interpolate(double a, double x, double b, double y)
    {
        a = a switch
        {
            < 0.0f => 0.0f,
            _ => a
        };
        b = b switch
        {
            < 0.0f => 0.0f,
            _ => b
        };
        return a <= b ? b == 0.0f ? (x + y) / 2.0f
            : x + (y - x) * (a / (a + b))
            : y + (x - y) * (b / (a + b));
    }

    private static void Swap(ref MeshUtils.Vertex a, ref MeshUtils.Vertex b)
    {
        (a, b) = (b, a);
    }

    /// <summary>
    /// Given edges (o1,d1) and (o2,d2), compute their point of intersection.
    /// The computed point is guaranteed to lie in the intersection of the
    /// bounding rectangles defined by each edge.
    /// </summary>
    public static void EdgeIntersect(MeshUtils.Vertex o1, MeshUtils.Vertex d1, MeshUtils.Vertex o2, MeshUtils.Vertex d2, MeshUtils.Vertex v)
    {
        // This is certainly not the most efficient way to find the intersection
        // of two line segments, but it is very numerically stable.
        // 
        // Strategy: find the two middle vertices in the VertLeq ordering,
        // and interpolate the intersection s-value from these.  Then repeat
        // using the TransLeq ordering to find the intersection t-value.

        if (!VertLeq(o1, d1)) { Swap(ref o1, ref d1); }
        if (!VertLeq(o2, d2)) { Swap(ref o2, ref d2); }
        if (!VertLeq(o1, o2)) { Swap(ref o1, ref o2); Swap(ref d1, ref d2); }

        if (!VertLeq(o2, d1))
        {
            // Technically, no intersection -- do our best
            v._s = (o2._s + d1._s) / 2.0f;
        }
        else if (VertLeq(d1, d2))
        {
            // Interpolate between o2 and d1
            double z1 = EdgeEval(o1, o2, d1);
            double z2 = EdgeEval(o2, d1, d2);
            switch (z1 + z2)
            {
                case < 0.0f:
                    z1 = -z1;
                    z2 = -z2;
                    break;
            }
            v._s = Interpolate(z1, o2._s, z2, d1._s);
        }
        else
        {
            // Interpolate between o2 and d2
            double z1 = EdgeSign(o1, o2, d1);
            double z2 = -EdgeSign(o1, d2, d1);
            switch (z1 + z2)
            {
                case < 0.0f:
                    z1 = -z1;
                    z2 = -z2;
                    break;
            }
            v._s = Interpolate(z1, o2._s, z2, d2._s);
        }

        // Now repeat the process for t

        if (!TransLeq(o1, d1)) { Swap(ref o1, ref d1); }
        if (!TransLeq(o2, d2)) { Swap(ref o2, ref d2); }
        if (!TransLeq(o1, o2)) { Swap(ref o1, ref o2); Swap(ref d1, ref d2); }

        if (!TransLeq(o2, d1))
        {
            // Technically, no intersection -- do our best
            v._t = (o2._t + d1._t) / 2.0f;
        }
        else if (TransLeq(d1, d2))
        {
            // Interpolate between o2 and d1
            double z1 = TransEval(o1, o2, d1);
            double z2 = TransEval(o2, d1, d2);
            switch (z1 + z2)
            {
                case < 0.0f:
                    z1 = -z1;
                    z2 = -z2;
                    break;
            }
            v._t = Interpolate(z1, o2._t, z2, d1._t);
        }
        else
        {
            // Interpolate between o2 and d2
            double z1 = TransSign(o1, o2, d1);
            double z2 = -TransSign(o1, d2, d1);
            switch (z1 + z2)
            {
                case < 0.0f:
                    z1 = -z1;
                    z2 = -z2;
                    break;
            }
            v._t = Interpolate(z1, o2._t, z2, d2._t);
        }
    }
}