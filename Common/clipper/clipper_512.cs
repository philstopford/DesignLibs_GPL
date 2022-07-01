/*******************************************************************************
*                                                                              *
* Author    :  Angus Johnson                                                   *
* Version   :  6.4.2                                                           *
* Date      :  27 February 2017                                                *
* Website   :  http://www.angusj.com                                           *
* Copyright :  Angus Johnson 2010-2017                                         *
*                                                                              *
* License:                                                                     *
* Use, modification & distribution is subject to Boost Software License Ver 1. *
* http://www.boost.org/LICENSE_1_0.txt                                         *
*                                                                              *
* Attributions:                                                                *
* The code in this library is an extension of Bala Vatti's clipping algorithm: *
* "A generic solution to polygon clipping"                                     *
* Communications of the ACM, Vol 35, Issue 7 (July 1992) pp 56-63.             *
* http://portal.acm.org/citation.cfm?id=129906                                 *
*                                                                              *
* Computer graphics and geometric modeling: implementation and algorithms      *
* By Max K. Agoston                                                            *
* Springer; 1 edition (January 4, 2005)                                        *
* http://books.google.com/books?q=vatti+clipping+agoston                       *
*                                                                              *
* See also:                                                                    *
* "Polygon Offsetting by Computing Winding Numbers"                            *
* Paper no. DETC2005-85513 pp. 565-575                                         *
* ASME 2005 International Design Engineering Technical Conferences             *
* and Computers and Information in Engineering Conference (IDETC/CIE2005)      *
* September 24-28, 2005 , Long Beach, California, USA                          *
* http://www.me.berkeley.edu/~mcmains/pubs/DAC05OffsetPolygon.pdf              *
*                                                                              *
*******************************************************************************/

/*******************************************************************************
*                                                                              *
* This is a translation of the Delphi Clipper library and the naming style     *
* used has retained a Delphi flavour.                                          *
*                                                                              *
*******************************************************************************/

//use_int32: When enabled 32bit ints are used instead of 64bit ints. This
//improve performance but coordinate values are limited to the range +/- 46340
//#define use_int32

//use_xyz: adds a Z member to IntPoint. Adds a minor cost to performance.

#define use_xyz

//use_lines: Enables open path clipping. Adds a very minor cost to performance.
#define use_lines


using System;
using System.Collections.Generic;
using System.Linq;

//using System.Text;          //for Int128.AsString() & StringBuilder
//using System.IO;            //debugging with streamReader & StreamWriter
//using System.Windows.Forms; //debugging to clipboard

namespace ClipperLib1;

#if use_int32
  using cInt = Int32;
#else
using cInt = Int64;
#endif

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

public struct DoublePoint
{
    public double X;
    public double Y;

    public DoublePoint(double x = 0, double y = 0)
    {
        X = x;
        Y = y;
    }

    public DoublePoint(DoublePoint dp)
    {
        X = dp.X;
        Y = dp.Y;
    }

    public DoublePoint(IntPoint ip)
    {
        X = ip.X;
        Y = ip.Y;
    }
}

//------------------------------------------------------------------------------
// PolyTree & PolyNode classes
//------------------------------------------------------------------------------

public class PolyTree : PolyNode
{
    internal List<PolyNode> m_AllPolys = new();

    //The GC probably handles this cleanup more efficiently ...
    //~PolyTree(){Clear();}

    public void Clear()
    {
        for (int i = 0; i < m_AllPolys.Count; i++)
            m_AllPolys[i] = null;
        m_AllPolys.Clear();
        m_Childs.Clear();
    }

    public PolyNode GetFirst()
    {
        return m_Childs.Count switch
        {
            > 0 => m_Childs[0],
            _ => null
        };
    }

    public int Total
    {
        get
        {
            int result = m_AllPolys.Count;
            switch (result)
            {
                //with negative offsets, ignore the hidden outer polygon ...
                case > 0 when m_Childs[0] != m_AllPolys[0]:
                    result--;
                    break;
            }

            return result;
        }
    }
}

public class PolyNode
{
    internal PolyNode m_Parent;
    internal Path m_polygon = new();
    internal int m_Index;
    internal JoinType m_jointype;
    internal EndType m_endtype;
    internal List<PolyNode> m_Childs = new();

    private bool IsHoleNode()
    {
        bool result = true;
        PolyNode node = m_Parent;
        while (node != null)
        {
            result = !result;
            node = node.m_Parent;
        }

        return result;
    }

    public int ChildCount => m_Childs.Count;

    public Path Contour => m_polygon;

    internal void AddChild(PolyNode Child)
    {
        int cnt = m_Childs.Count;
        m_Childs.Add(Child);
        Child.m_Parent = this;
        Child.m_Index = cnt;
    }

    public PolyNode GetNext()
    {
        return m_Childs.Count switch
        {
            > 0 => m_Childs[0],
            _ => GetNextSiblingUp()
        };
    }

    internal PolyNode GetNextSiblingUp()
    {
        switch (m_Parent)
        {
            case null:
                return null;
        }

        return m_Index == m_Parent.m_Childs.Count - 1 ? m_Parent.GetNextSiblingUp() : m_Parent.m_Childs[m_Index + 1];
    }

    public List<PolyNode> Childs => m_Childs;

    public PolyNode Parent => m_Parent;

    public bool IsHole => IsHoleNode();

    public bool IsOpen { get; set; }
}

//------------------------------------------------------------------------------
// Int128 struct (enables safe math on signed 64bit integers)
// eg Int128 val1((Int64)9223372036854775807); //ie 2^63 -1
//    Int128 val2((Int64)9223372036854775807);
//    Int128 val3 = val1 * val2;
//    val3.ToString => "85070591730234615847396907784232501249" (8.5e+37)
//------------------------------------------------------------------------------

internal struct Int128
{
    private long hi;
    private ulong lo;

    public Int128(long _lo)
    {
        lo = (ulong) _lo;
        hi = _lo switch
        {
            < 0 => -1,
            _ => 0
        };
    }

    public Int128(long _hi, ulong _lo)
    {
        lo = _lo;
        hi = _hi;
    }

    public Int128(Int128 val)
    {
        hi = val.hi;
        lo = val.lo;
    }

    public bool IsNegative()
    {
        return hi < 0;
    }

    public static bool operator ==(Int128 val1, Int128 val2)
    {
        if ((object) val1 == (object) val2)
        {
            return true;
        }

        if ((object) val1 == null || (object) val2 == null)
        {
            return false;
        }

        return val1.hi == val2.hi && val1.lo == val2.lo;
    }

    public static bool operator !=(Int128 val1, Int128 val2)
    {
        return !(val1 == val2);
    }

    public override bool Equals(object obj)
    {
        if (obj is not Int128 i128)
        {
            return false;
        }

        return i128.hi == hi && i128.lo == lo;
    }

    public override int GetHashCode()
    {
        return hi.GetHashCode() ^ lo.GetHashCode();
    }

    public static bool operator >(Int128 val1, Int128 val2)
    {
        if (val1.hi != val2.hi)
        {
            return val1.hi > val2.hi;
        }

        return val1.lo > val2.lo;
    }

    public static bool operator <(Int128 val1, Int128 val2)
    {
        if (val1.hi != val2.hi)
        {
            return val1.hi < val2.hi;
        }

        return val1.lo < val2.lo;
    }

    public static Int128 operator +(Int128 lhs, Int128 rhs)
    {
        lhs.hi += rhs.hi;
        lhs.lo += rhs.lo;
        if (lhs.lo < rhs.lo)
        {
            lhs.hi++;
        }

        return lhs;
    }

    public static Int128 operator -(Int128 lhs, Int128 rhs)
    {
        return lhs + -rhs;
    }

    public static Int128 operator -(Int128 val)
    {
        return val.lo switch
        {
            0 => new Int128(-val.hi, 0),
            _ => new Int128(~val.hi, ~val.lo + 1)
        };
    }

    public static explicit operator double(Int128 val)
    {
        const double shift64 = 18446744073709551616.0; //2^64
        return val.hi switch
        {
            < 0 when val.lo == 0 => val.hi * shift64,
            < 0 => -(~val.lo + ~val.hi * shift64),
            _ => val.lo + val.hi * shift64
        };
    }

    //nb: Constructing two new Int128 objects every time we want to multiply longs  
    //is slow. So, although calling the Int128Mul method doesn't look as clean, the 
    //code runs significantly faster than if we'd used the * operator.

    public static Int128 Int128Mul(long lhs, long rhs)
    {
        bool negate = lhs < 0 != rhs < 0;
        lhs = lhs switch
        {
            < 0 => -lhs,
            _ => lhs
        };

        rhs = rhs switch
        {
            < 0 => -rhs,
            _ => rhs
        };

        ulong int1Hi = (ulong) lhs >> 32;
        ulong int1Lo = (ulong) lhs & 0xFFFFFFFF;
        ulong int2Hi = (ulong) rhs >> 32;
        ulong int2Lo = (ulong) rhs & 0xFFFFFFFF;

        //nb: see comments in clipper.pas
        ulong a = int1Hi * int2Hi;
        ulong b = int1Lo * int2Lo;
        ulong c = int1Hi * int2Lo + int1Lo * int2Hi;

        ulong lo;
        long hi = (long) (a + (c >> 32));

        unchecked
        {
            lo = (c << 32) + b;
        }

        if (lo < b)
        {
            hi++;
        }

        Int128 result = new(hi, lo);
        return negate ? -result : result;
    }
}

//------------------------------------------------------------------------------

public struct IntPoint
{
    public long X;
    public long Y;
#if use_xyz
    public long Z;

    public IntPoint(long x, long y, long z = 0)
    {
        X = x;
        Y = y;
        Z = z;
    }

    public IntPoint(double x, double y, double z = 0)
    {
        X = (long) x;
        Y = (long) y;
        Z = (long) z;
    }

    public IntPoint(DoublePoint dp)
    {
        X = (long) dp.X;
        Y = (long) dp.Y;
        Z = 0;
    }

    public IntPoint(IntPoint pt)
    {
        X = pt.X;
        Y = pt.Y;
        Z = pt.Z;
    }
#else
    public IntPoint(cInt X, cInt Y)
    {
        this.X = X; this.Y = Y;
    }
    public IntPoint(double x, double y)
    {
      this.X = (cInt)x; this.Y = (cInt)y;
    }

    public IntPoint(IntPoint pt)
    {
        this.X = pt.X; this.Y = pt.Y;
    }
#endif

    public static bool operator ==(IntPoint a, IntPoint b)
    {
        return a.X == b.X && a.Y == b.Y;
    }

    public static bool operator !=(IntPoint a, IntPoint b)
    {
        return a.X != b.X || a.Y != b.Y;
    }

    public override bool Equals(object obj)
    {
        switch (obj)
        {
            case null:
                return false;
            case IntPoint point:
            {
                return X == point.X && Y == point.Y;
            }
            default:
                return false;
        }
    }

    public override int GetHashCode()
    {
        //simply prevents a compiler warning
        return base.GetHashCode();
    }
} // end struct IntPoint

public struct IntRect
{
    public long left;
    public long top;
    public long right;
    public long bottom;

    public IntRect(long l, long t, long r, long b)
    {
        left = l;
        top = t;
        right = r;
        bottom = b;
    }

    public IntRect(IntRect ir)
    {
        left = ir.left;
        top = ir.top;
        right = ir.right;
        bottom = ir.bottom;
    }
}

public enum ClipType
{
    ctIntersection,
    ctUnion,
    ctDifference,
    ctXor
}

public enum PolyType
{
    ptSubject,
    ptClip
}

//By far the most widely used winding rules for polygon filling are
//EvenOdd & NonZero (GDI, GDI+, XLib, OpenGL, Cairo, AGG, Quartz, SVG, Gr32)
//Others rules include Positive, Negative and ABS_GTR_EQ_TWO (only in OpenGL)
//see http://glprogramming.com/red/chapter11.html
public enum PolyFillType
{
    pftEvenOdd,
    pftNonZero,
    pftPositive,
    pftNegative
}

public enum JoinType
{
    jtSquare,
    jtRound,
    jtMiter
}

public enum EndType
{
    etClosedPolygon,
    etClosedLine,
    etOpenButt,
    etOpenSquare,
    etOpenRound
}

internal enum EdgeSide
{
    esLeft,
    esRight
}

internal enum Direction
{
    dRightToLeft,
    dLeftToRight
}

internal class TEdge
{
    internal IntPoint Bot;
    internal IntPoint Curr; //current (updated for every new scanbeam)
    internal IntPoint Top;
    internal IntPoint Delta;
    internal double Dx;
    internal PolyType PolyTyp;
    internal EdgeSide Side; //side only refers to current side of solution poly
    internal int WindDelta; //1 or -1 depending on winding direction
    internal int WindCnt;
    internal int WindCnt2; //winding count of the opposite polytype
    internal int OutIdx;
    internal TEdge Next;
    internal TEdge Prev;
    internal TEdge NextInLML;
    internal TEdge NextInAEL;
    internal TEdge PrevInAEL;
    internal TEdge NextInSEL;
    internal TEdge PrevInSEL;
}

public class IntersectNode
{
    internal TEdge Edge1;
    internal TEdge Edge2;
    internal IntPoint Pt;
}

public class MyIntersectNodeSort : IComparer<IntersectNode>
{
    public int Compare(IntersectNode node1, IntersectNode node2)
    {
        long i = node2.Pt.Y - node1.Pt.Y;
        return i switch
        {
            > 0 => 1,
            < 0 => -1,
            _ => 0
        };
    }
}

internal class LocalMinima
{
    internal long Y;
    internal TEdge LeftBound;
    internal TEdge RightBound;
    internal LocalMinima Next;
}

internal class Scanbeam
{
    internal long Y;
    internal Scanbeam Next;
}

internal class Maxima
{
    internal long X;
    internal Maxima Next;
    internal Maxima Prev;
}

//OutRec: contains a path in the clipping solution. Edges in the AEL will
//carry a pointer to an OutRec when they are part of the clipping solution.
internal class OutRec
{
    internal int Idx;
    internal bool IsHole;
    internal bool IsOpen;
    internal OutRec FirstLeft; //see comments in clipper.pas
    internal OutPt Pts;
    internal OutPt BottomPt;
    internal PolyNode PolyNode;
}

internal class OutPt
{
    internal int Idx;
    internal IntPoint Pt;
    internal OutPt Next;
    internal OutPt Prev;
}

internal class Join
{
    internal OutPt OutPt1;
    internal OutPt OutPt2;
    internal IntPoint OffPt;
}

public class ClipperBase
{
    internal const double horizontal = -3.4E+38;
    internal const int Skip = -2;
    internal const int Unassigned = -1;
    internal const double tolerance = 1.0E-20;

    internal static bool near_zero(double val)
    {
        return val is > -tolerance and < tolerance;
    }

#if use_int32
    public const cInt loRange = 0x7FFF;
    public const cInt hiRange = 0x7FFF;
#else
    public const long loRange = 0x3FFFFFFF;
    public const long hiRange = 0x3FFFFFFFFFFFFFFFL;
#endif

    internal LocalMinima m_MinimaList;
    internal LocalMinima m_CurrentLM;
    internal List<List<TEdge>> m_edges = new();
    internal Scanbeam m_Scanbeam;
    internal List<OutRec> m_PolyOuts;
    internal TEdge m_ActiveEdges;
    internal bool m_UseFullRange;
    internal bool m_HasOpenPaths;

    //------------------------------------------------------------------------------

    public bool PreserveCollinear { get; set; }
    //------------------------------------------------------------------------------

    public void Swap(ref long val1, ref long val2)
    {
        (val1, val2) = (val2, val1);
    }
    //------------------------------------------------------------------------------

    internal static bool IsHorizontal(TEdge e)
    {
        return e.Delta.Y == 0;
    }
    //------------------------------------------------------------------------------

    internal bool PointIsVertex(IntPoint pt, OutPt pp)
    {
        OutPt pp2 = pp;
        do
        {
            if (pp2.Pt == pt)
            {
                return true;
            }

            pp2 = pp2.Next;
        } while (pp2 != pp);

        return false;
    }
    //------------------------------------------------------------------------------

    internal bool PointOnLineSegment(IntPoint pt,
        IntPoint linePt1, IntPoint linePt2, bool UseFullRange)
    {
        return UseFullRange switch
        {
            true => pt.X == linePt1.X && pt.Y == linePt1.Y || pt.X == linePt2.X && pt.Y == linePt2.Y ||
                    pt.X > linePt1.X == pt.X < linePt2.X && pt.Y > linePt1.Y == pt.Y < linePt2.Y &&
                    Int128.Int128Mul(pt.X - linePt1.X, linePt2.Y - linePt1.Y) ==
                    Int128.Int128Mul(linePt2.X - linePt1.X, pt.Y - linePt1.Y),
            _ => pt.X == linePt1.X && pt.Y == linePt1.Y || pt.X == linePt2.X && pt.Y == linePt2.Y ||
                 pt.X > linePt1.X == pt.X < linePt2.X && pt.Y > linePt1.Y == pt.Y < linePt2.Y &&
                 (pt.X - linePt1.X) * (linePt2.Y - linePt1.Y) == (linePt2.X - linePt1.X) * (pt.Y - linePt1.Y)
        };
    }
    //------------------------------------------------------------------------------

    internal bool PointOnPolygon(IntPoint pt, OutPt pp, bool UseFullRange)
    {
        OutPt pp2 = pp;
        while (true)
        {
            if (PointOnLineSegment(pt, pp2.Pt, pp2.Next.Pt, UseFullRange))
            {
                return true;
            }

            pp2 = pp2.Next;
            if (pp2 == pp)
            {
                break;
            }
        }

        return false;
    }
    //------------------------------------------------------------------------------

    internal static bool SlopesEqual(TEdge e1, TEdge e2, bool UseFullRange)
    {
        return UseFullRange switch
        {
            true => Int128.Int128Mul(e1.Delta.Y, e2.Delta.X) == Int128.Int128Mul(e1.Delta.X, e2.Delta.Y),
            _ => e1.Delta.Y * e2.Delta.X == e1.Delta.X * e2.Delta.Y
        };
    }
    //------------------------------------------------------------------------------

    internal static bool SlopesEqual(IntPoint pt1, IntPoint pt2,
        IntPoint pt3, bool UseFullRange)
    {
        return UseFullRange switch
        {
            true => Int128.Int128Mul(pt1.Y - pt2.Y, pt2.X - pt3.X) == Int128.Int128Mul(pt1.X - pt2.X, pt2.Y - pt3.Y),
            _ => (pt1.Y - pt2.Y) * (pt2.X - pt3.X) - (pt1.X - pt2.X) * (pt2.Y - pt3.Y) == 0
        };
    }
    //------------------------------------------------------------------------------

    internal static bool SlopesEqual(IntPoint pt1, IntPoint pt2,
        IntPoint pt3, IntPoint pt4, bool UseFullRange)
    {
        return UseFullRange switch
        {
            true => Int128.Int128Mul(pt1.Y - pt2.Y, pt3.X - pt4.X) == Int128.Int128Mul(pt1.X - pt2.X, pt3.Y - pt4.Y),
            _ => (pt1.Y - pt2.Y) * (pt3.X - pt4.X) - (pt1.X - pt2.X) * (pt3.Y - pt4.Y) == 0
        };
    }
    //------------------------------------------------------------------------------

    internal ClipperBase() //constructor (nb: no external instantiation)
    {
        m_MinimaList = null;
        m_CurrentLM = null;
        m_UseFullRange = false;
        m_HasOpenPaths = false;
    }
    //------------------------------------------------------------------------------

    public virtual void Clear()
    {
        DisposeLocalMinimaList();
        foreach (List<TEdge> t in m_edges)
        {
            for (int j = 0; j < t.Count; ++j) t[j] = null;
            t.Clear();
        }

        m_edges.Clear();
        m_UseFullRange = false;
        m_HasOpenPaths = false;
    }
    //------------------------------------------------------------------------------

    private void DisposeLocalMinimaList()
    {
        while (m_MinimaList != null)
        {
            LocalMinima tmpLm = m_MinimaList.Next;
            m_MinimaList = null;
            m_MinimaList = tmpLm;
        }

        m_CurrentLM = null;
    }
    //------------------------------------------------------------------------------

    private void RangeTest(IntPoint Pt, ref bool useFullRange)
    {
        while (true)
        {
            switch (useFullRange)
            {
                case true:
                {
                    if (Pt.X > hiRange || Pt.Y > hiRange || -Pt.X > hiRange || -Pt.Y > hiRange)
                    {
                        throw new ClipperException("Coordinate outside allowed range");
                    }

                    break;
                }
                default:
                {
                    if (Pt.X > loRange || Pt.Y > loRange || -Pt.X > loRange || -Pt.Y > loRange)
                    {
                        useFullRange = true;
                        continue;
                    }

                    break;
                }
            }

            break;
        }
    }
    //------------------------------------------------------------------------------

    private void InitEdge(TEdge e, TEdge eNext,
        TEdge ePrev, IntPoint pt)
    {
        e.Next = eNext;
        e.Prev = ePrev;
        e.Curr = pt;
        e.OutIdx = Unassigned;
    }
    //------------------------------------------------------------------------------

    private void InitEdge2(TEdge e, PolyType polyType)
    {
        if (e.Curr.Y >= e.Next.Curr.Y)
        {
            e.Bot = e.Curr;
            e.Top = e.Next.Curr;
        }
        else
        {
            e.Top = e.Curr;
            e.Bot = e.Next.Curr;
        }

        SetDx(e);
        e.PolyTyp = polyType;
    }
    //------------------------------------------------------------------------------

    private TEdge FindNextLocMin(TEdge E)
    {
        for (;;)
        {
            while (E.Bot != E.Prev.Bot || E.Curr == E.Top)
            {
                E = E.Next;
            }

            if (Math.Abs(E.Dx - horizontal) > double.Epsilon && Math.Abs(E.Prev.Dx - horizontal) > double.Epsilon)
            {
                break;
            }

            while (Math.Abs(E.Prev.Dx - horizontal) <= double.Epsilon)
            {
                E = E.Prev;
            }

            TEdge E2 = E;
            while (Math.Abs(E.Dx - horizontal) <= double.Epsilon)
            {
                E = E.Next;
            }

            if (E.Top.Y == E.Prev.Bot.Y)
            {
                continue; //ie just an intermediate horz.
            }

            if (E2.Prev.Bot.X < E.Bot.X)
            {
                E = E2;
            }

            break;
        }

        return E;
    }
    //------------------------------------------------------------------------------

    private TEdge ProcessBound(TEdge E, bool LeftBoundIsForward)
    {
        TEdge EStart, Result = E;
        TEdge Horz;

        switch (Result.OutIdx)
        {
            case Skip:
            {
                //check if there are edges beyond the skip edge in the bound and if so
                //create another LocMin and calling ProcessBound once more ...
                E = Result;
                switch (LeftBoundIsForward)
                {
                    case true:
                    {
                        while (E.Top.Y == E.Next.Bot.Y)
                        {
                            E = E.Next;
                        }

                        while (E != Result && Math.Abs(E.Dx - horizontal) <= double.Epsilon)
                        {
                            E = E.Prev;
                        }

                        break;
                    }
                    default:
                    {
                        while (E.Top.Y == E.Prev.Bot.Y)
                        {
                            E = E.Prev;
                        }

                        while (E != Result && Math.Abs(E.Dx - horizontal) <= double.Epsilon)
                        {
                            E = E.Next;
                        }

                        break;
                    }
                }

                if (E == Result)
                {
                    Result = LeftBoundIsForward switch
                    {
                        true => E.Next,
                        _ => E.Prev
                    };
                }
                else
                {
                    E = LeftBoundIsForward switch
                    {
                        //there are more edges in the bound beyond result starting with E
                        true => Result.Next,
                        _ => Result.Prev
                    };

                    LocalMinima locMin = new()
                    {
                        Next = null,
                        Y = E.Bot.Y,
                        LeftBound = null,
                        RightBound = E
                    };
                    E.WindDelta = 0;
                    Result = ProcessBound(E, LeftBoundIsForward);
                    InsertLocalMinima(locMin);
                }

                return Result;
            }
        }

        switch (Math.Abs(E.Dx - horizontal))
        {
            case <= double.Epsilon:
            {
                EStart = LeftBoundIsForward switch
                {
                    //We need to be careful with open paths because this may not be a
                    //true local minima (ie E may be following a skip edge).
                    //Also, consecutive horz. edges may start heading left before going right.
                    true => E.Prev,
                    _ => E.Next
                };

                switch (Math.Abs(EStart.Dx - horizontal))
                {
                    //ie an adjoining horizontal skip edge
                    case <= double.Epsilon:
                    {
                        if (EStart.Bot.X != E.Bot.X && EStart.Top.X != E.Bot.X)
                        {
                            ReverseHorizontal(E);
                        }

                        break;
                    }
                    default:
                    {
                        if (EStart.Bot.X != E.Bot.X)
                        {
                            ReverseHorizontal(E);
                        }

                        break;
                    }
                }

                break;
            }
        }

        EStart = E;
        switch (LeftBoundIsForward)
        {
            case true:
            {
                while (Result.Top.Y == Result.Next.Bot.Y && Result.Next.OutIdx != Skip)
                {
                    Result = Result.Next;
                }

                switch (Math.Abs(Result.Dx - horizontal))
                {
                    case <= double.Epsilon when Result.Next.OutIdx != Skip:
                    {
                        //nb: at the top of a bound, horizontals are added to the bound
                        //only when the preceding edge attaches to the horizontal's left vertex
                        //unless a Skip edge is encountered when that becomes the top divide
                        Horz = Result;
                        while (Math.Abs(Horz.Prev.Dx - horizontal) <= double.Epsilon)
                        {
                            Horz = Horz.Prev;
                        }

                        if (Horz.Prev.Top.X > Result.Next.Top.X)
                        {
                            Result = Horz.Prev;
                        }

                        break;
                    }
                }

                while (E != Result)
                {
                    E.NextInLML = E.Next;
                    switch (Math.Abs(E.Dx - horizontal))
                    {
                        case <= double.Epsilon when E != EStart && E.Bot.X != E.Prev.Top.X:
                            ReverseHorizontal(E);
                            break;
                    }

                    E = E.Next;
                }

                switch (Math.Abs(E.Dx - horizontal))
                {
                    case <= double.Epsilon when E != EStart && E.Bot.X != E.Prev.Top.X:
                        ReverseHorizontal(E);
                        break;
                }

                Result = Result.Next; //move to the edge just beyond current bound
                break;
            }
            default:
            {
                while (Result.Top.Y == Result.Prev.Bot.Y && Result.Prev.OutIdx != Skip)
                {
                    Result = Result.Prev;
                }

                switch (Math.Abs(Result.Dx - horizontal))
                {
                    case <= double.Epsilon when Result.Prev.OutIdx != Skip:
                    {
                        Horz = Result;
                        while (Math.Abs(Horz.Next.Dx - horizontal) <= double.Epsilon)
                        {
                            Horz = Horz.Next;
                        }

                        if (Horz.Next.Top.X == Result.Prev.Top.X ||
                            Horz.Next.Top.X > Result.Prev.Top.X)
                        {
                            Result = Horz.Next;
                        }

                        break;
                    }
                }

                while (E != Result)
                {
                    E.NextInLML = E.Prev;
                    switch (Math.Abs(E.Dx - horizontal))
                    {
                        case <= double.Epsilon when E != EStart && E.Bot.X != E.Next.Top.X:
                            ReverseHorizontal(E);
                            break;
                    }

                    E = E.Prev;
                }

                switch (Math.Abs(E.Dx - horizontal))
                {
                    case <= double.Epsilon when E != EStart && E.Bot.X != E.Next.Top.X:
                        ReverseHorizontal(E);
                        break;
                }

                Result = Result.Prev; //move to the edge just beyond current bound
                break;
            }
        }

        return Result;
    }
    //------------------------------------------------------------------------------


    public bool AddPath(Path pg, PolyType polyType, bool Closed)
    {
#if use_lines
        switch (Closed)
        {
            case false when polyType == PolyType.ptClip:
                throw new ClipperException("AddPath: Open paths must be subject.");
        }
#else
      if (!Closed)
        throw new ClipperException("AddPath: Open paths have been disabled.");
#endif

        int highI = pg.Count - 1;
        switch (Closed)
        {
            case true:
            {
                while (highI > 0 && pg[highI] == pg[0])
                {
                    --highI;
                }

                break;
            }
        }

        while (highI > 0 && pg[highI] == pg[highI - 1])
        {
            --highI;
        }

        switch (Closed)
        {
            case true when highI < 2:
            case false when highI < 1:
                return false;
        }

        //create a new edge array ...
        List<TEdge> edges = new(highI + 1);
        for (int i = 0; i <= highI; i++) edges.Add(new TEdge());

        bool IsFlat = true;

        //1. Basic (first) edge initialization ...
        edges[1].Curr = pg[1];
        RangeTest(pg[0], ref m_UseFullRange);
        RangeTest(pg[highI], ref m_UseFullRange);
        InitEdge(edges[0], edges[1], edges[highI], pg[0]);
        InitEdge(edges[highI], edges[0], edges[highI - 1], pg[highI]);
        for (int i = highI - 1; i >= 1; --i)
        {
            RangeTest(pg[i], ref m_UseFullRange);
            InitEdge(edges[i], edges[i + 1], edges[i - 1], pg[i]);
        }

        TEdge eStart = edges[0];

        //2. Remove duplicate vertices, and (when closed) collinear edges ...
        TEdge E = eStart, eLoopStop = eStart;
        for (;;)
        {
            //nb: allows matching start and end points when not Closed ...
            if (E.Curr == E.Next.Curr && (Closed || E.Next != eStart))
            {
                if (E == E.Next)
                {
                    break;
                }

                if (E == eStart)
                {
                    eStart = E.Next;
                }

                E = RemoveEdge(E);
                eLoopStop = E;
                continue;
            }

            if (E.Prev == E.Next)
            {
                break; //only two vertices
            }

            switch (Closed)
            {
                case true when SlopesEqual(E.Prev.Curr, E.Curr, E.Next.Curr, m_UseFullRange) && (!PreserveCollinear ||
                    !Pt2IsBetweenPt1AndPt3(E.Prev.Curr, E.Curr, E.Next.Curr)):
                {
                    //Collinear edges are allowed for open paths but in closed paths
                    //the default is to merge adjacent collinear edges into a single edge.
                    //However, if the PreserveCollinear property is enabled, only overlapping
                    //collinear edges (ie spikes) will be removed from closed paths.
                    if (E == eStart)
                    {
                        eStart = E.Next;
                    }

                    E = RemoveEdge(E);
                    E = E.Prev;
                    eLoopStop = E;
                    continue;
                }
            }

            E = E.Next;
            if (E == eLoopStop || !Closed && E.Next == eStart)
            {
                break;
            }
        }

        switch (Closed)
        {
            case false when E == E.Next:
            case true when E.Prev == E.Next:
                return false;
            case false:
                m_HasOpenPaths = true;
                eStart.Prev.OutIdx = Skip;
                break;
        }

        //3. Do second stage of edge initialization ...
        E = eStart;
        do
        {
            InitEdge2(E, polyType);
            E = E.Next;
            IsFlat = IsFlat switch
            {
                true when E.Curr.Y != eStart.Curr.Y => false,
                _ => IsFlat
            };
        } while (E != eStart);

        switch (IsFlat)
        {
            //Totally flat paths must be handled differently when adding them
            //to LocalMinima list to avoid endless loops etc ...
            case true when Closed:
                return false;
            //4. Finally, add edge bounds to LocalMinima list ...
            case true:
            {
                E.Prev.OutIdx = Skip;
                LocalMinima locMin = new()
                {
                    Next = null,
                    Y = E.Bot.Y,
                    LeftBound = null,
                    RightBound = E
                };
                locMin.RightBound.Side = EdgeSide.esRight;
                locMin.RightBound.WindDelta = 0;
                for (;;)
                {
                    if (E.Bot.X != E.Prev.Top.X)
                    {
                        ReverseHorizontal(E);
                    }

                    if (E.Next.OutIdx == Skip)
                    {
                        break;
                    }

                    E.NextInLML = E.Next;
                    E = E.Next;
                }

                InsertLocalMinima(locMin);
                m_edges.Add(edges);
                return true;
            }
        }

        m_edges.Add(edges);
        bool leftBoundIsForward;
        TEdge EMin = null;

        //workaround to avoid an endless loop in the while loop below when
        //open paths have matching start and end points ...
        if (E.Prev.Bot == E.Prev.Top)
        {
            E = E.Next;
        }

        for (;;)
        {
            E = FindNextLocMin(E);
            if (E == EMin)
            {
                break;
            }

            EMin = EMin switch
            {
                null => E,
                _ => EMin
            };

            //E and E.Prev now share a local minima (left aligned if horizontal).
            //Compare their slopes to find which starts which bound ...
            LocalMinima locMin = new()
            {
                Next = null,
                Y = E.Bot.Y
            };
            if (E.Dx < E.Prev.Dx)
            {
                locMin.LeftBound = E.Prev;
                locMin.RightBound = E;
                leftBoundIsForward = false; //Q.nextInLML = Q.prev
            }
            else
            {
                locMin.LeftBound = E;
                locMin.RightBound = E.Prev;
                leftBoundIsForward = true; //Q.nextInLML = Q.next
            }

            locMin.LeftBound.Side = EdgeSide.esLeft;
            locMin.RightBound.Side = EdgeSide.esRight;

            switch (Closed)
            {
                case false:
                    locMin.LeftBound.WindDelta = 0;
                    break;
                default:
                {
                    if (locMin.LeftBound.Next == locMin.RightBound)
                    {
                        locMin.LeftBound.WindDelta = -1;
                    }
                    else
                    {
                        locMin.LeftBound.WindDelta = 1;
                    }

                    break;
                }
            }

            locMin.RightBound.WindDelta = -locMin.LeftBound.WindDelta;

            E = ProcessBound(locMin.LeftBound, leftBoundIsForward);
            switch (E.OutIdx)
            {
                case Skip:
                    E = ProcessBound(E, leftBoundIsForward);
                    break;
            }

            TEdge E2 = ProcessBound(locMin.RightBound, !leftBoundIsForward);
            switch (E2.OutIdx)
            {
                case Skip:
                    E2 = ProcessBound(E2, !leftBoundIsForward);
                    break;
            }

            switch (locMin.LeftBound.OutIdx)
            {
                case Skip:
                    locMin.LeftBound = null;
                    break;
                default:
                {
                    locMin.RightBound = locMin.RightBound.OutIdx switch
                    {
                        Skip => null,
                        _ => locMin.RightBound
                    };

                    break;
                }
            }

            InsertLocalMinima(locMin);
            E = leftBoundIsForward switch
            {
                false => E2,
                _ => E
            };
        }

        return true;
    }
    //------------------------------------------------------------------------------

    public bool AddPaths(Paths ppg, PolyType polyType, bool closed)
    {
        bool result = false;
        foreach (Path t in ppg.Where(t => AddPath(t, polyType, closed)))
            result = true;

        return result;
    }
    //------------------------------------------------------------------------------

    internal bool Pt2IsBetweenPt1AndPt3(IntPoint pt1, IntPoint pt2, IntPoint pt3)
    {
        if (pt1 == pt3 || pt1 == pt2 || pt3 == pt2)
        {
            return false;
        }

        if (pt1.X != pt3.X)
        {
            return pt2.X > pt1.X == pt2.X < pt3.X;
        }

        return pt2.Y > pt1.Y == pt2.Y < pt3.Y;
    }
    //------------------------------------------------------------------------------

    private TEdge RemoveEdge(TEdge e)
    {
        //removes e from double_linked_list (but without removing from memory)
        e.Prev.Next = e.Next;
        e.Next.Prev = e.Prev;
        TEdge result = e.Next;
        e.Prev = null; //flag as removed (see ClipperBase.Clear)
        return result;
    }
    //------------------------------------------------------------------------------

    private void SetDx(TEdge e)
    {
        e.Delta.X = e.Top.X - e.Bot.X;
        e.Delta.Y = e.Top.Y - e.Bot.Y;
        e.Dx = e.Delta.Y switch
        {
            0 => horizontal,
            _ => (double) e.Delta.X / e.Delta.Y
        };
    }
    //---------------------------------------------------------------------------

    private void InsertLocalMinima(LocalMinima newLm)
    {
        switch (m_MinimaList)
        {
            case null:
                m_MinimaList = newLm;
                break;
            default:
            {
                if (newLm.Y >= m_MinimaList.Y)
                {
                    newLm.Next = m_MinimaList;
                    m_MinimaList = newLm;
                }
                else
                {
                    LocalMinima tmpLm = m_MinimaList;
                    while (tmpLm.Next != null && newLm.Y < tmpLm.Next.Y)
                    {
                        tmpLm = tmpLm.Next;
                    }

                    newLm.Next = tmpLm.Next;
                    tmpLm.Next = newLm;
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    internal bool PopLocalMinima(long Y, out LocalMinima current)
    {
        current = m_CurrentLM;
        if (m_CurrentLM == null || m_CurrentLM.Y != Y)
        {
            return false;
        }

        m_CurrentLM = m_CurrentLM.Next;
        return true;
    }
    //------------------------------------------------------------------------------

    private void ReverseHorizontal(TEdge e)
    {
        //swap horizontal edges' top and bottom x's so they follow the natural
        //progression of the bounds - ie so their xbots will align with the
        //adjoining lower edge. [Helpful in the ProcessHorizontal() method.]
        Swap(ref e.Top.X, ref e.Bot.X);
#if use_xyz
        Swap(ref e.Top.Z, ref e.Bot.Z);
#endif
    }
    //------------------------------------------------------------------------------

    internal virtual void Reset()
    {
        m_CurrentLM = m_MinimaList;
        switch (m_CurrentLM)
        {
            case null:
                return; //ie nothing to process
        }

        //reset all edges ...
        m_Scanbeam = null;
        LocalMinima lm = m_MinimaList;
        while (lm != null)
        {
            InsertScanbeam(lm.Y);
            TEdge e = lm.LeftBound;
            if (e != null)
            {
                e.Curr = e.Bot;
                e.OutIdx = Unassigned;
            }

            e = lm.RightBound;
            if (e != null)
            {
                e.Curr = e.Bot;
                e.OutIdx = Unassigned;
            }

            lm = lm.Next;
        }

        m_ActiveEdges = null;
    }
    //------------------------------------------------------------------------------

    public static IntRect GetBounds(Paths paths)
    {
        int i = 0, cnt = paths.Count;
        while (i < cnt && paths[i].Count == 0)
        {
            i++;
        }

        if (i == cnt)
        {
            return new IntRect(0, 0, 0, 0);
        }

        IntRect result = new()
        {
            left = paths[i][0].X
        };
        result.right = result.left;
        result.top = paths[i][0].Y;
        result.bottom = result.top;
        for (; i < cnt; i++)
        for (int j = 0; j < paths[i].Count; j++)
        {
            if (paths[i][j].X < result.left)
            {
                result.left = paths[i][j].X;
            }
            else if (paths[i][j].X > result.right)
            {
                result.right = paths[i][j].X;
            }

            if (paths[i][j].Y < result.top)
            {
                result.top = paths[i][j].Y;
            }
            else if (paths[i][j].Y > result.bottom)
            {
                result.bottom = paths[i][j].Y;
            }
        }

        return result;
    }
    //------------------------------------------------------------------------------

    internal void InsertScanbeam(long Y)
    {
        switch (m_Scanbeam)
        {
            //single-linked list: sorted descending, ignoring dups.
            case null:
                m_Scanbeam = new Scanbeam
                {
                    Next = null,
                    Y = Y
                };
                break;
            default:
            {
                if (Y > m_Scanbeam.Y)
                {
                    Scanbeam newSb = new()
                    {
                        Y = Y,
                        Next = m_Scanbeam
                    };
                    m_Scanbeam = newSb;
                }
                else
                {
                    Scanbeam sb2 = m_Scanbeam;
                    while (sb2.Next != null && Y <= sb2.Next.Y)
                    {
                        sb2 = sb2.Next;
                    }

                    if (Y == sb2.Y)
                    {
                        return; //ie ignores duplicates
                    }

                    Scanbeam newSb = new()
                    {
                        Y = Y,
                        Next = sb2.Next
                    };
                    sb2.Next = newSb;
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    internal bool PopScanbeam(out long Y)
    {
        switch (m_Scanbeam)
        {
            case null:
                Y = 0;
                return false;
        }

        Y = m_Scanbeam.Y;
        m_Scanbeam = m_Scanbeam.Next;
        return true;
    }
    //------------------------------------------------------------------------------

    internal bool LocalMinimaPending()
    {
        return m_CurrentLM != null;
    }
    //------------------------------------------------------------------------------

    internal OutRec CreateOutRec()
    {
        OutRec result = new()
        {
            Idx = Unassigned,
            IsHole = false,
            IsOpen = false,
            FirstLeft = null,
            Pts = null,
            BottomPt = null,
            PolyNode = null
        };
        m_PolyOuts.Add(result);
        result.Idx = m_PolyOuts.Count - 1;
        return result;
    }
    //------------------------------------------------------------------------------

    internal void DisposeOutRec(int index)
    {
        OutRec outRec = m_PolyOuts[index];
        outRec.Pts = null;
        m_PolyOuts[index] = null;
    }
    //------------------------------------------------------------------------------

    internal void UpdateEdgeIntoAEL(ref TEdge e)
    {
        switch (e.NextInLML)
        {
            case null:
                throw new ClipperException("UpdateEdgeIntoAEL: invalid call");
        }

        TEdge AelPrev = e.PrevInAEL;
        TEdge AelNext = e.NextInAEL;
        e.NextInLML.OutIdx = e.OutIdx;
        if (AelPrev != null)
        {
            AelPrev.NextInAEL = e.NextInLML;
        }
        else
        {
            m_ActiveEdges = e.NextInLML;
        }

        if (AelNext != null)
        {
            AelNext.PrevInAEL = e.NextInLML;
        }

        e.NextInLML.Side = e.Side;
        e.NextInLML.WindDelta = e.WindDelta;
        e.NextInLML.WindCnt = e.WindCnt;
        e.NextInLML.WindCnt2 = e.WindCnt2;
        e = e.NextInLML;
        e.Curr = e.Bot;
        e.PrevInAEL = AelPrev;
        e.NextInAEL = AelNext;
        if (!IsHorizontal(e))
        {
            InsertScanbeam(e.Top.Y);
        }
    }
    //------------------------------------------------------------------------------

    internal void SwapPositionsInAEL(TEdge edge1, TEdge edge2)
    {
        //check that one or other edge hasn't already been removed from AEL ...
        if (edge1.NextInAEL == edge1.PrevInAEL ||
            edge2.NextInAEL == edge2.PrevInAEL)
        {
            return;
        }

        if (edge1.NextInAEL == edge2)
        {
            TEdge next = edge2.NextInAEL;
            if (next != null)
            {
                next.PrevInAEL = edge1;
            }

            TEdge prev = edge1.PrevInAEL;
            if (prev != null)
            {
                prev.NextInAEL = edge2;
            }

            edge2.PrevInAEL = prev;
            edge2.NextInAEL = edge1;
            edge1.PrevInAEL = edge2;
            edge1.NextInAEL = next;
        }
        else if (edge2.NextInAEL == edge1)
        {
            TEdge next = edge1.NextInAEL;
            if (next != null)
            {
                next.PrevInAEL = edge2;
            }

            TEdge prev = edge2.PrevInAEL;
            if (prev != null)
            {
                prev.NextInAEL = edge1;
            }

            edge1.PrevInAEL = prev;
            edge1.NextInAEL = edge2;
            edge2.PrevInAEL = edge1;
            edge2.NextInAEL = next;
        }
        else
        {
            TEdge next = edge1.NextInAEL;
            TEdge prev = edge1.PrevInAEL;
            edge1.NextInAEL = edge2.NextInAEL;
            if (edge1.NextInAEL != null)
            {
                edge1.NextInAEL.PrevInAEL = edge1;
            }

            edge1.PrevInAEL = edge2.PrevInAEL;
            if (edge1.PrevInAEL != null)
            {
                edge1.PrevInAEL.NextInAEL = edge1;
            }

            edge2.NextInAEL = next;
            if (edge2.NextInAEL != null)
            {
                edge2.NextInAEL.PrevInAEL = edge2;
            }

            edge2.PrevInAEL = prev;
            if (edge2.PrevInAEL != null)
            {
                edge2.PrevInAEL.NextInAEL = edge2;
            }
        }

        switch (edge1.PrevInAEL)
        {
            case null:
                m_ActiveEdges = edge1;
                break;
            default:
            {
                m_ActiveEdges = edge2.PrevInAEL switch
                {
                    null => edge2,
                    _ => m_ActiveEdges
                };

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    internal void DeleteFromAEL(TEdge e)
    {
        TEdge AelPrev = e.PrevInAEL;
        TEdge AelNext = e.NextInAEL;
        switch (AelPrev)
        {
            case null when AelNext == null && e != m_ActiveEdges:
                return; //already deleted
        }

        if (AelPrev != null)
        {
            AelPrev.NextInAEL = AelNext;
        }
        else
        {
            m_ActiveEdges = AelNext;
        }

        if (AelNext != null)
        {
            AelNext.PrevInAEL = AelPrev;
        }

        e.NextInAEL = null;
        e.PrevInAEL = null;
    }
    //------------------------------------------------------------------------------
} //end ClipperBase

public class Clipper : ClipperBase
{
    //InitOptions that can be passed to the constructor ...
    public const int ioReverseSolution = 1;
    public const int ioStrictlySimple = 2;
    public const int ioPreserveCollinear = 4;

    private ClipType m_ClipType;
    private Maxima m_Maxima;
    private TEdge m_SortedEdges;
    private List<IntersectNode> m_IntersectList;
    private IComparer<IntersectNode> m_IntersectNodeComparer;
    private bool m_ExecuteLocked;
    private PolyFillType m_ClipFillType;
    private PolyFillType m_SubjFillType;
    private List<Join> m_Joins;
    private List<Join> m_GhostJoins;
    private bool m_UsingPolyTree;
#if use_xyz
    public delegate void ZFillCallback(IntPoint bot1, IntPoint top1,
        IntPoint bot2, IntPoint top2, ref IntPoint pt);

    public ZFillCallback ZFillFunction { get; set; }
#endif
    public Clipper(int InitOptions = 0) //constructor
    {
        m_Scanbeam = null;
        m_Maxima = null;
        m_ActiveEdges = null;
        m_SortedEdges = null;
        m_IntersectList = new List<IntersectNode>();
        m_IntersectNodeComparer = new MyIntersectNodeSort();
        m_ExecuteLocked = false;
        m_UsingPolyTree = false;
        m_PolyOuts = new List<OutRec>();
        m_Joins = new List<Join>();
        m_GhostJoins = new List<Join>();
        ReverseSolution = (ioReverseSolution & InitOptions) != 0;
        StrictlySimple = (ioStrictlySimple & InitOptions) != 0;
        PreserveCollinear = (ioPreserveCollinear & InitOptions) != 0;
#if use_xyz
        ZFillFunction = null;
#endif
    }
    //------------------------------------------------------------------------------

    private void InsertMaxima(long X)
    {
        //double-linked list: sorted ascending, ignoring dups.
        Maxima newMax = new()
        {
            X = X
        };
        switch (m_Maxima)
        {
            case null:
                m_Maxima = newMax;
                m_Maxima.Next = null;
                m_Maxima.Prev = null;
                break;
            default:
            {
                if (X < m_Maxima.X)
                {
                    newMax.Next = m_Maxima;
                    newMax.Prev = null;
                    m_Maxima = newMax;
                }
                else
                {
                    Maxima m = m_Maxima;
                    while (m.Next != null && X >= m.Next.X)
                    {
                        m = m.Next;
                    }

                    if (X == m.X)
                    {
                        return; //ie ignores duplicates (& CG to clean up newMax)
                    }

                    //insert newMax between m and m.Next ...
                    newMax.Next = m.Next;
                    newMax.Prev = m;
                    if (m.Next != null)
                    {
                        m.Next.Prev = newMax;
                    }

                    m.Next = newMax;
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    public bool ReverseSolution { get; set; }
    //------------------------------------------------------------------------------

    public bool StrictlySimple { get; set; }
    //------------------------------------------------------------------------------

    public bool Execute(ClipType clipType, Paths solution,
        PolyFillType FillType = PolyFillType.pftEvenOdd)
    {
        return Execute(clipType, solution, FillType, FillType);
    }
    //------------------------------------------------------------------------------

    public bool Execute(ClipType clipType, PolyTree polytree,
        PolyFillType FillType = PolyFillType.pftEvenOdd)
    {
        return Execute(clipType, polytree, FillType, FillType);
    }
    //------------------------------------------------------------------------------

    public bool Execute(ClipType clipType, Paths solution,
        PolyFillType subjFillType, PolyFillType clipFillType)
    {
        switch (m_ExecuteLocked)
        {
            case true:
                return false;
        }

        switch (m_HasOpenPaths)
        {
            case true:
                throw
                    new ClipperException("Error: PolyTree struct is needed for open path clipping.");
        }

        m_ExecuteLocked = true;
        solution.Clear();
        m_SubjFillType = subjFillType;
        m_ClipFillType = clipFillType;
        m_ClipType = clipType;
        m_UsingPolyTree = false;
        bool succeeded;
        try
        {
            succeeded = ExecuteInternal();
            switch (succeeded)
            {
                //build the return polygons ...
                case true:
                    BuildResult(solution);
                    break;
            }
        }
        finally
        {
            DisposeAllPolyPts();
            m_ExecuteLocked = false;
        }

        return succeeded;
    }
    //------------------------------------------------------------------------------

    public bool Execute(ClipType clipType, PolyTree polytree,
        PolyFillType subjFillType, PolyFillType clipFillType)
    {
        switch (m_ExecuteLocked)
        {
            case true:
                return false;
        }

        m_ExecuteLocked = true;
        m_SubjFillType = subjFillType;
        m_ClipFillType = clipFillType;
        m_ClipType = clipType;
        m_UsingPolyTree = true;
        bool succeeded;
        try
        {
            succeeded = ExecuteInternal();
            switch (succeeded)
            {
                //build the return polygons ...
                case true:
                    BuildResult2(polytree);
                    break;
            }
        }
        finally
        {
            DisposeAllPolyPts();
            m_ExecuteLocked = false;
        }

        return succeeded;
    }
    //------------------------------------------------------------------------------

    internal void FixHoleLinkage(OutRec outRec)
    {
        //skip if an outermost polygon or
        //already already points to the correct FirstLeft ...
        if (outRec.FirstLeft == null ||
            outRec.IsHole != outRec.FirstLeft.IsHole &&
            outRec.FirstLeft.Pts != null)
        {
            return;
        }

        OutRec orfl = outRec.FirstLeft;
        while (orfl != null && (orfl.IsHole == outRec.IsHole || orfl.Pts == null))
        {
            orfl = orfl.FirstLeft;
        }

        outRec.FirstLeft = orfl;
    }
    //------------------------------------------------------------------------------

    private bool ExecuteInternal()
    {
        try
        {
            Reset();
            m_SortedEdges = null;
            m_Maxima = null;

            if (!PopScanbeam(out long botY))
            {
                return false;
            }

            InsertLocalMinimaIntoAEL(botY);
            while (PopScanbeam(out long topY) || LocalMinimaPending())
            {
                ProcessHorizontals();
                m_GhostJoins.Clear();
                if (!ProcessIntersections(topY))
                {
                    return false;
                }

                ProcessEdgesAtTopOfScanbeam(topY);
                botY = topY;
                InsertLocalMinimaIntoAEL(botY);
            }

            //fix orientations ...
            foreach (OutRec outRec in m_PolyOuts.Where(outRec => outRec.Pts != null && !outRec.IsOpen)
                         .Where(outRec => (outRec.IsHole ^ ReverseSolution) == Area(outRec) > 0))
            {
                ReversePolyPtLinks(outRec.Pts);
            }

            JoinCommonEdges();

            foreach (OutRec outRec in m_PolyOuts)
            {
                switch (outRec.Pts)
                {
                    case null:
                        continue;
                }

                switch (outRec.IsOpen)
                {
                    case true:
                        FixupOutPolyline(outRec);
                        break;
                    default:
                        FixupOutPolygon(outRec);
                        break;
                }
            }

            switch (StrictlySimple)
            {
                case true:
                    DoSimplePolygons();
                    break;
            }

            return true;
        }
        //catch { return false; }
        finally
        {
            m_Joins.Clear();
            m_GhostJoins.Clear();
        }
    }
    //------------------------------------------------------------------------------

    private void DisposeAllPolyPts()
    {
        for (int i = 0; i < m_PolyOuts.Count; ++i) DisposeOutRec(i);
        m_PolyOuts.Clear();
    }
    //------------------------------------------------------------------------------

    private void AddJoin(OutPt Op1, OutPt Op2, IntPoint OffPt)
    {
        Join j = new()
        {
            OutPt1 = Op1,
            OutPt2 = Op2,
            OffPt = OffPt
        };
        m_Joins.Add(j);
    }
    //------------------------------------------------------------------------------

    private void AddGhostJoin(OutPt Op, IntPoint OffPt)
    {
        Join j = new()
        {
            OutPt1 = Op,
            OffPt = OffPt
        };
        m_GhostJoins.Add(j);
    }
    //------------------------------------------------------------------------------

#if use_xyz
    internal void SetZ(ref IntPoint pt, TEdge e1, TEdge e2)
    {
        switch (ZFillFunction)
        {
            /*pt.Z != 0 ||*/
            case null:
                return;
            default:
                /*else if (pt == e1.Bot) pt.Z = e1.Bot.Z;
        else if (pt == e1.Top) pt.Z = e1.Top.Z;
        else if (pt == e2.Bot) pt.Z = e2.Bot.Z;
        else if (pt == e2.Top) pt.Z = e2.Top.Z;
        else*/
                ZFillFunction(e1.Bot, e1.Top, e2.Bot, e2.Top, ref pt);
                break;
        }
    }
    //------------------------------------------------------------------------------
#endif

    private void InsertLocalMinimaIntoAEL(long botY)
    {
        while (PopLocalMinima(botY, out LocalMinima lm))
        {
            TEdge lb = lm.LeftBound;
            TEdge rb = lm.RightBound;

            OutPt Op1 = null;
            switch (lb)
            {
                case null:
                {
                    InsertEdgeIntoAEL(rb, null);
                    SetWindingCount(rb);
                    if (IsContributing(rb))
                    {
                        Op1 = AddOutPt(rb, rb.Bot);
                    }

                    break;
                }
                default:
                {
                    switch (rb)
                    {
                        case null:
                        {
                            InsertEdgeIntoAEL(lb, null);
                            SetWindingCount(lb);
                            if (IsContributing(lb))
                            {
                                Op1 = AddOutPt(lb, lb.Bot);
                            }

                            InsertScanbeam(lb.Top.Y);
                            break;
                        }
                        default:
                        {
                            InsertEdgeIntoAEL(lb, null);
                            InsertEdgeIntoAEL(rb, lb);
                            SetWindingCount(lb);
                            rb.WindCnt = lb.WindCnt;
                            rb.WindCnt2 = lb.WindCnt2;
                            if (IsContributing(lb))
                            {
                                Op1 = AddLocalMinPoly(lb, rb, lb.Bot);
                            }

                            InsertScanbeam(lb.Top.Y);
                            break;
                        }
                    }

                    break;
                }
            }

            if (rb != null)
            {
                if (IsHorizontal(rb))
                {
                    if (rb.NextInLML != null)
                    {
                        InsertScanbeam(rb.NextInLML.Top.Y);
                    }

                    AddEdgeToSEL(rb);
                }
                else
                {
                    InsertScanbeam(rb.Top.Y);
                }
            }

            if (lb == null || rb == null)
            {
                continue;
            }

            //if output polygons share an Edge with a horizontal rb, they'll need joining later ...
            if (Op1 != null && IsHorizontal(rb) &&
                m_GhostJoins.Count > 0 && rb.WindDelta != 0)
            {
                foreach (Join j in m_GhostJoins.Where(j =>
                             HorzSegmentsOverlap(j.OutPt1.Pt.X, j.OffPt.X, rb.Bot.X, rb.Top.X)))
                {
                    AddJoin(j.OutPt1, Op1, j.OffPt);
                }
            }

            switch (lb.OutIdx)
            {
                case >= 0 when lb.PrevInAEL != null && lb.PrevInAEL.Curr.X == lb.Bot.X && lb.PrevInAEL.OutIdx >= 0 &&
                               SlopesEqual(lb.PrevInAEL.Curr, lb.PrevInAEL.Top, lb.Curr, lb.Top, m_UseFullRange) &&
                               lb.WindDelta != 0 && lb.PrevInAEL.WindDelta != 0:
                {
                    OutPt Op2 = AddOutPt(lb.PrevInAEL, lb.Bot);
                    AddJoin(Op1, Op2, lb.Top);
                    break;
                }
            }

            if (lb.NextInAEL == rb)
            {
                continue;
            }

            {
                switch (rb.OutIdx)
                {
                    case >= 0 when rb.PrevInAEL.OutIdx >= 0 &&
                                   SlopesEqual(rb.PrevInAEL.Curr, rb.PrevInAEL.Top, rb.Curr, rb.Top, m_UseFullRange) &&
                                   rb.WindDelta != 0 && rb.PrevInAEL.WindDelta != 0:
                    {
                        OutPt Op2 = AddOutPt(rb.PrevInAEL, rb.Bot);
                        AddJoin(Op1, Op2, rb.Top);
                        break;
                    }
                }

                TEdge e = lb.NextInAEL;
                if (e == null)
                {
                    continue;
                }

                while (e != rb)
                {
                    //nb: For calculating winding counts etc, IntersectEdges() assumes
                    //that param1 will be to the right of param2 ABOVE the intersection ...
                    IntersectEdges(rb, e, lb.Curr); //order important here
                    e = e.NextInAEL;
                }
            }
        }
    }
    //------------------------------------------------------------------------------

    private void InsertEdgeIntoAEL(TEdge edge, TEdge startEdge)
    {
        switch (m_ActiveEdges)
        {
            case null:
                edge.PrevInAEL = null;
                edge.NextInAEL = null;
                m_ActiveEdges = edge;
                break;
            default:
            {
                switch (startEdge)
                {
                    case null when E2InsertsBeforeE1(m_ActiveEdges, edge):
                        edge.PrevInAEL = null;
                        edge.NextInAEL = m_ActiveEdges;
                        m_ActiveEdges.PrevInAEL = edge;
                        m_ActiveEdges = edge;
                        break;
                    default:
                    {
                        startEdge = startEdge switch
                        {
                            null => m_ActiveEdges,
                            _ => startEdge
                        };

                        while (startEdge.NextInAEL != null &&
                               !E2InsertsBeforeE1(startEdge.NextInAEL, edge))
                        {
                            startEdge = startEdge.NextInAEL;
                        }

                        edge.NextInAEL = startEdge.NextInAEL;
                        if (startEdge.NextInAEL != null)
                        {
                            startEdge.NextInAEL.PrevInAEL = edge;
                        }

                        edge.PrevInAEL = startEdge;
                        startEdge.NextInAEL = edge;
                        break;
                    }
                }

                break;
            }
        }
    }
    //----------------------------------------------------------------------

    private bool E2InsertsBeforeE1(TEdge e1, TEdge e2)
    {
        if (e2.Curr.X != e1.Curr.X)
        {
            return e2.Curr.X < e1.Curr.X;
        }

        if (e2.Top.Y > e1.Top.Y)
        {
            return e2.Top.X < TopX(e1, e2.Top.Y);
        }

        return e1.Top.X > TopX(e2, e1.Top.Y);
    }
    //------------------------------------------------------------------------------

    private bool IsEvenOddFillType(TEdge edge)
    {
        return edge.PolyTyp switch
        {
            PolyType.ptSubject => m_SubjFillType == PolyFillType.pftEvenOdd,
            _ => m_ClipFillType == PolyFillType.pftEvenOdd
        };
    }
    //------------------------------------------------------------------------------

    private bool IsEvenOddAltFillType(TEdge edge)
    {
        return edge.PolyTyp switch
        {
            PolyType.ptSubject => m_ClipFillType == PolyFillType.pftEvenOdd,
            _ => m_SubjFillType == PolyFillType.pftEvenOdd
        };
    }
    //------------------------------------------------------------------------------

    private bool IsContributing(TEdge edge)
    {
        PolyFillType pft, pft2;
        switch (edge.PolyTyp)
        {
            case PolyType.ptSubject:
                pft = m_SubjFillType;
                pft2 = m_ClipFillType;
                break;
            default:
                pft = m_ClipFillType;
                pft2 = m_SubjFillType;
                break;
        }

        switch (pft)
        {
            case PolyFillType.pftEvenOdd:
                switch (edge.WindDelta)
                {
                    //return false if a subj line has been flagged as inside a subj polygon
                    case 0 when edge.WindCnt != 1:
                        return false;
                }

                break;
            case PolyFillType.pftNonZero:
                if (Math.Abs(edge.WindCnt) != 1)
                {
                    return false;
                }

                break;
            case PolyFillType.pftPositive:
                if (edge.WindCnt != 1)
                {
                    return false;
                }

                break;
            default: //PolyFillType.pftNegative
                if (edge.WindCnt != -1)
                {
                    return false;
                }

                break;
        }

        switch (m_ClipType)
        {
            case ClipType.ctIntersection:
                switch (pft2)
                {
                    case PolyFillType.pftEvenOdd:
                    case PolyFillType.pftNonZero:
                        return edge.WindCnt2 != 0;
                    case PolyFillType.pftPositive:
                        return edge.WindCnt2 > 0;
                    default:
                        return edge.WindCnt2 < 0;
                }
            case ClipType.ctUnion:
                switch (pft2)
                {
                    case PolyFillType.pftEvenOdd:
                    case PolyFillType.pftNonZero:
                        return edge.WindCnt2 == 0;
                    case PolyFillType.pftPositive:
                        return edge.WindCnt2 <= 0;
                    default:
                        return edge.WindCnt2 >= 0;
                }
            case ClipType.ctDifference:
                switch (edge.PolyTyp)
                {
                    case PolyType.ptSubject:
                        switch (pft2)
                        {
                            case PolyFillType.pftEvenOdd:
                            case PolyFillType.pftNonZero:
                                return edge.WindCnt2 == 0;
                            case PolyFillType.pftPositive:
                                return edge.WindCnt2 <= 0;
                            default:
                                return edge.WindCnt2 >= 0;
                        }

                        break;
                    default:
                        switch (pft2)
                        {
                            case PolyFillType.pftEvenOdd:
                            case PolyFillType.pftNonZero:
                                return edge.WindCnt2 != 0;
                            case PolyFillType.pftPositive:
                                return edge.WindCnt2 > 0;
                            default:
                                return edge.WindCnt2 < 0;
                        }

                        break;
                }
            case ClipType.ctXor:
                switch (edge.WindDelta)
                {
                    //XOr always contributing unless open
                    case 0:
                        switch (pft2)
                        {
                            case PolyFillType.pftEvenOdd:
                            case PolyFillType.pftNonZero:
                                return edge.WindCnt2 == 0;
                            case PolyFillType.pftPositive:
                                return edge.WindCnt2 <= 0;
                            default:
                                return edge.WindCnt2 >= 0;
                        }

                        break;
                    default:
                        return true;
                }
        }

        return true;
    }
    //------------------------------------------------------------------------------

    private void SetWindingCount(TEdge edge)
    {
        TEdge e = edge.PrevInAEL;
        //find the edge of the same polytype that immediately preceeds 'edge' in AEL
        while (e != null && (e.PolyTyp != edge.PolyTyp || e.WindDelta == 0))
        {
            e = e.PrevInAEL;
        }

        switch (e)
        {
            case null:
            {
                PolyFillType pft = edge.PolyTyp == PolyType.ptSubject ? m_SubjFillType : m_ClipFillType;
                edge.WindCnt = edge.WindDelta switch
                {
                    0 => pft == PolyFillType.pftNegative ? -1 : 1,
                    _ => edge.WindDelta
                };

                edge.WindCnt2 = 0;
                e = m_ActiveEdges; //ie get ready to calc WindCnt2
                break;
            }
            default:
            {
                switch (edge.WindDelta)
                {
                    case 0 when m_ClipType != ClipType.ctUnion:
                        edge.WindCnt = 1;
                        edge.WindCnt2 = e.WindCnt2;
                        e = e.NextInAEL; //ie get ready to calc WindCnt2
                        break;
                    default:
                    {
                        if (IsEvenOddFillType(edge))
                        {
                            switch (edge.WindDelta)
                            {
                                //EvenOdd filling ...
                                case 0:
                                {
                                    //are we inside a subj polygon ...
                                    bool Inside = true;
                                    TEdge e2 = e.PrevInAEL;
                                    while (e2 != null)
                                    {
                                        if (e2.PolyTyp == e.PolyTyp && e2.WindDelta != 0)
                                        {
                                            Inside = !Inside;
                                        }

                                        e2 = e2.PrevInAEL;
                                    }

                                    edge.WindCnt = Inside ? 0 : 1;
                                    break;
                                }
                                default:
                                    edge.WindCnt = edge.WindDelta;
                                    break;
                            }

                            edge.WindCnt2 = e.WindCnt2;
                            e = e.NextInAEL; //ie get ready to calc WindCnt2
                        }
                        else
                        {
                            switch (e.WindCnt * e.WindDelta)
                            {
                                //nonZero, Positive or Negative filling ...
                                //prev edge is 'decreasing' WindCount (WC) toward zero
                                //so we're outside the previous polygon ...
                                case < 0 when Math.Abs(e.WindCnt) > 1:
                                {
                                    edge.WindCnt = (e.WindDelta * edge.WindDelta) switch
                                    {
                                        //outside prev poly but still inside another.
                                        //when reversing direction of prev poly use the same WC 
                                        < 0 => e.WindCnt,
                                        _ => e.WindCnt + edge.WindDelta
                                    };

                                    break;
                                }
                                //now outside all polys of same polytype so set own WC ...
                                case < 0:
                                    edge.WindCnt = edge.WindDelta == 0 ? 1 : edge.WindDelta;
                                    break;
                                default:
                                {
                                    switch (edge.WindDelta)
                                    {
                                        //prev edge is 'increasing' WindCount (WC) away from zero
                                        //so we're inside the previous polygon ...
                                        case 0:
                                            edge.WindCnt = e.WindCnt < 0 ? e.WindCnt - 1 : e.WindCnt + 1;
                                            break;
                                        //if wind direction is reversing prev then use same WC
                                        default:
                                        {
                                            edge.WindCnt = (e.WindDelta * edge.WindDelta) switch
                                            {
                                                < 0 => e.WindCnt,
                                                _ => e.WindCnt + edge.WindDelta
                                            };

                                            break;
                                        }
                                    }

                                    break;
                                }
                            }

                            edge.WindCnt2 = e.WindCnt2;
                            e = e.NextInAEL; //ie get ready to calc WindCnt2
                        }

                        break;
                    }
                }

                break;
            }
        }

        //update WindCnt2 ...
        if (IsEvenOddAltFillType(edge))
        {
            //EvenOdd filling ...
            while (e != edge)
            {
                if (e.WindDelta != 0)
                {
                    edge.WindCnt2 = edge.WindCnt2 == 0 ? 1 : 0;
                }

                e = e.NextInAEL;
            }
        }
        else
        {
            //nonZero, Positive or Negative filling ...
            while (e != edge)
            {
                edge.WindCnt2 += e.WindDelta;
                e = e.NextInAEL;
            }
        }
    }
    //------------------------------------------------------------------------------

    private void AddEdgeToSEL(TEdge edge)
    {
        switch (m_SortedEdges)
        {
            //SEL pointers in PEdge are use to build transient lists of horizontal edges.
            //However, since we don't need to worry about processing order, all additions
            //are made to the front of the list ...
            case null:
                m_SortedEdges = edge;
                edge.PrevInSEL = null;
                edge.NextInSEL = null;
                break;
            default:
                edge.NextInSEL = m_SortedEdges;
                edge.PrevInSEL = null;
                m_SortedEdges.PrevInSEL = edge;
                m_SortedEdges = edge;
                break;
        }
    }
    //------------------------------------------------------------------------------

    internal bool PopEdgeFromSEL(out TEdge e)
    {
        //Pop edge from front of SEL (ie SEL is a FILO list)
        e = m_SortedEdges;
        switch (e)
        {
            case null:
                return false;
        }

        TEdge oldE = e;
        m_SortedEdges = e.NextInSEL;
        if (m_SortedEdges != null)
        {
            m_SortedEdges.PrevInSEL = null;
        }

        oldE.NextInSEL = null;
        oldE.PrevInSEL = null;
        return true;
    }
    //------------------------------------------------------------------------------

    private void CopyAELToSEL()
    {
        TEdge e = m_ActiveEdges;
        m_SortedEdges = e;
        while (e != null)
        {
            e.PrevInSEL = e.PrevInAEL;
            e.NextInSEL = e.NextInAEL;
            e = e.NextInAEL;
        }
    }
    //------------------------------------------------------------------------------

    private void SwapPositionsInSEL(TEdge edge1, TEdge edge2)
    {
        switch (edge1.NextInSEL)
        {
            case null when edge1.PrevInSEL == null:
                return;
        }

        switch (edge2.NextInSEL)
        {
            case null when edge2.PrevInSEL == null:
                return;
        }

        if (edge1.NextInSEL == edge2)
        {
            TEdge next = edge2.NextInSEL;
            if (next != null)
            {
                next.PrevInSEL = edge1;
            }

            TEdge prev = edge1.PrevInSEL;
            if (prev != null)
            {
                prev.NextInSEL = edge2;
            }

            edge2.PrevInSEL = prev;
            edge2.NextInSEL = edge1;
            edge1.PrevInSEL = edge2;
            edge1.NextInSEL = next;
        }
        else if (edge2.NextInSEL == edge1)
        {
            TEdge next = edge1.NextInSEL;
            if (next != null)
            {
                next.PrevInSEL = edge2;
            }

            TEdge prev = edge2.PrevInSEL;
            if (prev != null)
            {
                prev.NextInSEL = edge1;
            }

            edge1.PrevInSEL = prev;
            edge1.NextInSEL = edge2;
            edge2.PrevInSEL = edge1;
            edge2.NextInSEL = next;
        }
        else
        {
            TEdge next = edge1.NextInSEL;
            TEdge prev = edge1.PrevInSEL;
            edge1.NextInSEL = edge2.NextInSEL;
            if (edge1.NextInSEL != null)
            {
                edge1.NextInSEL.PrevInSEL = edge1;
            }

            edge1.PrevInSEL = edge2.PrevInSEL;
            if (edge1.PrevInSEL != null)
            {
                edge1.PrevInSEL.NextInSEL = edge1;
            }

            edge2.NextInSEL = next;
            if (edge2.NextInSEL != null)
            {
                edge2.NextInSEL.PrevInSEL = edge2;
            }

            edge2.PrevInSEL = prev;
            if (edge2.PrevInSEL != null)
            {
                edge2.PrevInSEL.NextInSEL = edge2;
            }
        }

        switch (edge1.PrevInSEL)
        {
            case null:
                m_SortedEdges = edge1;
                break;
            default:
            {
                m_SortedEdges = edge2.PrevInSEL switch
                {
                    null => edge2,
                    _ => m_SortedEdges
                };

                break;
            }
        }
    }
    //------------------------------------------------------------------------------


    private void AddLocalMaxPoly(TEdge e1, TEdge e2, IntPoint pt)
    {
        AddOutPt(e1, pt);
        switch (e2.WindDelta)
        {
            case 0:
                AddOutPt(e2, pt);
                break;
        }

        if (e1.OutIdx == e2.OutIdx)
        {
            e1.OutIdx = Unassigned;
            e2.OutIdx = Unassigned;
        }
        else if (e1.OutIdx < e2.OutIdx)
        {
            AppendPolygon(e1, e2);
        }
        else
        {
            AppendPolygon(e2, e1);
        }
    }
    //------------------------------------------------------------------------------

    private OutPt AddLocalMinPoly(TEdge e1, TEdge e2, IntPoint pt)
    {
        OutPt result;
        TEdge e, prevE;
        if (IsHorizontal(e2) || e1.Dx > e2.Dx)
        {
            result = AddOutPt(e1, pt);
            e2.OutIdx = e1.OutIdx;
            e1.Side = EdgeSide.esLeft;
            e2.Side = EdgeSide.esRight;
            e = e1;
            prevE = e.PrevInAEL == e2 ? e2.PrevInAEL : e.PrevInAEL;
        }
        else
        {
            result = AddOutPt(e2, pt);
            e1.OutIdx = e2.OutIdx;
            e1.Side = EdgeSide.esRight;
            e2.Side = EdgeSide.esLeft;
            e = e2;
            prevE = e.PrevInAEL == e1 ? e1.PrevInAEL : e.PrevInAEL;
        }

        if (prevE is not {OutIdx: >= 0} || prevE.Top.Y >= pt.Y || e.Top.Y >= pt.Y)
        {
            return result;
        }

        long xPrev = TopX(prevE, pt.Y);
        long xE = TopX(e, pt.Y);
        if (xPrev != xE || e.WindDelta == 0 || prevE.WindDelta == 0 || !SlopesEqual(new IntPoint(xPrev, pt.Y),
                prevE.Top, new IntPoint(xE, pt.Y), e.Top, m_UseFullRange))
        {
            return result;
        }

        OutPt outPt = AddOutPt(prevE, pt);
        AddJoin(result, outPt, e.Top);
        return result;
    }
    //------------------------------------------------------------------------------

    private OutPt AddOutPt(TEdge e, IntPoint pt)
    {
        switch (e.OutIdx)
        {
            case < 0:
            {
                OutRec outRec = CreateOutRec();
                outRec.IsOpen = e.WindDelta == 0;
                OutPt newOp = new();
                outRec.Pts = newOp;
                newOp.Idx = outRec.Idx;
                newOp.Pt = pt;
                newOp.Next = newOp;
                newOp.Prev = newOp;
                switch (outRec.IsOpen)
                {
                    case false:
                        SetHoleState(e, outRec);
                        break;
                }

                e.OutIdx = outRec.Idx; //nb: do this after SetZ !
                return newOp;
            }
            default:
            {
                OutRec outRec = m_PolyOuts[e.OutIdx];
                //OutRec.Pts is the 'Left-most' point & OutRec.Pts.Prev is the 'Right-most'
                OutPt op = outRec.Pts;
                bool ToFront = e.Side == EdgeSide.esLeft;
                switch (ToFront)
                {
                    case true when pt == op.Pt:
                        return op;
                    case false when pt == op.Prev.Pt:
                        return op.Prev;
                }

                OutPt newOp = new()
                {
                    Idx = outRec.Idx,
                    Pt = pt,
                    Next = op,
                    Prev = op.Prev
                };
                newOp.Prev.Next = newOp;
                op.Prev = newOp;
                outRec.Pts = ToFront switch
                {
                    true => newOp,
                    _ => outRec.Pts
                };

                return newOp;
            }
        }
    }
    //------------------------------------------------------------------------------

    private OutPt GetLastOutPt(TEdge e)
    {
        OutRec outRec = m_PolyOuts[e.OutIdx];
        return e.Side switch
        {
            EdgeSide.esLeft => outRec.Pts,
            _ => outRec.Pts.Prev
        };
    }
    //------------------------------------------------------------------------------

    internal void SwapPoints(ref IntPoint pt1, ref IntPoint pt2)
    {
        IntPoint tmp = new(pt1);
        pt1 = pt2;
        pt2 = tmp;
    }
    //------------------------------------------------------------------------------

    private bool HorzSegmentsOverlap(long seg1a, long seg1b, long seg2a, long seg2b)
    {
        if (seg1a > seg1b)
        {
            Swap(ref seg1a, ref seg1b);
        }

        if (seg2a > seg2b)
        {
            Swap(ref seg2a, ref seg2b);
        }

        return seg1a < seg2b && seg2a < seg1b;
    }
    //------------------------------------------------------------------------------

    private void SetHoleState(TEdge e, OutRec outRec)
    {
        TEdge e2 = e.PrevInAEL;
        TEdge eTmp = null;
        while (e2 != null)
        {
            switch (e2.OutIdx)
            {
                case >= 0 when e2.WindDelta != 0:
                {
                    switch (eTmp)
                    {
                        case null:
                            eTmp = e2;
                            break;
                        default:
                        {
                            if (eTmp.OutIdx == e2.OutIdx)
                            {
                                eTmp = null; //paired               
                            }

                            break;
                        }
                    }

                    break;
                }
            }

            e2 = e2.PrevInAEL;
        }

        switch (eTmp)
        {
            case null:
                outRec.FirstLeft = null;
                outRec.IsHole = false;
                break;
            default:
                outRec.FirstLeft = m_PolyOuts[eTmp.OutIdx];
                outRec.IsHole = !outRec.FirstLeft.IsHole;
                break;
        }
    }
    //------------------------------------------------------------------------------

    private double GetDx(IntPoint pt1, IntPoint pt2)
    {
        if (pt1.Y == pt2.Y)
        {
            return horizontal;
        }

        return (double) (pt2.X - pt1.X) / (pt2.Y - pt1.Y);
    }
    //---------------------------------------------------------------------------

    private bool FirstIsBottomPt(OutPt btmPt1, OutPt btmPt2)
    {
        OutPt p = btmPt1.Prev;
        while (p.Pt == btmPt1.Pt && p != btmPt1)
        {
            p = p.Prev;
        }

        double dx1p = Math.Abs(GetDx(btmPt1.Pt, p.Pt));
        p = btmPt1.Next;
        while (p.Pt == btmPt1.Pt && p != btmPt1)
        {
            p = p.Next;
        }

        double dx1n = Math.Abs(GetDx(btmPt1.Pt, p.Pt));

        p = btmPt2.Prev;
        while (p.Pt == btmPt2.Pt && p != btmPt2)
        {
            p = p.Prev;
        }

        double dx2p = Math.Abs(GetDx(btmPt2.Pt, p.Pt));
        p = btmPt2.Next;
        while (p.Pt == btmPt2.Pt && p != btmPt2)
        {
            p = p.Next;
        }

        double dx2n = Math.Abs(GetDx(btmPt2.Pt, p.Pt));

        return Math.Abs(Math.Max(dx1p, dx1n) - Math.Max(dx2p, dx2n)) switch
        {
            <= double.Epsilon when Math.Abs(Math.Min(dx1p, dx1n) - Math.Min(dx2p, dx2n)) <= double.Epsilon =>
                Area(btmPt1) > 0 //if otherwise identical use orientation
            ,
            _ => dx1p >= dx2p && dx1p >= dx2n || dx1n >= dx2p && dx1n >= dx2n
        };
    }
    //------------------------------------------------------------------------------

    private OutPt GetBottomPt(OutPt pp)
    {
        OutPt dups = null;
        OutPt p = pp.Next;
        while (p != pp)
        {
            if (p.Pt.Y > pp.Pt.Y)
            {
                pp = p;
                dups = null;
            }
            else if (p.Pt.Y == pp.Pt.Y && p.Pt.X <= pp.Pt.X)
            {
                if (p.Pt.X < pp.Pt.X)
                {
                    dups = null;
                    pp = p;
                }
                else
                {
                    if (p.Next != pp && p.Prev != pp)
                    {
                        dups = p;
                    }
                }
            }

            p = p.Next;
        }

        if (dups == null)
        {
            return pp;
        }

        //there appears to be at least 2 vertices at bottomPt so ...
        while (dups != p)
        {
            if (!FirstIsBottomPt(p, dups))
            {
                pp = dups;
            }

            dups = dups.Next;
            while (dups.Pt != pp.Pt)
            {
                dups = dups.Next;
            }
        }

        return pp;
    }
    //------------------------------------------------------------------------------

    private OutRec GetLowermostRec(OutRec outRec1, OutRec outRec2)
    {
        outRec1.BottomPt = outRec1.BottomPt switch
        {
            //work out which polygon fragment has the correct hole state ...
            null => GetBottomPt(outRec1.Pts),
            _ => outRec1.BottomPt
        };

        outRec2.BottomPt = outRec2.BottomPt switch
        {
            null => GetBottomPt(outRec2.Pts),
            _ => outRec2.BottomPt
        };

        OutPt bPt1 = outRec1.BottomPt;
        OutPt bPt2 = outRec2.BottomPt;
        if (bPt1.Pt.Y > bPt2.Pt.Y)
        {
            return outRec1;
        }

        if (bPt1.Pt.Y < bPt2.Pt.Y)
        {
            return outRec2;
        }

        if (bPt1.Pt.X < bPt2.Pt.X)
        {
            return outRec1;
        }

        if (bPt1.Pt.X > bPt2.Pt.X)
        {
            return outRec2;
        }

        if (bPt1.Next == bPt1)
        {
            return outRec2;
        }

        if (bPt2.Next == bPt2)
        {
            return outRec1;
        }

        return FirstIsBottomPt(bPt1, bPt2) ? outRec1 : outRec2;
    }
    //------------------------------------------------------------------------------

    private bool OutRec1RightOfOutRec2(OutRec outRec1, OutRec outRec2)
    {
        do
        {
            outRec1 = outRec1.FirstLeft;
            if (outRec1 == outRec2)
            {
                return true;
            }
        } while (outRec1 != null);

        return false;
    }
    //------------------------------------------------------------------------------

    private OutRec GetOutRec(int idx)
    {
        OutRec outrec = m_PolyOuts[idx];
        while (outrec != m_PolyOuts[outrec.Idx])
        {
            outrec = m_PolyOuts[outrec.Idx];
        }

        return outrec;
    }
    //------------------------------------------------------------------------------

    private void AppendPolygon(TEdge e1, TEdge e2)
    {
        OutRec outRec1 = m_PolyOuts[e1.OutIdx];
        OutRec outRec2 = m_PolyOuts[e2.OutIdx];

        OutRec holeStateRec;
        if (OutRec1RightOfOutRec2(outRec1, outRec2))
        {
            holeStateRec = outRec2;
        }
        else if (OutRec1RightOfOutRec2(outRec2, outRec1))
        {
            holeStateRec = outRec1;
        }
        else
        {
            holeStateRec = GetLowermostRec(outRec1, outRec2);
        }

        //get the start and ends of both output polygons and
        //join E2 poly onto E1 poly and delete pointers to E2 ...
        OutPt p1_lft = outRec1.Pts;
        OutPt p1_rt = p1_lft.Prev;
        OutPt p2_lft = outRec2.Pts;
        OutPt p2_rt = p2_lft.Prev;

        switch (e1.Side)
        {
            //join e2 poly onto e1 poly and delete pointers to e2 ...
            case EdgeSide.esLeft when e2.Side == EdgeSide.esLeft:
                //z y x a b c
                ReversePolyPtLinks(p2_lft);
                p2_lft.Next = p1_lft;
                p1_lft.Prev = p2_lft;
                p1_rt.Next = p2_rt;
                p2_rt.Prev = p1_rt;
                outRec1.Pts = p2_rt;
                break;
            case EdgeSide.esLeft:
                //x y z a b c
                p2_rt.Next = p1_lft;
                p1_lft.Prev = p2_rt;
                p2_lft.Prev = p1_rt;
                p1_rt.Next = p2_lft;
                outRec1.Pts = p2_lft;
                break;
            default:
            {
                switch (e2.Side)
                {
                    case EdgeSide.esRight:
                        //a b c z y x
                        ReversePolyPtLinks(p2_lft);
                        p1_rt.Next = p2_rt;
                        p2_rt.Prev = p1_rt;
                        p2_lft.Next = p1_lft;
                        p1_lft.Prev = p2_lft;
                        break;
                    default:
                        //a b c x y z
                        p1_rt.Next = p2_lft;
                        p2_lft.Prev = p1_rt;
                        p1_lft.Prev = p2_rt;
                        p2_rt.Next = p1_lft;
                        break;
                }

                break;
            }
        }

        outRec1.BottomPt = null;
        if (holeStateRec == outRec2)
        {
            if (outRec2.FirstLeft != outRec1)
            {
                outRec1.FirstLeft = outRec2.FirstLeft;
            }

            outRec1.IsHole = outRec2.IsHole;
        }

        outRec2.Pts = null;
        outRec2.BottomPt = null;

        outRec2.FirstLeft = outRec1;

        int OKIdx = e1.OutIdx;
        int ObsoleteIdx = e2.OutIdx;

        e1.OutIdx = Unassigned; //nb: safe because we only get here via AddLocalMaxPoly
        e2.OutIdx = Unassigned;

        TEdge e = m_ActiveEdges;
        while (e != null)
        {
            if (e.OutIdx == ObsoleteIdx)
            {
                e.OutIdx = OKIdx;
                e.Side = e1.Side;
                break;
            }

            e = e.NextInAEL;
        }

        outRec2.Idx = outRec1.Idx;
    }
    //------------------------------------------------------------------------------

    private void ReversePolyPtLinks(OutPt pp)
    {
        switch (pp)
        {
            case null:
                return;
        }

        OutPt pp1 = pp;
        do
        {
            OutPt pp2 = pp1.Next;
            pp1.Next = pp1.Prev;
            pp1.Prev = pp2;
            pp1 = pp2;
        } while (pp1 != pp);
    }
    //------------------------------------------------------------------------------

    private static void SwapSides(TEdge edge1, TEdge edge2)
    {
        (edge1.Side, edge2.Side) = (edge2.Side, edge1.Side);
    }
    //------------------------------------------------------------------------------

    private static void SwapPolyIndexes(TEdge edge1, TEdge edge2)
    {
        (edge1.OutIdx, edge2.OutIdx) = (edge2.OutIdx, edge1.OutIdx);
    }
    //------------------------------------------------------------------------------

    private void IntersectEdges(TEdge e1, TEdge e2, IntPoint pt)
    {
        //e1 will be to the left of e2 BELOW the intersection. Therefore e1 is before
        //e2 in AEL except when e1 is being inserted at the intersection point ...

        bool e1Contributing = e1.OutIdx >= 0;
        bool e2Contributing = e2.OutIdx >= 0;

#if use_xyz
        SetZ(ref pt, e1, e2);
#endif

#if use_lines
        //if either edge is on an OPEN path ...
        if (e1.WindDelta == 0 || e2.WindDelta == 0)
        {
            switch (e1.WindDelta)
            {
                //ignore subject-subject open path intersections UNLESS they
                //are both open paths, AND they are both 'contributing maximas' ...
                case 0 when e2.WindDelta == 0:
                    return;
            }
            //if intersecting a subj line with a subj poly ...

            if (e1.PolyTyp == e2.PolyTyp &&
                e1.WindDelta != e2.WindDelta && m_ClipType == ClipType.ctUnion)
            {
                switch (e1.WindDelta)
                {
                    case 0:
                    {
                        switch (e2Contributing)
                        {
                            case true:
                            {
                                AddOutPt(e1, pt);
                                e1.OutIdx = e1Contributing switch
                                {
                                    true => Unassigned,
                                    _ => e1.OutIdx
                                };

                                break;
                            }
                        }

                        break;
                    }
                    default:
                    {
                        switch (e1Contributing)
                        {
                            case true:
                            {
                                AddOutPt(e2, pt);
                                e2.OutIdx = e2Contributing switch
                                {
                                    true => Unassigned,
                                    _ => e2.OutIdx
                                };

                                break;
                            }
                        }

                        break;
                    }
                }
            }
            else if (e1.PolyTyp != e2.PolyTyp)
            {
                switch (e1.WindDelta)
                {
                    case 0 when Math.Abs(e2.WindCnt) == 1 && (m_ClipType != ClipType.ctUnion || e2.WindCnt2 == 0):
                    {
                        AddOutPt(e1, pt);
                        e1.OutIdx = e1Contributing switch
                        {
                            true => Unassigned,
                            _ => e1.OutIdx
                        };

                        break;
                    }
                    default:
                    {
                        switch (e2.WindDelta)
                        {
                            case 0 when Math.Abs(e1.WindCnt) == 1 &&
                                        (m_ClipType != ClipType.ctUnion || e1.WindCnt2 == 0):
                            {
                                AddOutPt(e2, pt);
                                e2.OutIdx = e2Contributing switch
                                {
                                    true => Unassigned,
                                    _ => e2.OutIdx
                                };

                                break;
                            }
                        }

                        break;
                    }
                }
            }

            return;
        }
#endif

        //update winding counts...
        //assumes that e1 will be to the Right of e2 ABOVE the intersection
        if (e1.PolyTyp == e2.PolyTyp)
        {
            if (IsEvenOddFillType(e1))
            {
                (e1.WindCnt, e2.WindCnt) = (e2.WindCnt, e1.WindCnt);
            }
            else
            {
                switch (e1.WindCnt + e2.WindDelta)
                {
                    case 0:
                        e1.WindCnt = -e1.WindCnt;
                        break;
                    default:
                        e1.WindCnt += e2.WindDelta;
                        break;
                }

                switch (e2.WindCnt - e1.WindDelta)
                {
                    case 0:
                        e2.WindCnt = -e2.WindCnt;
                        break;
                    default:
                        e2.WindCnt -= e1.WindDelta;
                        break;
                }
            }
        }
        else
        {
            if (!IsEvenOddFillType(e2))
            {
                e1.WindCnt2 += e2.WindDelta;
            }
            else
            {
                e1.WindCnt2 = e1.WindCnt2 == 0 ? 1 : 0;
            }

            if (!IsEvenOddFillType(e1))
            {
                e2.WindCnt2 -= e1.WindDelta;
            }
            else
            {
                e2.WindCnt2 = e2.WindCnt2 == 0 ? 1 : 0;
            }
        }

        PolyFillType e1FillType, e2FillType, e1FillType2, e2FillType2;
        switch (e1.PolyTyp)
        {
            case PolyType.ptSubject:
                e1FillType = m_SubjFillType;
                e1FillType2 = m_ClipFillType;
                break;
            default:
                e1FillType = m_ClipFillType;
                e1FillType2 = m_SubjFillType;
                break;
        }

        switch (e2.PolyTyp)
        {
            case PolyType.ptSubject:
                e2FillType = m_SubjFillType;
                e2FillType2 = m_ClipFillType;
                break;
            default:
                e2FillType = m_ClipFillType;
                e2FillType2 = m_SubjFillType;
                break;
        }

        int e1Wc = e1FillType switch
        {
            PolyFillType.pftPositive => e1.WindCnt,
            PolyFillType.pftNegative => -e1.WindCnt,
            _ => Math.Abs(e1.WindCnt)
        };
        int e2Wc = e2FillType switch
        {
            PolyFillType.pftPositive => e2.WindCnt,
            PolyFillType.pftNegative => -e2.WindCnt,
            _ => Math.Abs(e2.WindCnt)
        };

        switch (e1Contributing)
        {
            case true when e2Contributing:
            {
                if (e1Wc != 0 && e1Wc != 1 || e2Wc != 0 && e2Wc != 1 ||
                    e1.PolyTyp != e2.PolyTyp && m_ClipType != ClipType.ctXor)
                {
                    AddLocalMaxPoly(e1, e2, pt);
                }
                else
                {
                    AddOutPt(e1, pt);
                    AddOutPt(e2, pt);
                    SwapSides(e1, e2);
                    SwapPolyIndexes(e1, e2);
                }

                break;
            }
            case true:
            {
                switch (e2Wc)
                {
                    case 0:
                    case 1:
                        AddOutPt(e1, pt);
                        SwapSides(e1, e2);
                        SwapPolyIndexes(e1, e2);
                        break;
                }

                break;
            }
            default:
            {
                switch (e2Contributing)
                {
                    case true:
                    {
                        switch (e1Wc)
                        {
                            case 0:
                            case 1:
                                AddOutPt(e2, pt);
                                SwapSides(e1, e2);
                                SwapPolyIndexes(e1, e2);
                                break;
                        }

                        break;
                    }
                    default:
                    {
                        switch (e1Wc)
                        {
                            case 0 or 1 when e2Wc is 0 or 1:
                            {
                                //neither edge is currently contributing ...
                                long e1Wc2 = e1FillType2 switch
                                {
                                    PolyFillType.pftPositive => e1.WindCnt2,
                                    PolyFillType.pftNegative => -e1.WindCnt2,
                                    _ => Math.Abs(e1.WindCnt2)
                                };
                                long e2Wc2 = e2FillType2 switch
                                {
                                    PolyFillType.pftPositive => e2.WindCnt2,
                                    PolyFillType.pftNegative => -e2.WindCnt2,
                                    _ => Math.Abs(e2.WindCnt2)
                                };

                                if (e1.PolyTyp != e2.PolyTyp)
                                {
                                    AddLocalMinPoly(e1, e2, pt);
                                }
                                else
                                {
                                    switch (e1Wc)
                                    {
                                        case 1 when e2Wc == 1:
                                            switch (m_ClipType)
                                            {
                                                case ClipType.ctIntersection:
                                                    switch (e1Wc2)
                                                    {
                                                        case > 0 when e2Wc2 > 0:
                                                            AddLocalMinPoly(e1, e2, pt);
                                                            break;
                                                    }

                                                    break;
                                                case ClipType.ctUnion:
                                                    switch (e1Wc2)
                                                    {
                                                        case <= 0 when e2Wc2 <= 0:
                                                            AddLocalMinPoly(e1, e2, pt);
                                                            break;
                                                    }

                                                    break;
                                                case ClipType.ctDifference:
                                                    switch (e1.PolyTyp)
                                                    {
                                                        case PolyType.ptClip when e1Wc2 > 0 && e2Wc2 > 0:
                                                        case PolyType.ptSubject when e1Wc2 <= 0 && e2Wc2 <= 0:
                                                            AddLocalMinPoly(e1, e2, pt);
                                                            break;
                                                    }

                                                    break;
                                                case ClipType.ctXor:
                                                    AddLocalMinPoly(e1, e2, pt);
                                                    break;
                                            }

                                            break;
                                        default:
                                            SwapSides(e1, e2);
                                            break;
                                    }
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    //------------------------------------------------------------------------------

    private void ProcessHorizontals()
    {
        while (PopEdgeFromSEL(out TEdge horzEdge))
        {
            ProcessHorizontal(horzEdge);
        }
    }
    //------------------------------------------------------------------------------

    private void GetHorzDirection(TEdge HorzEdge, out Direction Dir, out long Left, out long Right)
    {
        if (HorzEdge.Bot.X < HorzEdge.Top.X)
        {
            Left = HorzEdge.Bot.X;
            Right = HorzEdge.Top.X;
            Dir = Direction.dLeftToRight;
        }
        else
        {
            Left = HorzEdge.Top.X;
            Right = HorzEdge.Bot.X;
            Dir = Direction.dRightToLeft;
        }
    }
    //------------------------------------------------------------------------

    private void ProcessHorizontal(TEdge horzEdge)
    {
        bool IsOpen = horzEdge.WindDelta == 0;

        GetHorzDirection(horzEdge, out Direction dir, out long horzLeft, out long horzRight);

        TEdge eLastHorz = horzEdge, eMaxPair = null;
        while (eLastHorz.NextInLML != null && IsHorizontal(eLastHorz.NextInLML))
        {
            eLastHorz = eLastHorz.NextInLML;
        }

        eMaxPair = eLastHorz.NextInLML switch
        {
            null => GetMaximaPair(eLastHorz),
            _ => eMaxPair
        };

        Maxima currMax = m_Maxima;
        if (currMax != null)
        {
            switch (dir)
            {
                //get the first maxima in range (X) ...
                case Direction.dLeftToRight:
                {
                    while (currMax != null && currMax.X <= horzEdge.Bot.X)
                    {
                        currMax = currMax.Next;
                    }

                    if (currMax != null && currMax.X >= eLastHorz.Top.X)
                    {
                        currMax = null;
                    }

                    break;
                }
                default:
                {
                    while (currMax.Next != null && currMax.Next.X < horzEdge.Bot.X)
                    {
                        currMax = currMax.Next;
                    }

                    if (currMax.X <= eLastHorz.Top.X)
                    {
                        currMax = null;
                    }

                    break;
                }
            }
        }

        OutPt op1 = null;
        for (;;) //loop through consec. horizontal edges
        {
            bool IsLastHorz = horzEdge == eLastHorz;
            TEdge e = GetNextInAEL(horzEdge, dir);
            while (e != null)
            {
                //this code block inserts extra coords into horizontal edges (in output
                //polygons) whereever maxima touch these horizontal edges. This helps
                //'simplifying' polygons (ie if the Simplify property is set).
                if (currMax != null)
                {
                    switch (dir)
                    {
                        case Direction.dLeftToRight:
                        {
                            while (currMax != null && currMax.X < e.Curr.X)
                            {
                                switch (horzEdge.OutIdx)
                                {
                                    case >= 0 when !IsOpen:
                                        AddOutPt(horzEdge, new IntPoint(currMax.X, horzEdge.Bot.Y));
                                        break;
                                }

                                currMax = currMax.Next;
                            }

                            break;
                        }
                        default:
                        {
                            while (currMax != null && currMax.X > e.Curr.X)
                            {
                                switch (horzEdge.OutIdx)
                                {
                                    case >= 0 when !IsOpen:
                                        AddOutPt(horzEdge, new IntPoint(currMax.X, horzEdge.Bot.Y));
                                        break;
                                }

                                currMax = currMax.Prev;
                            }

                            break;
                        }
                    }
                }

                if (dir == Direction.dLeftToRight && e.Curr.X > horzRight ||
                    dir == Direction.dRightToLeft && e.Curr.X < horzLeft)
                {
                    break;
                }

                //Also break if we've got to the end of an intermediate horizontal edge ...
                //nb: Smaller Dx's are to the right of larger Dx's ABOVE the horizontal.
                if (e.Curr.X == horzEdge.Top.X && horzEdge.NextInLML != null &&
                    e.Dx < horzEdge.NextInLML.Dx)
                {
                    break;
                }

                switch (horzEdge.OutIdx)
                {
                    //note: may be done multiple times
                    case >= 0 when !IsOpen:
                    {
#if use_xyz
                        switch (dir)
                        {
                            case Direction.dLeftToRight:
                                SetZ(ref e.Curr, horzEdge, e);
                                break;
                            default:
                                SetZ(ref e.Curr, e, horzEdge);
                                break;
                        }
#endif

                        op1 = AddOutPt(horzEdge, e.Curr);
                        TEdge eNextHorz = m_SortedEdges;
                        while (eNextHorz != null)
                        {
                            switch (eNextHorz.OutIdx)
                            {
                                case >= 0 when HorzSegmentsOverlap(horzEdge.Bot.X,
                                    horzEdge.Top.X, eNextHorz.Bot.X, eNextHorz.Top.X):
                                {
                                    OutPt op2 = GetLastOutPt(eNextHorz);
                                    AddJoin(op2, op1, eNextHorz.Top);
                                    break;
                                }
                            }

                            eNextHorz = eNextHorz.NextInSEL;
                        }

                        AddGhostJoin(op1, horzEdge.Bot);
                        break;
                    }
                }

                //OK, so far we're still in range of the horizontal Edge  but make sure
                //we're at the last of consec. horizontals when matching with eMaxPair
                if (e == eMaxPair && IsLastHorz)
                {
                    switch (horzEdge.OutIdx)
                    {
                        case >= 0:
                            AddLocalMaxPoly(horzEdge, eMaxPair, horzEdge.Top);
                            break;
                    }

                    DeleteFromAEL(horzEdge);
                    DeleteFromAEL(eMaxPair);
                    return;
                }

                switch (dir)
                {
                    case Direction.dLeftToRight:
                    {
                        IntPoint Pt = new(e.Curr.X, horzEdge.Curr.Y);
                        IntersectEdges(horzEdge, e, Pt);
                        break;
                    }
                    default:
                    {
                        IntPoint Pt = new(e.Curr.X, horzEdge.Curr.Y);
                        IntersectEdges(e, horzEdge, Pt);
                        break;
                    }
                }

                TEdge eNext = GetNextInAEL(e, dir);
                SwapPositionsInAEL(horzEdge, e);
                e = eNext;
            } //end while(e != null)

            //Break out of loop if HorzEdge.NextInLML is not also horizontal ...
            if (horzEdge.NextInLML == null || !IsHorizontal(horzEdge.NextInLML))
            {
                break;
            }

            UpdateEdgeIntoAEL(ref horzEdge);
            switch (horzEdge.OutIdx)
            {
                case >= 0:
                    AddOutPt(horzEdge, horzEdge.Bot);
                    break;
            }

            GetHorzDirection(horzEdge, out dir, out horzLeft, out horzRight);
        } //end for (;;)

        switch (horzEdge.OutIdx)
        {
            case >= 0 when op1 == null:
            {
                op1 = GetLastOutPt(horzEdge);
                TEdge eNextHorz = m_SortedEdges;
                while (eNextHorz != null)
                {
                    switch (eNextHorz.OutIdx)
                    {
                        case >= 0 when HorzSegmentsOverlap(horzEdge.Bot.X,
                            horzEdge.Top.X, eNextHorz.Bot.X, eNextHorz.Top.X):
                        {
                            OutPt op2 = GetLastOutPt(eNextHorz);
                            AddJoin(op2, op1, eNextHorz.Top);
                            break;
                        }
                    }

                    eNextHorz = eNextHorz.NextInSEL;
                }

                AddGhostJoin(op1, horzEdge.Top);
                break;
            }
        }

        if (horzEdge.NextInLML != null)
        {
            switch (horzEdge.OutIdx)
            {
                case >= 0:
                {
                    op1 = AddOutPt(horzEdge, horzEdge.Top);

                    UpdateEdgeIntoAEL(ref horzEdge);
                    switch (horzEdge.WindDelta)
                    {
                        case 0:
                            return;
                    }

                    //nb: HorzEdge is no longer horizontal here
                    TEdge ePrev = horzEdge.PrevInAEL;
                    TEdge eNext = horzEdge.NextInAEL;
                    if (ePrev != null && ePrev.Curr.X == horzEdge.Bot.X &&
                        ePrev.Curr.Y == horzEdge.Bot.Y && ePrev.WindDelta != 0 && ePrev.OutIdx >= 0 &&
                        ePrev.Curr.Y > ePrev.Top.Y && SlopesEqual(horzEdge, ePrev, m_UseFullRange))
                    {
                        OutPt op2 = AddOutPt(ePrev, horzEdge.Bot);
                        AddJoin(op1, op2, horzEdge.Top);
                    }
                    else if (eNext != null && eNext.Curr.X == horzEdge.Bot.X &&
                             eNext.Curr.Y == horzEdge.Bot.Y && eNext.WindDelta != 0 &&
                             eNext.OutIdx >= 0 && eNext.Curr.Y > eNext.Top.Y &&
                             SlopesEqual(horzEdge, eNext, m_UseFullRange))
                    {
                        OutPt op2 = AddOutPt(eNext, horzEdge.Bot);
                        AddJoin(op1, op2, horzEdge.Top);
                    }

                    break;
                }
                default:
                    UpdateEdgeIntoAEL(ref horzEdge);
                    break;
            }
        }
        else
        {
            switch (horzEdge.OutIdx)
            {
                case >= 0:
                    AddOutPt(horzEdge, horzEdge.Top);
                    break;
            }

            DeleteFromAEL(horzEdge);
        }
    }
    //------------------------------------------------------------------------------

    private TEdge GetNextInAEL(TEdge e, Direction Direction)
    {
        return Direction == Direction.dLeftToRight ? e.NextInAEL : e.PrevInAEL;
    }

    //------------------------------------------------------------------------------

    private bool IsMaxima(TEdge e, double Y)
    {
        return e != null && Math.Abs(e.Top.Y - Y) <= double.Epsilon && e.NextInLML == null;
    }
    //------------------------------------------------------------------------------

    private bool IsIntermediate(TEdge e, double Y)
    {
        return Math.Abs(e.Top.Y - Y) <= double.Epsilon && e.NextInLML != null;
    }
    //------------------------------------------------------------------------------

    internal TEdge GetMaximaPair(TEdge e)
    {
        if (e.Next.Top == e.Top && e.Next.NextInLML == null)
        {
            return e.Next;
        }

        if (e.Prev.Top == e.Top && e.Prev.NextInLML == null)
        {
            return e.Prev;
        }

        return null;
    }
    //------------------------------------------------------------------------------

    internal TEdge GetMaximaPairEx(TEdge e)
    {
        //as above but returns null if MaxPair isn't in AEL (unless it's horizontal)
        TEdge result = GetMaximaPair(e);
        if (result == null || result.OutIdx == Skip ||
            result.NextInAEL == result.PrevInAEL && !IsHorizontal(result))
        {
            return null;
        }

        return result;
    }
    //------------------------------------------------------------------------------

    private bool ProcessIntersections(long topY)
    {
        switch (m_ActiveEdges)
        {
            case null:
                return true;
        }

        try
        {
            BuildIntersectList(topY);
            switch (m_IntersectList.Count)
            {
                case 0:
                    return true;
            }

            if (m_IntersectList.Count == 1 || FixupIntersectionOrder())
            {
                ProcessIntersectList();
            }
            else
            {
                return false;
            }
        }
        catch
        {
            m_SortedEdges = null;
            m_IntersectList.Clear();
            throw new ClipperException("ProcessIntersections error");
        }

        m_SortedEdges = null;
        return true;
    }
    //------------------------------------------------------------------------------

    private void BuildIntersectList(long topY)
    {
        switch (m_ActiveEdges)
        {
            case null:
                return;
        }

        //prepare for sorting ...
        TEdge e = m_ActiveEdges;
        m_SortedEdges = e;
        while (e != null)
        {
            e.PrevInSEL = e.PrevInAEL;
            e.NextInSEL = e.NextInAEL;
            e.Curr.X = TopX(e, topY);
            e = e.NextInAEL;
        }

        //bubblesort ...
        bool isModified = true;
        while (isModified && m_SortedEdges != null)
        {
            isModified = false;
            e = m_SortedEdges;
            while (e.NextInSEL != null)
            {
                TEdge eNext = e.NextInSEL;
                if (e.Curr.X > eNext.Curr.X)
                {
                    IntersectPoint(e, eNext, out IntPoint pt);
                    if (pt.Y < topY)
                    {
                        pt = new IntPoint(TopX(e, topY), topY);
                    }

                    IntersectNode newNode = new()
                    {
                        Edge1 = e,
                        Edge2 = eNext,
                        Pt = pt
                    };
                    m_IntersectList.Add(newNode);

                    SwapPositionsInSEL(e, eNext);
                    isModified = true;
                }
                else
                {
                    e = eNext;
                }
            }

            if (e.PrevInSEL != null)
            {
                e.PrevInSEL.NextInSEL = null;
            }
            else
            {
                break;
            }
        }

        m_SortedEdges = null;
    }
    //------------------------------------------------------------------------------

    private bool EdgesAdjacent(IntersectNode inode)
    {
        return inode.Edge1.NextInSEL == inode.Edge2 ||
               inode.Edge1.PrevInSEL == inode.Edge2;
    }

    //------------------------------------------------------------------------------

    private bool FixupIntersectionOrder()
    {
        //pre-condition: intersections are sorted bottom-most first.
        //Now it's crucial that intersections are made only between adjacent edges,
        //so to ensure this the order of intersections may need adjusting ...
        m_IntersectList.Sort(m_IntersectNodeComparer);

        CopyAELToSEL();
        int cnt = m_IntersectList.Count;
        for (int i = 0; i < cnt; i++)
        {
            if (!EdgesAdjacent(m_IntersectList[i]))
            {
                int j = i + 1;
                while (j < cnt && !EdgesAdjacent(m_IntersectList[j]))
                {
                    j++;
                }

                if (j == cnt)
                {
                    return false;
                }

                (m_IntersectList[i], m_IntersectList[j]) = (m_IntersectList[j], m_IntersectList[i]);
            }

            SwapPositionsInSEL(m_IntersectList[i].Edge1, m_IntersectList[i].Edge2);
        }

        return true;
    }
    //------------------------------------------------------------------------------

    private void ProcessIntersectList()
    {
        foreach (IntersectNode iNode in m_IntersectList)
        {
            {
                IntersectEdges(iNode.Edge1, iNode.Edge2, iNode.Pt);
                SwapPositionsInAEL(iNode.Edge1, iNode.Edge2);
            }
        }

        m_IntersectList.Clear();
    }
    //------------------------------------------------------------------------------

    internal static long Round(double value)
    {
        return value < 0 ? (long) (value - 0.5) : (long) (value + 0.5);
    }
    //------------------------------------------------------------------------------

    private static long TopX(TEdge edge, long currentY)
    {
        if (currentY == edge.Top.Y)
        {
            return edge.Top.X;
        }

        return edge.Bot.X + Round(edge.Dx * (currentY - edge.Bot.Y));
    }
    //------------------------------------------------------------------------------

    private void IntersectPoint(TEdge edge1, TEdge edge2, out IntPoint ip)
    {
        ip = new IntPoint();
        double b2;
        switch (Math.Abs(edge1.Dx - edge2.Dx))
        {
            //nb: with very large coordinate values, it's possible for SlopesEqual() to 
            //return false but for the edge.Dx value be equal due to double precision rounding.
            case <= double.Epsilon:
                ip.Y = edge1.Curr.Y;
                ip.X = TopX(edge1, ip.Y);
                return;
        }

        switch (edge1.Delta.X)
        {
            case 0:
            {
                ip.X = edge1.Bot.X;
                if (IsHorizontal(edge2))
                {
                    ip.Y = edge2.Bot.Y;
                }
                else
                {
                    b2 = edge2.Bot.Y - edge2.Bot.X / edge2.Dx;
                    ip.Y = Round(ip.X / edge2.Dx + b2);
                }

                break;
            }
            default:
            {
                double b1;
                switch (edge2.Delta.X)
                {
                    case 0:
                    {
                        ip.X = edge2.Bot.X;
                        if (IsHorizontal(edge1))
                        {
                            ip.Y = edge1.Bot.Y;
                        }
                        else
                        {
                            b1 = edge1.Bot.Y - edge1.Bot.X / edge1.Dx;
                            ip.Y = Round(ip.X / edge1.Dx + b1);
                        }

                        break;
                    }
                    default:
                    {
                        b1 = edge1.Bot.X - edge1.Bot.Y * edge1.Dx;
                        b2 = edge2.Bot.X - edge2.Bot.Y * edge2.Dx;
                        double q = (b2 - b1) / (edge1.Dx - edge2.Dx);
                        ip.Y = Round(q);
                        ip.X = Math.Abs(edge1.Dx) < Math.Abs(edge2.Dx)
                            ? Round(edge1.Dx * q + b1)
                            : Round(edge2.Dx * q + b2);

                        break;
                    }
                }

                break;
            }
        }

        if (ip.Y < edge1.Top.Y || ip.Y < edge2.Top.Y)
        {
            ip.Y = edge1.Top.Y > edge2.Top.Y ? edge1.Top.Y : edge2.Top.Y;

            ip.X = TopX(Math.Abs(edge1.Dx) < Math.Abs(edge2.Dx) ? edge1 : edge2, ip.Y);
        }

        //finally, don't allow 'ip' to be BELOW curr.Y (ie bottom of scanbeam) ...
        if (ip.Y <= edge1.Curr.Y)
        {
            return;
        }

        ip.Y = edge1.Curr.Y;
        //better to use the more vertical edge to derive X ...
        ip.X = TopX(Math.Abs(edge1.Dx) > Math.Abs(edge2.Dx) ? edge2 : edge1, ip.Y);
    }
    //------------------------------------------------------------------------------

    private void ProcessEdgesAtTopOfScanbeam(long topY)
    {
        TEdge e = m_ActiveEdges;
        while (e != null)
        {
            //1. process maxima, treating them as if they're 'bent' horizontal edges,
            //   but exclude maxima with horizontal edges. nb: e can't be a horizontal.
            bool IsMaximaEdge = IsMaxima(e, topY);

            switch (IsMaximaEdge)
            {
                case true:
                {
                    TEdge eMaxPair = GetMaximaPairEx(e);
                    IsMaximaEdge = eMaxPair == null || !IsHorizontal(eMaxPair);
                    break;
                }
            }

            switch (IsMaximaEdge)
            {
                case true:
                {
                    switch (StrictlySimple)
                    {
                        case true:
                            InsertMaxima(e.Top.X);
                            break;
                    }

                    TEdge ePrev = e.PrevInAEL;
                    DoMaxima(e);
                    e = ePrev switch
                    {
                        null => m_ActiveEdges,
                        _ => ePrev.NextInAEL
                    };

                    break;
                }
                default:
                {
                    //2. promote horizontal edges, otherwise update Curr.X and Curr.Y ...
                    if (IsIntermediate(e, topY) && IsHorizontal(e.NextInLML))
                    {
                        UpdateEdgeIntoAEL(ref e);
                        switch (e.OutIdx)
                        {
                            case >= 0:
                                AddOutPt(e, e.Bot);
                                break;
                        }

                        AddEdgeToSEL(e);
                    }
                    else
                    {
                        e.Curr.X = TopX(e, topY);
                        e.Curr.Y = topY;
#if use_xyz
                        if (e.Top.Y == topY)
                        {
                            e.Curr.Z = e.Top.Z;
                        }
                        else if (e.Bot.Y == topY)
                        {
                            e.Curr.Z = e.Bot.Z;
                        }
                        else
                        {
                            e.Curr.Z = 0;
                        }
#endif
                    }

                    switch (StrictlySimple)
                    {
                        //When StrictlySimple and 'e' is being touched by another edge, then
                        //make sure both edges have a vertex here ...
                        case true:
                        {
                            TEdge ePrev = e.PrevInAEL;
                            switch (e.OutIdx)
                            {
                                case >= 0 when e.WindDelta != 0 && ePrev is {OutIdx: >= 0} &&
                                               ePrev.Curr.X == e.Curr.X && ePrev.WindDelta != 0:
                                {
                                    IntPoint ip = new(e.Curr);
#if use_xyz
                                    SetZ(ref ip, ePrev, e);
#endif
                                    OutPt op = AddOutPt(ePrev, ip);
                                    OutPt op2 = AddOutPt(e, ip);
                                    AddJoin(op, op2, ip); //StrictlySimple (type-3) join
                                    break;
                                }
                            }

                            break;
                        }
                    }

                    e = e.NextInAEL;
                    break;
                }
            }
        }

        //3. Process horizontals at the Top of the scanbeam ...
        ProcessHorizontals();
        m_Maxima = null;

        //4. Promote intermediate vertices ...
        e = m_ActiveEdges;
        while (e != null)
        {
            if (IsIntermediate(e, topY))
            {
                OutPt op = e.OutIdx switch
                {
                    >= 0 => AddOutPt(e, e.Top),
                    _ => null
                };

                UpdateEdgeIntoAEL(ref e);

                //if output polygons share an edge, they'll need joining later ...
                TEdge ePrev = e.PrevInAEL;
                TEdge eNext = e.NextInAEL;
                if (ePrev != null && ePrev.Curr.X == e.Bot.X &&
                    ePrev.Curr.Y == e.Bot.Y && op != null &&
                    ePrev.OutIdx >= 0 && ePrev.Curr.Y > ePrev.Top.Y &&
                    SlopesEqual(e.Curr, e.Top, ePrev.Curr, ePrev.Top, m_UseFullRange) &&
                    e.WindDelta != 0 && ePrev.WindDelta != 0)
                {
                    OutPt op2 = AddOutPt(ePrev, e.Bot);
                    AddJoin(op, op2, e.Top);
                }
                else if (eNext != null && eNext.Curr.X == e.Bot.X &&
                         eNext.Curr.Y == e.Bot.Y && op != null &&
                         eNext.OutIdx >= 0 && eNext.Curr.Y > eNext.Top.Y &&
                         SlopesEqual(e.Curr, e.Top, eNext.Curr, eNext.Top, m_UseFullRange) &&
                         e.WindDelta != 0 && eNext.WindDelta != 0)
                {
                    OutPt op2 = AddOutPt(eNext, e.Bot);
                    AddJoin(op, op2, e.Top);
                }
            }

            e = e.NextInAEL;
        }
    }
    //------------------------------------------------------------------------------

    private void DoMaxima(TEdge e)
    {
        TEdge eMaxPair = GetMaximaPairEx(e);
        switch (eMaxPair)
        {
            case null:
            {
                switch (e.OutIdx)
                {
                    case >= 0:
                        AddOutPt(e, e.Top);
                        break;
                }

                DeleteFromAEL(e);
                return;
            }
        }

        TEdge eNext = e.NextInAEL;
        while (eNext != null && eNext != eMaxPair)
        {
            IntersectEdges(e, eNext, e.Top);
            SwapPositionsInAEL(e, eNext);
            eNext = e.NextInAEL;
        }

        switch (e.OutIdx)
        {
            case Unassigned when eMaxPair.OutIdx == Unassigned:
                DeleteFromAEL(e);
                DeleteFromAEL(eMaxPair);
                break;
            case >= 0 when eMaxPair.OutIdx >= 0:
            {
                switch (e.OutIdx)
                {
                    case >= 0:
                        AddLocalMaxPoly(e, eMaxPair, e.Top);
                        break;
                }

                DeleteFromAEL(e);
                DeleteFromAEL(eMaxPair);
                break;
            }
            default:
            {
                switch (e.WindDelta)
                {
                    case 0:
                    {
                        switch (e.OutIdx)
                        {
                            case >= 0:
                                AddOutPt(e, e.Top);
                                e.OutIdx = Unassigned;
                                break;
                        }

                        DeleteFromAEL(e);

                        switch (eMaxPair.OutIdx)
                        {
                            case >= 0:
                                AddOutPt(eMaxPair, e.Top);
                                eMaxPair.OutIdx = Unassigned;
                                break;
                        }

                        DeleteFromAEL(eMaxPair);
                        break;
                    }
                    default:
                        throw new ClipperException("DoMaxima error");
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    public static void ReversePaths(Paths polys)
    {
        foreach (Path poly in polys)
        {
            poly.Reverse();
        }
    }
    //------------------------------------------------------------------------------

    public static bool Orientation(Path poly)
    {
        return Area(poly) >= 0;
    }
    //------------------------------------------------------------------------------

    private int PointCount(OutPt pts)
    {
        switch (pts)
        {
            case null:
                return 0;
        }

        int result = 0;
        OutPt p = pts;
        do
        {
            result++;
            p = p.Next;
        } while (p != pts);

        return result;
    }
    //------------------------------------------------------------------------------

    private void BuildResult(Paths polyg)
    {
        polyg.Clear();
        polyg.Capacity = m_PolyOuts.Count;
        foreach (OutRec outRec in m_PolyOuts)
        {
            switch (outRec.Pts)
            {
                case null:
                    continue;
            }

            OutPt p = outRec.Pts.Prev;
            int cnt = PointCount(p);
            switch (cnt)
            {
                case < 2:
                    continue;
            }

            Path pg = new(cnt);
            for (int j = 0; j < cnt; j++)
            {
                pg.Add(p.Pt);
                p = p.Prev;
            }

            polyg.Add(pg);
        }
    }
    //------------------------------------------------------------------------------

    private void BuildResult2(PolyTree polytree)
    {
        polytree.Clear();

        //add each output polygon/contour to polytree ...
        polytree.m_AllPolys.Capacity = m_PolyOuts.Count;
        foreach (OutRec outRec in m_PolyOuts)
        {
            int cnt = PointCount(outRec.Pts);
            switch (outRec.IsOpen)
            {
                case true when cnt < 2:
                case false when cnt < 3:
                    continue;
            }

            FixHoleLinkage(outRec);
            PolyNode pn = new();
            polytree.m_AllPolys.Add(pn);
            outRec.PolyNode = pn;
            pn.m_polygon.Capacity = cnt;
            OutPt op = outRec.Pts.Prev;
            for (int j = 0; j < cnt; j++)
            {
                pn.m_polygon.Add(op.Pt);
                op = op.Prev;
            }
        }

        //fixup PolyNode links etc ...
        polytree.m_Childs.Capacity = m_PolyOuts.Count;
        foreach (OutRec outRec in m_PolyOuts)
        {
            switch (outRec.PolyNode)
            {
                case null:
                    continue;
            }

            switch (outRec.IsOpen)
            {
                case true:
                    outRec.PolyNode.IsOpen = true;
                    polytree.AddChild(outRec.PolyNode);
                    break;
                default:
                {
                    if (outRec.FirstLeft is {PolyNode: { }})
                    {
                        outRec.FirstLeft.PolyNode.AddChild(outRec.PolyNode);
                    }
                    else
                    {
                        polytree.AddChild(outRec.PolyNode);
                    }

                    break;
                }
            }
        }
    }
    //------------------------------------------------------------------------------

    private void FixupOutPolyline(OutRec outrec)
    {
        OutPt pp = outrec.Pts;
        OutPt lastPP = pp.Prev;
        while (pp != lastPP)
        {
            pp = pp.Next;
            if (pp.Pt != pp.Prev.Pt)
            {
                continue;
            }

            if (pp == lastPP)
            {
                lastPP = pp.Prev;
            }

            OutPt tmpPP = pp.Prev;
            tmpPP.Next = pp.Next;
            pp.Next.Prev = tmpPP;
            pp = tmpPP;
        }

        if (pp == pp.Prev)
        {
            outrec.Pts = null;
        }
    }
    //------------------------------------------------------------------------------

    private void FixupOutPolygon(OutRec outRec)
    {
        //FixupOutPolygon() - removes duplicate points and simplifies consecutive
        //parallel edges by removing the middle vertex.
        OutPt lastOK = null;
        outRec.BottomPt = null;
        OutPt pp = outRec.Pts;
        bool preserveCol = PreserveCollinear || StrictlySimple;
        for (;;)
        {
            if (pp.Prev == pp || pp.Prev == pp.Next)
            {
                outRec.Pts = null;
                return;
            }

            //test for duplicate points and collinear edges ...
            if (pp.Pt == pp.Next.Pt || pp.Pt == pp.Prev.Pt ||
                SlopesEqual(pp.Prev.Pt, pp.Pt, pp.Next.Pt, m_UseFullRange) &&
                (!preserveCol || !Pt2IsBetweenPt1AndPt3(pp.Prev.Pt, pp.Pt, pp.Next.Pt)))
            {
                lastOK = null;
                pp.Prev.Next = pp.Next;
                pp.Next.Prev = pp.Prev;
                pp = pp.Prev;
            }
            else if (pp == lastOK)
            {
                break;
            }
            else
            {
                lastOK = lastOK switch
                {
                    null => pp,
                    _ => lastOK
                };

                pp = pp.Next;
            }
        }

        outRec.Pts = pp;
    }
    //------------------------------------------------------------------------------

    private OutPt DupOutPt(OutPt outPt, bool InsertAfter)
    {
        OutPt result = new()
        {
            Pt = outPt.Pt,
            Idx = outPt.Idx
        };
        switch (InsertAfter)
        {
            case true:
                result.Next = outPt.Next;
                result.Prev = outPt;
                outPt.Next.Prev = result;
                outPt.Next = result;
                break;
            default:
                result.Prev = outPt.Prev;
                result.Next = outPt;
                outPt.Prev.Next = result;
                outPt.Prev = result;
                break;
        }

        return result;
    }
    //------------------------------------------------------------------------------

    private bool GetOverlap(long a1, long a2, long b1, long b2, out long Left, out long Right)
    {
        if (a1 < a2)
        {
            if (b1 < b2)
            {
                Left = Math.Max(a1, b1);
                Right = Math.Min(a2, b2);
            }
            else
            {
                Left = Math.Max(a1, b2);
                Right = Math.Min(a2, b1);
            }
        }
        else
        {
            if (b1 < b2)
            {
                Left = Math.Max(a2, b1);
                Right = Math.Min(a1, b2);
            }
            else
            {
                Left = Math.Max(a2, b2);
                Right = Math.Min(a1, b1);
            }
        }

        return Left < Right;
    }
    //------------------------------------------------------------------------------

    private bool JoinHorz(OutPt op1, OutPt op1b, OutPt op2, OutPt op2b,
        IntPoint Pt, bool DiscardLeft)
    {
        Direction Dir1 = op1.Pt.X > op1b.Pt.X ? Direction.dRightToLeft : Direction.dLeftToRight;
        Direction Dir2 = op2.Pt.X > op2b.Pt.X ? Direction.dRightToLeft : Direction.dLeftToRight;
        if (Dir1 == Dir2)
        {
            return false;
        }

        switch (Dir1)
        {
            //When DiscardLeft, we want Op1b to be on the Left of Op1, otherwise we
            //want Op1b to be on the Right. (And likewise with Op2 and Op2b.)
            //So, to facilitate this while inserting Op1b and Op2b ...
            //when DiscardLeft, make sure we're AT or RIGHT of Pt before adding Op1b,
            //otherwise make sure we're AT or LEFT of Pt. (Likewise with Op2b.)
            case Direction.dLeftToRight:
            {
                while (op1.Next.Pt.X <= Pt.X &&
                       op1.Next.Pt.X >= op1.Pt.X && op1.Next.Pt.Y == Pt.Y)
                {
                    op1 = op1.Next;
                }

                op1 = DiscardLeft switch
                {
                    true when op1.Pt.X != Pt.X => op1.Next,
                    _ => op1
                };

                op1b = DupOutPt(op1, !DiscardLeft);
                if (op1b.Pt != Pt)
                {
                    op1 = op1b;
                    op1.Pt = Pt;
                    op1b = DupOutPt(op1, !DiscardLeft);
                }

                break;
            }
            default:
            {
                while (op1.Next.Pt.X >= Pt.X &&
                       op1.Next.Pt.X <= op1.Pt.X && op1.Next.Pt.Y == Pt.Y)
                {
                    op1 = op1.Next;
                }

                op1 = DiscardLeft switch
                {
                    false when op1.Pt.X != Pt.X => op1.Next,
                    _ => op1
                };

                op1b = DupOutPt(op1, DiscardLeft);
                if (op1b.Pt != Pt)
                {
                    op1 = op1b;
                    op1.Pt = Pt;
                    op1b = DupOutPt(op1, DiscardLeft);
                }

                break;
            }
        }

        switch (Dir2)
        {
            case Direction.dLeftToRight:
            {
                while (op2.Next.Pt.X <= Pt.X &&
                       op2.Next.Pt.X >= op2.Pt.X && op2.Next.Pt.Y == Pt.Y)
                {
                    op2 = op2.Next;
                }

                op2 = DiscardLeft switch
                {
                    true when op2.Pt.X != Pt.X => op2.Next,
                    _ => op2
                };

                op2b = DupOutPt(op2, !DiscardLeft);
                if (op2b.Pt != Pt)
                {
                    op2 = op2b;
                    op2.Pt = Pt;
                    op2b = DupOutPt(op2, !DiscardLeft);
                }

                break;
            }
            default:
            {
                while (op2.Next.Pt.X >= Pt.X &&
                       op2.Next.Pt.X <= op2.Pt.X && op2.Next.Pt.Y == Pt.Y)
                {
                    op2 = op2.Next;
                }

                op2 = DiscardLeft switch
                {
                    false when op2.Pt.X != Pt.X => op2.Next,
                    _ => op2
                };

                op2b = DupOutPt(op2, DiscardLeft);
                if (op2b.Pt != Pt)
                {
                    op2 = op2b;
                    op2.Pt = Pt;
                    op2b = DupOutPt(op2, DiscardLeft);
                }

                break;
            }
        }

        if (Dir1 == Direction.dLeftToRight == DiscardLeft)
        {
            op1.Prev = op2;
            op2.Next = op1;
            op1b.Next = op2b;
            op2b.Prev = op1b;
        }
        else
        {
            op1.Next = op2;
            op2.Prev = op1;
            op1b.Prev = op2b;
            op2b.Next = op1b;
        }

        return true;
    }
    //------------------------------------------------------------------------------

    private bool JoinPoints(Join j, OutRec outRec1, OutRec outRec2)
    {
        OutPt op1 = j.OutPt1, op1b;
        OutPt op2 = j.OutPt2, op2b;

        //There are 3 kinds of joins for output polygons ...
        //1. Horizontal joins where Join.OutPt1 & Join.OutPt2 are vertices anywhere
        //along (horizontal) collinear edges (& Join.OffPt is on the same horizontal).
        //2. Non-horizontal joins where Join.OutPt1 & Join.OutPt2 are at the same
        //location at the Bottom of the overlapping segment (& Join.OffPt is above).
        //3. StrictlySimple joins where edges touch but are not collinear and where
        //Join.OutPt1, Join.OutPt2 & Join.OffPt all share the same point.
        bool isHorizontal = j.OutPt1.Pt.Y == j.OffPt.Y;

        switch (isHorizontal)
        {
            case true when j.OffPt == j.OutPt1.Pt && j.OffPt == j.OutPt2.Pt:
            {
                //Strictly Simple join ...
                if (outRec1 != outRec2)
                {
                    return false;
                }

                op1b = j.OutPt1.Next;
                while (op1b != op1 && op1b.Pt == j.OffPt)
                {
                    op1b = op1b.Next;
                }

                bool reverse1 = op1b.Pt.Y > j.OffPt.Y;
                op2b = j.OutPt2.Next;
                while (op2b != op2 && op2b.Pt == j.OffPt)
                {
                    op2b = op2b.Next;
                }

                bool reverse2 = op2b.Pt.Y > j.OffPt.Y;
                if (reverse1 == reverse2)
                {
                    return false;
                }

                switch (reverse1)
                {
                    case true:
                        op1b = DupOutPt(op1, false);
                        op2b = DupOutPt(op2, true);
                        op1.Prev = op2;
                        op2.Next = op1;
                        op1b.Next = op2b;
                        op2b.Prev = op1b;
                        j.OutPt1 = op1;
                        j.OutPt2 = op1b;
                        return true;
                }

                op1b = DupOutPt(op1, true);
                op2b = DupOutPt(op2, false);
                op1.Next = op2;
                op2.Prev = op1;
                op1b.Prev = op2b;
                op2b.Next = op1b;
                j.OutPt1 = op1;
                j.OutPt2 = op1b;
                return true;
            }
            case true:
            {
                //treat horizontal joins differently to non-horizontal joins since with
                //them we're not yet sure where the overlapping is. OutPt1.Pt & OutPt2.Pt
                //may be anywhere along the horizontal edge.
                op1b = op1;
                while (op1.Prev.Pt.Y == op1.Pt.Y && op1.Prev != op1b && op1.Prev != op2)
                {
                    op1 = op1.Prev;
                }

                while (op1b.Next.Pt.Y == op1b.Pt.Y && op1b.Next != op1 && op1b.Next != op2)
                {
                    op1b = op1b.Next;
                }

                if (op1b.Next == op1 || op1b.Next == op2)
                {
                    return false; //a flat 'polygon'
                }

                op2b = op2;
                while (op2.Prev.Pt.Y == op2.Pt.Y && op2.Prev != op2b && op2.Prev != op1b)
                {
                    op2 = op2.Prev;
                }

                while (op2b.Next.Pt.Y == op2b.Pt.Y && op2b.Next != op2 && op2b.Next != op1)
                {
                    op2b = op2b.Next;
                }

                if (op2b.Next == op2 || op2b.Next == op1)
                {
                    return false; //a flat 'polygon'
                }

                //Op1 -. Op1b & Op2 -. Op2b are the extremites of the horizontal edges
                if (!GetOverlap(op1.Pt.X, op1b.Pt.X, op2.Pt.X, op2b.Pt.X, out long Left, out long Right))
                {
                    return false;
                }

                //DiscardLeftSide: when overlapping edges are joined, a spike will created
                //which needs to be cleaned up. However, we don't want Op1 or Op2 caught up
                //on the discard Side as either may still be needed for other joins ...
                IntPoint Pt;
                bool DiscardLeftSide;
                if (op1.Pt.X >= Left && op1.Pt.X <= Right)
                {
                    Pt = op1.Pt;
                    DiscardLeftSide = op1.Pt.X > op1b.Pt.X;
                }
                else if (op2.Pt.X >= Left && op2.Pt.X <= Right)
                {
                    Pt = op2.Pt;
                    DiscardLeftSide = op2.Pt.X > op2b.Pt.X;
                }
                else if (op1b.Pt.X >= Left && op1b.Pt.X <= Right)
                {
                    Pt = op1b.Pt;
                    DiscardLeftSide = op1b.Pt.X > op1.Pt.X;
                }
                else
                {
                    Pt = op2b.Pt;
                    DiscardLeftSide = op2b.Pt.X > op2.Pt.X;
                }

                j.OutPt1 = op1;
                j.OutPt2 = op2;
                return JoinHorz(op1, op1b, op2, op2b, Pt, DiscardLeftSide);
            }
            default:
            {
                //nb: For non-horizontal joins ...
                //    1. Jr.OutPt1.Pt.Y == Jr.OutPt2.Pt.Y
                //    2. Jr.OutPt1.Pt > Jr.OffPt.Y

                //make sure the polygons are correctly oriented ...
                op1b = op1.Next;
                while (op1b.Pt == op1.Pt && op1b != op1)
                {
                    op1b = op1b.Next;
                }

                bool Reverse1 = op1b.Pt.Y > op1.Pt.Y ||
                                !SlopesEqual(op1.Pt, op1b.Pt, j.OffPt, m_UseFullRange);
                switch (Reverse1)
                {
                    case true:
                    {
                        op1b = op1.Prev;
                        while (op1b.Pt == op1.Pt && op1b != op1)
                        {
                            op1b = op1b.Prev;
                        }

                        if (op1b.Pt.Y > op1.Pt.Y ||
                            !SlopesEqual(op1.Pt, op1b.Pt, j.OffPt, m_UseFullRange))
                        {
                            return false;
                        }

                        break;
                    }
                }

                op2b = op2.Next;
                while (op2b.Pt == op2.Pt && op2b != op2)
                {
                    op2b = op2b.Next;
                }

                bool Reverse2 = op2b.Pt.Y > op2.Pt.Y ||
                                !SlopesEqual(op2.Pt, op2b.Pt, j.OffPt, m_UseFullRange);
                switch (Reverse2)
                {
                    case true:
                    {
                        op2b = op2.Prev;
                        while (op2b.Pt == op2.Pt && op2b != op2)
                        {
                            op2b = op2b.Prev;
                        }

                        if (op2b.Pt.Y > op2.Pt.Y ||
                            !SlopesEqual(op2.Pt, op2b.Pt, j.OffPt, m_UseFullRange))
                        {
                            return false;
                        }

                        break;
                    }
                }

                if (op1b == op1 || op2b == op2 || op1b == op2b ||
                    outRec1 == outRec2 && Reverse1 == Reverse2)
                {
                    return false;
                }

                switch (Reverse1)
                {
                    case true:
                        op1b = DupOutPt(op1, false);
                        op2b = DupOutPt(op2, true);
                        op1.Prev = op2;
                        op2.Next = op1;
                        op1b.Next = op2b;
                        op2b.Prev = op1b;
                        j.OutPt1 = op1;
                        j.OutPt2 = op1b;
                        return true;
                }

                op1b = DupOutPt(op1, true);
                op2b = DupOutPt(op2, false);
                op1.Next = op2;
                op2.Prev = op1;
                op1b.Prev = op2b;
                op2b.Next = op1b;
                j.OutPt1 = op1;
                j.OutPt2 = op1b;
                return true;
            }
        }
    }
    //----------------------------------------------------------------------

    public static int PointInPolygon(IntPoint pt, Path path)
    {
        //returns 0 if false, +1 if true, -1 if pt ON polygon boundary
        //See "The Point in Polygon Problem for Arbitrary Polygons" by Hormann & Agathos
        //http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.88.5498&rep=rep1&type=pdf
        int result = 0, cnt = path.Count;
        switch (cnt)
        {
            case < 3:
                return 0;
        }

        IntPoint ip = path[0];
        for (int i = 1; i <= cnt; ++i)
        {
            IntPoint ipNext = i == cnt ? path[0] : path[i];
            if (ipNext.Y == pt.Y)
            {
                if (ipNext.X == pt.X || ip.Y == pt.Y &&
                    ipNext.X > pt.X == ip.X < pt.X)
                {
                    return -1;
                }
            }

            if (ip.Y < pt.Y != ipNext.Y < pt.Y)
            {
                if (ip.X >= pt.X)
                {
                    if (ipNext.X > pt.X)
                    {
                        result = 1 - result;
                    }
                    else
                    {
                        double d = (double) (ip.X - pt.X) * (ipNext.Y - pt.Y) -
                                   (double) (ipNext.X - pt.X) * (ip.Y - pt.Y);
                        switch (d)
                        {
                            case 0:
                                return -1;
                        }

                        if (d > 0 == ipNext.Y > ip.Y)
                        {
                            result = 1 - result;
                        }
                    }
                }
                else
                {
                    if (ipNext.X > pt.X)
                    {
                        double d = (double) (ip.X - pt.X) * (ipNext.Y - pt.Y) -
                                   (double) (ipNext.X - pt.X) * (ip.Y - pt.Y);
                        switch (d)
                        {
                            case 0:
                                return -1;
                        }

                        if (d > 0 == ipNext.Y > ip.Y)
                        {
                            result = 1 - result;
                        }
                    }
                }
            }

            ip = ipNext;
        }

        return result;
    }
    //------------------------------------------------------------------------------

    //See "The Point in Polygon Problem for Arbitrary Polygons" by Hormann & Agathos
    //http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.88.5498&rep=rep1&type=pdf
    private static int PointInPolygon(IntPoint pt, OutPt op)
    {
        //returns 0 if false, +1 if true, -1 if pt ON polygon boundary
        int result = 0;
        OutPt startOp = op;
        long ptx = pt.X, pty = pt.Y;
        long poly0x = op.Pt.X, poly0y = op.Pt.Y;
        do
        {
            op = op.Next;
            long poly1x = op.Pt.X, poly1y = op.Pt.Y;

            if (poly1y == pty)
            {
                if (poly1x == ptx || poly0y == pty &&
                    poly1x > ptx == poly0x < ptx)
                {
                    return -1;
                }
            }

            if (poly0y < pty != poly1y < pty)
            {
                if (poly0x >= ptx)
                {
                    if (poly1x > ptx)
                    {
                        result = 1 - result;
                    }
                    else
                    {
                        double d = (double) (poly0x - ptx) * (poly1y - pty) -
                                   (double) (poly1x - ptx) * (poly0y - pty);
                        switch (d)
                        {
                            case 0:
                                return -1;
                        }

                        if (d > 0 == poly1y > poly0y)
                        {
                            result = 1 - result;
                        }
                    }
                }
                else
                {
                    if (poly1x > ptx)
                    {
                        double d = (double) (poly0x - ptx) * (poly1y - pty) -
                                   (double) (poly1x - ptx) * (poly0y - pty);
                        switch (d)
                        {
                            case 0:
                                return -1;
                        }

                        if (d > 0 == poly1y > poly0y)
                        {
                            result = 1 - result;
                        }
                    }
                }
            }

            poly0x = poly1x;
            poly0y = poly1y;
        } while (startOp != op);

        return result;
    }
    //------------------------------------------------------------------------------

    private static bool Poly2ContainsPoly1(OutPt outPt1, OutPt outPt2)
    {
        OutPt op = outPt1;
        do
        {
            //nb: PointInPolygon returns 0 if false, +1 if true, -1 if pt on polygon
            int res = PointInPolygon(op.Pt, outPt2);
            switch (res)
            {
                case >= 0:
                    return res > 0;
                default:
                    op = op.Next;
                    break;
            }
        } while (op != outPt1);

        return true;
    }
    //----------------------------------------------------------------------

    private void FixupFirstLefts1(OutRec OldOutRec, OutRec NewOutRec)
    {
        foreach (OutRec outRec in from outRec in m_PolyOuts
                 let firstLeft = ParseFirstLeft(outRec.FirstLeft)
                 where outRec.Pts != null && firstLeft == OldOutRec
                 where Poly2ContainsPoly1(outRec.Pts, NewOutRec.Pts)
                 select outRec)
        {
            outRec.FirstLeft = NewOutRec;
        }
    }
    //----------------------------------------------------------------------

    private void FixupFirstLefts2(OutRec innerOutRec, OutRec outerOutRec)
    {
        //A polygon has split into two such that one is now the inner of the other.
        //It's possible that these polygons now wrap around other polygons, so check
        //every polygon that's also contained by OuterOutRec's FirstLeft container
        //(including nil) to see if they've become inner to the new inner polygon ...
        OutRec orfl = outerOutRec.FirstLeft;
        foreach (OutRec outRec in from outRec in m_PolyOuts
                 where outRec.Pts != null && outRec != outerOutRec && outRec != innerOutRec
                 let firstLeft = ParseFirstLeft(outRec.FirstLeft)
                 where firstLeft == orfl || firstLeft == innerOutRec || firstLeft == outerOutRec
                 select outRec)
        {
            if (Poly2ContainsPoly1(outRec.Pts, innerOutRec.Pts))
            {
                outRec.FirstLeft = innerOutRec;
            }
            else if (Poly2ContainsPoly1(outRec.Pts, outerOutRec.Pts))
            {
                outRec.FirstLeft = outerOutRec;
            }
            else if (outRec.FirstLeft == innerOutRec || outRec.FirstLeft == outerOutRec)
            {
                outRec.FirstLeft = orfl;
            }
        }
    }
    //----------------------------------------------------------------------

    private void FixupFirstLefts3(OutRec OldOutRec, OutRec NewOutRec)
    {
        //same as FixupFirstLefts1 but doesn't call Poly2ContainsPoly1()
        foreach (OutRec outRec in from outRec in m_PolyOuts
                 let firstLeft = ParseFirstLeft(outRec.FirstLeft)
                 where outRec.Pts != null && firstLeft == OldOutRec
                 select outRec)
        {
            outRec.FirstLeft = NewOutRec;
        }
    }
    //----------------------------------------------------------------------

    private static OutRec ParseFirstLeft(OutRec FirstLeft)
    {
        while (FirstLeft is {Pts: null})
        {
            FirstLeft = FirstLeft.FirstLeft;
        }

        return FirstLeft;
    }
    //------------------------------------------------------------------------------

    private void JoinCommonEdges()
    {
        foreach (Join join in m_Joins)
        {
            OutRec outRec1 = GetOutRec(join.OutPt1.Idx);
            OutRec outRec2 = GetOutRec(join.OutPt2.Idx);

            if (outRec1.Pts == null || outRec2.Pts == null)
            {
                continue;
            }

            if (outRec1.IsOpen || outRec2.IsOpen)
            {
                continue;
            }

            //get the polygon fragment with the correct hole state (FirstLeft)
            //before calling JoinPoints() ...
            OutRec holeStateRec;
            if (outRec1 == outRec2)
            {
                holeStateRec = outRec1;
            }
            else if (OutRec1RightOfOutRec2(outRec1, outRec2))
            {
                holeStateRec = outRec2;
            }
            else if (OutRec1RightOfOutRec2(outRec2, outRec1))
            {
                holeStateRec = outRec1;
            }
            else
            {
                holeStateRec = GetLowermostRec(outRec1, outRec2);
            }

            if (!JoinPoints(join, outRec1, outRec2))
            {
                continue;
            }

            if (outRec1 == outRec2)
            {
                //instead of joining two polygons, we've just created a new one by
                //splitting one polygon into two.
                outRec1.Pts = join.OutPt1;
                outRec1.BottomPt = null;
                outRec2 = CreateOutRec();
                outRec2.Pts = join.OutPt2;

                //update all OutRec2.Pts Idx's ...
                UpdateOutPtIdxs(outRec2);

                if (Poly2ContainsPoly1(outRec2.Pts, outRec1.Pts))
                {
                    //outRec1 contains outRec2 ...
                    outRec2.IsHole = !outRec1.IsHole;
                    outRec2.FirstLeft = outRec1;

                    switch (m_UsingPolyTree)
                    {
                        case true:
                            FixupFirstLefts2(outRec2, outRec1);
                            break;
                    }

                    if ((outRec2.IsHole ^ ReverseSolution) == Area(outRec2) > 0)
                    {
                        ReversePolyPtLinks(outRec2.Pts);
                    }
                }
                else if (Poly2ContainsPoly1(outRec1.Pts, outRec2.Pts))
                {
                    //outRec2 contains outRec1 ...
                    outRec2.IsHole = outRec1.IsHole;
                    outRec1.IsHole = !outRec2.IsHole;
                    outRec2.FirstLeft = outRec1.FirstLeft;
                    outRec1.FirstLeft = outRec2;

                    switch (m_UsingPolyTree)
                    {
                        case true:
                            FixupFirstLefts2(outRec1, outRec2);
                            break;
                    }

                    if ((outRec1.IsHole ^ ReverseSolution) == Area(outRec1) > 0)
                    {
                        ReversePolyPtLinks(outRec1.Pts);
                    }
                }
                else
                {
                    //the 2 polygons are completely separate ...
                    outRec2.IsHole = outRec1.IsHole;
                    outRec2.FirstLeft = outRec1.FirstLeft;

                    switch (m_UsingPolyTree)
                    {
                        //fixup FirstLeft pointers that may need reassigning to OutRec2
                        case true:
                            FixupFirstLefts1(outRec1, outRec2);
                            break;
                    }
                }
            }
            else
            {
                //joined 2 polygons together ...

                outRec2.Pts = null;
                outRec2.BottomPt = null;
                outRec2.Idx = outRec1.Idx;

                outRec1.IsHole = holeStateRec.IsHole;
                if (holeStateRec == outRec2)
                {
                    outRec1.FirstLeft = outRec2.FirstLeft;
                }

                outRec2.FirstLeft = outRec1;

                switch (m_UsingPolyTree)
                {
                    //fixup FirstLeft pointers that may need reassigning to OutRec1
                    case true:
                        FixupFirstLefts3(outRec2, outRec1);
                        break;
                }
            }
        }
    }
    //------------------------------------------------------------------------------

    private void UpdateOutPtIdxs(OutRec outrec)
    {
        OutPt op = outrec.Pts;
        do
        {
            op.Idx = outrec.Idx;
            op = op.Prev;
        } while (op != outrec.Pts);
    }
    //------------------------------------------------------------------------------

    private void DoSimplePolygons()
    {
        int i = 0;
        while (i < m_PolyOuts.Count)
        {
            OutRec outrec = m_PolyOuts[i++];
            OutPt op = outrec.Pts;
            if (op == null || outrec.IsOpen)
            {
                continue;
            }

            do //for each Pt in Polygon until duplicate found do ...
            {
                OutPt op2 = op.Next;
                while (op2 != outrec.Pts)
                {
                    if (op.Pt == op2.Pt && op2.Next != op && op2.Prev != op)
                    {
                        //split the polygon into two ...
                        OutPt op3 = op.Prev;
                        OutPt op4 = op2.Prev;
                        op.Prev = op4;
                        op4.Next = op;
                        op2.Prev = op3;
                        op3.Next = op2;

                        outrec.Pts = op;
                        OutRec outrec2 = CreateOutRec();
                        outrec2.Pts = op2;
                        UpdateOutPtIdxs(outrec2);
                        if (Poly2ContainsPoly1(outrec2.Pts, outrec.Pts))
                        {
                            //OutRec2 is contained by OutRec1 ...
                            outrec2.IsHole = !outrec.IsHole;
                            outrec2.FirstLeft = outrec;
                            switch (m_UsingPolyTree)
                            {
                                case true:
                                    FixupFirstLefts2(outrec2, outrec);
                                    break;
                            }
                        }
                        else if (Poly2ContainsPoly1(outrec.Pts, outrec2.Pts))
                        {
                            //OutRec1 is contained by OutRec2 ...
                            outrec2.IsHole = outrec.IsHole;
                            outrec.IsHole = !outrec2.IsHole;
                            outrec2.FirstLeft = outrec.FirstLeft;
                            outrec.FirstLeft = outrec2;
                            switch (m_UsingPolyTree)
                            {
                                case true:
                                    FixupFirstLefts2(outrec, outrec2);
                                    break;
                            }
                        }
                        else
                        {
                            //the 2 polygons are separate ...
                            outrec2.IsHole = outrec.IsHole;
                            outrec2.FirstLeft = outrec.FirstLeft;
                            switch (m_UsingPolyTree)
                            {
                                case true:
                                    FixupFirstLefts1(outrec, outrec2);
                                    break;
                            }
                        }

                        op2 = op; //ie get ready for the next iteration
                    }

                    op2 = op2.Next;
                }

                op = op.Next;
            } while (op != outrec.Pts);
        }
    }
    //------------------------------------------------------------------------------

    public static double Area(Path poly)
    {
        int cnt = poly.Count;
        switch (cnt)
        {
            case < 3:
                return 0;
        }

        double a = 0;
        for (int i = 0, j = cnt - 1; i < cnt; ++i)
        {
            a += ((double) poly[j].X + poly[i].X) * ((double) poly[j].Y - poly[i].Y);
            j = i;
        }

        return -a * 0.5;
    }
    //------------------------------------------------------------------------------

    internal double Area(OutRec outRec)
    {
        return Area(outRec.Pts);
    }
    //------------------------------------------------------------------------------

    internal double Area(OutPt op)
    {
        OutPt opFirst = op;
        switch (op)
        {
            case null:
                return 0;
        }

        double a = 0;
        do
        {
            a += (op.Prev.Pt.X + op.Pt.X) * (double) (op.Prev.Pt.Y - op.Pt.Y);
            op = op.Next;
        } while (op != opFirst);

        return a * 0.5;
    }

    //------------------------------------------------------------------------------
    // SimplifyPolygon functions ...
    // Convert self-intersecting polygons into simple polygons
    //------------------------------------------------------------------------------

    public static Paths SimplifyPolygon(Path poly,
        PolyFillType fillType = PolyFillType.pftEvenOdd, bool preserveColinear = false)
    {
        Paths result = new();
        Clipper c = new()
        {
            StrictlySimple = true,
            PreserveCollinear = preserveColinear
        };
        c.AddPath(poly, PolyType.ptSubject, true);
        c.Execute(ClipType.ctUnion, result, fillType, fillType);
        return result;
    }
    //------------------------------------------------------------------------------

    public static Paths SimplifyPolygons(Paths polys,
        PolyFillType fillType = PolyFillType.pftEvenOdd, bool preserveColinear = false)
    {
        Paths result = new();
        Clipper c = new()
        {
            StrictlySimple = true,
            PreserveCollinear = preserveColinear
        };
        c.AddPaths(polys, PolyType.ptSubject, true);
        c.Execute(ClipType.ctUnion, result, fillType, fillType);
        return result;
    }

    //------------------------------------------------------------------------------

    private static double DistanceFromLineSqrd(IntPoint pt, IntPoint ln1, IntPoint ln2)
    {
        //The equation of a line in general form (Ax + By + C = 0)
        //given 2 points (x¹,y¹) & (x²,y²) is ...
        //(y¹ - y²)x + (x² - x¹)y + (y² - y¹)x¹ - (x² - x¹)y¹ = 0
        //A = (y¹ - y²); B = (x² - x¹); C = (y² - y¹)x¹ - (x² - x¹)y¹
        //perpendicular distance of point (x³,y³) = (Ax³ + By³ + C)/Sqrt(A² + B²)
        //see http://en.wikipedia.org/wiki/Perpendicular_distance
        double A = ln1.Y - ln2.Y;
        double B = ln2.X - ln1.X;
        double C = A * ln1.X + B * ln1.Y;
        C = A * pt.X + B * pt.Y - C;
        return C * C / (A * A + B * B);
    }
    //---------------------------------------------------------------------------

    private static bool SlopesNearCollinear(IntPoint pt1,
        IntPoint pt2, IntPoint pt3, double distSqrd)
    {
        //this function is more accurate when the point that's GEOMETRICALLY 
        //between the other 2 points is the one that's tested for distance.  
        //nb: with 'spikes', either pt1 or pt3 is geometrically between the other pts                    
        if (Math.Abs(pt1.X - pt2.X) > Math.Abs(pt1.Y - pt2.Y))
        {
            if (pt1.X > pt2.X == pt1.X < pt3.X)
            {
                return DistanceFromLineSqrd(pt1, pt2, pt3) < distSqrd;
            }

            if (pt2.X > pt1.X == pt2.X < pt3.X)
            {
                return DistanceFromLineSqrd(pt2, pt1, pt3) < distSqrd;
            }

            return DistanceFromLineSqrd(pt3, pt1, pt2) < distSqrd;
        }

        if (pt1.Y > pt2.Y == pt1.Y < pt3.Y)
        {
            return DistanceFromLineSqrd(pt1, pt2, pt3) < distSqrd;
        }

        if (pt2.Y > pt1.Y == pt2.Y < pt3.Y)
        {
            return DistanceFromLineSqrd(pt2, pt1, pt3) < distSqrd;
        }

        return DistanceFromLineSqrd(pt3, pt1, pt2) < distSqrd;
    }
    //------------------------------------------------------------------------------

    private static bool PointsAreClose(IntPoint pt1, IntPoint pt2, double distSqrd)
    {
        double dx = (double) pt1.X - pt2.X;
        double dy = (double) pt1.Y - pt2.Y;
        return dx * dx + dy * dy <= distSqrd;
    }
    //------------------------------------------------------------------------------

    private static OutPt ExcludeOp(OutPt op)
    {
        OutPt result = op.Prev;
        result.Next = op.Next;
        op.Next.Prev = result;
        result.Idx = 0;
        return result;
    }
    //------------------------------------------------------------------------------

    public static Path CleanPolygon(Path path, double distance = 1.415)
    {
        //distance = proximity in units/pixels below which vertices will be stripped. 
        //Default ~= sqrt(2) so when adjacent vertices or semi-adjacent vertices have 
        //both x & y coords within 1 unit, then the second vertex will be stripped.

        int cnt = path.Count;

        switch (cnt)
        {
            case 0:
                return new Path();
        }

        OutPt[] outPts = new OutPt[cnt];
        for (int i = 0; i < cnt; ++i) outPts[i] = new OutPt();

        for (int i = 0; i < cnt; ++i)
        {
            outPts[i].Pt = path[i];
            outPts[i].Next = outPts[(i + 1) % cnt];
            outPts[i].Next.Prev = outPts[i];
            outPts[i].Idx = 0;
        }

        double distSqrd = distance * distance;
        OutPt op = outPts[0];
        while (op.Idx == 0 && op.Next != op.Prev)
        {
            if (PointsAreClose(op.Pt, op.Prev.Pt, distSqrd))
            {
                op = ExcludeOp(op);
                cnt--;
            }
            else if (PointsAreClose(op.Prev.Pt, op.Next.Pt, distSqrd))
            {
                ExcludeOp(op.Next);
                op = ExcludeOp(op);
                cnt -= 2;
            }
            else if (SlopesNearCollinear(op.Prev.Pt, op.Pt, op.Next.Pt, distSqrd))
            {
                op = ExcludeOp(op);
                cnt--;
            }
            else
            {
                op.Idx = 1;
                op = op.Next;
            }
        }

        cnt = cnt switch
        {
            < 3 => 0,
            _ => cnt
        };

        Path result = new(cnt);
        for (int i = 0; i < cnt; ++i)
        {
            result.Add(op.Pt);
            op = op.Next;
        }

        return result;
    }
    //------------------------------------------------------------------------------

    public static Paths CleanPolygons(Paths polys,
        double distance = 1.415)
    {
        Paths result = new(polys.Count);
        result.AddRange(polys.Select(t => CleanPolygon(t, distance)));

        return result;
    }
    //------------------------------------------------------------------------------

    internal static Paths Minkowski(Path pattern, Path path, bool IsSum, bool IsClosed)
    {
        int delta = IsClosed ? 1 : 0;
        int polyCnt = pattern.Count;
        int pathCnt = path.Count;
        Paths result = new(pathCnt);
        switch (IsSum)
        {
            case true:
            {
                for (int i = 0; i < pathCnt; i++)
                {
                    Path p = new(polyCnt);
                    p.AddRange(pattern.Select(ip => new IntPoint(path[i].X + ip.X, path[i].Y + ip.Y)));
                    result.Add(p);
                }

                break;
            }
            default:
            {
                for (int i = 0; i < pathCnt; i++)
                {
                    Path p = new(polyCnt);
                    p.AddRange(pattern.Select(ip => new IntPoint(path[i].X - ip.X, path[i].Y - ip.Y)));
                    result.Add(p);
                }

                break;
            }
        }

        Paths quads = new((pathCnt + delta) * (polyCnt + 1));
        for (int i = 0; i < pathCnt - 1 + delta; i++)
        for (int j = 0; j < polyCnt; j++)
        {
            Path quad = new(4)
            {
                result[i % pathCnt][j % polyCnt],
                result[(i + 1) % pathCnt][j % polyCnt],
                result[(i + 1) % pathCnt][(j + 1) % polyCnt],
                result[i % pathCnt][(j + 1) % polyCnt]
            };
            if (!Orientation(quad))
            {
                quad.Reverse();
            }

            quads.Add(quad);
        }

        return quads;
    }
    //------------------------------------------------------------------------------

    public static Paths MinkowskiSum(Path pattern, Path path, bool pathIsClosed)
    {
        Paths paths = Minkowski(pattern, path, true, pathIsClosed);
        Clipper c = new();
        c.AddPaths(paths, PolyType.ptSubject, true);
        c.Execute(ClipType.ctUnion, paths, PolyFillType.pftNonZero, PolyFillType.pftNonZero);
        return paths;
    }
    //------------------------------------------------------------------------------

    private static Path TranslatePath(Path path, IntPoint delta)
    {
        Path outPath = new(path.Count);
        for (int i = 0; i < path.Count; i++)
            outPath.Add(new IntPoint(path[i].X + delta.X, path[i].Y + delta.Y));
        return outPath;
    }
    //------------------------------------------------------------------------------

    public static Paths MinkowskiSum(Path pattern, Paths paths, bool pathIsClosed)
    {
        Paths solution = new();
        Clipper c = new();
        foreach (Path t in paths)
        {
            Paths tmp = Minkowski(pattern, t, true, pathIsClosed);
            c.AddPaths(tmp, PolyType.ptSubject, true);
            switch (pathIsClosed)
            {
                case true:
                {
                    Path path = TranslatePath(t, pattern[0]);
                    c.AddPath(path, PolyType.ptClip, true);
                    break;
                }
            }
        }

        c.Execute(ClipType.ctUnion, solution,
            PolyFillType.pftNonZero, PolyFillType.pftNonZero);
        return solution;
    }
    //------------------------------------------------------------------------------

    public static Paths MinkowskiDiff(Path poly1, Path poly2)
    {
        Paths paths = Minkowski(poly1, poly2, false, true);
        Clipper c = new();
        c.AddPaths(paths, PolyType.ptSubject, true);
        c.Execute(ClipType.ctUnion, paths, PolyFillType.pftNonZero, PolyFillType.pftNonZero);
        return paths;
    }
    //------------------------------------------------------------------------------

    internal enum NodeType
    {
        ntAny,
        ntOpen,
        ntClosed
    }

    public static Paths PolyTreeToPaths(PolyTree polytree)
    {
        Paths result = new()
        {
            Capacity = polytree.Total
        };
        AddPolyNodeToPaths(polytree, NodeType.ntAny, result);
        return result;
    }
    //------------------------------------------------------------------------------

    internal static void AddPolyNodeToPaths(PolyNode polynode, NodeType nt, Paths paths)
    {
        bool match = true;
        switch (nt)
        {
            case NodeType.ntOpen: return;
            case NodeType.ntClosed:
                match = !polynode.IsOpen;
                break;
        }

        switch (polynode.m_polygon.Count)
        {
            case > 0 when match:
                paths.Add(polynode.m_polygon);
                break;
        }

        foreach (PolyNode pn in polynode.Childs)
            AddPolyNodeToPaths(pn, nt, paths);
    }
    //------------------------------------------------------------------------------

    public static Paths OpenPathsFromPolyTree(PolyTree polytree)
    {
        Paths result = new()
        {
            Capacity = polytree.ChildCount
        };
        for (int i = 0; i < polytree.ChildCount; i++)
            switch (polytree.Childs[i].IsOpen)
            {
                case true:
                    result.Add(polytree.Childs[i].m_polygon);
                    break;
            }

        return result;
    }
    //------------------------------------------------------------------------------

    public static Paths ClosedPathsFromPolyTree(PolyTree polytree)
    {
        Paths result = new()
        {
            Capacity = polytree.Total
        };
        AddPolyNodeToPaths(polytree, NodeType.ntClosed, result);
        return result;
    }
    //------------------------------------------------------------------------------
} //end Clipper

public class ClipperOffset
{
    private Paths m_destPolys;
    private Path m_srcPoly;
    private Path m_destPoly;
    private List<DoublePoint> m_normals = new();
    private double m_delta, m_sinA, m_sin, m_cos;
    private double m_miterLim, m_StepsPerRad;

    private IntPoint m_lowest;
    private PolyNode m_polyNodes = new();

    public double ArcTolerance { get; set; }
    public double MiterLimit { get; set; }

    public bool PreserveCollinear { get; set; }

    private const double two_pi = Math.PI * 2;
    private const double def_arc_tolerance = 0.25;

    public ClipperOffset(
        double miterLimit = 2.0, double arcTolerance = def_arc_tolerance, bool def_preserveColinear = false)
    {
        MiterLimit = miterLimit;
        ArcTolerance = arcTolerance;
        m_lowest.X = -1;
        PreserveCollinear = def_preserveColinear;
    }
    //------------------------------------------------------------------------------

    public void Clear()
    {
        m_polyNodes.Childs.Clear();
        m_lowest.X = -1;
    }
    //------------------------------------------------------------------------------

    internal static long Round(double value)
    {
        return value < 0 ? (long) (value - 0.5) : (long) (value + 0.5);
    }
    //------------------------------------------------------------------------------

    public void AddPath(Path path, JoinType joinType, EndType endType)
    {
        int highI = path.Count - 1;
        switch (highI)
        {
            case < 0:
                return;
        }

        PolyNode newNode = new()
        {
            m_jointype = joinType,
            m_endtype = endType
        };

        switch (endType)
        {
            //strip duplicate points from path and also get index to the lowest point ...
            case EndType.etClosedLine:
            case EndType.etClosedPolygon:
            {
                while (highI > 0 && path[0] == path[highI])
                {
                    highI--;
                }

                break;
            }
        }

        newNode.m_polygon.Capacity = highI + 1;
        newNode.m_polygon.Add(path[0]);
        int j = 0, k = 0;
        for (int i = 1; i <= highI; i++)
            if ( /*(ioPreserveCollinear != 0) || */newNode.m_polygon[j] != path[i])
            {
                j++;
                newNode.m_polygon.Add(path[i]);
                if (path[i].Y > newNode.m_polygon[k].Y ||
                    path[i].Y == newNode.m_polygon[k].Y &&
                    path[i].X < newNode.m_polygon[k].X)
                {
                    k = j;
                }
            }

        switch (endType)
        {
            case EndType.etClosedPolygon when j < 2:
                return;
        }

        m_polyNodes.AddChild(newNode);

        //if this path's lowest pt is lower than all the others then update m_lowest
        if (endType != EndType.etClosedPolygon)
        {
            return;
        }

        switch (m_lowest.X)
        {
            case < 0:
                m_lowest = new IntPoint(m_polyNodes.ChildCount - 1, k);
                break;
            default:
            {
                IntPoint ip = m_polyNodes.Childs[(int) m_lowest.X].m_polygon[(int) m_lowest.Y];
                if (newNode.m_polygon[k].Y > ip.Y ||
                    newNode.m_polygon[k].Y == ip.Y &&
                    newNode.m_polygon[k].X < ip.X)
                {
                    m_lowest = new IntPoint(m_polyNodes.ChildCount - 1, k);
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    public void AddPaths(Paths paths, JoinType joinType, EndType endType)
    {
        foreach (Path p in paths)
            AddPath(p, joinType, endType);
    }
    //------------------------------------------------------------------------------

    private void FixOrientations()
    {
        switch (m_lowest.X)
        {
            //fixup orientations of all closed paths if the orientation of the
            //closed path with the lowermost vertex is wrong ...
            case >= 0 when !Clipper.Orientation(m_polyNodes.Childs[(int) m_lowest.X].m_polygon):
            {
                for (int i = 0; i < m_polyNodes.ChildCount; i++)
                {
                    PolyNode node = m_polyNodes.Childs[i];
                    switch (node.m_endtype)
                    {
                        case EndType.etClosedPolygon:
                        case EndType.etClosedLine when Clipper.Orientation(node.m_polygon):
                            node.m_polygon.Reverse();
                            break;
                    }
                }

                break;
            }
            default:
            {
                for (int i = 0; i < m_polyNodes.ChildCount; i++)
                {
                    PolyNode node = m_polyNodes.Childs[i];
                    switch (node.m_endtype)
                    {
                        case EndType.etClosedLine when !Clipper.Orientation(node.m_polygon):
                            node.m_polygon.Reverse();
                            break;
                    }
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    internal static DoublePoint GetUnitNormal(IntPoint pt1, IntPoint pt2)
    {
        double dx = pt2.X - pt1.X;
        double dy = pt2.Y - pt1.Y;
        switch (dx)
        {
            case 0 when dy == 0:
                return new DoublePoint();
        }

        double f = 1 * 1.0 / Math.Sqrt(dx * dx + dy * dy);
        dx *= f;
        dy *= f;

        return new DoublePoint(dy, -dx);
    }
    //------------------------------------------------------------------------------

    private void DoOffset(double delta)
    {
        m_destPolys = new Paths();
        m_delta = delta;

        //if Zero offset, just copy any CLOSED polygons to m_p and return ...
        if (ClipperBase.near_zero(delta))
        {
            m_destPolys.Capacity = m_polyNodes.ChildCount;
            for (int i = 0; i < m_polyNodes.ChildCount; i++)
            {
                PolyNode node = m_polyNodes.Childs[i];
                switch (node.m_endtype)
                {
                    case EndType.etClosedPolygon:
                        m_destPolys.Add(node.m_polygon);
                        break;
                }
            }

            return;
        }

        m_miterLim = MiterLimit switch
        {
            //see offset_triginometry3.svg in the documentation folder ...
            > 2 => 2 / (MiterLimit * MiterLimit),
            _ => 0.5
        };

        double y;
        switch (ArcTolerance)
        {
            case <= 0.0:
                y = def_arc_tolerance;
                break;
            default:
            {
                if (ArcTolerance > Math.Abs(delta) * def_arc_tolerance)
                {
                    y = Math.Abs(delta) * def_arc_tolerance;
                }
                else
                {
                    y = ArcTolerance;
                }

                break;
            }
        }

        //see offset_triginometry2.svg in the documentation folder ...
        double steps = Math.PI / Math.Acos(1 - y / Math.Abs(delta));
        m_sin = Math.Sin(two_pi / steps);
        m_cos = Math.Cos(two_pi / steps);
        m_StepsPerRad = steps / two_pi;
        m_sin = delta switch
        {
            < 0.0 => -m_sin,
            _ => m_sin
        };

        m_destPolys.Capacity = m_polyNodes.ChildCount * 2;
        for (int i = 0; i < m_polyNodes.ChildCount; i++)
        {
            PolyNode node = m_polyNodes.Childs[i];
            m_srcPoly = node.m_polygon;

            int len = m_srcPoly.Count;

            if (len == 0 || delta <= 0 && (len < 3 ||
                                           node.m_endtype != EndType.etClosedPolygon))
            {
                continue;
            }

            m_destPoly = new Path();

            switch (len)
            {
                case 1:
                {
                    switch (node.m_jointype)
                    {
                        case JoinType.jtRound:
                        {
                            double X = 1.0, Y = 0.0;
                            for (int j = 1; j <= steps; j++)
                            {
                                m_destPoly.Add(new IntPoint(
                                    Round(m_srcPoly[0].X + X * delta),
                                    Round(m_srcPoly[0].Y + Y * delta)));
                                double X2 = X;
                                X = X * m_cos - m_sin * Y;
                                Y = X2 * m_sin + Y * m_cos;
                            }

                            break;
                        }
                        default:
                        {
                            double X = -1.0, Y = -1.0;
                            for (int j = 0; j < 4; ++j)
                            {
                                m_destPoly.Add(new IntPoint(
                                    Round(m_srcPoly[0].X + X * delta),
                                    Round(m_srcPoly[0].Y + Y * delta)));
                                switch (X)
                                {
                                    case < 0:
                                        X = 1;
                                        break;
                                    default:
                                    {
                                        switch (Y)
                                        {
                                            case < 0:
                                                Y = 1;
                                                break;
                                            default:
                                                X = -1;
                                                break;
                                        }

                                        break;
                                    }
                                }
                            }

                            break;
                        }
                    }

                    m_destPolys.Add(m_destPoly);
                    continue;
                }
            }

            //build m_normals ...
            m_normals.Clear();
            m_normals.Capacity = len;
            for (int j = 0; j < len - 1; j++)
                m_normals.Add(GetUnitNormal(m_srcPoly[j], m_srcPoly[j + 1]));
            switch (node.m_endtype)
            {
                case EndType.etClosedLine:
                case EndType.etClosedPolygon:
                    m_normals.Add(GetUnitNormal(m_srcPoly[len - 1], m_srcPoly[0]));
                    break;
                default:
                    m_normals.Add(new DoublePoint(m_normals[len - 2]));
                    break;
            }

            switch (node.m_endtype)
            {
                case EndType.etClosedPolygon:
                {
                    int k = len - 1;
                    for (int j = 0; j < len; j++)
                        OffsetPoint(j, ref k, node.m_jointype);
                    m_destPolys.Add(m_destPoly);
                    break;
                }
                case EndType.etClosedLine:
                {
                    int k = len - 1;
                    for (int j = 0; j < len; j++)
                        OffsetPoint(j, ref k, node.m_jointype);
                    m_destPolys.Add(m_destPoly);
                    m_destPoly = new Path();
                    //re-build m_normals ...
                    DoublePoint n = m_normals[len - 1];
                    for (int j = len - 1; j > 0; j--)
                        m_normals[j] = new DoublePoint(-m_normals[j - 1].X, -m_normals[j - 1].Y);
                    m_normals[0] = new DoublePoint(-n.X, -n.Y);
                    k = 0;
                    for (int j = len - 1; j >= 0; j--)
                        OffsetPoint(j, ref k, node.m_jointype);
                    m_destPolys.Add(m_destPoly);
                    break;
                }
                default:
                {
                    int k = 0;
                    for (int j = 1; j < len - 1; ++j)
                        OffsetPoint(j, ref k, node.m_jointype);

                    IntPoint pt1;
                    switch (node.m_endtype)
                    {
                        case EndType.etOpenButt:
                        {
                            int j = len - 1;
                            pt1 = new IntPoint(Round(m_srcPoly[j].X + m_normals[j].X *
                                delta), Round(m_srcPoly[j].Y + m_normals[j].Y * delta));
                            m_destPoly.Add(pt1);
                            pt1 = new IntPoint(Round(m_srcPoly[j].X - m_normals[j].X *
                                delta), Round(m_srcPoly[j].Y - m_normals[j].Y * delta));
                            m_destPoly.Add(pt1);
                            break;
                        }
                        default:
                        {
                            int j = len - 1;
                            k = len - 2;
                            m_sinA = 0;
                            m_normals[j] = new DoublePoint(-m_normals[j].X, -m_normals[j].Y);
                            switch (node.m_endtype)
                            {
                                case EndType.etOpenSquare:
                                    DoSquare(j, k);
                                    break;
                                default:
                                    DoRound(j, k);
                                    break;
                            }

                            break;
                        }
                    }

                    //re-build m_normals ...
                    for (int j = len - 1; j > 0; j--)
                        m_normals[j] = new DoublePoint(-m_normals[j - 1].X, -m_normals[j - 1].Y);

                    m_normals[0] = new DoublePoint(-m_normals[1].X, -m_normals[1].Y);

                    k = len - 1;
                    for (int j = k - 1; j > 0; --j)
                        OffsetPoint(j, ref k, node.m_jointype);

                    switch (node.m_endtype)
                    {
                        case EndType.etOpenButt:
                            pt1 = new IntPoint(Round(m_srcPoly[0].X - m_normals[0].X * delta),
                                Round(m_srcPoly[0].Y - m_normals[0].Y * delta));
                            m_destPoly.Add(pt1);
                            pt1 = new IntPoint(Round(m_srcPoly[0].X + m_normals[0].X * delta),
                                Round(m_srcPoly[0].Y + m_normals[0].Y * delta));
                            m_destPoly.Add(pt1);
                            break;
                        default:
                        {
                            m_sinA = 0;
                            switch (node.m_endtype)
                            {
                                case EndType.etOpenSquare:
                                    DoSquare(0, 1);
                                    break;
                                default:
                                    DoRound(0, 1);
                                    break;
                            }

                            break;
                        }
                    }

                    m_destPolys.Add(m_destPoly);
                    break;
                }
            }
        }
    }
    //------------------------------------------------------------------------------

    public void Execute(ref Paths solution, double delta)
    {
        solution.Clear();
        FixOrientations();
        DoOffset(delta);
        //now clean up 'corners' ...
        Clipper clpr = new()
        {
            PreserveCollinear = PreserveCollinear
        };
        clpr.AddPaths(m_destPolys, PolyType.ptSubject, true);
        switch (delta)
        {
            case > 0:
                clpr.Execute(ClipType.ctUnion, solution,
                    PolyFillType.pftPositive, PolyFillType.pftPositive);
                break;
            default:
            {
                IntRect r = ClipperBase.GetBounds(m_destPolys);
                Path outer = new(4)
                {
                    new IntPoint(r.left - 10, r.bottom + 10),
                    new IntPoint(r.right + 10, r.bottom + 10),
                    new IntPoint(r.right + 10, r.top - 10),
                    new IntPoint(r.left - 10, r.top - 10)
                };

                clpr.AddPath(outer, PolyType.ptSubject, true);
                clpr.ReverseSolution = true;
                clpr.Execute(ClipType.ctUnion, solution, PolyFillType.pftNegative, PolyFillType.pftNegative);
                switch (solution.Count)
                {
                    case > 0:
                        solution.RemoveAt(0);
                        break;
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    public void Execute(ref PolyTree solution, double delta)
    {
        solution.Clear();
        FixOrientations();
        DoOffset(delta);

        //now clean up 'corners' ...
        Clipper clpr = new()
        {
            PreserveCollinear = PreserveCollinear
        };
        clpr.AddPaths(m_destPolys, PolyType.ptSubject, true);
        switch (delta)
        {
            case > 0:
                clpr.Execute(ClipType.ctUnion, solution,
                    PolyFillType.pftPositive, PolyFillType.pftPositive);
                break;
            default:
            {
                IntRect r = ClipperBase.GetBounds(m_destPolys);
                Path outer = new(4)
                {
                    new IntPoint(r.left - 10, r.bottom + 10),
                    new IntPoint(r.right + 10, r.bottom + 10),
                    new IntPoint(r.right + 10, r.top - 10),
                    new IntPoint(r.left - 10, r.top - 10)
                };

                clpr.AddPath(outer, PolyType.ptSubject, true);
                clpr.ReverseSolution = true;
                clpr.Execute(ClipType.ctUnion, solution, PolyFillType.pftNegative, PolyFillType.pftNegative);
                switch (solution.ChildCount)
                {
                    //remove the outer PolyNode rectangle ...
                    case 1 when solution.Childs[0].ChildCount > 0:
                    {
                        PolyNode outerNode = solution.Childs[0];
                        solution.Childs.Capacity = outerNode.ChildCount;
                        solution.Childs[0] = outerNode.Childs[0];
                        solution.Childs[0].m_Parent = solution;
                        for (int i = 1; i < outerNode.ChildCount; i++)
                            solution.AddChild(outerNode.Childs[i]);
                        break;
                    }
                    default:
                        solution.Clear();
                        break;
                }

                break;
            }
        }
    }
    //------------------------------------------------------------------------------

    private void OffsetPoint(int j, ref int k, JoinType jointype)
    {
        //cross product ...
        m_sinA = m_normals[k].X * m_normals[j].Y - m_normals[j].X * m_normals[k].Y;

        switch (Math.Abs(m_sinA * m_delta))
        {
            case < 1.0:
            {
                //dot product ...
                double cosA = m_normals[k].X * m_normals[j].X + m_normals[j].Y * m_normals[k].Y;
                switch (cosA)
                {
                    // angle ==> 0 degrees
                    case > 0:
                        m_destPoly.Add(new IntPoint(Round(m_srcPoly[j].X + m_normals[k].X * m_delta),
                            Round(m_srcPoly[j].Y + m_normals[k].Y * m_delta)));
                        return;
                }

                //else angle ==> 180 degrees   
                break;
            }
            default:
                m_sinA = m_sinA switch
                {
                    > 1.0 => 1.0,
                    < -1.0 => -1.0,
                    _ => m_sinA
                };

                break;
        }

        switch (m_sinA * m_delta)
        {
            case < 0:
                m_destPoly.Add(new IntPoint(Round(m_srcPoly[j].X + m_normals[k].X * m_delta),
                    Round(m_srcPoly[j].Y + m_normals[k].Y * m_delta)));
                m_destPoly.Add(m_srcPoly[j]);
                m_destPoly.Add(new IntPoint(Round(m_srcPoly[j].X + m_normals[j].X * m_delta),
                    Round(m_srcPoly[j].Y + m_normals[j].Y * m_delta)));
                break;
            default:
                switch (jointype)
                {
                    case JoinType.jtMiter:
                    {
                        double r = 1 + (m_normals[j].X * m_normals[k].X +
                                        m_normals[j].Y * m_normals[k].Y);
                        if (r >= m_miterLim)
                        {
                            DoMiter(j, k, r);
                        }
                        else
                        {
                            DoSquare(j, k);
                        }

                        break;
                    }
                    case JoinType.jtSquare:
                        DoSquare(j, k);
                        break;
                    case JoinType.jtRound:
                        DoRound(j, k);
                        break;
                }

                break;
        }

        k = j;
    }
    //------------------------------------------------------------------------------

    internal void DoSquare(int j, int k)
    {
        double dx = Math.Tan(Math.Atan2(m_sinA,
            m_normals[k].X * m_normals[j].X + m_normals[k].Y * m_normals[j].Y) / 4);
        m_destPoly.Add(new IntPoint(
            Round(m_srcPoly[j].X + m_delta * (m_normals[k].X - m_normals[k].Y * dx)),
            Round(m_srcPoly[j].Y + m_delta * (m_normals[k].Y + m_normals[k].X * dx))));
        m_destPoly.Add(new IntPoint(
            Round(m_srcPoly[j].X + m_delta * (m_normals[j].X + m_normals[j].Y * dx)),
            Round(m_srcPoly[j].Y + m_delta * (m_normals[j].Y - m_normals[j].X * dx))));
    }
    //------------------------------------------------------------------------------

    internal void DoMiter(int j, int k, double r)
    {
        double q = m_delta / r;
        m_destPoly.Add(new IntPoint(Round(m_srcPoly[j].X + (m_normals[k].X + m_normals[j].X) * q),
            Round(m_srcPoly[j].Y + (m_normals[k].Y + m_normals[j].Y) * q)));
    }
    //------------------------------------------------------------------------------

    internal void DoRound(int j, int k)
    {
        double a = Math.Atan2(m_sinA,
            m_normals[k].X * m_normals[j].X + m_normals[k].Y * m_normals[j].Y);
        int steps = Math.Max((int) Round(m_StepsPerRad * Math.Abs(a)), 1);

        double X = m_normals[k].X, Y = m_normals[k].Y, X2;
        for (int i = 0; i < steps; ++i)
        {
            m_destPoly.Add(new IntPoint(
                Round(m_srcPoly[j].X + X * m_delta),
                Round(m_srcPoly[j].Y + Y * m_delta)));
            X2 = X;
            X = X * m_cos - m_sin * Y;
            Y = X2 * m_sin + Y * m_cos;
        }

        m_destPoly.Add(new IntPoint(
            Round(m_srcPoly[j].X + m_normals[j].X * m_delta),
            Round(m_srcPoly[j].Y + m_normals[j].Y * m_delta)));
    }
    //------------------------------------------------------------------------------
}

internal class ClipperException : Exception
{
    public ClipperException(string description) : base(description)
    {
    }
}
//------------------------------------------------------------------------------

//end ClipperLib namespace