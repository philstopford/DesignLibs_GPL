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
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
namespace LibTessDotNet.Double;

public struct Vec3
{
    public static readonly Vec3 Zero = new();

    public double X, Y, Z;

    public double this[int index]
    {
        get
        {
            return index switch
            {
                0 => X,
                1 => Y,
                2 => Z,
                _ => throw new IndexOutOfRangeException()
            };
        }
        set
        {
            switch (index)
            {
                case 0:
                    X = value;
                    break;
                case 1:
                    Y = value;
                    break;
                case 2:
                    Z = value;
                    break;
                default:
                    throw new IndexOutOfRangeException();
            }
        }
    }

    public Vec3(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    public static void Sub(ref Vec3 lhs, ref Vec3 rhs, out Vec3 result)
    {
        result.X = lhs.X - rhs.X;
        result.Y = lhs.Y - rhs.Y;
        result.Z = lhs.Z - rhs.Z;
    }

    public static void Neg(ref Vec3 v)
    {
        v.X = -v.X;
        v.Y = -v.Y;
        v.Z = -v.Z;
    }

    public static void Dot(ref Vec3 u, ref Vec3 v, out double dot)
    {
        dot = u.X * v.X + u.Y * v.Y + u.Z * v.Z;
    }

    public static void Normalize(ref Vec3 v)
    {
        double len = v.X * v.X + v.Y * v.Y + v.Z * v.Z;
        Debug.Assert(len >= 0.0f);
        len = 1.0f / Math.Sqrt(len);
        v.X *= len;
        v.Y *= len;
        v.Z *= len;
    }

    public static int LongAxis(ref Vec3 v)
    {
        int i = 0;
        if (Math.Abs(v.Y) > Math.Abs(v.X))
        {
            i = 1;
        }

        if (Math.Abs(v.Z) > Math.Abs(i == 0 ? v.X : v.Y))
        {
            i = 2;
        }

        return i;
    }

    public override string ToString()
    {
        return $"{X}, {Y}, {Z}";
    }
}

public interface ITypePool
{
    object Get();
    void Return(object obj);
}

public class DefaultTypePool<T> : ITypePool where T : class, Pooled<T>, new()
{
    private Queue<T> _pool = new();

    public object Get()
    {
        lock (_pool)
        {
            switch (_pool.Count)
            {
                case > 0:
                    return _pool.Dequeue();
            }
        }
        return new T();
    }

    public void Return(object obj)
    {
        lock (_pool)
        {
#if DEBUG
            if (_pool.Any(other => other == obj))
            {
                throw new InvalidOperationException("object already pooled");
            }
#endif
            _pool.Enqueue(obj as T);
        }
    }
}

public abstract class IPool
{
    public IPool()
    {
        Register<MeshUtils.Vertex>(new DefaultTypePool<MeshUtils.Vertex>());
        Register<MeshUtils.Face>(new DefaultTypePool<MeshUtils.Face>());
        Register<MeshUtils.Edge>(new DefaultTypePool<MeshUtils.Edge>());
        Register<Tess.ActiveRegion>(new DefaultTypePool<Tess.ActiveRegion>());
    }
    public abstract void Register<T>(ITypePool typePool) where T : class, Pooled<T>, new();
    public abstract T Get<T>() where T : class, Pooled<T>, new();
    public abstract void Return<T>(T obj) where T : class, Pooled<T>, new();
}

public class NullPool : IPool
{
    public override T Get<T>()
    {
        T obj = new();
        obj.Init(this);
        return obj;
    }

    public override void Register<T>(ITypePool typePool)
    {
    }

    public override void Return<T>(T obj)
    {
    }
}

public class DefaultPool : IPool
{
    private IDictionary<Type, ITypePool> _register;

    public override void Register<T>(ITypePool typePool)
    {
        _register = _register switch
        {
            null =>
                // can support multiple readers as long as it's not modified
                new Dictionary<Type, ITypePool>(),
            _ => _register
        };

        _register[typeof(T)] = typePool;
    }

    public override T Get<T>()
    {
        T obj = null;
        if (_register.TryGetValue(typeof(T), out ITypePool typePool))
        {
            obj = typePool.Get() as T;
        }

        obj = obj switch
        {
            null => new T(),
            _ => obj
        };
        obj.Init(this);
        return obj;
    }

    public override void Return<T>(T obj)
    {
        switch (obj)
        {
            case null:
                return;
        }
        obj.Reset(this);
        if (_register.TryGetValue(typeof(T), out ITypePool typePool))
        {
            typePool.Return(obj);
        }
    }
}

public interface Pooled<T> where T : class, Pooled<T>, new()
{
    void Init(IPool pool);
    void Reset(IPool pool);
}

internal static class MeshUtils
{
    internal class Vertex : Pooled<Vertex>
    {
        internal Vertex _prev, _next;
        internal Edge _anEdge;

        internal Vec3 _coords;
        internal double _s, _t;
        internal PQHandle _pqHandle;
        internal int _n;
        internal object _data;

        public void Init(IPool pool)
        {
        }

        public void Reset(IPool pool)
        {
            _prev = _next = null;
            _anEdge = null;
            _coords = Vec3.Zero;
            _s = 0;
            _t = 0;
            _pqHandle = new PQHandle();
            _n = 0;
            _data = null;
        }
    }

    internal class Face : Pooled<Face>
    {
        internal Face _prev, _next;
        internal Edge _anEdge;

        internal Face _trail;
        internal int _n;
        internal bool _marked, _inside;

        internal int VertsCount
        {
            get
            {
                int n = 0;
                Edge eCur = _anEdge;
                do
                {
                    n++;
                    eCur = eCur._Lnext;
                } while (eCur != _anEdge);
                return n;
            }
        }

        public void Init(IPool pool)
        {
        }

        public void Reset(IPool pool)
        {
            _prev = _next = null;
            _anEdge = null;
            _trail = null;
            _n = 0;
            _marked = false;
            _inside = false;
        }
    }

    internal struct EdgePair
    {
        internal Edge _e, _eSym;

        public static EdgePair Create(IPool pool)
        {
            Edge e = pool.Get<Edge>();
            Edge eSym = pool.Get<Edge>();

            e._pair._e = e;
            e._pair._eSym = eSym;
            eSym._pair = e._pair;
            return e._pair;
        }

        public void Reset(IPool pool)
        {
            _e = _eSym = null;
        }
    }

    internal class Edge : Pooled<Edge>
    {
        internal EdgePair _pair;
        internal Edge _next, _Sym, _Onext, _Lnext;
        internal Vertex _Org;
        internal Face _Lface;
        internal Tess.ActiveRegion _activeRegion;
        internal int _winding;

        internal Face _Rface { get => _Sym._Lface;
            set => _Sym._Lface = value;
        }
        internal Vertex _Dst { get => _Sym._Org;
            set => _Sym._Org = value;
        }

        internal Edge _Oprev { get => _Sym._Lnext;
            set => _Sym._Lnext = value;
        }
        internal Edge _Lprev { get => _Onext._Sym;
            set => _Onext._Sym = value;
        }
        internal Edge _Dprev { get => _Lnext._Sym;
            set => _Lnext._Sym = value;
        }
        internal Edge _Rprev { get => _Sym._Onext;
            set => _Sym._Onext = value;
        }
        internal Edge _Dnext { get => _Rprev._Sym;
            set => _Rprev._Sym = value;
        }
        internal Edge _Rnext { get => _Oprev._Sym;
            set => _Oprev._Sym = value;
        }

        internal static void EnsureFirst(ref Edge e)
        {
            if (e == e._pair._eSym)
            {
                e = e._Sym;
            }
        }

        public void Init(IPool pool)
        {
        }

        public void Reset(IPool pool)
        {
            _pair.Reset(pool);
            _next = _Sym = _Onext = _Lnext = null;
            _Org = null;
            _Lface = null;
            _activeRegion = null;
            _winding = 0;
        }
    }

    /// <summary>
    /// Splice( a, b ) is best described by the Guibas/Stolfi paper or the
    /// CS348a notes (see Mesh.cs). Basically it modifies the mesh so that
    /// a->Onext and b->Onext are exchanged. This can have various effects
    /// depending on whether a and b belong to different face or vertex rings.
    /// For more explanation see Mesh.Splice().
    /// </summary>
    public static void Splice(Edge a, Edge b)
    {
        Edge aOnext = a._Onext;
        Edge bOnext = b._Onext;

        aOnext._Sym._Lnext = b;
        bOnext._Sym._Lnext = a;
        a._Onext = bOnext;
        b._Onext = aOnext;
    }

    /// <summary>
    /// MakeVertex( eOrig, vNext ) attaches a new vertex and makes it the
    /// origin of all edges in the vertex loop to which eOrig belongs. "vNext" gives
    /// a place to insert the new vertex in the global vertex list. We insert
    /// the new vertex *before* vNext so that algorithms which walk the vertex
    /// list will not see the newly created vertices.
    /// </summary>
    public static void MakeVertex(IPool pool, Edge eOrig, Vertex vNext)
    {
        Vertex vNew = pool.Get<Vertex>();

        // insert in circular doubly-linked list before vNext
        Vertex vPrev = vNext._prev;
        vNew._prev = vPrev;
        vPrev._next = vNew;
        vNew._next = vNext;
        vNext._prev = vNew;

        vNew._anEdge = eOrig;
        // leave coords, s, t undefined

        // fix other edges on this vertex loop
        Edge e = eOrig;
        do
        {
            e._Org = vNew;
            e = e._Onext;
        } while (e != eOrig);
    }

    /// <summary>
    /// MakeFace( eOrig, fNext ) attaches a new face and makes it the left
    /// face of all edges in the face loop to which eOrig belongs. "fNext" gives
    /// a place to insert the new face in the global face list. We insert
    /// the new face *before* fNext so that algorithms which walk the face
    /// list will not see the newly created faces.
    /// </summary>
    public static void MakeFace(IPool pool, Edge eOrig, Face fNext)
    {
        Face fNew = pool.Get<Face>();

        // insert in circular doubly-linked list before fNext
        Face fPrev = fNext._prev;
        fNew._prev = fPrev;
        fPrev._next = fNew;
        fNew._next = fNext;
        fNext._prev = fNew;

        fNew._anEdge = eOrig;
        fNew._trail = null;
        fNew._marked = false;

        // The new face is marked "inside" if the old one was. This is a
        // convenience for the common case where a face has been split in two.
        fNew._inside = fNext._inside;

        // fix other edges on this face loop
        Edge e = eOrig;
        do
        {
            e._Lface = fNew;
            e = e._Lnext;
        } while (e != eOrig);
    }

    /// <summary>
    /// MakeEdge creates a new pair of half-edges which form their own loop.
    /// No vertex or face structures are allocated, but these must be assigned
    /// before the current edge operation is completed.
    /// </summary>
    public static Edge MakeEdge(IPool pool, Edge eNext)
    {
        Debug.Assert(eNext != null);

        EdgePair pair = EdgePair.Create(pool);
        Edge e = pair._e;
        Edge eSym = pair._eSym;

        // Make sure eNext points to the first edge of the edge pair
        Edge.EnsureFirst(ref eNext);

        // Insert in circular doubly-linked list before eNext.
        // Note that the prev pointer is stored in Sym->next.
        Edge ePrev = eNext._Sym._next;
        eSym._next = ePrev;
        ePrev._Sym._next = e;
        e._next = eNext;
        eNext._Sym._next = eSym;

        e._Sym = eSym;
        e._Onext = e;
        e._Lnext = eSym;
        e._Org = null;
        e._Lface = null;
        e._winding = 0;
        e._activeRegion = null;

        eSym._Sym = e;
        eSym._Onext = eSym;
        eSym._Lnext = e;
        eSym._Org = null;
        eSym._Lface = null;
        eSym._winding = 0;
        eSym._activeRegion = null;

        return e;
    }

    /// <summary>
    /// KillEdge( eDel ) destroys an edge (the half-edges eDel and eDel->Sym),
    /// and removes from the global edge list.
    /// </summary>
    public static void KillEdge(IPool pool, Edge eDel)
    {
        // Half-edges are allocated in pairs, see EdgePair above
        Edge.EnsureFirst(ref eDel);

        // delete from circular doubly-linked list
        Edge eNext = eDel._next;
        Edge ePrev = eDel._Sym._next;
        eNext._Sym._next = ePrev;
        ePrev._Sym._next = eNext;

        pool.Return(eDel._Sym);
        pool.Return(eDel);
    }

    /// <summary>
    /// KillVertex( vDel ) destroys a vertex and removes it from the global
    /// vertex list. It updates the vertex loop to point to a given new vertex.
    /// </summary>
    public static void KillVertex(IPool pool, Vertex vDel, Vertex newOrg)
    {
        Edge eStart = vDel._anEdge;

        // change the origin of all affected edges
        Edge e = eStart;
        do
        {
            e._Org = newOrg;
            e = e._Onext;
        } while (e != eStart);

        // delete from circular doubly-linked list
        Vertex vPrev = vDel._prev;
        Vertex vNext = vDel._next;
        vNext._prev = vPrev;
        vPrev._next = vNext;

        pool.Return(vDel);
    }

    /// <summary>
    /// KillFace( fDel ) destroys a face and removes it from the global face
    /// list. It updates the face loop to point to a given new face.
    /// </summary>
    public static void KillFace(IPool pool, Face fDel, Face newLFace)
    {
        Edge eStart = fDel._anEdge;

        // change the left face of all affected edges
        Edge e = eStart;
        do
        {
            e._Lface = newLFace;
            e = e._Lnext;
        } while (e != eStart);

        // delete from circular doubly-linked list
        Face fPrev = fDel._prev;
        Face fNext = fDel._next;
        fNext._prev = fPrev;
        fPrev._next = fNext;

        pool.Return(fDel);
    }

    /// <summary>
    /// Return signed area of face.
    /// </summary>
    public static double FaceArea(Face f)
    {
        double area = 0;
        Edge e = f._anEdge;
        do
        {
            area += (e._Org._s - e._Dst._s) * (e._Org._t + e._Dst._t);
            e = e._Lnext;
        } while (e != f._anEdge);
        return area;
    }
}