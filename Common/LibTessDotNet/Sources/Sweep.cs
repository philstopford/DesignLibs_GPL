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

public partial class Tess
{
    internal class ActiveRegion : Pooled<ActiveRegion>
    {
        internal MeshUtils.Edge _eUp;
        internal Dict<ActiveRegion>.Node _nodeUp;
        internal int _windingNumber;
        internal bool _inside, _sentinel, _dirty, _fixUpperEdge;

        public void Init(IPool pool)
        {
        }

        public void Reset(IPool pool)
        {
            _eUp = null;
            _nodeUp = null;
            _windingNumber = 0;
            _inside = false;
            _sentinel = false;
            _dirty = false;
            _fixUpperEdge = false;
        }
    }

    private static ActiveRegion RegionBelow(ActiveRegion reg)
    {
        return reg._nodeUp._prev._key;
    }

    private static ActiveRegion RegionAbove(ActiveRegion reg)
    {
        return reg._nodeUp._next._key;
    }

    /// <summary>
    /// Both edges must be directed from right to left (this is the canonical
    /// direction for the upper edge of each region).
    /// 
    /// The strategy is to evaluate a "t" value for each edge at the
    /// current sweep line position, given by tess->event. The calculations
    /// are designed to be very stable, but of course they are not perfect.
    /// 
    /// Special case: if both edge destinations are at the sweep event,
    /// we sort the edges by slope (they would otherwise compare equally).
    /// </summary>
    private bool EdgeLeq(ActiveRegion reg1, ActiveRegion reg2)
    {
        MeshUtils.Edge e1 = reg1._eUp;
        MeshUtils.Edge e2 = reg2._eUp;

        if (e1._Dst == _event)
        {
            if (e2._Dst != _event)
            {
                return Geom.EdgeSign(e2._Dst, _event, e2._Org) <= 0.0f;
            }

            // Two edges right of the sweep line which meet at the sweep event.
            // Sort them by slope.
            if (Geom.VertLeq(e1._Org, e2._Org))
            {
                return Geom.EdgeSign(e2._Dst, e1._Org, e2._Org) <= 0.0f;
            }
            return Geom.EdgeSign(e1._Dst, e2._Org, e1._Org) >= 0.0f;
        }
        if (e2._Dst == _event)
        {
            return Geom.EdgeSign(e1._Dst, _event, e1._Org) >= 0.0f;
        }

        // General case - compute signed distance *from* e1, e2 to event
        double t1 = Geom.EdgeEval(e1._Dst, _event, e1._Org);
        double t2 = Geom.EdgeEval(e2._Dst, _event, e2._Org);
        return t1 >= t2;
    }

    private void DeleteRegion(ActiveRegion reg)
    {
        switch (reg._fixUpperEdge)
        {
            case true:
                // It was created with zero winding number, so it better be
                // deleted with zero winding number (ie. it better not get merged
                // with a real edge).
                Debug.Assert(reg._eUp._winding == 0);
                break;
        }
        reg._eUp._activeRegion = null;
        _dict.Remove(reg._nodeUp);
        _pool.Return(reg);
    }

    /// <summary>
    /// Replace an upper edge which needs fixing (see ConnectRightVertex).
    /// </summary>
    private void FixUpperEdge(ActiveRegion reg, MeshUtils.Edge newEdge)
    {
        Debug.Assert(reg._fixUpperEdge);
        Mesh.Delete(_pool, reg._eUp);
        reg._fixUpperEdge = false;
        reg._eUp = newEdge;
        newEdge._activeRegion = reg;
    }

    private ActiveRegion TopLeftRegion(ActiveRegion reg)
    {
        MeshUtils.Vertex org = reg._eUp._Org;

        // Find the region above the uppermost edge with the same origin
        do
        {
            reg = RegionAbove(reg);
        } while (reg._eUp._Org == org);

        switch (reg._fixUpperEdge)
        {
            // If the edge above was a temporary edge introduced by ConnectRightVertex,
            // now is the time to fix it.
            case true:
            {
                MeshUtils.Edge e = Mesh.Connect(_pool, RegionBelow(reg)._eUp._Sym, reg._eUp._Lnext);
                FixUpperEdge(reg, e);
                reg = RegionAbove(reg);
                break;
            }
        }

        return reg;
    }

    private ActiveRegion TopRightRegion(ActiveRegion reg)
    {
        MeshUtils.Vertex dst = reg._eUp._Dst;

        // Find the region above the uppermost edge with the same destination
        do
        {
            reg = RegionAbove(reg);
        } while (reg._eUp._Dst == dst);

        return reg;
    }

    /// <summary>
    /// Add a new active region to the sweep line, *somewhere* below "regAbove"
    /// (according to where the new edge belongs in the sweep-line dictionary).
    /// The upper edge of the new region will be "eNewUp".
    /// Winding number and "inside" flag are not updated.
    /// </summary>
    private ActiveRegion AddRegionBelow(ActiveRegion regAbove, MeshUtils.Edge eNewUp)
    {
        ActiveRegion regNew = _pool.Get<ActiveRegion>();

        regNew._eUp = eNewUp;
        regNew._nodeUp = _dict.InsertBefore(regAbove._nodeUp, regNew);
        regNew._fixUpperEdge = false;
        regNew._sentinel = false;
        regNew._dirty = false;

        eNewUp._activeRegion = regNew;

        return regNew;
    }

    private void ComputeWinding(ActiveRegion reg)
    {
        reg._windingNumber = RegionAbove(reg)._windingNumber + reg._eUp._winding;
        reg._inside = Geom.IsWindingInside(_windingRule, reg._windingNumber);
    }

    /// <summary>
    /// Delete a region from the sweep line. This happens when the upper
    /// and lower chains of a region meet (at a vertex on the sweep line).
    /// The "inside" flag is copied to the appropriate mesh face (we could
    /// not do this before -- since the structure of the mesh is always
    /// changing, this face may not have even existed until now).
    /// </summary>
    private void FinishRegion(ActiveRegion reg)
    {
        MeshUtils.Edge e = reg._eUp;
        MeshUtils.Face f = e._Lface;

        f._inside = reg._inside;
        f._anEdge = e;
        DeleteRegion(reg);
    }

    /// <summary>
    /// We are given a vertex with one or more left-going edges.  All affected
    /// edges should be in the edge dictionary.  Starting at regFirst->eUp,
    /// we walk down deleting all regions where both edges have the same
    /// origin vOrg.  At the same time we copy the "inside" flag from the
    /// active region to the face, since at this point each face will belong
    /// to at most one region (this was not necessarily true until this point
    /// in the sweep).  The walk stops at the region above regLast; if regLast
    /// is null we walk as far as possible.  At the same time we relink the
    /// mesh if necessary, so that the ordering of edges around vOrg is the
    /// same as in the dictionary.
    /// </summary>
    private MeshUtils.Edge FinishLeftRegions(ActiveRegion regFirst, ActiveRegion regLast)
    {
        ActiveRegion regPrev = regFirst;
        MeshUtils.Edge ePrev = regFirst._eUp;

        while (regPrev != regLast)
        {
            regPrev._fixUpperEdge = false;	// placement was OK
            ActiveRegion reg = RegionBelow(regPrev);
            MeshUtils.Edge e = reg._eUp;
            if (e._Org != ePrev._Org)
            {
                if (!reg._fixUpperEdge)
                {
                    // Remove the last left-going edge.  Even though there are no further
                    // edges in the dictionary with this origin, there may be further
                    // such edges in the mesh (if we are adding left edges to a vertex
                    // that has already been processed).  Thus it is important to call
                    // FinishRegion rather than just DeleteRegion.
                    FinishRegion(regPrev);
                    break;
                }
                // If the edge below was a temporary edge introduced by
                // ConnectRightVertex, now is the time to fix it.
                e = Mesh.Connect(_pool, ePrev._Lprev, e._Sym);
                FixUpperEdge(reg, e);
            }

            // Relink edges so that ePrev.Onext == e
            if (ePrev._Onext != e)
            {
                Mesh.Splice(_pool, e._Oprev, e);
                Mesh.Splice(_pool, ePrev, e);
            }
            FinishRegion(regPrev); // may change reg.eUp
            ePrev = reg._eUp;
            regPrev = reg;
        }

        return ePrev;
    }

    /// <summary>
    /// Purpose: insert right-going edges into the edge dictionary, and update
    /// winding numbers and mesh connectivity appropriately.  All right-going
    /// edges share a common origin vOrg.  Edges are inserted CCW starting at
    /// eFirst; the last edge inserted is eLast.Oprev.  If vOrg has any
    /// left-going edges already processed, then eTopLeft must be the edge
    /// such that an imaginary upward vertical segment from vOrg would be
    /// contained between eTopLeft.Oprev and eTopLeft; otherwise eTopLeft
    /// should be null.
    /// </summary>
    private void AddRightEdges(ActiveRegion regUp, MeshUtils.Edge eFirst, MeshUtils.Edge eLast, MeshUtils.Edge eTopLeft, bool cleanUp)
    {
        bool firstTime = true;

        MeshUtils.Edge e = eFirst; do
        {
            Debug.Assert(Geom.VertLeq(e._Org, e._Dst));
            AddRegionBelow(regUp, e._Sym);
            e = e._Onext;
        } while (e != eLast);

        eTopLeft = eTopLeft switch
        {
            // Walk *all* right-going edges from e.Org, in the dictionary order,
            // updating the winding numbers of each region, and re-linking the mesh
            // edges to match the dictionary ordering (if necessary).
            null => RegionBelow(regUp)._eUp._Rprev,
            _ => eTopLeft
        };

        ActiveRegion regPrev = regUp, reg;
        MeshUtils.Edge ePrev = eTopLeft;
        while (true)
        {
            reg = RegionBelow(regPrev);
            e = reg._eUp._Sym;
            if (e._Org != ePrev._Org)
            {
                break;
            }

            if (e._Onext != ePrev)
            {
                // Unlink e from its current position, and relink below ePrev
                Mesh.Splice(_pool, e._Oprev, e);
                Mesh.Splice(_pool, ePrev._Oprev, e);
            }
            // Compute the winding number and "inside" flag for the new regions
            reg._windingNumber = regPrev._windingNumber - e._winding;
            reg._inside = Geom.IsWindingInside(_windingRule, reg._windingNumber);

            // Check for two outgoing edges with same slope -- process these
            // before any intersection tests (see example in tessComputeInterior).
            regPrev._dirty = true;
            switch (firstTime)
            {
                case false when CheckForRightSplice(regPrev):
                    Geom.AddWinding(e, ePrev);
                    DeleteRegion(regPrev);
                    Mesh.Delete(_pool, ePrev);
                    break;
            }
            firstTime = false;
            regPrev = reg;
            ePrev = e;
        }
        regPrev._dirty = true;
        Debug.Assert(regPrev._windingNumber - e._winding == reg._windingNumber);

        switch (cleanUp)
        {
            case true:
                // Check for intersections between newly adjacent edges.
                WalkDirtyRegions(regPrev);
                break;
        }
    }

    /// <summary>
    /// Two vertices with identical coordinates are combined into one.
    /// e1.Org is kept, while e2.Org is discarded.
    /// </summary>
    private void SpliceMergeVertices(MeshUtils.Edge e1, MeshUtils.Edge e2)
    {
        Mesh.Splice(_pool, e1, e2);
    }

    /// <summary>
    /// Find some weights which describe how the intersection vertex is
    /// a linear combination of "org" and "dest".  Each of the two edges
    /// which generated "isect" is allocated 50% of the weight; each edge
    /// splits the weight between its org and dst according to the
    /// relative distance to "isect".
    /// </summary>
    private static void VertexWeights(MeshUtils.Vertex isect, MeshUtils.Vertex org, MeshUtils.Vertex dst, out double w0, out double w1)
    {
        double t1 = Geom.VertL1dist(org, isect);
        double t2 = Geom.VertL1dist(dst, isect);

        w0 = t2 / (t1 + t2) / 2.0f;
        w1 = t1 / (t1 + t2) / 2.0f;

        isect._coords.X += w0 * org._coords.X + w1 * dst._coords.X;
        isect._coords.Y += w0 * org._coords.Y + w1 * dst._coords.Y;
        isect._coords.Z += w0 * org._coords.Z + w1 * dst._coords.Z;
    }

    /// <summary>
    /// We've computed a new intersection point, now we need a "data" pointer
    /// from the user so that we can refer to this new vertex in the
    /// rendering callbacks.
    /// </summary>
    private void GetIntersectData(MeshUtils.Vertex isect, MeshUtils.Vertex orgUp, MeshUtils.Vertex dstUp, MeshUtils.Vertex orgLo, MeshUtils.Vertex dstLo)
    {
        isect._coords = Vec3.Zero;
        VertexWeights(isect, orgUp, dstUp, out double w0, out double w1);
        VertexWeights(isect, orgLo, dstLo, out double w2, out double w3);

        if (_combineCallback != null)
        {
            isect._data = _combineCallback(
                isect._coords,
                [orgUp._data, dstUp._data, orgLo._data, dstLo._data],
                [w0, w1, w2, w3]
            );
        }
    }

    /// <summary>
    /// Check the upper and lower edge of "regUp", to make sure that the
    /// eUp->Org is above eLo, or eLo->Org is below eUp (depending on which
    /// origin is leftmost).
    /// 
    /// The main purpose is to splice right-going edges with the same
    /// dest vertex and nearly identical slopes (ie. we can't distinguish
    /// the slopes numerically).  However the splicing can also help us
    /// to recover from numerical errors.  For example, suppose at one
    /// point we checked eUp and eLo, and decided that eUp->Org is barely
    /// above eLo.  Then later, we split eLo into two edges (eg. from
    /// a splice operation like this one).  This can change the result of
    /// our test so that now eUp->Org is incident to eLo, or barely below it.
    /// We must correct this condition to maintain the dictionary invariants.
    /// 
    /// One possibility is to check these edges for intersection again
    /// (ie. CheckForIntersect).  This is what we do if possible.  However
    /// CheckForIntersect requires that tess->event lies between eUp and eLo,
    /// so that it has something to fall back on when the intersection
    /// calculation gives us an unusable answer.  So, for those cases where
    /// we can't check for intersection, this routine fixes the problem
    /// by just splicing the offending vertex into the other edge.
    /// This is a guaranteed solution, no matter how degenerate things get.
    /// Basically this is a combinatorial solution to a numerical problem.
    /// </summary>
    private bool CheckForRightSplice(ActiveRegion regUp)
    {
        ActiveRegion regLo = RegionBelow(regUp);
        MeshUtils.Edge eUp = regUp._eUp;
        MeshUtils.Edge eLo = regLo._eUp;

        if (Geom.VertLeq(eUp._Org, eLo._Org))
        {
            if (Geom.EdgeSign(eLo._Dst, eUp._Org, eLo._Org) > 0.0f)
            {
                return false;
            }

            // eUp.Org appears to be below eLo
            if (!Geom.VertEq(eUp._Org, eLo._Org))
            {
                // Splice eUp._Org into eLo
                _mesh.SplitEdge(_pool, eLo._Sym);
                Mesh.Splice(_pool, eUp, eLo._Oprev);
                regUp._dirty = regLo._dirty = true;
            }
            else if (eUp._Org != eLo._Org)
            {
                // merge the two vertices, discarding eUp.Org
                _pq.Remove(eUp._Org._pqHandle);
                SpliceMergeVertices(eLo._Oprev, eUp);
            }
        }
        else
        {
            if (Geom.EdgeSign(eUp._Dst, eLo._Org, eUp._Org) < 0.0f)
            {
                return false;
            }

            // eLo.Org appears to be above eUp, so splice eLo.Org into eUp
            RegionAbove(regUp)._dirty = regUp._dirty = true;
            _mesh.SplitEdge(_pool, eUp._Sym);
            Mesh.Splice(_pool, eLo._Oprev, eUp);
        }
        return true;
    }

    /// <summary>
    /// Check the upper and lower edge of "regUp", to make sure that the
    /// eUp->Dst is above eLo, or eLo->Dst is below eUp (depending on which
    /// destination is rightmost).
    /// 
    /// Theoretically, this should always be true.  However, splitting an edge
    /// into two pieces can change the results of previous tests.  For example,
    /// suppose at one point we checked eUp and eLo, and decided that eUp->Dst
    /// is barely above eLo.  Then later, we split eLo into two edges (eg. from
    /// a splice operation like this one).  This can change the result of
    /// the test so that now eUp->Dst is incident to eLo, or barely below it.
    /// We must correct this condition to maintain the dictionary invariants
    /// (otherwise new edges might get inserted in the wrong place in the
    /// dictionary, and bad stuff will happen).
    /// 
    /// We fix the problem by just splicing the offending vertex into the
    /// other edge.
    /// </summary>
    private bool CheckForLeftSplice(ActiveRegion regUp)
    {
        ActiveRegion regLo = RegionBelow(regUp);
        MeshUtils.Edge eUp = regUp._eUp;
        MeshUtils.Edge eLo = regLo._eUp;

        Debug.Assert(!Geom.VertEq(eUp._Dst, eLo._Dst));

        if (Geom.VertLeq(eUp._Dst, eLo._Dst))
        {
            if (Geom.EdgeSign(eUp._Dst, eLo._Dst, eUp._Org) < 0.0f)
            {
                return false;
            }

            // eLo.Dst is above eUp, so splice eLo.Dst into eUp
            RegionAbove(regUp)._dirty = regUp._dirty = true;
            MeshUtils.Edge e = _mesh.SplitEdge(_pool, eUp);
            Mesh.Splice(_pool, eLo._Sym, e);
            e._Lface._inside = regUp._inside;
        }
        else
        {
            if (Geom.EdgeSign(eLo._Dst, eUp._Dst, eLo._Org) > 0.0f)
            {
                return false;
            }

            // eUp.Dst is below eLo, so splice eUp.Dst into eLo
            regUp._dirty = regLo._dirty = true;
            MeshUtils.Edge e = _mesh.SplitEdge(_pool, eLo);
            Mesh.Splice(_pool, eUp._Lnext, eLo._Sym);
            e._Rface._inside = regUp._inside;
        }
        return true;
    }

    /// <summary>
    /// Check the upper and lower edges of the given region to see if
    /// they intersect.  If so, create the intersection and add it
    /// to the data structures.
    /// 
    /// Returns TRUE if adding the new intersection resulted in a recursive
    /// call to AddRightEdges(); in this case all "dirty" regions have been
    /// checked for intersections, and possibly regUp has been deleted.
    /// </summary>
    private bool CheckForIntersect(ActiveRegion regUp)
    {
        ActiveRegion regLo = RegionBelow(regUp);
        MeshUtils.Edge eUp = regUp._eUp;
        MeshUtils.Edge eLo = regLo._eUp;
        MeshUtils.Vertex orgUp = eUp._Org;
        MeshUtils.Vertex orgLo = eLo._Org;
        MeshUtils.Vertex dstUp = eUp._Dst;
        MeshUtils.Vertex dstLo = eLo._Dst;

        Debug.Assert(!Geom.VertEq(dstLo, dstUp));
        Debug.Assert(Geom.EdgeSign(dstUp, _event, orgUp) <= 0.0f);
        Debug.Assert(Geom.EdgeSign(dstLo, _event, orgLo) >= 0.0f);
        Debug.Assert(orgUp != _event && orgLo != _event);
        Debug.Assert(!regUp._fixUpperEdge && !regLo._fixUpperEdge);

        if (orgUp == orgLo)
        {
            // right endpoints are the same
            return false;
        }

        double tMinUp = Math.Min(orgUp._t, dstUp._t);
        double tMaxLo = Math.Max(orgLo._t, dstLo._t);
        if (tMinUp > tMaxLo)
        {
            // t ranges do not overlap
            return false;
        }

        if (Geom.VertLeq(orgUp, orgLo))
        {
            if (Geom.EdgeSign(dstLo, orgUp, orgLo) > 0.0f)
            {
                return false;
            }
        }
        else
        {
            if (Geom.EdgeSign(dstUp, orgLo, orgUp) < 0.0f)
            {
                return false;
            }
        }

        // At this point the edges intersect, at least marginally

        MeshUtils.Vertex isect = _pool.Get<MeshUtils.Vertex>();
        Geom.EdgeIntersect(dstUp, orgUp, dstLo, orgLo, isect);
        // The following properties are guaranteed:
        Debug.Assert(Math.Min(orgUp._t, dstUp._t) <= isect._t);
        Debug.Assert(isect._t <= Math.Max(orgLo._t, dstLo._t));
        Debug.Assert(Math.Min(dstLo._s, dstUp._s) <= isect._s);
        Debug.Assert(isect._s <= Math.Max(orgLo._s, orgUp._s));

        if (Geom.VertLeq(isect, _event))
        {
            // The intersection point lies slightly to the left of the sweep line,
            // so move it until it's slightly to the right of the sweep line.
            // (If we had perfect numerical precision, this would never happen
            // in the first place). The easiest and safest thing to do is
            // replace the intersection by tess._event.
            isect._s = _event._s;
            isect._t = _event._t;
        }
        // Similarly, if the computed intersection lies to the right of the
        // rightmost origin (which should rarely happen), it can cause
        // unbelievable inefficiency on sufficiently degenerate inputs.
        // (If you have the test program, try running test54.d with the
        // "X zoom" option turned on).
        MeshUtils.Vertex orgMin = Geom.VertLeq(orgUp, orgLo) ? orgUp : orgLo;
        if (Geom.VertLeq(orgMin, isect))
        {
            isect._s = orgMin._s;
            isect._t = orgMin._t;
        }

        if (Geom.VertEq(isect, orgUp) || Geom.VertEq(isect, orgLo))
        {
            // Easy case -- intersection at one of the right endpoints
            CheckForRightSplice(regUp);
            _pool.Return(isect);
            return false;
        }

        if (!Geom.VertEq(dstUp, _event)
            && Geom.EdgeSign(dstUp, _event, isect) >= 0.0f
            || !Geom.VertEq(dstLo, _event)
            && Geom.EdgeSign(dstLo, _event, isect) <= 0.0f)
        {
            // Very unusual -- the new upper or lower edge would pass on the
            // wrong side of the sweep event, or through it. This can happen
            // due to very small numerical errors in the intersection calculation.
            if (dstLo == _event)
            {
                // Splice dstLo into eUp, and process the new region(s)
                _mesh.SplitEdge(_pool, eUp._Sym);
                Mesh.Splice(_pool, eLo._Sym, eUp);
                regUp = TopLeftRegion(regUp);
                eUp = RegionBelow(regUp)._eUp;
                FinishLeftRegions(RegionBelow(regUp), regLo);
                AddRightEdges(regUp, eUp._Oprev, eUp, eUp, true);
                _pool.Return(isect);
                return true;
            }
            if (dstUp == _event)
            {
                /* Splice dstUp into eLo, and process the new region(s) */
                _mesh.SplitEdge(_pool, eLo._Sym);
                Mesh.Splice(_pool, eUp._Lnext, eLo._Oprev);
                regLo = regUp;
                regUp = TopRightRegion(regUp);
                MeshUtils.Edge e = RegionBelow(regUp)._eUp._Rprev;
                regLo._eUp = eLo._Oprev;
                eLo = FinishLeftRegions(regLo, null);
                AddRightEdges(regUp, eLo._Onext, eUp._Rprev, e, true);
                _pool.Return(isect);
                return true;
            }
            // Special case: called from ConnectRightVertex. If either
            // edge passes on the wrong side of tess._event, split it
            // (and wait for ConnectRightVertex to splice it appropriately).
            if (Geom.EdgeSign(dstUp, _event, isect) >= 0.0f)
            {
                RegionAbove(regUp)._dirty = regUp._dirty = true;
                _mesh.SplitEdge(_pool, eUp._Sym);
                eUp._Org._s = _event._s;
                eUp._Org._t = _event._t;
            }
            if (Geom.EdgeSign(dstLo, _event, isect) <= 0.0f)
            {
                regUp._dirty = regLo._dirty = true;
                _mesh.SplitEdge(_pool, eLo._Sym);
                eLo._Org._s = _event._s;
                eLo._Org._t = _event._t;
            }
            // leave the rest for ConnectRightVertex
            _pool.Return(isect);
            return false;
        }

        // General case -- split both edges, splice into new vertex.
        // When we do the splice operation, the order of the arguments is
        // arbitrary as far as correctness goes. However, when the operation
        // creates a new face, the work done is proportional to the size of
        // the new face.  We expect the faces in the processed part of
        // the mesh (ie. eUp._Lface) to be smaller than the faces in the
        // unprocessed original contours (which will be eLo._Oprev._Lface).
        _mesh.SplitEdge(_pool, eUp._Sym);
        _mesh.SplitEdge(_pool, eLo._Sym);
        Mesh.Splice(_pool, eLo._Oprev, eUp);
        eUp._Org._s = isect._s;
        eUp._Org._t = isect._t;
        _pool.Return(isect);
        isect = null;
        eUp._Org._pqHandle = _pq.Insert(eUp._Org);
        if (eUp._Org._pqHandle._handle == PQHandle.Invalid)
        {
            throw new InvalidOperationException("PQHandle should not be invalid");
        }
        GetIntersectData(eUp._Org, orgUp, dstUp, orgLo, dstLo);
        RegionAbove(regUp)._dirty = regUp._dirty = regLo._dirty = true;
        return false;
    }

    /// <summary>
    /// When the upper or lower edge of any region changes, the region is
    /// marked "dirty".  This routine walks through all the dirty regions
    /// and makes sure that the dictionary invariants are satisfied
    /// (see the comments at the beginning of this file).  Of course
    /// new dirty regions can be created as we make changes to restore
    /// the invariants.
    /// </summary>
    private void WalkDirtyRegions(ActiveRegion regUp)
    {
        ActiveRegion regLo = RegionBelow(regUp);

        while (true)
        {
            // Find the lowest dirty region (we walk from the bottom up).
            while (regLo._dirty)
            {
                regUp = regLo;
                regLo = RegionBelow(regLo);
            }
            switch (regUp._dirty)
            {
                case false:
                {
                    regLo = regUp;
                    regUp = RegionAbove(regUp);
                    if (regUp is not {_dirty: true})
                    {
                        // We've walked all the dirty regions
                        return;
                    }

                    break;
                }
            }
            regUp._dirty = false;
            MeshUtils.Edge eUp = regUp._eUp;
            MeshUtils.Edge eLo = regLo._eUp;

            if (eUp._Dst != eLo._Dst)
            {
                // Check that the edge ordering is obeyed at the Dst vertices.
                if (CheckForLeftSplice(regUp))
                {
                    switch (regLo._fixUpperEdge)
                    {
                        // If the upper or lower edge was marked fixUpperEdge, then
                        // we no longer need it (since these edges are needed only for
                        // vertices which otherwise have no right-going edges).
                        case true:
                            DeleteRegion(regLo);
                            Mesh.Delete(_pool, eLo);
                            regLo = RegionBelow(regUp);
                            eLo = regLo._eUp;
                            break;
                        default:
                        {
                            switch (regUp._fixUpperEdge)
                            {
                                case true:
                                    DeleteRegion(regUp);
                                    Mesh.Delete(_pool, eUp);
                                    regUp = RegionAbove(regLo);
                                    eUp = regUp._eUp;
                                    break;
                            }

                            break;
                        }
                    }
                }
            }
            if (eUp._Org != eLo._Org)
            {
                if (eUp._Dst != eLo._Dst
                    && !regUp._fixUpperEdge && !regLo._fixUpperEdge
                    && (eUp._Dst == _event || eLo._Dst == _event))
                {
                    // When all else fails in CheckForIntersect(), it uses tess._event
                    // as the intersection location. To make this possible, it requires
                    // that tess._event lie between the upper and lower edges, and also
                    // that neither of these is marked fixUpperEdge (since in the worst
                    // case it might splice one of these edges into tess.event, and
                    // violate the invariant that fixable edges are the only right-going
                    // edge from their associated vertex).
                    if (CheckForIntersect(regUp))
                    {
                        // WalkDirtyRegions() was called recursively; we're done
                        return;
                    }
                }
                else
                {
                    // Even though we can't use CheckForIntersect(), the Org vertices
                    // may violate the dictionary edge ordering. Check and correct this.
                    CheckForRightSplice(regUp);
                }
            }

            if (eUp._Org != eLo._Org || eUp._Dst != eLo._Dst)
            {
                continue;
            }

            // A degenerate loop consisting of only two edges -- delete it.
            Geom.AddWinding(eLo, eUp);
            DeleteRegion(regUp);
            Mesh.Delete(_pool, eUp);
            regUp = RegionAbove(regLo);
        }
    }

    /// <summary>
    /// Purpose: connect a "right" vertex vEvent (one where all edges go left)
    /// to the unprocessed portion of the mesh.  Since there are no right-going
    /// edges, two regions (one above vEvent and one below) are being merged
    /// into one.  "regUp" is the upper of these two regions.
    /// 
    /// There are two reasons for doing this (adding a right-going edge):
    ///  - if the two regions being merged are "inside", we must add an edge
    ///    to keep them separated (the combined region would not be monotone).
    ///  - in any case, we must leave some record of vEvent in the dictionary,
    ///    so that we can merge vEvent with features that we have not seen yet.
    ///    For example, maybe there is a vertical edge which passes just to
    ///    the right of vEvent; we would like to splice vEvent into this edge.
    /// 
    /// However, we don't want to connect vEvent to just any vertex.  We don't
    /// want the new edge to cross any other edges; otherwise we will create
    /// intersection vertices even when the input data had no self-intersections.
    /// (This is a bad thing; if the user's input data has no intersections,
    /// we don't want to generate any false intersections ourselves.)
    /// 
    /// Our eventual goal is to connect vEvent to the leftmost unprocessed
    /// vertex of the combined region (the union of regUp and regLo).
    /// But because of unseen vertices with all right-going edges, and also
    /// new vertices which may be created by edge intersections, we don't
    /// know where that leftmost unprocessed vertex is.  In the meantime, we
    /// connect vEvent to the closest vertex of either chain, and mark the region
    /// as "fixUpperEdge".  This flag says to delete and reconnect this edge
    /// to the next processed vertex on the boundary of the combined region.
    /// Quite possibly the vertex we connected to will turn out to be the
    /// closest one, in which case we won't need to make any changes.
    /// </summary>
    private void ConnectRightVertex(ActiveRegion regUp, MeshUtils.Edge eBottomLeft)
    {
        MeshUtils.Edge eTopLeft = eBottomLeft._Onext;
        ActiveRegion regLo = RegionBelow(regUp);
        MeshUtils.Edge eUp = regUp._eUp;
        MeshUtils.Edge eLo = regLo._eUp;
        bool degenerate = false;

        if (eUp._Dst != eLo._Dst)
        {
            CheckForIntersect(regUp);
        }

        // Possible new degeneracies: upper or lower edge of regUp may pass
        // through vEvent, or may coincide with new intersection vertex
        if (Geom.VertEq(eUp._Org, _event))
        {
            Mesh.Splice(_pool, eTopLeft._Oprev, eUp);
            regUp = TopLeftRegion(regUp);
            eTopLeft = RegionBelow(regUp)._eUp;
            FinishLeftRegions(RegionBelow(regUp), regLo);
            degenerate = true;
        }
        if (Geom.VertEq(eLo._Org, _event))
        {
            Mesh.Splice(_pool, eBottomLeft, eLo._Oprev);
            eBottomLeft = FinishLeftRegions(regLo, null);
            degenerate = true;
        }
        switch (degenerate)
        {
            case true:
                AddRightEdges(regUp, eBottomLeft._Onext, eTopLeft, eTopLeft, true);
                return;
        }

        // Non-degenerate situation -- need to add a temporary, fixable edge.
        // Connect to the closer of eLo.Org, eUp.Org.
        MeshUtils.Edge eNew = Geom.VertLeq(eLo._Org, eUp._Org) ? eLo._Oprev : eUp;
        eNew = Mesh.Connect(_pool, eBottomLeft._Lprev, eNew);

        // Prevent cleanup, otherwise eNew might disappear before we've even
        // had a chance to mark it as a temporary edge.
        AddRightEdges(regUp, eNew, eNew._Onext, eNew._Onext, false);
        eNew._Sym._activeRegion._fixUpperEdge = true;
        WalkDirtyRegions(regUp);
    }

    /// <summary>
    /// The event vertex lies exactly on an already-processed edge or vertex.
    /// Adding the new vertex involves splicing it into the already-processed
    /// part of the mesh.
    /// </summary>
    private void ConnectLeftDegenerate(ActiveRegion regUp, MeshUtils.Vertex vEvent)
    {
        MeshUtils.Edge e = regUp._eUp;
        if (Geom.VertEq(e._Org, vEvent))
        {
            // e.Org is an unprocessed vertex - just combine them, and wait
            // for e.Org to be pulled from the queue
            // C# : in the C version, there is a flag but it was never implemented
            // the vertices are before beginning the tessellation
            throw new InvalidOperationException("Vertices should have been merged before");
        }

        if (Geom.VertEq(e._Dst, vEvent))
        {
            throw new InvalidOperationException("Vertices should have been merged before");
        }

        // General case -- splice vEvent into edge e which passes through it
        _mesh.SplitEdge(_pool, e._Sym);
        switch (regUp._fixUpperEdge)
        {
            case true:
                // This edge was fixable -- delete unused portion of original edge
                Mesh.Delete(_pool, e._Onext);
                regUp._fixUpperEdge = false;
                break;
        }
        Mesh.Splice(_pool, vEvent._anEdge, e);
        SweepEvent(vEvent);	// recurse

        // See above
    }

    /// <summary>
    /// Purpose: connect a "left" vertex (one where both edges go right)
    /// to the processed portion of the mesh.  Let R be the active region
    /// containing vEvent, and let U and L be the upper and lower edge
    /// chains of R.  There are two possibilities:
    /// 
    /// - the normal case: split R into two regions, by connecting vEvent to
    ///   the rightmost vertex of U or L lying to the left of the sweep line
    /// 
    /// - the degenerate case: if vEvent is close enough to U or L, we
    ///   merge vEvent into that edge chain.  The subcases are:
    ///     - merging with the rightmost vertex of U or L
    ///     - merging with the active edge of U or L
    ///     - merging with an already-processed portion of U or L
    /// </summary>
    private void ConnectLeftVertex(MeshUtils.Vertex vEvent)
    {
        ActiveRegion tmp = _pool.Get<ActiveRegion>();

        // Get a pointer to the active region containing vEvent
        tmp._eUp = vEvent._anEdge._Sym;
        ActiveRegion regUp = _dict.Find(tmp).Key;
        _pool.Return(tmp);
        ActiveRegion regLo = RegionBelow(regUp);
        switch (regLo)
        {
            case null:
                // This may happen if the input polygon is coplanar.
                return;
        }
        MeshUtils.Edge eUp = regUp._eUp;
        MeshUtils.Edge eLo = regLo._eUp;

        // Try merging with U or L first
        if (Geom.EdgeSign(eUp._Dst, vEvent, eUp._Org) == 0.0f)
        {
            ConnectLeftDegenerate(regUp, vEvent);
            return;
        }

        // Connect vEvent to rightmost processed vertex of either chain.
        // e._Dst is the vertex that we will connect to vEvent.
        ActiveRegion reg = Geom.VertLeq(eLo._Dst, eUp._Dst) ? regUp : regLo;

        if (regUp._inside || reg._fixUpperEdge)
        {
            MeshUtils.Edge eNew = reg == regUp ? Mesh.Connect(_pool, vEvent._anEdge._Sym, eUp._Lnext) : Mesh.Connect(_pool, eLo._Dnext, vEvent._anEdge)._Sym;
            switch (reg._fixUpperEdge)
            {
                case true:
                    FixUpperEdge(reg, eNew);
                    break;
                default:
                    ComputeWinding(AddRegionBelow(regUp, eNew));
                    break;
            }
            SweepEvent(vEvent);
        }
        else
        {
            // The new vertex is in a region which does not belong to the polygon.
            // We don't need to connect this vertex to the rest of the mesh.
            AddRightEdges(regUp, vEvent._anEdge, vEvent._anEdge, null, true);
        }
    }

    /// <summary>
    /// Does everything necessary when the sweep line crosses a vertex.
    /// Updates the mesh and the edge dictionary.
    /// </summary>
    private void SweepEvent(MeshUtils.Vertex vEvent)
    {
        _event = vEvent;

        // Check if this vertex is the right endpoint of an edge that is
        // already in the dictionary. In this case we don't need to waste
        // time searching for the location to insert new edges.
        MeshUtils.Edge e = vEvent._anEdge;
        while (e._activeRegion == null)
        {
            e = e._Onext;
            if (e != vEvent._anEdge)
            {
                continue;
            }

            // All edges go right -- not incident to any processed edges
            ConnectLeftVertex(vEvent);
            return;
        }

        // Processing consists of two phases: first we "finish" all the
        // active regions where both the upper and lower edges terminate
        // at vEvent (ie. vEvent is closing off these regions).
        // We mark these faces "inside" or "outside" the polygon according
        // to their winding number, and delete the edges from the dictionary.
        // This takes care of all the left-going edges from vEvent.
        ActiveRegion regUp = TopLeftRegion(e._activeRegion);
        ActiveRegion reg = RegionBelow(regUp);
        MeshUtils.Edge eTopLeft = reg._eUp;
        MeshUtils.Edge eBottomLeft = FinishLeftRegions(reg, null);

        // Next we process all the right-going edges from vEvent. This
        // involves adding the edges to the dictionary, and creating the
        // associated "active regions" which record information about the
        // regions between adjacent dictionary edges.
        if (eBottomLeft._Onext == eTopLeft)
        {
            // No right-going edges -- add a temporary "fixable" edge
            ConnectRightVertex(regUp, eBottomLeft);
        }
        else
        {
            AddRightEdges(regUp, eBottomLeft._Onext, eTopLeft, eTopLeft, true);
        }
    }

    /// <summary>
    /// Make the sentinel coordinates big enough that they will never be
    /// merged with real input features.
    /// 
    /// We add two sentinel edges above and below all other edges,
    /// to avoid special cases at the top and bottom.
    /// </summary>
    private void AddSentinel(double smin, double smax, double t)
    {
        MeshUtils.Edge e = _mesh.MakeEdge(_pool);
        e._Org._s = smax;
        e._Org._t = t;
        e._Dst._s = smin;
        e._Dst._t = t;
        _event = e._Dst; // initialize it

        ActiveRegion reg = _pool.Get<ActiveRegion>();
        reg._eUp = e;
        reg._windingNumber = 0;
        reg._inside = false;
        reg._fixUpperEdge = false;
        reg._sentinel = true;
        reg._dirty = false;
        reg._nodeUp = _dict.Insert(reg);
    }

    /// <summary>
    /// We maintain an ordering of edge intersections with the sweep line.
    /// This order is maintained in a dynamic dictionary.
    /// </summary>
    private void InitEdgeDict()
    {
        _dict = new Dict<ActiveRegion>(EdgeLeq);

        AddSentinel(-SentinelCoord, SentinelCoord, -SentinelCoord);
        AddSentinel(-SentinelCoord, SentinelCoord, +SentinelCoord);
    }

    private void DoneEdgeDict()
    {
        while (_dict.Min().Key is { } reg)
        {
            switch (reg._sentinel)
            {
                // At the end of all processing, the dictionary should contain
                // only the two sentinel edges, plus at most one "fixable" edge
                // created by ConnectRightVertex().
                case false:
                    Debug.Assert(reg._fixUpperEdge);
                    int fixedEdges = 0;
                    Debug.Assert(++fixedEdges == 1);
                    break;
            }
            Debug.Assert(reg._windingNumber == 0);
            DeleteRegion(reg);
        }

        _dict = null;
    }

    /// <summary>
    /// Remove zero-length edges, and contours with fewer than 3 vertices.
    /// </summary>
    private void RemoveDegenerateEdges()
    {
        MeshUtils.Edge eHead = _mesh._eHead, e, eNext;

        for (e = eHead._next; e != eHead; e = eNext)
        {
            eNext = e._next;
            MeshUtils.Edge eLnext = e._Lnext;

            if (Geom.VertEq(e._Org, e._Dst) && e._Lnext._Lnext != e)
            {
                // Zero-length edge, contour has at least 3 edges

                SpliceMergeVertices(eLnext, e);	// deletes e.Org
                Mesh.Delete(_pool, e); // e is a self-loop
                e = eLnext;
                eLnext = e._Lnext;
            }

            if (eLnext._Lnext != e)
            {
                continue;
            }
            // Degenerate contour (one or two edges)

            if (eLnext != e)
            {
                if (eLnext == eNext || eLnext == eNext._Sym)
                {
                    eNext = eNext._next;
                }
                Mesh.Delete(_pool, eLnext);
            }
            if (e == eNext || e == eNext._Sym)
            {
                eNext = eNext._next;
            }
            Mesh.Delete(_pool, e);
        }
    }

    /// <summary>
    /// Insert all vertices into the priority queue which determines the
    /// order in which vertices cross the sweep line.
    /// </summary>
    private void InitPriorityQ()
    {
        MeshUtils.Vertex vHead = _mesh._vHead, v;
        int vertexCount = 0;

        for (v = vHead._next; v != vHead; v = v._next)
        {
            vertexCount++;
        }
        // Make sure there is enough space for sentinels.
        vertexCount += 8;

        _pq = new PriorityQueue<MeshUtils.Vertex>(vertexCount, Geom.VertLeq);

        vHead = _mesh._vHead;
        for (v = vHead._next; v != vHead; v = v._next)
        {
            v._pqHandle = _pq.Insert(v);
            if (v._pqHandle._handle == PQHandle.Invalid)
            {
                throw new InvalidOperationException("PQHandle should not be invalid");
            }
        }
        _pq.Init();
    }

    private void DonePriorityQ()
    {
        _pq = null;
    }

    /// <summary>
    /// Delete any degenerate faces with only two edges.  WalkDirtyRegions()
    /// will catch almost all of these, but it won't catch degenerate faces
    /// produced by splice operations on already-processed edges.
    /// The two places this can happen are in FinishLeftRegions(), when
    /// we splice in a "temporary" edge produced by ConnectRightVertex(),
    /// and in CheckForLeftSplice(), where we splice already-processed
    /// edges to ensure that our dictionary invariants are not violated
    /// by numerical errors.
    /// 
    /// In both these cases it is *very* dangerous to delete the offending
    /// edge at the time, since one of the routines further up the stack
    /// will sometimes be keeping a pointer to that edge.
    /// </summary>
    private void RemoveDegenerateFaces()
    {
        MeshUtils.Face f, fNext;

        for (f = _mesh._fHead._next; f != _mesh._fHead; f = fNext)
        {
            fNext = f._next;
            MeshUtils.Edge e = f._anEdge;
            Debug.Assert(e._Lnext != e);

            if (e._Lnext._Lnext != e)
            {
                continue;
            }

            // A face with only two edges
            Geom.AddWinding(e._Onext, e);
            Mesh.Delete(_pool, e);
        }
    }

    /// <summary>
    /// ComputeInterior computes the planar arrangement specified
    /// by the given contours, and further subdivides this arrangement
    /// into regions.  Each region is marked "inside" if it belongs
    /// to the polygon, according to the rule given by windingRule.
    /// Each interior region is guaranteed to be monotone.
    /// </summary>
    protected void ComputeInterior()
    {
        // Each vertex defines an event for our sweep line. Start by inserting
        // all the vertices in a priority queue. Events are processed in
        // lexicographic order, ie.
        // 
        // e1 < e2  iff  e1.x < e2.x || (e1.x == e2.x && e1.y < e2.y)
        RemoveDegenerateEdges();
        InitPriorityQ();
        RemoveDegenerateFaces();
        InitEdgeDict();

        while (_pq.ExtractMin() is { } v)
        {
            while (true)
            {
                MeshUtils.Vertex vNext = _pq.Minimum();
                if (vNext == null || !Geom.VertEq(vNext, v))
                {
                    break;
                }

                // Merge together all vertices at exactly the same location.
                // This is more efficient than processing them one at a time,
                // simplifies the code (see ConnectLeftDegenerate), and is also
                // important for correct handling of certain degenerate cases.
                // For example, suppose there are two identical edges A and B
                // that belong to different contours (so without this code they would
                // be processed by separate sweep events). Suppose another edge C
                // crosses A and B from above. When A is processed, we split it
                // at its intersection point with C. However this also splits C,
                // so when we insert B we may compute a slightly different
                // intersection point. This might leave two edges with a small
                // gap between them. This kind of error is especially obvious
                // when using boundary extraction (BoundaryOnly).
                vNext = _pq.ExtractMin();
                SpliceMergeVertices(v._anEdge, vNext._anEdge);
            }
            SweepEvent(v);
        }

        DoneEdgeDict();
        DonePriorityQ();

        RemoveDegenerateFaces();
        _mesh.Check();
    }
}