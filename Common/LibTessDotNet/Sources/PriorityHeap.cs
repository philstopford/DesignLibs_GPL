/*
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

internal struct PQHandle
{
    public const int Invalid = 0x0fffffff;
    internal int _handle;
}

internal class PriorityHeap<TValue> where TValue : class
{
    public delegate bool LessOrEqual(TValue lhs, TValue rhs);

    protected class HandleElem
    {
        internal TValue _key;
        internal int _node;
    }

    private LessOrEqual _leq;
    private int[] _nodes;
    private HandleElem[] _handles;
    private int _size, _max;
    private int _freeList;
    private bool _initialized;

    public bool Empty => _size == 0;

    public PriorityHeap(int initialSize, LessOrEqual leq)
    {
        _leq = leq;

        _nodes = new int[initialSize + 1];
        _handles = new HandleElem[initialSize + 1];

        _size = 0;
        _max = initialSize;
        _freeList = 0;
        _initialized = false;

        _nodes[1] = 1;
        _handles[1] = new HandleElem { _key = null };
    }

    private void FloatDown(int curr)
    {
        int hCurr = _nodes[curr];
        int size = _size; // Cache field access
        HandleElem[] handles = _handles; // Cache field access
        int[] nodes = _nodes; // Cache field access
        
        while (true)
        {
            int child = curr << 1;
            if (child < size && _leq(handles[nodes[child + 1]]._key, handles[nodes[child]]._key))
            {
                ++child;
            }

            Debug.Assert(child <= _max);

            int hChild = nodes[child];
            if (child > size || _leq(handles[hCurr]._key, handles[hChild]._key))
            {
                nodes[curr] = hCurr;
                handles[hCurr]._node = curr;
                break;
            }

            nodes[curr] = hChild;
            handles[hChild]._node = curr;
            curr = child;
        }
    }

    private void FloatUp(int curr)
    {
        int hCurr = _nodes[curr];
        HandleElem[] handles = _handles; // Cache field access
        int[] nodes = _nodes; // Cache field access
        
        while (true)
        {
            int parent = curr >> 1;
            int hParent = nodes[parent];
            if (parent == 0 || _leq(handles[hParent]._key, handles[hCurr]._key))
            {
                nodes[curr] = hCurr;
                handles[hCurr]._node = curr;
                break;
            }
            nodes[curr] = hParent;
            handles[hParent]._node = curr;
            curr = parent;
        }
    }

    public void Init()
    {
        for (int i = _size; i >= 1; --i)
        {
            FloatDown(i);
        }
        _initialized = true;
    }

    public PQHandle Insert(TValue value)
    {
        int curr = ++_size;
        if (curr * 2 > _max)
        {
            _max <<= 1;
            Array.Resize(ref _nodes, _max + 1);
            Array.Resize(ref _handles, _max + 1);
        }

        int free;
        switch (_freeList)
        {
            case 0:
                free = curr;
                break;
            default:
                free = _freeList;
                _freeList = _handles[free]._node;
                break;
        }

        _nodes[curr] = free;
        switch (_handles[free])
        {
            case null:
                _handles[free] = new HandleElem { _key = value, _node = curr };
                break;
            default:
                _handles[free]._node = curr;
                _handles[free]._key = value;
                break;
        }

        switch (_initialized)
        {
            case true:
                FloatUp(curr);
                break;
        }

        Debug.Assert(free != PQHandle.Invalid);
        return new PQHandle { _handle = free };
    }

    public TValue ExtractMin()
    {
        Debug.Assert(_initialized);

        int hMin = _nodes[1];
        TValue min = _handles[hMin]._key;

        switch (_size)
        {
            case > 0:
            {
                _nodes[1] = _nodes[_size];
                _handles[_nodes[1]]._node = 1;

                _handles[hMin]._key = null;
                _handles[hMin]._node = _freeList;
                _freeList = hMin;

                if (--_size > 0)
                {
                    FloatDown(1);
                }

                break;
            }
        }

        return min;
    }

    public TValue Minimum()
    {
        Debug.Assert(_initialized);
        return _handles[_nodes[1]]._key;
    }

    public void Remove(PQHandle handle)
    {
        Debug.Assert(_initialized);

        int hCurr = handle._handle;
        Debug.Assert(hCurr >= 1 && hCurr <= _max && _handles[hCurr]._key != null);

        int curr = _handles[hCurr]._node;
        _nodes[curr] = _nodes[_size];
        _handles[_nodes[curr]]._node = curr;

        if (curr <= --_size)
        {
            if (curr <= 1 || _leq(_handles[_nodes[curr >> 1]]._key, _handles[_nodes[curr]]._key))
            {
                FloatDown(curr);
            }
            else
            {
                FloatUp(curr);
            }
        }

        _handles[hCurr]._key = null;
        _handles[hCurr]._node = _freeList;
        _freeList = hCurr;
    }
}