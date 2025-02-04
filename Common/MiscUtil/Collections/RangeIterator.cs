﻿using MiscUtil.Collections.Extensions;
using MiscUtil.Extensions;
using System;
using System.Collections;
using System.Collections.Generic;

namespace MiscUtil.Collections;

/// <summary>
/// Iterates over a range. Despite its name, this implements IEnumerable{T} rather than
/// IEnumerator{T} - it just sounds better, frankly.
/// </summary>
public class RangeIterator<T> : IEnumerable<T>
{
    /// <summary>
    /// Returns the range this object iterates over
    /// </summary>
    public Range<T> Range { get; }

    /// <summary>
    /// Returns the step function used for this range
    /// </summary>
    public Func<T, T> Step { get; }

    /// <summary>
    /// Returns whether or not this iterator works up from the start point (ascending)
    /// or down from the end point (descending)
    /// </summary>
    public bool Ascending { get; }

    /// <summary>
    /// Creates an ascending iterator over the given range with the given step function
    /// </summary>
    public RangeIterator(Range<T> range, Func<T, T> step)
        : this(range, step, true)
    {
    }

    /// <summary>
    /// Creates an iterator over the given range with the given step function,
    /// with the specified direction.
    /// </summary>
    public RangeIterator(Range<T> range, Func<T, T> step, bool ascending)
    {
        step.ThrowIfNull("step");

        switch (ascending)
        {
            case true when range.Comparer.Compare(range.Start, step(range.Start)) >= 0:
            case false when range.Comparer.Compare(range.End, step(range.End)) <= 0:
                throw new ArgumentException("step does nothing, or progresses the wrong way");
        }
        Ascending = ascending;
        Range = range;
        Step = step;
    }

    /// <summary>
    /// Returns an IEnumerator{T} running over the range.
    /// </summary>
    public IEnumerator<T> GetEnumerator()
    {
        // A descending range effectively has the start and end points (and inclusions)
        // reversed, and a reverse comparer.
        bool includesStart = Ascending ? Range.IncludesStart : Range.IncludesEnd;
        bool includesEnd = Ascending ? Range.IncludesEnd : Range.IncludesStart;
        T start = Ascending ? Range.Start : Range.End;
        T end = Ascending ? Range.End : Range.Start;
        IComparer<T> comparer = Ascending ? Range.Comparer : Range.Comparer.Reverse();

        // Now we can use our local version of the range variables to iterate

        T value = start;

        switch (includesStart)
        {
            case true:
            {
                // Deal with possibility that start point = end point
                if (includesEnd || comparer.Compare(value, end) < 0)
                {
                    yield return value;
                }

                break;
            }
        }
        value = Step(value);

        while (comparer.Compare(value, end) < 0)
        {
            yield return value;
            value = Step(value);
        }

        switch (includesEnd)
        {
            // We've already performed a step, therefore we can't
            // still be at the start point
            case true when comparer.Compare(value, end) == 0:
                yield return value;
                break;
        }
    }

    /// <summary>
    /// Returns an IEnumerator running over the range.
    /// </summary>
    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
}