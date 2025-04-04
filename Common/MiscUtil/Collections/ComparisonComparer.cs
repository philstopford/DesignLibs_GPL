using System;
using System.Collections.Generic;

namespace MiscUtil.Collections;

/// <summary>
/// Utility to build an IComparer implementation from a Comparison delegate,
/// and a static method to do the reverse.
/// </summary>
public sealed class ComparisonComparer<T> : IComparer<T>
{
    private readonly Comparison<T> comparison;

    /// <summary>
    /// Creates a new instance which will proxy to the given Comparison
    /// delegate when called.
    /// </summary>
    /// <param name="comparison">Comparison delegate to proxy to. Must not be null.</param>
    public ComparisonComparer(Comparison<T> comparison)
    {
        this.comparison = comparison switch
        {
            null => throw new ArgumentNullException(nameof(comparison)),
            _ => comparison
        };
    }

    /// <summary>
    /// Implementation of IComparer.Compare which simply proxies
    /// to the originally specified Comparison delegate.
    /// </summary>
    public int Compare(T x, T y)
    {
        return comparison(x, y);
    }

    /// <summary>
    /// Creates a Comparison delegate from the given Comparer.
    /// </summary>
    /// <param name="comparer">Comparer to use when the returned delegate is called. Must not be null.</param>
    /// <returns>A Comparison delegate which proxies to the given Comparer.</returns>
    public static Comparison<T> CreateComparison(IComparer<T> comparer)
    {
        return comparer switch
        {
            null => throw new ArgumentNullException(nameof(comparer)),
            _ => delegate(T x, T y) { return comparer.Compare(x, y); }
        };
    }
}