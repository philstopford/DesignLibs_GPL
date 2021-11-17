using MiscUtil.Collections;
using MiscUtil.Collections.Extensions;
using System;
using System.Collections.Generic;
namespace MiscUtil.Linq.Extensions;

/// <summary>
/// Provides extension methods to List&lt;T&gt;
/// </summary>
public static class ListExt
{
    /// <summary>
    /// Sorts the elements in the entire System.Collections.Generic.List{T} using
    /// a projection.
    /// </summary>
    /// <param name="source">Data source</param>
    /// <param name="selector">The projection to use to obtain values for comparison</param>
    /// <param name="comparer">The comparer to use to compare projected values (on null to use the default comparer)</param>
    /// <param name="descending">Should the list be sorted ascending or descending?</param>
    public static void Sort<T, TValue>(this List<T> source, Func<T, TValue> selector, IComparer<TValue> comparer, bool descending)
    {
        switch (source)
        {
            case null:
                throw new ArgumentNullException("source");
        }

        comparer = comparer switch
        {
            null => Comparer<TValue>.Default,
            _ => comparer
        };
        IComparer<T> itemComparer = new ProjectionComparer<T, TValue>(selector, comparer);
        switch (@descending)
        {
            case true:
                itemComparer = itemComparer.Reverse();
                break;
        }
        source.Sort(itemComparer);
    }

    /// <summary>
    /// Sorts the elements in the entire System.Collections.Generic.List{T} using
    /// a projection.
    /// </summary>
    /// <param name="source">Data source</param>
    /// <param name="selector">The projection to use to obtain values for comparison</param>
    public static void Sort<T, TValue>(this List<T> source, Func<T, TValue> selector)
    {
        Sort<T, TValue>(source, selector, null, false);
    }
}