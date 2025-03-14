﻿
namespace MiscUtil;

/// <summary>
/// Generic event argument taking a single value
/// </summary>
public class EventArgs<T> : System.EventArgs
{
    /// <summary>
    /// The typed value of the EventArgs&lt;T&gt;
    /// </summary>
    public T Value { get; }

    /// <summary>
    /// Creates a new EventArgs&lt;T&gt; with the specified value.
    /// </summary>
    /// <param name="value">The Value of the EventArgs&lt;T&gt; instance.</param>
    public EventArgs(T value)
    {
        Value = value;
    }
}