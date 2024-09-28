namespace MiscUtil.Linq;

/// <summary>
/// Generic tuple for a key and a single value
/// </summary>
/// <typeparam name="TKey">The Type of the key</typeparam>
/// <typeparam name="T">The Type of the value</typeparam>
public struct KeyValueTuple<TKey, T>
{
    /// <summary>
    /// The key for the tuple
    /// </summary>
    public TKey Key { get; }

    /// <summary>
    /// The value for the tuple
    /// </summary>
    public T Value { get; }

    /// <summary>
    /// Creates a new tuple with the given key and value
    /// </summary>
    public KeyValueTuple(TKey key, T value)
    {
        this.Key = key;
        this.Value = value;
    }
}
/// <summary>
/// Generic tuple for a key and a pair of values
/// </summary>
/// <typeparam name="TKey">The Type of the key</typeparam>
/// <typeparam name="T1">The Type of the first value</typeparam>
/// <typeparam name="T2">The Type of the second value</typeparam>    
public struct KeyValueTuple<TKey, T1, T2>
{
    /// <summary>
    /// The key for the tuple
    /// </summary>
    public TKey Key { get; }

    /// <summary>
    /// The first value
    /// </summary>
    public T1 Value1 { get; }

    /// <summary>
    /// The second value
    /// </summary>
    public T2 Value2 { get; }

    /// <summary>
    /// Creates a new tuple with the given key and values
    /// </summary>
    public KeyValueTuple(TKey key, T1 value1, T2 value2)
    {
        this.Key = key;
        this.Value1 = value1;
        this.Value2 = value2;
    }
}
/// <summary>
/// Generic tuple for a key and a trio of values
/// </summary>
/// <typeparam name="TKey">The Type of the key</typeparam>
/// <typeparam name="T1">The Type of the first value</typeparam>
/// <typeparam name="T2">The Type of the second value</typeparam>
/// <typeparam name="T3">The Type of the third value</typeparam>
public struct KeyValueTuple<TKey, T1, T2, T3>
{
    /// <summary>
    /// The key for the tuple
    /// </summary>
    public TKey Key { get; }

    /// <summary>
    /// The first value
    /// </summary>
    public T1 Value1 { get; }

    /// <summary>
    /// The second value
    /// </summary>
    public T2 Value2 { get; }

    /// <summary>
    /// The third value
    /// </summary>
    public T3 Value3 { get; }

    /// <summary>
    /// Creates a new tuple with the given key and values
    /// </summary>
    public KeyValueTuple(TKey key, T1 value1, T2 value2, T3 value3)
    {
        this.Key = key;
        this.Value1 = value1;
        this.Value2 = value2;
        this.Value3 = value3;
    }
}
/// <summary>
/// Generic tuple for a key and a quartet of values
/// </summary>
/// <typeparam name="TKey">The Type of the key</typeparam>
/// <typeparam name="T1">The Type of the first value</typeparam>
/// <typeparam name="T2">The Type of the second value</typeparam>
/// <typeparam name="T3">The Type of the third value</typeparam>
/// <typeparam name="T4">The Type of the fourth value</typeparam>
public struct KeyValueTuple<TKey, T1, T2, T3, T4>
{
    /// <summary>
    /// The key for the tuple
    /// </summary>
    public TKey Key { get; }

    /// <summary>
    /// The first value
    /// </summary>
    public T1 Value1 { get; }

    /// <summary>
    /// The second value
    /// </summary>
    public T2 Value2 { get; }

    /// <summary>
    /// The third value
    /// </summary>
    public T3 Value3 { get; }

    /// <summary>
    /// The fourth value
    /// </summary>
    public T4 Value4 { get; }

    /// <summary>
    /// Creates a new tuple with the given key and values
    /// </summary>
    public KeyValueTuple(TKey key, T1 value1, T2 value2, T3 value3, T4 value4)
    {
        this.Key = key;
        this.Value1 = value1;
        this.Value2 = value2;
        this.Value3 = value3;
        this.Value4 = value4;
    }
}