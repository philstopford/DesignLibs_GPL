/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/
using System.Collections.Generic;
using System.Linq;

namespace info.lundin.math;

/// <summary>
/// Dictionary of values
/// </summary>
public class ValuesDictionary : IDictionary<string, Value>
{
    private IDictionary<string, Value> dictionary;

    public ValuesDictionary()
    {
        dictionary = new Dictionary<string, Value>();
    }

    /// <summary>
    /// Adds a variable and string value 
    /// </summary>
    /// <param name="variable">variable</param>
    /// <param name="value">string value</param>
    public void Add(string variable, string value)
    {
        Add(variable, new StringValue { Value = value });
    }

    /// <summary>
    /// Adds a variable and a double value
    /// </summary>
    /// <param name="variable">variable</param>
    /// <param name="value">double value</param>
    public void Add(string variable, double value)
    {
        Add(variable, new DoubleValue { Value = value });
    }

    /// <summary>
    /// Adds a variable and a Value instance as value
    /// </summary>
    /// <param name="variable">variable</param>
    /// <param name="value">value instance</param>
    public void Add(string variable, Value value)
    {
        dictionary.Add(variable, value);
    }

    #region Standard Dictionary Method Implementations

    public bool ContainsKey(string key)
    {
        return dictionary.ContainsKey(key);
    }

    public ICollection<string> Keys
    {
        get { return dictionary.Keys; }
    }

    public bool Remove(string key)
    {
        Value value = null;

        dictionary.TryGetValue(key, out value);
        bool removed = dictionary.Remove(key);

        value.Type = removed switch
        {
            true when value != null => ValueType.Invalid,
            _ => value.Type
        };

        return removed;
    }

    public bool TryGetValue(string key, out Value value)
    {
        return dictionary.TryGetValue(key, out value);
    }

    public ICollection<Value> Values
    {
        get { return dictionary.Values; }
    }

    public Value this[string key]
    {
        get
        {
            return dictionary[key];
        }
        set
        {
            switch (value)
            {
                case null:
                    Remove(key);
                    return;
                default:
                    dictionary[key] = value;
                    break;
            }
        }
    }

    public void Add(KeyValuePair<string, Value> item)
    {
        dictionary.Add(item);
    }

    public void Clear()
    {
        List<KeyValuePair<string, Value>> items = dictionary.ToList();
        foreach (KeyValuePair<string, Value> item in items) Remove(item);
    }

    public bool Contains(KeyValuePair<string, Value> item)
    {
        return dictionary.Contains(item);
    }

    public void CopyTo(KeyValuePair<string, Value>[] array, int arrayIndex)
    {
        dictionary.CopyTo(array, arrayIndex);
    }

    public int Count
    {
        get { return dictionary.Count; }
    }

    public bool IsReadOnly
    {
        get { return dictionary.IsReadOnly; }
    }

    public bool Remove(KeyValuePair<string, Value> item)
    {
        Value value = item.Value;

        bool removed = dictionary.Remove(item);

        value.Type = removed switch
        {
            true when value != null => ValueType.Invalid,
            _ => value.Type
        };

        return removed;
    }

    public IEnumerator<KeyValuePair<string, Value>> GetEnumerator()
    {
        return dictionary.GetEnumerator();
    }

    System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
    {
        return dictionary.GetEnumerator();
    }

    #endregion
}