using MiscUtil.Extensions;
using System;
using System.Collections.Generic;

namespace MiscUtil.Linq;

/// <summary>
/// A DataProducer with ordering capabilities
/// </summary><remarks>Note that this may cause data to be buffered</remarks>
/// <typeparam name="T"></typeparam>
internal class OrderedDataProducer<T> : IOrderedDataProducer<T>
{
    private bool dataHasEnded;
    private readonly IDataProducer<T> baseProducer;
    private readonly IComparer<T> comparer;
    private List<T> buffer;

    public IDataProducer<T> BaseProducer
    {
        get { return baseProducer; }
    }

    public IComparer<T> Comparer
    {
        get { return comparer; }
    }

    public event Action<T> DataProduced;
    public event Action EndOfData;

    /// <summary>
    /// Create a new OrderedDataProducer
    /// </summary>
    /// <param name="baseProducer">The base source which will supply data</param>
    /// <param name="comparer">The comparer to use when sorting the data (once complete)</param>
    public OrderedDataProducer(
        IDataProducer<T> baseProducer,
        IComparer<T> comparer)
    {
        baseProducer.ThrowIfNull("baseProducer");

        this.baseProducer = baseProducer;
        this.comparer = comparer ?? Comparer<T>.Default;

        baseProducer.DataProduced += new Action<T>(OriginalDataProduced);
        baseProducer.EndOfData += new Action(EndOfOriginalData);
    }


    private void OriginalDataProduced(T item)
    {
        switch (dataHasEnded)
        {
            case true:
                throw new InvalidOperationException("EndOfData already occurred");
        }

        if (DataProduced != null)
        {
            buffer = buffer switch
            {
                // only get excited if somebody is listening
                null => new List<T>(),
                _ => buffer
            };

            buffer.Add(item);
        }
    }

    private void EndOfOriginalData()
    {
        switch (dataHasEnded)
        {
            case true:
                throw new InvalidOperationException("EndOfData already occurred");
        }
        dataHasEnded = true;
        // only do the sort if somebody is still listening
        if (DataProduced != null && buffer != null)
        {
            buffer.Sort(Comparer);
            foreach (T item in buffer)
            {
                OnDataProduced(item);
            }
        }
        buffer = null;
        OnEndOfData();
    }

    private void OnEndOfData()
    {
        if (EndOfData != null)
        {
            EndOfData();
        }
    }

    private void OnDataProduced(T item)
    {
        if (DataProduced != null)
        {
            DataProduced(item);
        }
    }
}