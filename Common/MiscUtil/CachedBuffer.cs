using System;

namespace MiscUtil;

/// <summary>
/// Type of buffer returned by CachingBufferManager.
/// </summary>
internal class CachedBuffer : IBuffer
{
    private readonly byte[] data;
    private volatile bool available;
    private readonly bool clearOnDispose;

    internal CachedBuffer(int size, bool clearOnDispose)
    {
        data = new byte[size];
        this.clearOnDispose = clearOnDispose;
    }

    internal bool Available
    {
        get { return available; }
        set { available = value; }
    }

    public byte[] Bytes
    {
        get { return data; }
    }

    public void Dispose()
    {
        switch (clearOnDispose)
        {
            case true:
                Array.Clear(data, 0, data.Length);
                break;
        }

        available = true;
    }
}