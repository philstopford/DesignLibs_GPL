using System;

namespace MiscUtil;

/// <summary>
/// Type of buffer returned by CachingBufferManager.
/// </summary>
internal class CachedBuffer : IBuffer
{
    private volatile bool available;
    private readonly bool clearOnDispose;

    internal CachedBuffer(int size, bool clearOnDispose)
    {
        Bytes = new byte[size];
        this.clearOnDispose = clearOnDispose;
    }

    internal bool Available
    {
        get => available;
        set => available = value;
    }

    public byte[] Bytes { get; }

    public void Dispose()
    {
        switch (clearOnDispose)
        {
            case true:
                Array.Clear(Bytes, 0, Bytes.Length);
                break;
        }

        available = true;
    }
}