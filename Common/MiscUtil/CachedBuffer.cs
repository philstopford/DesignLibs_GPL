using System;
using System.Runtime.CompilerServices;

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

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Dispose()
    {
        if (clearOnDispose)
        {
            Array.Clear(Bytes, 0, Bytes.Length);
        }
        available = true;
    }
}