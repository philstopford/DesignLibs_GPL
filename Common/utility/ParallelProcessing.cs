using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace utility;

/// <summary>
/// High-performance parallel processing utilities optimized for numerical computations
/// </summary>
public static class ParallelProcessing
{
    // Performance optimization: Cache CPU capabilities
    private static readonly bool IsAvx2Available = Avx2.IsSupported;
    private static readonly bool IsAvxAvailable = Avx.IsSupported;
    private static readonly bool IsSse2Available = Sse2.IsSupported;
    
    /// <summary>
    /// High-performance parallel execution optimized for numerical operations with better work distribution
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void OptimizedParallelFor<T>(T[] array, Func<int, T> generator)
    {
        if (array.Length <= 1000)
        {
            // For small arrays, sequential processing is often faster due to overhead
            for (int i = 0; i < array.Length; i++)
            {
                array[i] = generator(i);
            }
            return;
        }

        int processorCount = Math.Min(Environment.ProcessorCount, array.Length / 100);
        int chunkSize = array.Length / processorCount;
        
        // Performance optimization: Use Parallel.For for better thread management
        Parallel.For(0, processorCount, new ParallelOptions 
        { 
            MaxDegreeOfParallelism = processorCount 
        }, p =>
        {
            int start = p * chunkSize;
            int end = p == processorCount - 1 ? array.Length : start + chunkSize;
            
            for (int i = start; i < end; i++)
            {
                array[i] = generator(i);
            }
        });
    }

    /// <summary>
    /// High-performance parallel execution with SIMD-friendly batching
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void VectorizedParallelFor<T>(T[] array, Func<int, T> generator, int batchSize = 16)
    {
        if (array.Length <= batchSize * 4)
        {
            OptimizedParallelFor(array, generator);
            return;
        }

        // Performance optimization: Adjust batch size based on available SIMD capabilities
        int optimizedBatchSize = batchSize;
        if (IsAvx2Available)
            optimizedBatchSize = Math.Max(batchSize, 32); // AVX2 works best with 256-bit (32 byte) alignment
        else if (IsAvxAvailable)
            optimizedBatchSize = Math.Max(batchSize, 16); // AVX works best with 128-bit (16 byte) alignment

        int processorCount = Environment.ProcessorCount;
        int alignedBatchCount = (array.Length + optimizedBatchSize - 1) / optimizedBatchSize;
        int batchesPerProcessor = Math.Max(1, alignedBatchCount / processorCount);
        
        Parallel.For(0, processorCount, new ParallelOptions 
        { 
            MaxDegreeOfParallelism = processorCount 
        }, p =>
        {
            int startBatch = p * batchesPerProcessor;
            int endBatch = p == processorCount - 1 ? alignedBatchCount : startBatch + batchesPerProcessor;
            
            for (int batch = startBatch; batch < endBatch; batch++)
            {
                int start = batch * optimizedBatchSize;
                int end = Math.Min(start + optimizedBatchSize, array.Length);
                
                // Process batch with better cache locality
                for (int i = start; i < end; i++)
                {
                    array[i] = generator(i);
                }
            }
        });
    }
    
        /// <summary>
    /// Optimized parallel execution for operations without return values.
    /// Uses smart work distribution to minimize thread overhead.
    /// </summary>
    /// <param name="fromInclusive">The start index, inclusive</param>
    /// <param name="toExclusive">The end index, exclusive</param>
    /// <param name="body">The delegate that is invoked for each iteration</param>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void OptimizedParallelFor(int fromInclusive, int toExclusive, Action<int> body)
    {
        int length = toExclusive - fromInclusive;

        if (length <= 1000)
        {
            // For small ranges, sequential processing is often faster
            for (int i = fromInclusive; i < toExclusive; i++)
            {
                body(i);
            }
            return;
        }

        // Performance optimization: Use optimal parallelism for the workload
        var options = new ParallelOptions 
        { 
            MaxDegreeOfParallelism = Math.Min(Environment.ProcessorCount, Math.Max(1, length / 1000))
        };
        Parallel.For(fromInclusive, toExclusive, options, body);
    }

    /// <summary>
    /// Optimized parallel execution with ParallelOptions for maximum control.
    /// Automatically optimizes degree of parallelism based on workload size.
    /// </summary>
    /// <param name="fromInclusive">The start index, inclusive</param>
    /// <param name="toExclusive">The end index, exclusive</param>
    /// <param name="body">The delegate that is invoked for each iteration</param>
    /// <param name="options">Additional parallel options (will be optimized internally)</param>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void OptimizedParallelFor(int fromInclusive, int toExclusive, ParallelOptions options, Action<int> body)
    {
        int length = toExclusive - fromInclusive;

        if (length <= 1000)
        {
            // For small ranges, sequential processing is often faster
            for (int i = fromInclusive; i < toExclusive; i++)
            {
                body(i);
            }
            return;
        }

        // Optimize the degree of parallelism based on workload
        var optimizedOptions = new ParallelOptions
        {
            MaxDegreeOfParallelism = Math.Min(
                options.MaxDegreeOfParallelism > 0 ? options.MaxDegreeOfParallelism : Environment.ProcessorCount,
                Math.Max(1, length / 100)
            ),
            CancellationToken = options.CancellationToken,
            TaskScheduler = options.TaskScheduler
        };

        Parallel.For(fromInclusive, toExclusive, optimizedOptions, body);
    }
}