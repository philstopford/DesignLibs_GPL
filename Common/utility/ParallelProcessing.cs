using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace utility;

/// <summary>
/// High-performance parallel processing utilities optimized for numerical computations
/// </summary>
public static class ParallelProcessing
{
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
        Parallel.For(0, processorCount, p =>
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

        int processorCount = Environment.ProcessorCount;
        int alignedBatchCount = (array.Length + batchSize - 1) / batchSize;
        int batchesPerProcessor = alignedBatchCount / processorCount;
        
        Parallel.For(0, processorCount, p =>
        {
            int startBatch = p * batchesPerProcessor;
            int endBatch = p == processorCount - 1 ? alignedBatchCount : startBatch + batchesPerProcessor;
            
            for (int batch = startBatch; batch < endBatch; batch++)
            {
                int start = batch * batchSize;
                int end = Math.Min(start + batchSize, array.Length);
                
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

        var options = new ParallelOptions { MaxDegreeOfParallelism = Environment.ProcessorCount };
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