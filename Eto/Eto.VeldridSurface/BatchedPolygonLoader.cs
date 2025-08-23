using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace VeldridEto;

public class BatchedPolygonLoader
{
	private readonly OVPSettings _settings;
	private volatile bool _isLoading = false;
	private volatile bool _shouldCancel = false;
	private CancellationTokenSource? _cancellationTokenSource;
	
	// Progress tracking
	public event Action<int, int>? ProgressChanged; // (processedCount, totalCount)
	public event Action<TimeSpan>? LoadingCompleted; // elapsed time
	public event Action<string>? LoadingStatusChanged; // status message
	
	// Batch state tracking
	private int _processedForegroundCount = 0;
	private int _processedBackgroundCount = 0;
	private int _processedTessellatedCount = 0;
	private int _totalForegroundCount = 0;
	private int _totalBackgroundCount = 0;
	private int _totalTessellatedCount = 0;
	
	public bool IsLoading => _isLoading;
	public float Progress => GetOverallProgress();
	
	public BatchedPolygonLoader(OVPSettings settings)
	{
		_settings = settings ?? throw new ArgumentNullException(nameof(settings));
	}
	
	private float GetOverallProgress()
	{
		int totalProcessed = _processedForegroundCount + _processedBackgroundCount + _processedTessellatedCount;
		int totalCount = _totalForegroundCount + _totalBackgroundCount + _totalTessellatedCount;
		
		return totalCount > 0 ? (float)totalProcessed / totalCount : 0f;
	}
	
	public void CancelLoading()
	{
		_shouldCancel = true;
		_cancellationTokenSource?.Cancel();
	}
	
	public async Task<BatchLoadResult> LoadPolygonsProgressivelyAsync(
		Func<PolygonBatch, Task> onBatchReady,
		IProgress<BatchLoadProgress>? progress = null)
	{
		if (_isLoading)
		{
			return new BatchLoadResult { Success = false, ErrorMessage = "Loading already in progress" };
		}
		
		_isLoading = true;
		_shouldCancel = false;
		_cancellationTokenSource = new CancellationTokenSource();
		
		var stopwatch = Stopwatch.StartNew();
		
		try
		{
			LoadingStatusChanged?.Invoke("Initializing progressive loading...");
			
			// Count total polygons
			_totalForegroundCount = _settings.polyList?.Count ?? 0;
			_totalBackgroundCount = _settings.bgPolyList?.Count ?? 0;
			_totalTessellatedCount = _settings.tessPolyList?.Count ?? 0;
			
			_processedForegroundCount = 0;
			_processedBackgroundCount = 0;
			_processedTessellatedCount = 0;
			
			int batchSize = _settings.getBatchSize();
			var totalCount = _totalForegroundCount + _totalBackgroundCount + _totalTessellatedCount;
			
			LoadingStatusChanged?.Invoke($"Loading {totalCount} polygons in batches of {batchSize}...");
			
			// Process foreground polygons in batches
			if (_totalForegroundCount > 0 && _settings.polyList != null)
			{
				await ProcessPolygonListInBatches(
					_settings.polyList,
					PolygonType.Foreground,
					batchSize,
					onBatchReady,
					progress,
					_cancellationTokenSource.Token);
			}
			
			// Process background polygons in batches
			if (_totalBackgroundCount > 0 && !_shouldCancel && _settings.bgPolyList != null)
			{
				await ProcessPolygonListInBatches(
					_settings.bgPolyList,
					PolygonType.Background,
					batchSize,
					onBatchReady,
					progress,
					_cancellationTokenSource.Token);
			}
			
			// Process tessellated polygons in batches
			if (_totalTessellatedCount > 0 && !_shouldCancel && _settings.tessPolyList != null)
			{
				await ProcessPolygonListInBatches(
					_settings.tessPolyList,
					PolygonType.Tessellated,
					batchSize,
					onBatchReady,
					progress,
					_cancellationTokenSource.Token);
			}
			
			stopwatch.Stop();
			
			if (_shouldCancel)
			{
				LoadingStatusChanged?.Invoke("Loading cancelled");
				return new BatchLoadResult { Success = false, ErrorMessage = "Loading was cancelled" };
			}
			
			LoadingStatusChanged?.Invoke($"Loading completed in {stopwatch.ElapsedMilliseconds}ms");
			LoadingCompleted?.Invoke(stopwatch.Elapsed);
			
			return new BatchLoadResult
			{
				Success = true,
				ElapsedTime = stopwatch.Elapsed,
				ProcessedCount = _processedForegroundCount + _processedBackgroundCount + _processedTessellatedCount,
				TotalCount = totalCount
			};
		}
		catch (OperationCanceledException)
		{
			LoadingStatusChanged?.Invoke("Loading cancelled");
			return new BatchLoadResult { Success = false, ErrorMessage = "Loading was cancelled" };
		}
		catch (Exception ex)
		{
			LoadingStatusChanged?.Invoke($"Loading failed: {ex.Message}");
			return new BatchLoadResult { Success = false, ErrorMessage = ex.Message };
		}
		finally
		{
			_isLoading = false;
			_cancellationTokenSource?.Dispose();
			_cancellationTokenSource = null;
		}
	}
	
	private async Task ProcessPolygonListInBatches(
		List<ovp_Poly> polygons,
		PolygonType polygonType,
		int batchSize,
		Func<PolygonBatch, Task> onBatchReady,
		IProgress<BatchLoadProgress>? progress,
		CancellationToken cancellationToken)
	{
		var totalCount = polygons.Count;
		var processedCount = 0;
		
		for (int i = 0; i < totalCount; i += batchSize)
		{
			if (_shouldCancel || cancellationToken.IsCancellationRequested)
			{
				break;
			}
			
			var batchStopwatch = Stopwatch.StartNew();
			
			// Create batch
			var endIndex = Math.Min(i + batchSize, totalCount);
			var batchPolygons = polygons.Skip(i).Take(endIndex - i).ToList();
			
			var batch = new PolygonBatch
			{
				Polygons = batchPolygons,
				PolygonType = polygonType,
				BatchIndex = i / batchSize,
				StartIndex = i,
				EndIndex = endIndex - 1,
				IsLastBatch = endIndex >= totalCount
			};
			
			// Process batch
			try
			{
				await onBatchReady(batch);
			}
			catch (Exception ex)
			{
				Console.WriteLine($"Error processing batch {batch.BatchIndex}: {ex.Message}");
				// Continue processing other batches instead of failing completely
				continue;
			}
			
			// Update progress
			processedCount = endIndex;
			
			switch (polygonType)
			{
				case PolygonType.Foreground:
					_processedForegroundCount = processedCount;
					break;
				case PolygonType.Background:
					_processedBackgroundCount = processedCount;
					break;
				case PolygonType.Tessellated:
					_processedTessellatedCount = processedCount;
					break;
			}
			
			var overallProgress = GetOverallProgress();
			var totalProcessed = _processedForegroundCount + _processedBackgroundCount + _processedTessellatedCount;
			var grandTotal = _totalForegroundCount + _totalBackgroundCount + _totalTessellatedCount;
			
			ProgressChanged?.Invoke(totalProcessed, grandTotal);
			
			progress?.Report(new BatchLoadProgress
			{
				ProcessedCount = totalProcessed,
				TotalCount = grandTotal,
				CurrentBatch = batch.BatchIndex + 1,
				PolygonType = polygonType,
				BatchElapsedMs = (int)batchStopwatch.ElapsedMilliseconds,
				OverallProgress = overallProgress
			});
			
			batchStopwatch.Stop();
			
			// Throttle processing to maintain UI responsiveness
			var maxProcessingTime = _settings.getMaxBatchProcessingTimeMs();
			if (batchStopwatch.ElapsedMilliseconds < maxProcessingTime)
			{
				var remainingTime = maxProcessingTime - (int)batchStopwatch.ElapsedMilliseconds;
				if (remainingTime > 0)
				{
					await Task.Delay(Math.Min(remainingTime, 5), cancellationToken); // Max 5ms delay
				}
			}
		}
	}
}

public enum PolygonType
{
	Foreground,
	Background,
	Tessellated
}

public class PolygonBatch
{
	public List<ovp_Poly> Polygons { get; set; } = new();
	public PolygonType PolygonType { get; set; }
	public int BatchIndex { get; set; }
	public int StartIndex { get; set; }
	public int EndIndex { get; set; }
	public bool IsLastBatch { get; set; }
}

public class BatchLoadProgress
{
	public int ProcessedCount { get; set; }
	public int TotalCount { get; set; }
	public int CurrentBatch { get; set; }
	public PolygonType PolygonType { get; set; }
	public int BatchElapsedMs { get; set; }
	public float OverallProgress { get; set; }
}

public class BatchLoadResult
{
	public bool Success { get; set; }
	public string? ErrorMessage { get; set; }
	public TimeSpan ElapsedTime { get; set; }
	public int ProcessedCount { get; set; }
	public int TotalCount { get; set; }
}