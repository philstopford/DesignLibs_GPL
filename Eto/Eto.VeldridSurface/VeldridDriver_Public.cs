using Eto;
using Eto.Drawing;
using Eto.Forms;
using VeldridEto;
using System.Diagnostics;

namespace VeldridEto;

public partial class VeldridDriver
{
	// This gets connected to the reference passed in, to allow easier updates from the client.
	public OVPSettings ovpSettings;

	// Delegates to allow client to work with updates/selections.
	public delegate void updateHost();

	public updateHost? updateHostFunc { get; set; }

	public delegate void updateHostSelection(int index);

	public updateHostSelection? updateHostSelectionFunc { get; set; }

	// Clock is exposed to allow client side to enable automatic time-based updates
	public UITimer Clock { get; } = new();

	// Is saved location valid?
	public bool savedLocation_valid { get; private set; }

	// Diagnostic flag to enable/disable performance logging
	private bool _enableDiagnostics = false;
	public bool EnableDiagnostics 
	{ 
		get => _enableDiagnostics || Environment.GetEnvironmentVariable("VIEWPORT_DIAGNOSTICS") == "1";
		set => _enableDiagnostics = value;
	}

	public void updateViewport()
	{
		var totalTimer = Stopwatch.StartNew();
		
		if (EnableDiagnostics)
		{
			Console.WriteLine($"[VIEWPORT DIAG] updateViewport() called on {(Platform.Instance.IsGtk ? "GTK" : Platform.Instance.IsWpf ? "WPF" : "Other")} platform");
			Console.WriteLine($"[VIEWPORT DIAG] ovpSettings.changed={ovpSettings.changed}, drawing={drawing}, done_drawing={done_drawing}");
		}
		
		var asyncTimer = Stopwatch.StartNew();
		pUpdateViewportAsync().GetAwaiter().GetResult();
		asyncTimer.Stop();
		
		if (EnableDiagnostics)
		{
			Console.WriteLine($"[VIEWPORT DIAG] pUpdateViewportAsync completed in {asyncTimer.ElapsedMilliseconds}ms");
			Console.WriteLine($"[VIEWPORT DIAG] After async: done_drawing={done_drawing}");
		}
		
		if (done_drawing)
		{
			var callbackTimer = Stopwatch.StartNew();
			updateHostFunc?.Invoke();
			callbackTimer.Stop();
			
			if (EnableDiagnostics && callbackTimer.ElapsedMilliseconds > 0)
			{
				Console.WriteLine($"[VIEWPORT DIAG] updateHostFunc callback took {callbackTimer.ElapsedMilliseconds}ms");
			}
			
			var invalidateTimer = Stopwatch.StartNew();
			Surface!.Invalidate();
			invalidateTimer.Stop();
			
			if (EnableDiagnostics)
			{
				Console.WriteLine($"[VIEWPORT DIAG] Surface.Invalidate() took {invalidateTimer.ElapsedMilliseconds}ms");
			}
			
			ovpSettings.changed = false;
			drawing = false;
			done_drawing = false;
		}
		
		totalTimer.Stop();
		
		if (EnableDiagnostics)
		{
			Console.WriteLine($"[VIEWPORT DIAG] Total updateViewport() time: {totalTimer.ElapsedMilliseconds}ms");
			Console.WriteLine($"[VIEWPORT DIAG] ----------------------------------------");
		}
	}

	public void setContextMenu(ref ContextMenu menu_)
	{
		menu = menu_;
	}

	public void changeSettingsRef(ref OVPSettings newSettings)
	{
		ovpSettings = newSettings;
		updateViewport();
	}

	public void saveLocation()
	{
		savedLocation = new PointF(ovpSettings.getCameraX(), ovpSettings.getCameraY());
		savedLocation_valid = true;
	}

	public void zoomExtents(int index)
	{
		getExtents(index);

		if (ovpSettings.polyList!.Count == 0 && ovpSettings.lineList!.Count == 0 ||
		    ovpSettings is { minX: 0, maxX: 0 } ||
		    ovpSettings is { minY: 0, maxY: 0 })
		{
			reset();
			return;
		}

		// Locate camera at center of the polygon field.
		float dX = ovpSettings.maxX - ovpSettings.minX;
		float dY = ovpSettings.maxY - ovpSettings.minY;
		float cX = dX / 2.0f + ovpSettings.minX;
		float cY = dY / 2.0f + ovpSettings.minY;

		// Now need to get the zoom level organized.
		float zoomLevel_x = dX / Surface!.Width;
		float zoomLevel_y = dY / Surface!.Height;

		if (zoomLevel_x > zoomLevel_y)
		{
			ovpSettings.setZoomFactor(zoomLevel_x / ovpSettings.getBaseZoom());
		}
		else
		{
			ovpSettings.setZoomFactor(zoomLevel_y / ovpSettings.getBaseZoom());
		}

		goToLocation(cX, cY);
	}

	public void loadLocation()
	{
		switch (savedLocation_valid)
		{
			case true:
				ovpSettings.setCameraPos(savedLocation.X, savedLocation.Y);
				updateViewport();
				break;
		}
	}

	public void goToLocation(float x, float y)
	{
		ovpSettings.setCameraPos(x, y);
		updateViewport();
	}

	public void freeze_thaw()
	{
		ovpSettings.lockVP(!ovpSettings.isLocked());
		updateHostFunc?.Invoke();
	}

	public void zoomIn(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() + ovpSettings.getZoomStep() * 0.01f * delta);
		updateHostFunc?.Invoke();
	}

	public void zoomOut(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() - ovpSettings.getZoomStep() * 0.01f * delta);
		updateHostFunc?.Invoke();
	}

	public void fastZoomIn(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() * calcZoom(delta));
		updateHostFunc?.Invoke();
	}

	public void fastZoomOut(float delta)
	{
		ovpSettings.setZoomFactor(ovpSettings.getZoomFactor() / calcZoom(delta));
		updateHostFunc?.Invoke();
	}

	public void reset()
	{
		ovpSettings.resetCamera();
	}
}
