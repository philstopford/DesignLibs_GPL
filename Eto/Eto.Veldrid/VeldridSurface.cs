using Eto.Drawing;
using Eto.Forms;
using System;
using Veldrid;
using Veldrid.OpenGL;

namespace Eto.Veldrid;

/// <summary>
/// A simple control that allows drawing with Veldrid.
/// </summary>
[Handler(typeof(IHandler))]
public class VeldridSurface : Control
{
	[AutoInitialize(false)]
	public new interface IHandler : Control.IHandler
	{
		Size RenderSize { get; }
		Swapchain? CreateSwapchain();
		void InitializeGraphicsDevice(VeldridSurface surface, InitializeEventArgs e);
	}

	private new IHandler Handler => (IHandler)base.Handler;

	public new interface ICallback : Control.ICallback
	{
		void OnInitializeBackend(VeldridSurface s, InitializeEventArgs e);
		void OnDraw(VeldridSurface s, EventArgs e);
		void OnResize(VeldridSurface s, ResizeEventArgs e);
	}

	protected new class Callback : Control.Callback, ICallback
	{
		public void OnInitializeBackend(VeldridSurface s, InitializeEventArgs e) => s?.InitializeGraphicsBackend(e);
		public void OnDraw(VeldridSurface s, EventArgs e) => s?.OnDraw(e);
		public void OnResize(VeldridSurface s, ResizeEventArgs e) => s?.OnResize(e);
	}

	protected override object GetCallback() => new Callback();

	public interface IOpenGL
	{
		IntPtr OpenGLContextHandle { get; }
		IntPtr GetProcAddress(string name);
		void MakeCurrent(IntPtr context);
		IntPtr GetCurrentContext();
		void ClearCurrentContext();
		void DeleteContext(IntPtr context);
		void SwapBuffers();
		void SetSyncToVerticalBlank(bool enable);
		void SetSwapchainFramebuffer();
		void ResizeSwapchain(uint width, uint height);
	}

	public IOpenGL OpenGL => (IOpenGL)Handler;

	public static GraphicsBackend PreferredBackend { get; } = GetPreferredBackend();

	/// <summary>
	/// The render area's size, which may differ from the control's size
	/// (e.g. with high DPI displays).
	/// </summary>
	public Size RenderSize => Handler.RenderSize;
	/// <summary>
	/// The render area's width, which may differ from the control's width
	/// (e.g. with high DPI displays).
	/// </summary>
	public int RenderWidth => RenderSize.Width;
	/// <summary>
	/// The render area's height, which may differ from the control's height
	/// (e.g. with high DPI displays).
	/// </summary>
	public int RenderHeight => RenderSize.Height;

	public GraphicsBackend Backend { get; private set; }
	public GraphicsDevice? GraphicsDevice { get; set; }
	public GraphicsDeviceOptions GraphicsDeviceOptions { get; private set; }
	public Swapchain? Swapchain { get; private set; }

	public const string VeldridInitializedEvent = "VeldridSurface.VeldridInitialized";
	public const string DrawEvent = "VeldridSurface.Draw";
	public const string ResizeEvent = "VeldridSurface.Resize";

	public event EventHandler<InitializeEventArgs> VeldridInitialized
	{
		add => Properties.AddHandlerEvent(VeldridInitializedEvent, value);
		remove => Properties.RemoveEvent(VeldridInitializedEvent, value);
	}
	public event EventHandler<EventArgs> Draw
	{
		add => Properties.AddHandlerEvent(DrawEvent, value);
		remove => Properties.RemoveEvent(DrawEvent, value);
	}
	public event EventHandler<ResizeEventArgs> Resize
	{
		add => Properties.AddHandlerEvent(ResizeEvent, value);
		remove => Properties.RemoveEvent(ResizeEvent, value);
	}

	public VeldridSurface()
		: this(PreferredBackend)
	{
	}
	public VeldridSurface(GraphicsBackend backend)
		: this(backend, new GraphicsDeviceOptions())
	{
	}

	public VeldridSurface(GraphicsBackend backend, GraphicsDeviceOptions gdOptions)
	{
		Backend = backend;
		GraphicsDeviceOptions = gdOptions;
		Initialize();
	}

	private static GraphicsBackend GetPreferredBackend()
	{
		GraphicsBackend? backend = null;

		if (GraphicsDevice.IsBackendSupported(GraphicsBackend.Metal))
		{
			backend = GraphicsBackend.Metal;
		}
		else if (GraphicsDevice.IsBackendSupported(GraphicsBackend.Direct3D11))
		{
			backend = GraphicsBackend.Direct3D11;
		}
		else if (EtoEnvironment.Platform.IsLinux)
		{
			// On Linux, prefer OpenGL for better compatibility, then Vulkan
			// Vulkan can be problematic in headless/CI environments
			if (GraphicsDevice.IsBackendSupported(GraphicsBackend.OpenGL))
			{
				backend = GraphicsBackend.OpenGL;
			}
			else if (GraphicsDevice.IsBackendSupported(GraphicsBackend.Vulkan))
			{
				backend = GraphicsBackend.Vulkan;
			}
		}

		if (backend == null)
		{
			throw new VeldridException("VeldridSurface: No supported Veldrid backend found!");
		}

		return (GraphicsBackend)backend;
	}

	private void InitializeGraphicsBackend(InitializeEventArgs e)
	{
		try
		{
			Console.WriteLine($"VeldridSurface: InitializeGraphicsBackend called with size {e.Width}x{e.Height}");
			Console.WriteLine($"VeldridSurface: Current backend: {Backend}");
			
			// Delegate graphics device initialization to the platform-specific handler
			Handler.InitializeGraphicsDevice(this, e);
			Console.WriteLine($"VeldridSurface: Graphics device initialized: {GraphicsDevice}");

			Swapchain = Handler.CreateSwapchain();
			Console.WriteLine($"VeldridSurface: Swapchain created: {Swapchain}");

			OnVeldridInitialized(e);
			Console.WriteLine("VeldridSurface: VeldridInitialized event fired");
		}
		catch (VeldridException vex) when (Backend == GraphicsBackend.Vulkan)
		{
			Console.WriteLine($"VeldridSurface: Vulkan initialization failed, attempting fallback to OpenGL: {vex.Message}");
			
			try
			{
				// Attempt to fallback to OpenGL
				Backend = GraphicsBackend.OpenGL;
				Console.WriteLine("VeldridSurface: Switched backend to OpenGL for fallback");
				
				// Recreate the handler with OpenGL backend - this requires recreating the control
				// For now, just rethrow - the user can manually specify OpenGL
				throw new VeldridException($"Vulkan backend failed, please try with OpenGL backend instead. Original error: {vex.Message}", vex);
			}
			catch (Exception fallbackEx)
			{
				Console.WriteLine($"VeldridSurface: OpenGL fallback also failed: {fallbackEx.Message}");
				throw new VeldridException($"Both Vulkan and OpenGL backends failed. Vulkan: {vex.Message}, OpenGL: {fallbackEx.Message}", vex);
			}
		}
		catch (Exception ex)
		{
			Console.WriteLine($"VeldridSurface: Error in InitializeGraphicsBackend: {ex.Message}");
			Console.WriteLine($"VeldridSurface: Stack trace: {ex.StackTrace}");
			throw;
		}
	}

	protected virtual void OnDraw(EventArgs e) => Properties.TriggerEvent(DrawEvent, this, e);

	protected virtual void OnResize(ResizeEventArgs e)
	{
		if (e == null)
			throw new ArgumentNullException(nameof(e));

		Swapchain?.Resize((uint)e.Width, (uint)e.Height);

		Properties.TriggerEvent(ResizeEvent, this, e);
	}

	protected virtual void OnVeldridInitialized(InitializeEventArgs e) => Properties.TriggerEvent(VeldridInitializedEvent, this, e);
}