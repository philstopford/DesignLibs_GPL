using Eto.Drawing;
using Eto.Veldrid;
using Eto.Veldrid.Wpf;
using System;
using System.Runtime.InteropServices;
using Veldrid;
using Veldrid.OpenGL;

[assembly: Eto.ExportHandler(typeof(VeldridSurface), typeof(WpfVeldridSurfaceHandler))]

namespace Eto.Veldrid.Wpf
{
	public class WpfVeldridSurfaceHandler : Eto.Wpf.Forms.ManualBubbleWindowsFormsHostHandler<WinFormsVeldridUserControl, VeldridSurface, VeldridSurface.ICallback>, VeldridSurface.IHandler
	{
		public Size RenderSize => Size.Round((SizeF)Widget.Size * Scale);

		float Scale => Widget.ParentWindow?.LogicalPixelSize ?? 1;

		public WpfVeldridSurfaceHandler() : base(new WinFormsVeldridUserControl())
		{
			Control.Loaded += Control_Loaded;
		}

		public void InitializeGraphicsDevice(VeldridSurface surface, InitializeEventArgs e)
		{
			switch (surface.Backend)
			{
				case GraphicsBackend.Metal:
					surface.GraphicsDevice = GraphicsDevice.CreateMetal(surface.GraphicsDeviceOptions);
					break;
				case GraphicsBackend.Direct3D11:
					surface.GraphicsDevice = GraphicsDevice.CreateD3D11(surface.GraphicsDeviceOptions);
					break;
				case GraphicsBackend.Vulkan:
					surface.GraphicsDevice = GraphicsDevice.CreateVulkan(surface.GraphicsDeviceOptions);
					break;
				case GraphicsBackend.OpenGL:
					surface.GraphicsDevice = GraphicsDevice.CreateOpenGL(
						surface.GraphicsDeviceOptions,
						new OpenGLPlatformInfo(
							surface.OpenGL.OpenGLContextHandle,
							surface.OpenGL.GetProcAddress,
							surface.OpenGL.MakeCurrent,
							surface.OpenGL.GetCurrentContext,
							surface.OpenGL.ClearCurrentContext,
							surface.OpenGL.DeleteContext,
							surface.OpenGL.SwapBuffers,
							surface.OpenGL.SetSyncToVerticalBlank,
							surface.OpenGL.SetSwapchainFramebuffer,
							surface.OpenGL.ResizeSwapchain),
						(uint)e.Width,
						(uint)e.Height);
					break;
				default:
					string message;
					if (!Enum.IsDefined(typeof(GraphicsBackend), surface.Backend))
					{
						message = "Unrecognized backend!";
					}
					else
					{
						message = "Specified backend not supported on this platform!";
					}

					throw new ArgumentException(message);
			}
		}

		public Swapchain? CreateSwapchain()
		{
			Swapchain? swapchain;

			if (Widget.Backend == GraphicsBackend.OpenGL)
			{
				swapchain = Widget.GraphicsDevice?.MainSwapchain;
			}
			else
			{
				// To embed Veldrid in an Eto control, these platform-specific
				// versions of CreateSwapchain use the technique outlined here:
				//
				//   https://github.com/mellinoe/veldrid/issues/155
				//
				var source = SwapchainSource.CreateWin32(
					WinFormsControl.Handle,
					Marshal.GetHINSTANCE(typeof(VeldridSurface).Module));

				var renderSize = RenderSize;
				swapchain = Widget.GraphicsDevice?.ResourceFactory.CreateSwapchain(
					new SwapchainDescription(
						source,
						(uint)renderSize.Width,
						(uint)renderSize.Height,
						Widget.GraphicsDeviceOptions.SwapchainDepthFormat,
						Widget.GraphicsDeviceOptions.SyncToVerticalBlank,
						Widget.GraphicsDeviceOptions.SwapchainSrgbFormat));
			}

			return swapchain;
		}

		private void Control_Loaded(object sender, System.Windows.RoutedEventArgs e)
		{
			Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));

			Control.Loaded -= Control_Loaded;
			Widget.SizeChanged += Widget_SizeChanged;
		}

		private void Widget_SizeChanged(object? sender, EventArgs e)
		{
			Callback.OnResize(Widget, new ResizeEventArgs(RenderSize));
		}

		public override void AttachEvent(string id)
		{
			switch (id)
			{
				case VeldridSurface.DrawEvent:
					WinFormsControl.Paint += (sender, e) => Callback.OnDraw(Widget, EventArgs.Empty);
					break;
				default:
					base.AttachEvent(id);
					break;
			}
		}
	}
}
