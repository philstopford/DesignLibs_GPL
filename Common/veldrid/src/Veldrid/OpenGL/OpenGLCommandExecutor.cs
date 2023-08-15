﻿using System;
using static Veldrid.OpenGLBinding.OpenGLNative;
using static Veldrid.OpenGL.OpenGLUtil;
using Veldrid.OpenGLBinding;
using System.Runtime.CompilerServices;
using System.Diagnostics;

namespace Veldrid.OpenGL
{
    internal sealed unsafe class OpenGLCommandExecutor
    {
        private readonly OpenGLGraphicsDevice _gd;
        private readonly GraphicsBackend _backend;
        private readonly OpenGLTextureSamplerManager _textureSamplerManager;
        private readonly StagingMemoryPool _stagingMemoryPool;
        private readonly OpenGLExtensions _extensions;
        private readonly OpenGLPlatformInfo _platformInfo;
        private readonly GraphicsDeviceFeatures _features;

        private Framebuffer? _fb;
        private bool _isSwapchainFB;
        private OpenGLPipeline? _graphicsPipeline;
        private BoundResourceSetInfo[] _graphicsResourceSets = Array.Empty<BoundResourceSetInfo>();
        private bool _graphicsResourcesFlushed;
        private bool[] _newGraphicsResourceSets = Array.Empty<bool>();
        private uint[] _vertexBuffers = Array.Empty<uint>();
        private nuint[] _vbOffsets = Array.Empty<nuint>();
        private bool[] _newVertexBuffers = Array.Empty<bool>();
        private DrawElementsType _drawElementsType;
        private uint _ibOffset;
        private PrimitiveType _primitiveType;
        private OpenGLBuffer? _indexBuffer;

        private OpenGLPipeline? _computePipeline;
        private BoundResourceSetInfo[] _computeResourceSets = Array.Empty<BoundResourceSetInfo>();
        private bool _computeResourcesFlushed;
        private bool[] _newComputeResourceSets = Array.Empty<bool>();

        private bool _graphicsPipelineActive;
        private bool _computePipelineActive;
        private bool _vertexLayoutFlushed;
        private bool _indexBufferBound;

        public OpenGLCommandExecutor(OpenGLGraphicsDevice gd, OpenGLPlatformInfo platformInfo)
        {
            _gd = gd;
            _backend = gd.BackendType;
            _extensions = gd.Extensions;
            _textureSamplerManager = gd.TextureSamplerManager;
            _stagingMemoryPool = gd.StagingMemoryPool;
            _platformInfo = platformInfo;
            _features = gd.Features;
        }

        public void Begin()
        {
        }

        public void ClearColorTarget(uint index, RgbaFloat clearColor)
        {
            if (!_isSwapchainFB)
            {
                DrawBuffersEnum bufs = (DrawBuffersEnum)((uint)DrawBuffersEnum.ColorAttachment0 + index);
                glDrawBuffers(1, &bufs);
                CheckLastError();
            }

            RgbaFloat color = clearColor;
            glClearColor(color.R, color.G, color.B, color.A);
            CheckLastError();

            if (_graphicsPipeline != null && _graphicsPipeline.RasterizerState.ScissorTestEnabled)
            {
                glDisable(EnableCap.ScissorTest);
                CheckLastError();
            }

            glClear(ClearBufferMask.ColorBufferBit);
            CheckLastError();

            if (_graphicsPipeline != null && _graphicsPipeline.RasterizerState.ScissorTestEnabled)
            {
                glEnable(EnableCap.ScissorTest);
            }

            if (!_isSwapchainFB)
            {
                int colorCount = _fb!.ColorTargets.Length;
                DrawBuffersEnum* bufs = stackalloc DrawBuffersEnum[colorCount];
                for (int i = 0; i < colorCount; i++)
                {
                    bufs[i] = DrawBuffersEnum.ColorAttachment0 + i;
                }
                glDrawBuffers((uint)colorCount, bufs);
                CheckLastError();
            }
        }

        public void ClearDepthStencil(float depth, byte stencil)
        {
            glClearDepth_Compat(depth);
            CheckLastError();

            glStencilMask(~0u);
            CheckLastError();

            glClearStencil(stencil);
            CheckLastError();

            if (_graphicsPipeline != null && _graphicsPipeline.RasterizerState.ScissorTestEnabled)
            {
                glDisable(EnableCap.ScissorTest);
                CheckLastError();
            }

            glDepthMask(true);
            glClear(ClearBufferMask.DepthBufferBit | ClearBufferMask.StencilBufferBit);
            CheckLastError2();

            if (_graphicsPipeline != null && _graphicsPipeline.RasterizerState.ScissorTestEnabled)
            {
                glEnable(EnableCap.ScissorTest);
            }
        }

        public void Draw(uint vertexCount, uint instanceCount, uint vertexStart, uint instanceStart)
        {
            PreDrawCommand();

            if (instanceCount == 1 && instanceStart == 0)
            {
                glDrawArrays(_primitiveType, (int)vertexStart, vertexCount);
            }
            else
            {
                if (instanceStart == 0)
                {
                    glDrawArraysInstanced(_primitiveType, (int)vertexStart, vertexCount, instanceCount);
                }
                else
                {
                    glDrawArraysInstancedBaseInstance(_primitiveType, (int)vertexStart, vertexCount, instanceCount, instanceStart);
                }
            }
            CheckLastError();
        }

        public void DrawIndexed(uint indexCount, uint instanceCount, uint indexStart, int vertexOffset, uint instanceStart)
        {
            PreDrawCommand();
            PreDrawIndexedCommand();

            uint indexSize = _drawElementsType == DrawElementsType.UnsignedShort ? 2u : 4u;
            void* indices = (void*)((indexStart * indexSize) + _ibOffset);

            if (instanceCount == 1 && instanceStart == 0)
            {
                if (vertexOffset == 0)
                {
                    glDrawElements(_primitiveType, indexCount, _drawElementsType, indices);
                }
                else
                {
                    glDrawElementsBaseVertex(_primitiveType, indexCount, _drawElementsType, indices, vertexOffset);
                }
            }
            else
            {
                if (instanceStart > 0)
                {
                    glDrawElementsInstancedBaseVertexBaseInstance(
                        _primitiveType,
                        indexCount,
                        _drawElementsType,
                        indices,
                        instanceCount,
                        vertexOffset,
                        instanceStart);
                }
                else if (vertexOffset == 0)
                {
                    glDrawElementsInstanced(_primitiveType, indexCount, _drawElementsType, indices, instanceCount);
                }
                else
                {
                    glDrawElementsInstancedBaseVertex(
                        _primitiveType,
                        indexCount,
                        _drawElementsType,
                        indices,
                        instanceCount,
                        vertexOffset);
                }
            }
            CheckLastError();
        }

        public void DrawIndirect(DeviceBuffer indirectBuffer, uint offset, uint drawCount, uint stride)
        {
            PreDrawCommand();

            OpenGLBuffer glBuffer = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(indirectBuffer);
            glBindBuffer(BufferTarget.DrawIndirectBuffer, glBuffer.Buffer);
            CheckLastError();

            if (_extensions.MultiDrawIndirect)
            {
                glMultiDrawArraysIndirect(_primitiveType, (IntPtr)offset, drawCount, stride);
                CheckLastError();
            }
            else
            {
                uint indirect = offset;
                for (uint i = 0; i < drawCount; i++)
                {
                    glDrawArraysIndirect(_primitiveType, (IntPtr)indirect);
                    CheckLastError();

                    indirect += stride;
                }
            }
        }

        public void DrawIndexedIndirect(DeviceBuffer indirectBuffer, uint offset, uint drawCount, uint stride)
        {
            PreDrawCommand();
            PreDrawIndexedCommand();

            OpenGLBuffer glBuffer = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(indirectBuffer);
            glBindBuffer(BufferTarget.DrawIndirectBuffer, glBuffer.Buffer);
            CheckLastError();

            if (_extensions.MultiDrawIndirect)
            {
                glMultiDrawElementsIndirect(_primitiveType, _drawElementsType, (IntPtr)offset, drawCount, stride);
                CheckLastError();
            }
            else
            {
                uint indirect = offset;
                for (uint i = 0; i < drawCount; i++)
                {
                    glDrawElementsIndirect(_primitiveType, _drawElementsType, (IntPtr)indirect);
                    CheckLastError();

                    indirect += stride;
                }
            }
        }

        private void PreDrawCommand()
        {
            if (!_graphicsPipelineActive)
            {
                ActivateGraphicsPipeline();
            }

            if (!_graphicsResourcesFlushed)
            {
                FlushResourceSets(graphics: true);
                _graphicsResourcesFlushed = true;
            }

            if (!_vertexLayoutFlushed)
            {
                FlushVertexLayouts();
                _vertexLayoutFlushed = true;
            }
        }

        private void PreDrawIndexedCommand()
        {
            if (!_indexBufferBound)
            {
                BindIndexBuffer();
                _indexBufferBound = true;
            }
        }

        private void BindIndexBuffer()
        {
            Debug.Assert(_indexBuffer != null);

            glBindBuffer(BufferTarget.ElementArrayBuffer, _indexBuffer.Buffer);
            CheckLastError();
        }

        private void FlushResourceSets(bool graphics)
        {
            int setCount = graphics
                ? _graphicsPipeline!.ResourceLayouts.Length
                : _computePipeline!.ResourceLayouts.Length;

            Span<BoundResourceSetInfo> sets = (graphics ? _graphicsResourceSets : _computeResourceSets).AsSpan(0, setCount);
            Span<bool> newSets = (graphics ? _newGraphicsResourceSets : _newComputeResourceSets).AsSpan(0, setCount);

            for (int slot = 0; slot < setCount; slot++)
            {
                ref BoundResourceSetInfo brsi = ref sets[slot];
                OpenGLResourceSet glSet = Util.AssertSubtype<ResourceSet, OpenGLResourceSet>(brsi.Set);
                ResourceLayoutElementDescription[] layoutElements = glSet.Layout.Elements;
                bool isNew = newSets[slot];
                newSets[slot] = false;

                ActivateResourceSet((uint)slot, graphics, ref brsi, layoutElements, isNew);
            }
        }

        private void FlushVertexLayouts()
        {
            ReadOnlySpan<uint> strides = _graphicsPipeline!.VertexStrides;
            ReadOnlySpan<nuint> offsets = _vbOffsets;
            ReadOnlySpan<uint> buffers = _vertexBuffers;
            bool[] newVertexBuffers = _newVertexBuffers;

            uint totalSlotsBound = 0;
            VertexLayoutDescription[] layouts = _graphicsPipeline.VertexLayouts;
            bool separateBinding = _gd.Extensions.ARB_vertex_attrib_binding;

            for (int i = 0; i < layouts.Length; i++)
            {
                VertexLayoutDescription input = layouts[i];

                if (newVertexBuffers[i])
                {
                    if (separateBinding)
                    {
                        glBindVertexBuffer((uint)i, buffers[i], (nint)offsets[i], strides[i]);
                        CheckLastError();
                    }
                    else
                    {
                        BindVertexAttribPointers(input.Elements, totalSlotsBound, buffers[i], offsets[i], strides[i]);
                    }
                    newVertexBuffers[i] = false;
                }

                totalSlotsBound += (uint)input.Elements.Length;
            }
        }

        private static void BindVertexAttribPointers(
            VertexElementDescription[] elements, uint totalSlotsBound, uint vbBuffer, nuint vbOffset, uint stride)
        {
            glBindBuffer(BufferTarget.ArrayBuffer, vbBuffer);
            CheckLastError();

            uint offset = 0;

            for (uint slot = 0; slot < elements.Length; slot++)
            {
                ref readonly VertexElementDescription element = ref elements[slot];
                uint actualSlot = totalSlotsBound + slot;

                int elementCount = FormatHelpers.GetElementCount(element.Format);
                VertexAttribPointerType type = OpenGLFormats.VdToGLVertexAttribPointerType(
                    element.Format,
                    out bool normalized,
                    out bool isInteger);

                nuint actualOffset = element.Offset != 0 ? element.Offset : offset;
                actualOffset += vbOffset;

                if (isInteger && !normalized)
                {
                    glVertexAttribIPointer(
                        actualSlot,
                        elementCount,
                        type,
                        stride,
                        (void*)actualOffset);
                }
                else
                {
                    glVertexAttribPointer(
                        actualSlot,
                        elementCount,
                        type,
                        normalized,
                        stride,
                        (void*)actualOffset);
                }
                CheckLastError();

                offset += FormatSizeHelpers.GetSizeInBytes(element.Format);
            }
        }

        internal void Dispatch(uint groupCountX, uint groupCountY, uint groupCountZ)
        {
            PreDispatchCommand();

            glDispatchCompute(groupCountX, groupCountY, groupCountZ);
            CheckLastError();

            PostDispatchCommand();
        }

        public void DispatchIndirect(DeviceBuffer indirectBuffer, uint offset)
        {
            PreDispatchCommand();

            OpenGLBuffer glBuffer = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(indirectBuffer);
            glBindBuffer(BufferTarget.DrawIndirectBuffer, glBuffer.Buffer);
            CheckLastError();

            glDispatchComputeIndirect((IntPtr)offset);
            CheckLastError();

            PostDispatchCommand();
        }

        private void PreDispatchCommand()
        {
            if (!_computePipelineActive)
            {
                ActivateComputePipeline();
            }

            if (!_computeResourcesFlushed)
            {
                FlushResourceSets(graphics: false);
                _computeResourcesFlushed = true;
            }
        }

        private static void PostDispatchCommand()
        {
            // TODO: Smart barriers?
            glMemoryBarrier(MemoryBarrierFlags.AllBarrierBits);
            CheckLastError();
        }

        public void End()
        {
            if (_graphicsPipelineActive)
            {
                UnbindGraphicsPipeline();
            }

            if (_computePipelineActive)
            {
                UnbindComputePipeline();
            }
        }

        private void UnbindGraphicsPipeline()
        {
            // TODO: unbind buffers/resources?

            glBindVertexArray(0);
            CheckLastError2();

            _vertexLayoutFlushed = false;
            _indexBufferBound = false;
            _graphicsResourcesFlushed = false;
            _graphicsPipelineActive = false;
        }

        private void UnbindComputePipeline()
        {
            // TODO: unbind buffers/resources?

            _computeResourcesFlushed = false;
            _computePipelineActive = false;
        }

        public void SetFramebuffer(Framebuffer fb)
        {
            if (fb is OpenGLFramebuffer glFB)
            {
                if (_backend == GraphicsBackend.OpenGL || _extensions.EXT_sRGBWriteControl)
                {
                    glEnable(EnableCap.FramebufferSrgb);
                    CheckLastError();
                }

                glFB.EnsureResourcesCreated();

                glBindFramebuffer(FramebufferTarget.Framebuffer, glFB.Framebuffer);
                CheckLastError();
                _isSwapchainFB = false;
            }
            else if (fb is OpenGLSwapchainFramebuffer swapchainFB)
            {
                if ((_backend == GraphicsBackend.OpenGL || _extensions.EXT_sRGBWriteControl))
                {
                    if (swapchainFB.DisableSrgbConversion)
                    {
                        glDisable(EnableCap.FramebufferSrgb);
                    }
                    else
                    {
                        glEnable(EnableCap.FramebufferSrgb);
                    }
                    CheckLastError();
                }

                if (_platformInfo.SetSwapchainFramebuffer != null)
                {
                    _platformInfo.SetSwapchainFramebuffer();
                }
                else
                {
                    glBindFramebuffer(FramebufferTarget.Framebuffer, 0);
                    CheckLastError();
                }

                _isSwapchainFB = true;
            }
            else
            {
                throw new VeldridException("Invalid Framebuffer type: " + fb.GetType().Name);
            }

            _fb = fb;
        }

        public void SetIndexBuffer(DeviceBuffer ib, IndexFormat format, uint offset)
        {
            OpenGLBuffer glIB = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(ib);
            glIB.EnsureResourcesCreated();

            if (_gd.IsDebug)
            {
                _gd.ThrowIfMapped(glIB, 0);
            }

            _drawElementsType = OpenGLFormats.VdToGLDrawElementsType(format);

            if (_indexBuffer != glIB || _ibOffset != offset)
            {
                _indexBuffer = glIB;
                _ibOffset = offset;
                _indexBufferBound = false;
            }
        }

        public void SetPipeline(Pipeline pipeline)
        {
            if (!pipeline.IsComputePipeline && _graphicsPipeline != pipeline)
            {
                _graphicsPipeline = Util.AssertSubtype<Pipeline, OpenGLPipeline>(pipeline);
                ActivateGraphicsPipeline();
            }
            else if (pipeline.IsComputePipeline && _computePipeline != pipeline)
            {
                _computePipeline = Util.AssertSubtype<Pipeline, OpenGLPipeline>(pipeline);
                ActivateComputePipeline();
            }
        }

        private void ActivateGraphicsPipeline()
        {
            if (_computePipelineActive)
            {
                UnbindComputePipeline();
            }

            _graphicsPipeline!.EnsureResourcesCreated();

            Util.EnsureArrayMinimumSize(ref _graphicsResourceSets, (uint)_graphicsPipeline.ResourceLayouts.Length);
            Util.EnsureArrayMinimumSize(ref _newGraphicsResourceSets, (uint)_graphicsPipeline.ResourceLayouts.Length);

            // Force ResourceSets to be re-bound.
            for (int i = 0; i < _graphicsPipeline.ResourceLayouts.Length; i++)
            {
                _newGraphicsResourceSets[i] = true;
            }
            _graphicsResourcesFlushed = false;

            // Blend State

            BlendStateDescription blendState = _graphicsPipeline.BlendState;
            glBlendColor(blendState.BlendFactor.R, blendState.BlendFactor.G, blendState.BlendFactor.B, blendState.BlendFactor.A);
            CheckLastError();

            if (blendState.AlphaToCoverageEnabled)
            {
                glEnable(EnableCap.SampleAlphaToCoverage);
            }
            else
            {
                glDisable(EnableCap.SampleAlphaToCoverage);
            }
            CheckLastError();

            if (_features.IndependentBlend)
            {
                for (uint i = 0; i < blendState.AttachmentStates.Length; i++)
                {
                    BlendAttachmentDescription attachment = blendState.AttachmentStates[i];
                    ColorWriteMask colorMask = attachment.ColorWriteMask.GetOrDefault();

                    glColorMaski(
                        i,
                        (colorMask & ColorWriteMask.Red) == ColorWriteMask.Red,
                        (colorMask & ColorWriteMask.Green) == ColorWriteMask.Green,
                        (colorMask & ColorWriteMask.Blue) == ColorWriteMask.Blue,
                        (colorMask & ColorWriteMask.Alpha) == ColorWriteMask.Alpha);
                    CheckLastError();

                    if (!attachment.BlendEnabled)
                    {
                        glDisablei(EnableCap.Blend, i);
                        CheckLastError();
                    }
                    else
                    {
                        glEnablei(EnableCap.Blend, i);
                        CheckLastError();

                        glBlendFuncSeparatei(
                            i,
                            OpenGLFormats.VdToGLBlendFactorSrc(attachment.SourceColorFactor),
                            OpenGLFormats.VdToGLBlendFactorDest(attachment.DestinationColorFactor),
                            OpenGLFormats.VdToGLBlendFactorSrc(attachment.SourceAlphaFactor),
                            OpenGLFormats.VdToGLBlendFactorDest(attachment.DestinationAlphaFactor));
                        CheckLastError();

                        glBlendEquationSeparatei(
                            i,
                            OpenGLFormats.VdToGLBlendEquationMode(attachment.ColorFunction),
                            OpenGLFormats.VdToGLBlendEquationMode(attachment.AlphaFunction));
                        CheckLastError();
                    }
                }
            }
            else if (blendState.AttachmentStates.Length > 0)
            {
                BlendAttachmentDescription attachment = blendState.AttachmentStates[0];
                ColorWriteMask colorMask = attachment.ColorWriteMask.GetOrDefault();

                glColorMask(
                    (colorMask & ColorWriteMask.Red) == ColorWriteMask.Red,
                    (colorMask & ColorWriteMask.Green) == ColorWriteMask.Green,
                    (colorMask & ColorWriteMask.Blue) == ColorWriteMask.Blue,
                    (colorMask & ColorWriteMask.Alpha) == ColorWriteMask.Alpha);
                CheckLastError();

                if (!attachment.BlendEnabled)
                {
                    glDisable(EnableCap.Blend);
                    CheckLastError();
                }
                else
                {
                    glEnable(EnableCap.Blend);
                    CheckLastError();

                    glBlendFuncSeparate(
                        OpenGLFormats.VdToGLBlendFactorSrc(attachment.SourceColorFactor),
                        OpenGLFormats.VdToGLBlendFactorDest(attachment.DestinationColorFactor),
                        OpenGLFormats.VdToGLBlendFactorSrc(attachment.SourceAlphaFactor),
                        OpenGLFormats.VdToGLBlendFactorDest(attachment.DestinationAlphaFactor));
                    CheckLastError();

                    glBlendEquationSeparate(
                        OpenGLFormats.VdToGLBlendEquationMode(attachment.ColorFunction),
                        OpenGLFormats.VdToGLBlendEquationMode(attachment.AlphaFunction));
                    CheckLastError();
                }
            }

            // Depth Stencil State

            DepthStencilStateDescription dss = _graphicsPipeline.DepthStencilState;
            if (!dss.DepthTestEnabled)
            {
                glDisable(EnableCap.DepthTest);
                CheckLastError();
            }
            else
            {
                glEnable(EnableCap.DepthTest);
                CheckLastError();

                glDepthFunc(OpenGLFormats.VdToGLDepthFunction(dss.DepthComparison));
                CheckLastError();
            }

            glDepthMask(dss.DepthWriteEnabled);
            CheckLastError();

            if (dss.StencilTestEnabled)
            {
                glEnable(EnableCap.StencilTest);
                CheckLastError();

                glStencilFuncSeparate(
                    CullFaceMode.Front,
                    OpenGLFormats.VdToGLStencilFunction(dss.StencilFront.Comparison),
                    (int)dss.StencilReference,
                    dss.StencilReadMask);
                CheckLastError();

                glStencilOpSeparate(
                    CullFaceMode.Front,
                    OpenGLFormats.VdToGLStencilOp(dss.StencilFront.Fail),
                    OpenGLFormats.VdToGLStencilOp(dss.StencilFront.DepthFail),
                    OpenGLFormats.VdToGLStencilOp(dss.StencilFront.Pass));
                CheckLastError();

                glStencilFuncSeparate(
                    CullFaceMode.Back,
                    OpenGLFormats.VdToGLStencilFunction(dss.StencilBack.Comparison),
                    (int)dss.StencilReference,
                    dss.StencilReadMask);
                CheckLastError();

                glStencilOpSeparate(
                    CullFaceMode.Back,
                    OpenGLFormats.VdToGLStencilOp(dss.StencilBack.Fail),
                    OpenGLFormats.VdToGLStencilOp(dss.StencilBack.DepthFail),
                    OpenGLFormats.VdToGLStencilOp(dss.StencilBack.Pass));
                CheckLastError();

                glStencilMask(dss.StencilWriteMask);
                CheckLastError();
            }
            else
            {
                glDisable(EnableCap.StencilTest);
                CheckLastError();
            }

            // Rasterizer State

            RasterizerStateDescription rs = _graphicsPipeline.RasterizerState;
            if (rs.CullMode == FaceCullMode.None)
            {
                glDisable(EnableCap.CullFace);
                CheckLastError();
            }
            else
            {
                glEnable(EnableCap.CullFace);
                CheckLastError();

                glCullFace(OpenGLFormats.VdToGLCullFaceMode(rs.CullMode));
                CheckLastError();
            }

            if (_backend == GraphicsBackend.OpenGL)
            {
                glPolygonMode(MaterialFace.FrontAndBack, OpenGLFormats.VdToGLPolygonMode(rs.FillMode));
                CheckLastError();
            }

            if (!rs.ScissorTestEnabled)
            {
                glDisable(EnableCap.ScissorTest);
            }
            else
            {
                glEnable(EnableCap.ScissorTest);
            }
            CheckLastError();

            if (_backend == GraphicsBackend.OpenGL)
            {
                if (!rs.DepthClipEnabled)
                {
                    glEnable(EnableCap.DepthClamp);
                }
                else
                {
                    glDisable(EnableCap.DepthClamp);
                }
                CheckLastError();
            }

            glFrontFace(OpenGLFormats.VdToGLFrontFaceDirection(rs.FrontFace));
            CheckLastError();

            // Primitive Topology
            _primitiveType = OpenGLFormats.VdToGLPrimitiveType(_graphicsPipeline.PrimitiveTopology);

            // Shader Set
            glUseProgram(_graphicsPipeline.Program);
            CheckLastError();

            // Force vertex buffers to be re-bound.
            for (int i = 0; i < _newVertexBuffers.Length; i++)
            {
                _newVertexBuffers[i] = true;
            }

            int vertexStridesCount = _graphicsPipeline.VertexStrides.Length;
            Util.EnsureArrayMinimumSize(ref _vertexBuffers, (uint)vertexStridesCount);

            glBindVertexArray(_graphicsPipeline.Vao);
            CheckLastError();

            _vertexLayoutFlushed = false;
            _indexBufferBound = false;
            _graphicsPipelineActive = true;
        }

        public void GenerateMipmaps(Texture texture)
        {
            OpenGLTexture glTex = Util.AssertSubtype<Texture, OpenGLTexture>(texture);
            glTex.EnsureResourcesCreated();
            if (_extensions.ARB_DirectStateAccess)
            {
                glGenerateTextureMipmap(glTex.Texture);
            }
            else
            {
                TextureTarget target = glTex.TextureTarget;
                _textureSamplerManager.SetTextureTransient(target, glTex.Texture);
                glGenerateMipmap(target);
            }
            CheckLastError();
        }

        [SkipLocalsInit]
        public void PushDebugGroup(ReadOnlySpan<char> name)
        {
            Span<byte> byteBuffer = stackalloc byte[1024];
            bool khr_debug = _extensions.KHR_Debug;
            bool ext_debugMarker = _extensions.EXT_DebugMarker;

            if (!khr_debug && !ext_debugMarker)
                return;

            int length = Util.GetNullTerminatedUtf8(name, ref byteBuffer);
            fixed (byte* utf8Ptr = byteBuffer)
            {
                if (khr_debug)
                {
                    glPushDebugGroup(DebugSource.DebugSourceApplication, 0, (uint)length, utf8Ptr);
                    CheckLastError();
                }
                else if (ext_debugMarker)
                {
                    glPushGroupMarker((uint)length, utf8Ptr);
                    CheckLastError();
                }
            }
        }

        public void PopDebugGroup()
        {
            if (_extensions.KHR_Debug)
            {
                glPopDebugGroup();
                CheckLastError();
            }
            else if (_extensions.EXT_DebugMarker)
            {
                glPopGroupMarker();
                CheckLastError();
            }
        }

        [SkipLocalsInit]
        public void InsertDebugMarker(string name)
        {
            Span<byte> byteBuffer = stackalloc byte[1024];
            bool khr_debug = _extensions.KHR_Debug;
            bool ext_debugMarker = _extensions.EXT_DebugMarker;

            if (!khr_debug && !ext_debugMarker)
                return;

            int length = Util.GetNullTerminatedUtf8(name, ref byteBuffer);
            fixed (byte* utf8Ptr = byteBuffer)
            {
                if (khr_debug)
                {
                    glDebugMessageInsert(
                        DebugSource.DebugSourceApplication,
                        DebugType.DebugTypeMarker,
                        0,
                        DebugSeverity.DebugSeverityNotification,
                        (uint)length,
                        utf8Ptr);
                    CheckLastError();
                }
                else if (ext_debugMarker)
                {
                    glInsertEventMarker((uint)length, utf8Ptr);
                    CheckLastError();
                }
            }
        }

        public void InsertFence(OpenGLFence fence)
        {
            fence.Set();
        }

        private void ActivateComputePipeline()
        {
            if (_graphicsPipelineActive)
            {
                UnbindGraphicsPipeline();
            }

            _computePipeline!.EnsureResourcesCreated();
            Util.EnsureArrayMinimumSize(ref _computeResourceSets, (uint)_computePipeline.ResourceLayouts.Length);
            Util.EnsureArrayMinimumSize(ref _newComputeResourceSets, (uint)_computePipeline.ResourceLayouts.Length);

            // Force ResourceSets to be re-bound.
            for (int i = 0; i < _computePipeline.ResourceLayouts.Length; i++)
            {
                _newComputeResourceSets[i] = true;
            }
            _computeResourcesFlushed = false;

            // Shader Set
            glUseProgram(_computePipeline.Program);
            CheckLastError();

            _computePipelineActive = true;
        }

        public void SetGraphicsResourceSet(uint slot, ResourceSet rs, ReadOnlySpan<uint> dynamicOffsets)
        {
            ref BoundResourceSetInfo set = ref _graphicsResourceSets[slot];
            if (!set.Equals(rs, dynamicOffsets))
            {
                set.Offsets.Dispose();
                set = new BoundResourceSetInfo(rs, dynamicOffsets);
                _newGraphicsResourceSets[slot] = true;
                _graphicsResourcesFlushed = false;
            }
        }

        public void SetComputeResourceSet(uint slot, ResourceSet rs, ReadOnlySpan<uint> dynamicOffsets)
        {
            ref BoundResourceSetInfo set = ref _computeResourceSets[slot];
            if (!set.Equals(rs, dynamicOffsets))
            {
                set.Offsets.Dispose();
                set = new BoundResourceSetInfo(rs, dynamicOffsets);
                _newComputeResourceSets[slot] = true;
                _computeResourcesFlushed = false;
            }
        }

        private void ActivateResourceSet(
            uint slot,
            bool graphics,
            ref BoundResourceSetInfo brsi,
            ResourceLayoutElementDescription[] layoutElements,
            bool isNew)
        {
            OpenGLResourceSet glResourceSet = Util.AssertSubtype<ResourceSet, OpenGLResourceSet>(brsi.Set);
            OpenGLPipeline pipeline = graphics ? _graphicsPipeline! : _computePipeline!;
            uint ubBaseIndex = GetUniformBaseIndex(slot, graphics);
            uint ssboBaseIndex = GetShaderStorageBaseIndex(slot, graphics);

            uint ubOffset = 0;
            uint ssboOffset = 0;
            uint dynamicOffsetIndex = 0;
            for (uint element = 0; element < glResourceSet.Resources.Length; element++)
            {
                ResourceKind kind = layoutElements[element].Kind;
                BindableResource resource = glResourceSet.Resources[(int)element];

                uint bufferOffset = 0;
                if (glResourceSet.Layout.IsDynamicBuffer(element))
                {
                    bufferOffset = brsi.Offsets.Get(dynamicOffsetIndex);
                    dynamicOffsetIndex += 1;
                }

                switch (kind)
                {
                    case ResourceKind.UniformBuffer:
                    {
                        if (!isNew)
                        {
                            continue;
                        }

                        DeviceBufferRange range = Util.GetBufferRange(resource, bufferOffset);
                        OpenGLBuffer glUB = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(range.Buffer);

                        glUB.EnsureResourcesCreated();
                        if (pipeline.GetUniformBindingForSlot(slot, element, out OpenGLUniformBinding uniformBindingInfo))
                        {
                            if (range.SizeInBytes < uniformBindingInfo.BlockSize)
                            {
                                string name = glResourceSet.Layout.Elements[element].Name;
                                throw new VeldridException(
                                    $"Not enough data in uniform buffer \"{name}\" (slot {slot}, element {element}). Shader expects at least {uniformBindingInfo.BlockSize} bytes, but buffer only contains {range.SizeInBytes} bytes");
                            }
                            glUniformBlockBinding(pipeline.Program, uniformBindingInfo.BlockLocation, ubBaseIndex + ubOffset);
                            CheckLastError();

                            glBindBufferRange(
                                BufferRangeTarget.UniformBuffer,
                                ubBaseIndex + ubOffset,
                                glUB.Buffer,
                                (IntPtr)range.Offset,
                                (UIntPtr)range.SizeInBytes);
                            CheckLastError();

                            ubOffset += 1;
                        }
                        break;
                    }

                    case ResourceKind.StructuredBufferReadWrite:
                    case ResourceKind.StructuredBufferReadOnly:
                    {
                        if (!isNew)
                        {
                            continue;
                        }

                        DeviceBufferRange range = Util.GetBufferRange(resource, bufferOffset);
                        OpenGLBuffer glBuffer = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(range.Buffer);

                        glBuffer.EnsureResourcesCreated();
                        if (pipeline.GetStorageBufferBindingForSlot(slot, element, out OpenGLShaderStorageBinding shaderStorageBinding))
                        {
                            if (_backend == GraphicsBackend.OpenGL)
                            {
                                glShaderStorageBlockBinding(
                                    pipeline.Program,
                                    shaderStorageBinding.StorageBlockBinding,
                                    ssboBaseIndex + ssboOffset);
                                CheckLastError();

                                glBindBufferRange(
                                    BufferRangeTarget.ShaderStorageBuffer,
                                    ssboBaseIndex + ssboOffset,
                                    glBuffer.Buffer,
                                    (IntPtr)range.Offset,
                                    (UIntPtr)range.SizeInBytes);
                            }
                            else
                            {
                                glBindBufferRange(
                                    BufferRangeTarget.ShaderStorageBuffer,
                                    shaderStorageBinding.StorageBlockBinding,
                                    glBuffer.Buffer,
                                    (IntPtr)range.Offset,
                                    (UIntPtr)range.SizeInBytes);
                            }
                            CheckLastError();
                            ssboOffset += 1;
                        }
                        break;
                    }

                    case ResourceKind.TextureReadOnly:
                        TextureView texView = Util.GetTextureView(_gd, resource);
                        OpenGLTextureView glTexView = Util.AssertSubtype<TextureView, OpenGLTextureView>(texView);
                        glTexView.EnsureResourcesCreated();
                        if (pipeline.GetTextureBindingInfo(slot, element, out OpenGLTextureBindingSlotInfo textureBindingInfo))
                        {
                            _textureSamplerManager.SetTexture((uint)textureBindingInfo.RelativeIndex, glTexView);
                            glUniform1i(textureBindingInfo.UniformLocation, textureBindingInfo.RelativeIndex);
                            CheckLastError();
                        }
                        break;

                    case ResourceKind.TextureReadWrite:
                        TextureView texViewRW = Util.GetTextureView(_gd, resource);
                        OpenGLTextureView glTexViewRW = Util.AssertSubtype<TextureView, OpenGLTextureView>(texViewRW);
                        glTexViewRW.EnsureResourcesCreated();
                        if (pipeline.GetTextureBindingInfo(slot, element, out OpenGLTextureBindingSlotInfo imageBindingInfo))
                        {
                            bool layered = (texViewRW.Target.Usage & TextureUsage.Cubemap) != 0 || texViewRW.ArrayLayers > 1;

                            if (layered && (texViewRW.BaseArrayLayer > 0
                                || (texViewRW.ArrayLayers > 1 && texViewRW.ArrayLayers < texViewRW.Target.ArrayLayers)))
                            {
                                throw new VeldridException(
                                    "Cannot bind texture with BaseArrayLayer > 0 and ArrayLayers > 1, or with an incomplete set of array layers (cubemaps have ArrayLayers == 6 implicitly).");
                            }

                            if (_backend == GraphicsBackend.OpenGL)
                            {
                                glBindImageTexture(
                                    (uint)imageBindingInfo.RelativeIndex,
                                    glTexViewRW.Target.Texture,
                                    (int)texViewRW.BaseMipLevel,
                                    layered,
                                    (int)texViewRW.BaseArrayLayer,
                                    TextureAccess.ReadWrite,
                                    glTexViewRW.GetReadWriteSizedInternalFormat());
                                CheckLastError();

                                glUniform1i(imageBindingInfo.UniformLocation, imageBindingInfo.RelativeIndex);
                            }
                            else
                            {
                                glBindImageTexture(
                                    (uint)imageBindingInfo.RelativeIndex,
                                    glTexViewRW.Target.Texture,
                                    (int)texViewRW.BaseMipLevel,
                                    layered,
                                    (int)texViewRW.BaseArrayLayer,
                                    TextureAccess.ReadWrite,
                                    glTexViewRW.GetReadWriteSizedInternalFormat());
                            }
                            CheckLastError();
                        }
                        break;

                    case ResourceKind.Sampler:
                        OpenGLSampler glSampler = Util.AssertSubtype<BindableResource, OpenGLSampler>(resource);
                        glSampler.EnsureResourcesCreated();

                        if (pipeline.GetSamplerBindingInfo(slot, element, out OpenGLSamplerBindingSlotInfo samplerBindingInfo))
                        {
                            foreach (int index in samplerBindingInfo.RelativeIndices)
                            {
                                _textureSamplerManager.SetSampler((uint)index, glSampler);
                            }
                        }
                        break;

                    default:
                        throw Illegal.Value<ResourceKind>();
                }
            }
        }

        public void ResolveTexture(Texture source, Texture destination)
        {
            OpenGLTexture glSourceTex = Util.AssertSubtype<Texture, OpenGLTexture>(source);
            OpenGLTexture glDestinationTex = Util.AssertSubtype<Texture, OpenGLTexture>(destination);
            glSourceTex.EnsureResourcesCreated();
            glDestinationTex.EnsureResourcesCreated();

            uint sourceFramebuffer = glSourceTex.GetFramebuffer(0, 0);
            uint destinationFramebuffer = glDestinationTex.GetFramebuffer(0, 0);

            glBindFramebuffer(FramebufferTarget.ReadFramebuffer, sourceFramebuffer);
            CheckLastError();

            glBindFramebuffer(FramebufferTarget.DrawFramebuffer, destinationFramebuffer);
            CheckLastError();

            glDisable(EnableCap.ScissorTest);
            CheckLastError();

            glBlitFramebuffer(
                0,
                0,
                (int)source.Width,
                (int)source.Height,
                0,
                0,
                (int)destination.Width,
                (int)destination.Height,
                ClearBufferMask.ColorBufferBit,
                BlitFramebufferFilter.Nearest);
            CheckLastError();
        }

        private uint GetUniformBaseIndex(uint slot, bool graphics)
        {
            OpenGLPipeline pipeline = graphics ? _graphicsPipeline! : _computePipeline!;
            uint ret = 0;
            for (uint i = 0; i < slot; i++)
            {
                ret += pipeline.GetUniformBufferCount(i);
            }

            return ret;
        }

        private uint GetShaderStorageBaseIndex(uint slot, bool graphics)
        {
            OpenGLPipeline pipeline = graphics ? _graphicsPipeline! : _computePipeline!;
            uint ret = 0;
            for (uint i = 0; i < slot; i++)
            {
                ret += pipeline.GetShaderStorageBufferCount(i);
            }

            return ret;
        }

        public void SetScissorRect(uint index, uint x, uint y, uint width, uint height)
        {
            if (_backend == GraphicsBackend.OpenGL)
            {
                glScissorIndexed(
                    index,
                    (int)x,
                    (int)(_fb!.Height - (int)height - y),
                    width,
                    height);
                CheckLastError();
            }
            else
            {
                if (index == 0)
                {
                    glScissor(
                        (int)x,
                        (int)(_fb!.Height - (int)height - y),
                        width,
                        height);
                    CheckLastError();
                }
            }
        }

        public void SetVertexBuffer(uint index, DeviceBuffer vb, uint offset)
        {
            OpenGLBuffer glVB = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(vb);
            glVB.EnsureResourcesCreated();

            if (_gd.IsDebug)
            {
                _gd.ThrowIfMapped(glVB, 0);
            }

            Util.EnsureArrayMinimumSize(ref _vertexBuffers, index + 1);
            Util.EnsureArrayMinimumSize(ref _vbOffsets, index + 1);
            Util.EnsureArrayMinimumSize(ref _newVertexBuffers, index + 1);

            uint buffer = glVB.Buffer;
            if (_vertexBuffers[index] != buffer || _vbOffsets[index] != offset)
            {
                _vertexBuffers[index] = buffer;
                _vbOffsets[index] = offset;
                _newVertexBuffers[index] = true;
                _vertexLayoutFlushed = false;
            }
        }

        public void SetViewport(uint index, ref Viewport viewport)
        {
            if (_backend == GraphicsBackend.OpenGL)
            {
                float left = viewport.X;
                float bottom = _fb!.Height - (viewport.Y + viewport.Height);

                glViewportIndexed(index, left, bottom, viewport.Width, viewport.Height);
                CheckLastError();

                glDepthRangeIndexed(index, viewport.MinDepth, viewport.MaxDepth);
                CheckLastError();
            }
            else
            {
                if (index == 0)
                {
                    glViewport((int)viewport.X, (int)viewport.Y, (uint)viewport.Width, (uint)viewport.Height);
                    CheckLastError();

                    glDepthRangef(viewport.MinDepth, viewport.MaxDepth);
                    CheckLastError();
                }
            }
        }

        public void UpdateBuffer(DeviceBuffer buffer, uint bufferOffsetInBytes, IntPtr dataPtr, uint sizeInBytes)
        {
            OpenGLBuffer glBuffer = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(buffer);
            glBuffer.EnsureResourcesCreated();

            if (glBuffer.CanBufferSubData)
            {
                if (_extensions.ARB_DirectStateAccess)
                {
                    glNamedBufferSubData(
                        glBuffer.Buffer,
                        (IntPtr)bufferOffsetInBytes,
                        sizeInBytes,
                        dataPtr.ToPointer());
                }
                else
                {
                    BufferTarget bufferTarget = BufferTarget.CopyWriteBuffer;
                    glBindBuffer(bufferTarget, glBuffer.Buffer);
                    CheckLastError();

                    glBufferSubData(
                        bufferTarget,
                        (IntPtr)bufferOffsetInBytes,
                        (UIntPtr)sizeInBytes,
                        dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else
            {
                uint tmpBuffer = 0;
                if (_extensions.ARB_DirectStateAccess)
                {
                    glCreateBuffers(1, &tmpBuffer);
                    CheckLastError();

                    glNamedBufferData(
                        tmpBuffer,
                        sizeInBytes,
                        dataPtr.ToPointer(),
                        BufferUsageHint.StreamCopy);
                }
                else
                {
                    glGenBuffers(1, &tmpBuffer);
                    CheckLastError();

                    BufferTarget bufferTarget = BufferTarget.CopyWriteBuffer;
                    glBindBuffer(bufferTarget, tmpBuffer);
                    CheckLastError();

                    glBufferData(
                        bufferTarget,
                        sizeInBytes,
                        dataPtr.ToPointer(),
                        BufferUsageHint.StreamCopy);
                }
                CheckLastError();

                CopyBufferCore(tmpBuffer, 0, glBuffer.Buffer, bufferOffsetInBytes, sizeInBytes);

                glDeleteBuffers(1, &tmpBuffer);
                CheckLastError();
            }
        }

        public void UpdateTexture(
            Texture texture,
            IntPtr dataPtr,
            uint x,
            uint y,
            uint z,
            uint width,
            uint height,
            uint depth,
            uint mipLevel,
            uint arrayLayer)
        {
            if (width == 0 || height == 0 || depth == 0)
            {
                return;
            }

            OpenGLTexture glTex = Util.AssertSubtype<Texture, OpenGLTexture>(texture);
            glTex.EnsureResourcesCreated();

            TextureTarget texTarget = glTex.TextureTarget;

            _textureSamplerManager.SetTextureTransient(texTarget, glTex.Texture);
            CheckLastError();

            bool isCompressed = FormatHelpers.IsCompressedFormat(texture.Format);
            uint blockSize = isCompressed ? 4u : 1u;

            uint blockAlignedWidth = Math.Max(width, blockSize);
            uint blockAlignedHeight = Math.Max(height, blockSize);

            uint rowPitch = FormatHelpers.GetRowPitch(blockAlignedWidth, texture.Format);
            uint depthPitch = FormatHelpers.GetDepthPitch(rowPitch, blockAlignedHeight, texture.Format);

            // Compressed textures can specify regions that are larger than the dimensions.
            // We should only pass up to the dimensions to OpenGL, though.
            Util.GetMipDimensions(glTex, mipLevel, out uint mipWidth, out uint mipHeight);
            width = Math.Min(width, mipWidth);
            height = Math.Min(height, mipHeight);

            uint unpackAlignment = 4;
            if (!isCompressed)
            {
                unpackAlignment = FormatSizeHelpers.GetSizeInBytes(glTex.Format);
            }
            if (unpackAlignment < 4)
            {
                glPixelStorei(PixelStoreParameter.UnpackAlignment, (int)unpackAlignment);
                CheckLastError();
            }

            if (texTarget == TextureTarget.Texture1D)
            {
                if (isCompressed)
                {
                    glCompressedTexSubImage1D(
                        TextureTarget.Texture1D,
                        (int)mipLevel,
                        (int)x,
                        width,
                        glTex.GLInternalFormat,
                        rowPitch,
                        dataPtr.ToPointer());
                }
                else
                {
                    glTexSubImage1D(
                        TextureTarget.Texture1D,
                        (int)mipLevel,
                        (int)x,
                        width,
                        glTex.GLPixelFormat,
                        glTex.GLPixelType,
                        dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else if (texTarget == TextureTarget.Texture1DArray)
            {
                if (isCompressed)
                {
                    glCompressedTexSubImage2D(
                        TextureTarget.Texture1DArray,
                        (int)mipLevel,
                        (int)x,
                        (int)arrayLayer,
                        width,
                        1,
                        glTex.GLInternalFormat,
                        rowPitch,
                        dataPtr.ToPointer());
                }
                else
                {
                    glTexSubImage2D(
                    TextureTarget.Texture1DArray,
                    (int)mipLevel,
                    (int)x,
                    (int)arrayLayer,
                    width,
                    1,
                    glTex.GLPixelFormat,
                    glTex.GLPixelType,
                    dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else if (texTarget == TextureTarget.Texture2D)
            {
                if (isCompressed)
                {
                    glCompressedTexSubImage2D(
                        TextureTarget.Texture2D,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        width,
                        height,
                        glTex.GLInternalFormat,
                        depthPitch,
                        dataPtr.ToPointer());
                }
                else
                {
                    glTexSubImage2D(
                        TextureTarget.Texture2D,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        width,
                        height,
                        glTex.GLPixelFormat,
                        glTex.GLPixelType,
                        dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else if (texTarget == TextureTarget.Texture2DArray)
            {
                if (isCompressed)
                {
                    glCompressedTexSubImage3D(
                        TextureTarget.Texture2DArray,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        (int)arrayLayer,
                        width,
                        height,
                        1,
                        glTex.GLInternalFormat,
                        depthPitch,
                        dataPtr.ToPointer());
                }
                else
                {
                    glTexSubImage3D(
                        TextureTarget.Texture2DArray,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        (int)arrayLayer,
                        width,
                        height,
                        1,
                        glTex.GLPixelFormat,
                        glTex.GLPixelType,
                        dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else if (texTarget == TextureTarget.Texture3D)
            {
                if (isCompressed)
                {
                    glCompressedTexSubImage3D(
                        TextureTarget.Texture3D,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        (int)z,
                        width,
                        height,
                        depth,
                        glTex.GLInternalFormat,
                        depthPitch * depth,
                        dataPtr.ToPointer());
                }
                else
                {
                    glTexSubImage3D(
                        TextureTarget.Texture3D,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        (int)z,
                        width,
                        height,
                        depth,
                        glTex.GLPixelFormat,
                        glTex.GLPixelType,
                        dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else if (texTarget == TextureTarget.TextureCubeMap)
            {
                TextureTarget cubeTarget = GetCubeTarget(arrayLayer);
                if (isCompressed)
                {
                    glCompressedTexSubImage2D(
                        cubeTarget,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        width,
                        height,
                        glTex.GLInternalFormat,
                        depthPitch,
                        dataPtr.ToPointer());
                }
                else
                {
                    glTexSubImage2D(
                        cubeTarget,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        width,
                        height,
                        glTex.GLPixelFormat,
                        glTex.GLPixelType,
                        dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else if (texTarget == TextureTarget.TextureCubeMapArray)
            {
                if (isCompressed)
                {
                    glCompressedTexSubImage3D(
                        TextureTarget.TextureCubeMapArray,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        (int)arrayLayer,
                        width,
                        height,
                        1,
                        glTex.GLInternalFormat,
                        depthPitch,
                        dataPtr.ToPointer());
                }
                else
                {
                    glTexSubImage3D(
                        TextureTarget.TextureCubeMapArray,
                        (int)mipLevel,
                        (int)x,
                        (int)y,
                        (int)arrayLayer,
                        width,
                        height,
                        1,
                        glTex.GLPixelFormat,
                        glTex.GLPixelType,
                        dataPtr.ToPointer());
                }
                CheckLastError();
            }
            else
            {
                throw new VeldridException($"Invalid OpenGL TextureTarget encountered: {glTex.TextureTarget}.");
            }

            if (unpackAlignment < 4)
            {
                glPixelStorei(PixelStoreParameter.UnpackAlignment, 4);
                CheckLastError();
            }
        }

        private static TextureTarget GetCubeTarget(uint arrayLayer)
        {
            return arrayLayer switch
            {
                0 => TextureTarget.TextureCubeMapPositiveX,
                1 => TextureTarget.TextureCubeMapNegativeX,
                2 => TextureTarget.TextureCubeMapPositiveY,
                3 => TextureTarget.TextureCubeMapNegativeY,
                4 => TextureTarget.TextureCubeMapPositiveZ,
                5 => TextureTarget.TextureCubeMapNegativeZ,
                _ => throw new VeldridException("Unexpected array layer in UpdateTexture called on a cubemap texture."),
            };
        }

        public void CopyBuffer(DeviceBuffer source, uint sourceOffset, DeviceBuffer destination, uint destinationOffset, uint sizeInBytes)
        {
            OpenGLBuffer srcGLBuffer = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(source);
            OpenGLBuffer dstGLBuffer = Util.AssertSubtype<DeviceBuffer, OpenGLBuffer>(destination);

            srcGLBuffer.EnsureResourcesCreated();
            dstGLBuffer.EnsureResourcesCreated();

            CopyBufferCore(srcGLBuffer.Buffer, sourceOffset, dstGLBuffer.Buffer, destinationOffset, sizeInBytes);
        }

        private void CopyBufferCore(uint srcBuffer, uint sourceOffset, uint dstBuffer, uint destinationOffset, uint sizeInBytes)
        {
            if (_extensions.ARB_DirectStateAccess)
            {
                glCopyNamedBufferSubData(
                    srcBuffer,
                    dstBuffer,
                    (IntPtr)sourceOffset,
                    (IntPtr)destinationOffset,
                    sizeInBytes);
            }
            else
            {
                glBindBuffer(BufferTarget.CopyReadBuffer, srcBuffer);
                CheckLastError();

                glBindBuffer(BufferTarget.CopyWriteBuffer, dstBuffer);
                CheckLastError();

                glCopyBufferSubData(
                    BufferTarget.CopyReadBuffer,
                    BufferTarget.CopyWriteBuffer,
                    (IntPtr)sourceOffset,
                    (IntPtr)destinationOffset,
                    (IntPtr)sizeInBytes);
            }
            CheckLastError();
        }

        public void CopyTexture(
            Texture source,
            uint srcX, uint srcY, uint srcZ,
            uint srcMipLevel,
            uint srcBaseArrayLayer,
            Texture destination,
            uint dstX, uint dstY, uint dstZ,
            uint dstMipLevel,
            uint dstBaseArrayLayer,
            uint width, uint height, uint depth,
            uint layerCount)
        {
            OpenGLTexture srcGLTexture = Util.AssertSubtype<Texture, OpenGLTexture>(source);
            OpenGLTexture dstGLTexture = Util.AssertSubtype<Texture, OpenGLTexture>(destination);

            srcGLTexture.EnsureResourcesCreated();
            dstGLTexture.EnsureResourcesCreated();

            if (_extensions.CopyImage && depth == 1)
            {
                // glCopyImageSubData does not work properly when depth > 1, so use the awful roundabout copy.
                uint srcZOrLayer = Math.Max(srcBaseArrayLayer, srcZ);
                uint dstZOrLayer = Math.Max(dstBaseArrayLayer, dstZ);
                uint depthOrLayerCount = Math.Max(depth, layerCount);
                // Copy width and height are allowed to be a full compressed block size, even if the mip level only contains a
                // region smaller than the block size.
                Util.GetMipDimensions(source, srcMipLevel, out uint mipWidth, out uint mipHeight);
                width = Math.Min(width, mipWidth);
                height = Math.Min(height, mipHeight);

                glCopyImageSubData(
                    srcGLTexture.Texture, srcGLTexture.TextureTarget, (int)srcMipLevel, (int)srcX, (int)srcY, (int)srcZOrLayer,
                    dstGLTexture.Texture, dstGLTexture.TextureTarget, (int)dstMipLevel, (int)dstX, (int)dstY, (int)dstZOrLayer,
                    width, height, depthOrLayerCount);
                CheckLastError();
            }
            else
            {
                for (uint layer = 0; layer < layerCount; layer++)
                {
                    uint srcLayer = layer + srcBaseArrayLayer;
                    uint dstLayer = layer + dstBaseArrayLayer;
                    CopyRoundabout(
                        srcGLTexture, dstGLTexture,
                        srcX, srcY, srcZ, srcMipLevel, srcLayer,
                        dstX, dstY, dstZ, dstMipLevel, dstLayer,
                        width, height, depth);
                }
            }
        }

        private void CopyRoundabout(
            OpenGLTexture srcGLTexture, OpenGLTexture dstGLTexture,
            uint srcX, uint srcY, uint srcZ, uint srcMipLevel, uint srcLayer,
            uint dstX, uint dstY, uint dstZ, uint dstMipLevel, uint dstLayer,
            uint width, uint height, uint depth)
        {
            bool isCompressed = FormatHelpers.IsCompressedFormat(srcGLTexture.Format);
            if (srcGLTexture.Format != dstGLTexture.Format)
            {
                throw new VeldridException("Copying to/from Textures with different formats is not supported.");
            }

            uint packAlignment = 4;
            uint depthSliceSize = 0;
            uint sizeInBytes;
            TextureTarget srcTarget = srcGLTexture.TextureTarget;
            if (isCompressed)
            {
                _textureSamplerManager.SetTextureTransient(srcTarget, srcGLTexture.Texture);
                CheckLastError();

                int compressedSize;
                glGetTexLevelParameteriv(
                    srcTarget,
                    (int)srcMipLevel,
                    GetTextureParameter.TextureCompressedImageSize,
                    &compressedSize);
                CheckLastError();
                sizeInBytes = (uint)compressedSize;
            }
            else
            {
                uint pixelSize = FormatSizeHelpers.GetSizeInBytes(srcGLTexture.Format);
                packAlignment = pixelSize;
                depthSliceSize = width * height * pixelSize;
                sizeInBytes = depthSliceSize * depth;
            }

            StagingBlock block = _stagingMemoryPool.GetStagingBlock(sizeInBytes);

            if (packAlignment < 4)
            {
                glPixelStorei(PixelStoreParameter.PackAlignment, (int)packAlignment);
                CheckLastError();
            }

            if (isCompressed)
            {
                if (_extensions.ARB_DirectStateAccess)
                {
                    glGetCompressedTextureImage(
                        srcGLTexture.Texture,
                        (int)srcMipLevel,
                        block.SizeInBytes,
                        block.Data);
                }
                else
                {
                    _textureSamplerManager.SetTextureTransient(srcTarget, srcGLTexture.Texture);
                    CheckLastError();

                    glGetCompressedTexImage(srcTarget, (int)srcMipLevel, block.Data);
                }
                CheckLastError();

                TextureTarget dstTarget = dstGLTexture.TextureTarget;
                _textureSamplerManager.SetTextureTransient(dstTarget, dstGLTexture.Texture);
                CheckLastError();

                Util.GetMipDimensions(srcGLTexture, srcMipLevel, out uint mipWidth, out uint mipHeight, out _);
                uint fullRowPitch = FormatHelpers.GetRowPitch(mipWidth, srcGLTexture.Format);
                uint fullDepthPitch = FormatHelpers.GetDepthPitch(
                    fullRowPitch,
                    mipHeight,
                    srcGLTexture.Format);

                uint denseRowPitch = FormatHelpers.GetRowPitch(width, srcGLTexture.Format);
                uint denseDepthPitch = FormatHelpers.GetDepthPitch(denseRowPitch, height, srcGLTexture.Format);
                uint numRows = FormatHelpers.GetNumRows(height, srcGLTexture.Format);
                uint trueCopySize = denseRowPitch * numRows;
                StagingBlock trueCopySrc = _stagingMemoryPool.GetStagingBlock(trueCopySize);

                uint layerStartOffset = denseDepthPitch * srcLayer;

                Util.CopyTextureRegion(
                    (byte*)block.Data + layerStartOffset,
                    srcX, srcY, srcZ,
                    fullRowPitch, fullDepthPitch,
                    trueCopySrc.Data,
                    0, 0, 0,
                    denseRowPitch,
                    denseDepthPitch,
                    width, height, depth,
                    srcGLTexture.Format);

                UpdateTexture(
                    dstGLTexture,
                    (IntPtr)trueCopySrc.Data,
                    dstX, dstY, dstZ,
                    width, height, 1,
                    dstMipLevel, dstLayer);

                _stagingMemoryPool.Free(trueCopySrc);
            }
            else // !isCompressed
            {
                if (_extensions.ARB_DirectStateAccess)
                {
                    glGetTextureSubImage(
                        srcGLTexture.Texture, (int)srcMipLevel, (int)srcX, (int)srcY, (int)srcZ,
                        width, height, depth,
                        srcGLTexture.GLPixelFormat, srcGLTexture.GLPixelType, block.SizeInBytes, block.Data);
                    CheckLastError();
                }
                else
                {
                    for (uint layer = 0; layer < depth; layer++)
                    {
                        uint curLayer = srcZ + srcLayer + layer;
                        uint curOffset = depthSliceSize * layer;
                        uint readFB;
                        glGenFramebuffers(1, &readFB);
                        CheckLastError();
                        glBindFramebuffer(FramebufferTarget.ReadFramebuffer, readFB);
                        CheckLastError();

                        if (srcGLTexture.ArrayLayers > 1 || srcGLTexture.Type == TextureType.Texture3D
                            || (srcGLTexture.Usage & TextureUsage.Cubemap) != 0)
                        {
                            glFramebufferTextureLayer(
                                FramebufferTarget.ReadFramebuffer,
                                GLFramebufferAttachment.ColorAttachment0,
                                srcGLTexture.Texture,
                                (int)srcMipLevel,
                                (int)curLayer);
                            CheckLastError();
                        }
                        else if (srcGLTexture.Type == TextureType.Texture1D)
                        {
                            glFramebufferTexture1D(
                                FramebufferTarget.ReadFramebuffer,
                                GLFramebufferAttachment.ColorAttachment0,
                                TextureTarget.Texture1D,
                                srcGLTexture.Texture,
                                (int)srcMipLevel);
                            CheckLastError();
                        }
                        else
                        {
                            glFramebufferTexture2D(
                                FramebufferTarget.ReadFramebuffer,
                                GLFramebufferAttachment.ColorAttachment0,
                                TextureTarget.Texture2D,
                                srcGLTexture.Texture,
                                (int)srcMipLevel);
                            CheckLastError();
                        }

                        CheckLastError();
                        glReadPixels(
                            (int)srcX, (int)srcY,
                            width, height,
                            srcGLTexture.GLPixelFormat,
                            srcGLTexture.GLPixelType,
                            (byte*)block.Data + curOffset);
                        CheckLastError();
                        glDeleteFramebuffers(1, &readFB);
                        CheckLastError();
                    }
                }

                UpdateTexture(
                    dstGLTexture,
                    (IntPtr)block.Data,
                    dstX, dstY, dstZ,
                    width, height, depth, dstMipLevel, dstLayer);
            }

            if (packAlignment < 4)
            {
                glPixelStorei(PixelStoreParameter.PackAlignment, 4);
                CheckLastError();
            }

            _stagingMemoryPool.Free(block);
        }

        private static void CopyWithFBO(
            OpenGLTextureSamplerManager textureSamplerManager,
            OpenGLTexture srcGLTexture, OpenGLTexture dstGLTexture,
            uint srcX, uint srcY, uint srcZ, uint srcMipLevel, uint srcBaseArrayLayer,
            uint dstX, uint dstY, uint dstZ, uint dstMipLevel, uint dstBaseArrayLayer,
            uint width, uint height, uint depth, uint layerCount, uint layer)
        {
            TextureTarget dstTarget = dstGLTexture.TextureTarget;
            if (dstTarget == TextureTarget.Texture2D)
            {
                glBindFramebuffer(
                    FramebufferTarget.ReadFramebuffer,
                    srcGLTexture.GetFramebuffer(srcMipLevel, srcBaseArrayLayer + layer));
                CheckLastError();

                textureSamplerManager.SetTextureTransient(TextureTarget.Texture2D, dstGLTexture.Texture);

                glCopyTexSubImage2D(
                    TextureTarget.Texture2D,
                    (int)dstMipLevel,
                    (int)dstX, (int)dstY,
                    (int)srcX, (int)srcY,
                    width, height);
                CheckLastError();
            }
            else if (dstTarget == TextureTarget.Texture2DArray)
            {
                glBindFramebuffer(
                    FramebufferTarget.ReadFramebuffer,
                    srcGLTexture.GetFramebuffer(srcMipLevel, srcBaseArrayLayer + layerCount));
                CheckLastError();

                textureSamplerManager.SetTextureTransient(TextureTarget.Texture2DArray, dstGLTexture.Texture);

                glCopyTexSubImage3D(
                    TextureTarget.Texture2DArray,
                    (int)dstMipLevel,
                    (int)dstX,
                    (int)dstY,
                    (int)(dstBaseArrayLayer + layer),
                    (int)srcX,
                    (int)srcY,
                    width,
                    height);
                CheckLastError();
            }
            else if (dstTarget == TextureTarget.Texture3D)
            {
                textureSamplerManager.SetTextureTransient(TextureTarget.Texture3D, dstGLTexture.Texture);

                for (uint i = srcZ; i < srcZ + depth; i++)
                {
                    glCopyTexSubImage3D(
                        TextureTarget.Texture3D,
                        (int)dstMipLevel,
                        (int)dstX,
                        (int)dstY,
                        (int)dstZ,
                        (int)srcX,
                        (int)srcY,
                        width,
                        height);
                    CheckLastError();
                }
            }
        }
    }
}
