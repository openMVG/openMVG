require 'cunn'
require 'inn'

local ffi = require 'ffi'

ffi.cdef[[
void cunnrelease_SpatialConvolution(THCState *state,
    THCudaTensor *input,
    THCudaTensor *weight,
    THCudaTensor *bias,
    THCudaTensor *columns,
    THCudaTensor *ones,
    THCudaTensor *output,
    int nInputPlane, int nOutputPlane, int kW, int kH, int dW, int dH, int padding);

void cunnrelease_SpatialMaxPooling(THCState* state,
    THCudaTensor* input, 
    THCudaTensor* output,
    int kW, int kH, int dW, int dH, bool is_ceil);

void cunnrelease_SpatialAveragePooling(THCState* state,
    THCudaTensor* input,
    THCudaTensor* output,
    int kW, int kH, int dW, int dH, bool is_ceil);

void cunnrelease_Linear(THCState *state,
    THCudaTensor *input,
    THCudaTensor *output,
    THCudaTensor *weight,
    THCudaTensor *bias,
    THCudaTensor *buffer);

void cunnrelease_ReLU(THCState *state,
    THCudaTensor *input);
]]

local C = ffi.load'./build/libcunnrelease.dylib'

local mytester = torch.Tester()

local cunnreleasetest = {}

precision_forward = 0

function cunnreleasetest.SpatialMaxPooling_floor()
  local from = math.random(1,5)
  local ki = math.random(1,4)
  local kj = math.random(1,4)
  local si = math.random(1,3)
  local sj = math.random(1,3)
  local outi = math.random(4,5)
  local outj = math.random(4,5)
  local ini = (outi-1)*si+ki
  local inj = (outj-1)*sj+kj

  local module = nn.SpatialMaxPooling(ki,kj,si,sj):cuda()
  local input = torch.rand(from,ini,inj):cuda()

  local ref_output = module:forward(input)

  local act_output = torch.CudaTensor()
  C.cunnrelease_SpatialMaxPooling(cutorch.getState(), input:cdata(), act_output:cdata(), ki, kj, si, sj, false)

  mytester:asserteq((ref_output - act_output):abs():max(), precision_forward, 'SpatialMaxPooling_floor')
end

function cunnreleasetest.SpatialMaxPooling_ceil()
  local from = math.random(1,5)
  local ki = math.random(1,4)
  local kj = math.random(1,4)
  local si = math.random(1,3)
  local sj = math.random(1,3)
  local outi = math.random(4,5)
  local outj = math.random(4,5)
  local ini = (outi-1)*si+ki
  local inj = (outj-1)*sj+kj

  local module = inn.SpatialMaxPooling(ki,kj,si,sj):cuda()
  local input = torch.rand(from,ini,inj):cuda()

  local ref_output = module:forward(input)

  local act_output = torch.CudaTensor()

  C.cunnrelease_SpatialMaxPooling(cutorch.getState(), input:cdata(), act_output:cdata(), ki, kj, si, sj, true)

  mytester:asserteq((ref_output - act_output):abs():max(), precision_forward, 'SpatialMaxPooling_ceil')
end

function cunnreleasetest.SpatialConvolution_forward_single()
   local from = math.random(1,32)
   local to = math.random(1,8) * 8
   local ki = math.random(3,15)
   local kj = math.random(3,15)
   local si = 1 -- not supported by CPU version yet
   local sj = si
   local outi = math.random(1,64)
   local outj = math.random(1,64)
   local ini = (outi-1)*si+ki
   local inj = (outj-1)*sj+kj

   local input = torch.randn(from,inj,ini):cuda()
   local sconv = nn.SpatialConvolutionMM(from,to,ki,kj,si,sj):cuda()
   local groundtruth = sconv:forward(input)

   local finput = torch.CudaTensor()
   local fgradinput = torch.CudaTensor()
   local output = torch.CudaTensor()

   C.cunnrelease_SpatialConvolution(cutorch.getState(),
   	input:cdata(),
	sconv.weight:cdata(),
	sconv.bias:cdata(),
	finput:cdata(),
	fgradinput:cdata(),
	output:cdata(),
	from, to, ki, kj, si, sj, 0)

   local error = output - groundtruth
   mytester:asserteq(error:abs():max(), precision_forward, 'error on state (forward) ')
end

function cunnreleasetest.Linear()
  local from = math.random(1,32)
  local to = math.random(1,32)
  local bs = math.random(2,32)

  local module = nn.Linear(from, to):cuda()
  local input = torch.rand(bs, from):cuda()
  local groundtruth = module:forward(input)

  local output = torch.CudaTensor()
  local buffer = torch.CudaTensor(bs):fill(1)
  C.cunnrelease_Linear(cutorch.getState(),
  	input:cdata(),
	output:cdata(),
	module.weight:cdata(),
	module.bias:cdata(),
	buffer:cdata())

  local error = output - groundtruth
  mytester:asserteq(error:abs():max(), precision_forward, 'error on state (forward) ')
end

mytester:add(cunnreleasetest)
mytester:run()
