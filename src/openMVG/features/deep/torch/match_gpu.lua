require 'cunn'
-- example of matching descriptors on GPU

N = 76  -- the number of patches to match
patches = torch.rand(N,2,64,64):cuda()

-- load the network and move to cuda
net = torch.load'../networks/2ch/2ch_liberty_nn.t7':cuda()

-- in place mean subtraction
local p = patches:view(N,2,64*64)
p:add(-p:mean(3):expandAs(p))

-- get the output similarities
output = net:forward(patches)


-- OR if there is a lot of patches to match
-- it is better to split them to smaller batches
batch_size = 128

-- preallocate the output tensor
similarities = torch.CudaTensor(N)
similarities_split = similarities:split(batch_size)

for i,v in ipairs(patches:split(batch_size)) do
  similarities_split[i]:copy(net:forward(v))
end
