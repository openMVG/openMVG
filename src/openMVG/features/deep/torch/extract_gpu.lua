require 'cunn'
-- example of extracting descriptors on GPU

N = 76  -- the number of patches to match
patches = torch.rand(N,1,64,64):cuda()

-- load the network
net = torch.load'../networks/siam2stream/siam2stream_liberty_nn.t7':cuda()
print(net)

-- to extract the descriptors we need a branch of the first parallel module
net = net:get(1):get(1)
print(net)

-- in place mean subtraction
local p = patches:view(N,1,64*64)
p:add(-p:mean(3):expandAs(p))

-- get the output descriptors
output = net:forward(patches)

-- print the size of the output tensor
print(#output)

-- OR if there is a lot of patches to extract descriptors
-- it is better to split them to smaller batches
batch_size = 128

-- preallocate the output tensor, here we know the dimensionality of descriptor
descriptors = torch.CudaTensor(N,512)
descriptors_split = descriptors:split(batch_size)

for i,v in ipairs(patches:split(batch_size)) do
  descriptors_split[i]:copy(net:forward(v))
end
