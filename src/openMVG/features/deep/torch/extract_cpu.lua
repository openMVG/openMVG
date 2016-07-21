require 'nn'
-- example of extracting descriptors on CPU

N = 76  -- the number of patches to match
patches = torch.rand(N,1,64,64):float()

-- load the network
net = torch.load'../networks/siam2stream/siam2stream_liberty_nn.t7'
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
