require 'nn'

N = 76  -- the number of patches to match
patches = torch.rand(N,2,64,64):float()

-- load the network
net = torch.load'../networks/2ch/2ch_liberty_nn.t7'

-- in place mean subtraction
local p = patches:view(N,2,64*64)
p:add(-p:mean(3):expandAs(p))

-- get the output similarities
output = net:forward(patches)
