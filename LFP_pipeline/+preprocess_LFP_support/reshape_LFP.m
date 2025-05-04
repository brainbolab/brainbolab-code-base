function lfp_struct = reshape_LFP( ...
  lfp_struct, ...
  window_sz, ...
  fs)
%%
% This function takes an `LFP.mat` data file's fields saved in a variable
% eg. from calling `lfp_struct = open('LFP.mat')` and organizes the data into
% an MxNxW tensor where M is the number of data points in a window, N is 
% the number of channels and W is the number of windows
%
% Parameters:
% 
%   lfp_struct:  The resulting structure from calling `open()` on one 
%                of the `LFP.mat` data files created from calling  
%                `NSX2MAT()` on an `.nsx` file. The fields of this 
%                struct are brain region channel labels eg. `BLA_01`
%                or `NAc_02` and their values are vectors of signal
%                data
%
%   window_size: In seconds, the length of time for a single window
% 
%   fs:          Sampling frequency of the data in Hz
%
% Return values:
%
%   Three dimensional MxNxW tensor from the input LFP structure.
%%
  
  lfp_struct = double(cell2mat(struct2cell(lfp_struct) ));
  
  M_dim = window_sz * fs;
  N_dim = size(lfp_struct, 1);
  W_dim = floor(size(lfp_struct, 2) / M_dim);

  lfp_struct = lfp_struct(:,1:(M_dim*W_dim) );
  
  lfp_struct = reshape(lfp_struct, [N_dim, M_dim, W_dim]);
  lfp_struct = permute(lfp_struct, [2, 1, 3]);
end