function code = train(traindir, n)
% Speaker Recognition: Training Stage
%
% Input:
%       traindir : string name of directory contains all train sound files
%       n        : number of train files in traindir
%
% Output:
%       code     : trained VQ codebooks, code{i} for i-th speaker
%
% Note:
%       Sound files in traindir is supposed to be: 
%                       s1.wav, s2.wav, ..., sn.wav
% Example:
%       >> code = train('C:\data\train\', 8);

k = 8;                         % number of centroids required

%code = zeros(8, 20, 16);
for i = 1:n                     % train a VQ codebook for each speaker
    file = sprintf('%ss%d.mp3', traindir, i);           
    disp(file);
   
    [s, fs] = audioread(file);
    s = isolate(s, fs);
    v = mfcc1(s, fs);            % Compute MFCC's
   
    code{i} = vqlbg(v, k);      % Train VQ codebook
end
