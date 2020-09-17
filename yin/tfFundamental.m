%% Yartsev Lab
% Bats and Ages Script
% Vidush Mukund vidush_mukund@berkeley.edu
% March 13, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function builds acts as a wrapper to call the YIN Algorithm package;
% based on the script written of Yosef Prat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F0, ap] = tfFundamental(x, fs, wind, overlap,N, P)

% if the number of arguments is less than 6 then build the P object
if nargin < 6
    P =[];
end

%   P.minf0:    Hz - minimum expected F0 (default: 30 Hz)
%   P.maxf0:    Hz - maximum expected F0 (default: SR/(4*dsratio))
%   P.thresh:   threshold (default: 0.1)
%   P.relfag:   if ~0, thresh is relative to min of difference function (default: 1)
%   P.hop:      samples - interval between estimates (default: 32)
%   P.range:    samples - range of samples ([start stop]) to process
%   P.bufsize:  samples - size of computation buffer (default: 10000)
%	P.sr:		Hz - sampling rate (usually taken from file header)
%	P.wsize:	samples - integration window size (defaut: SR/minf0)
%	P.lpf:		Hz - intial low-pass filtering (default: SR/4)
%	P.shift		0: shift symmetric, 1: shift right, -1: shift left (default: 0)
if ~isfield(P,'thresh'); P.thresh = 0.1; end%?
if ~isfield(P,'sr'); P.sr = fs; end
if ~isfield(P,'wsize'); P.wsize = min(wind, floor(length(x)/2)-2);end
if ~isfield(P,'hop'); P.hop = wind-overlap; end
if ~isfield(P,'range'); P.range = [1 length(x)]; end% -overlap];
if ~isfield(P,'bufsize'); P.bufsize = length(x)+2; end
if ~isfield(P,'APthresh'); P.APthresh = 0.25; end

if ~isfield(P,'maxf0'); P.maxf0 = 20000; end
if ~isfield(P,'minf0'); P.minf0 = max(100,ceil(1.5*fs/length(x))); end

% call the yin function from the yin package
R = yin(x,P);

F0 = R.f0(round(P.wsize/2/P.hop)+1:end);
ap = R.ap0(round(P.wsize/2/P.hop)+1:end);
ap(ap ~= real(ap)) = abs(ap(ap ~= real(ap)));
ap = smooth(ap,5);

smoo = max(ceil(P.sr/P.hop /1000),5);
F0(ap > P.APthresh) = NaN;

% F0 = smooth(F0,smoo,'rlowess');

% any frequency values that are below the accepted range of fundamentals or
% above the maximum fundamental frequency are replced with NaN
F0(F0<=P.minf0 | F0>=P.maxf0 | F0<= 110) = NaN; 

% drop all NaN values (i.e. f0 values outside of the accepted range of
% possible frequency values)
ap = ap(~isnan(F0));
F0 = F0(~isnan(F0));

% NaN interpolation
% if length(F0) == 1
%     F0 = [NaN];
% else
%     F0 = inpaint_nans(F0);
% end

% F0 = ExInfNan(F0,'linear');
%F0 = smooth(F0,smoo);

if nargin < 5
    N = length(F0);
end

% perform a 1-D interpolation of the fundamenetal frequencies
t=(1:length(F0)).*P.hop/fs + 0.5/P.minf0;
tt=linspace(t(1),length(x)/fs,2*N+1);
tt = tt(2:2:end-1);
F0 = interp1(t,F0,tt);
ap = interp1(t,ap,tt);

% drop any NaN values
ap = ap(~isnan(F0));
F0 = F0(~isnan(F0));
ap = ap.^0.5;


% NaN interpolation used by Yoseef Prat
% F0 = ExInfNan(F0,'linear');
% ap = ExInfNan(ap,'linear').^0.5;

% NaN interpolation
% F0 = inpaint_nans(F0);
% ap = inpaint_nans(ap).^0.5;

end
