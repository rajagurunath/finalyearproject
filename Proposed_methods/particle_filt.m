function varargout = pkalman(x, varargin)
% PKALMAN Power Spectral Density (PSD) estimate via Kalman filtering and smoothing.
% INPUT
%   x: the input signal, which is a column vector
%
% The following optional input arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%   Ns: the frequency resolution. If it is in the input argument list, Ns must 
%       show before Fs. [Ns = 64]
%   Fs: Fs can be either a scalar or a row vector.
%      (scalar) the sampling frequency Fs in Hz, or twice of the freqeunce range of a 
%       given signal. [Fs=6 (Hz)]
%      (vector) the given frequencies bases for spectrum estimation
%   t: the sampling time.   [(1:n)/(Fs*2)]
%   sigma: the process noise variance    [1e1*eye(2*length(Fs)+1)]
%   lambda: the observation noise variance   [1e-2]
%   m0: the mean of the prior distribution    [ones(2*length(Fs)+1,1)]
%   v0: the variance of the prior distribution  [1e2*eye(2*length(Fs)+1,1)]
%   disp:  0-no display; 1-displaying the result
%   fig_sz: display size in inch [matlab default size]
%   color: 0-gray display 1-color display [1];
%   log_scale: 0-not plotting the image in log scale 1-plotting the image in log scale. 
%   ts:  the index of one time slice. If it is a positive integer, then plot the 
%        spectrum at that time slice. if it is an empty matrix, then plot the 
%        whole spectrogram. [ [] ]
%   tp: the start point in time for displaying a spectrogram 
%   te: the end point in time for displaying a spectrogram 
%   axis: matlab display argument [[]]   
%   XTICK: matlab display argument [[]]   
%   YTICK: matlab display argument [[]]   
%   FontSize: matlab display argument [[]]   
%
% OUTPUT
% varargout = [amp, phs, vtar]  
%   amp: estimated frequency amplitudes from time 1 to n
%   phs: estimated frequency phases from time 1 to n
%   vtar: the variances from last block's smoothing; it's an empty 
%        matrix when we only do Kalman filtering
%
% EXAMPLE 1: Estimate the spectrum of a signal with fast decaying amplitude and 10 percent missing data
%   fs = 500;
%   tt = 0:(1/fs):2;
%   perc = 0.9;
%   tr = randperm(length(tt));
%   tc = sort(tr(1:floor(length(tt)*perc))); % truncated index of tt
%   t = tt(tc);
%   w = 50;
%   aplmd = exp(-t*1);
%   x = (sin(2*pi*t*w).*aplmd)';
%   m = pkalman(x,'t',t,'Ns',50,'Fs', 200, 'smoothing',1, 'LogScale',0,'disp',1);
% EXAMPLE 2: Estimate the spectrum of a signal consisting of three frequency components which are close to each other
%   fs = 50; 
%   t = 0:(1/fs):3;
%   y = sin(2*pi*21*t+ 20) + 1*sin(2*pi*20*t + 20)+.5*sin(2*pi*19*t + 20);
%   nslevel = sqrt(1e-1); 
%   x=y' + nslevel*randn(size(y))';
%   pkalman(x,'Ns',64,'Fs',fs,'ts',length(x),'axis',...
%   [0 25 -55 10],'XTICK',[0 5 10 15 19 20 21 25], 'YTICK',[-40 -20 0],'FontSize',6)
%
% References
%     [1] Bayesian Spectrum Estimation of Unevenly Sampled Nonstationary Data", 
%      Yuan Qi, Thomas P. Minka, and Rosalind W. Picard, MIT Media Lab 
%      Technical Report Vismod-TR-556, 2002
%      http://web.media.mit.edu/~yuanqi/papers/qi-spec-02.pdf
%     [2] A shorter version in ICASSP 02, May 2002, Orlando, Florida
%      http://web.media.mit.edu/~yuanqi/papers/qi-icassp02.pdf
%
% Author: Yuan Qi
%        MIT Media lab
% Date: April 12, 2003

x = x(:);
n = length(x);

% parsing the parameters
[t, Fs, Ns, smoothing, n_block, n_overlap, sigma, lambda, m0,v0, disp_flg, sparse] = parse_par(n,varargin{:});

% estimate the posterior mean by Kalman filtering/smoothing
[m, vtar] =  kalman_est(x, t, Fs, Ns, smoothing, n_block, n_overlap, sigma, lambda, m0,v0, sparse);

ds = length(Fs);
amp = [];
for j = 1:ds
  rowwd = [j j+ds];
  % amp: amplitude estimate
  amp(j,:) = sqrt(sum(m(rowwd, :).^2));
end
amp = ([abs(m(end,:)); amp]);
if nargout > 0
  varargout{1} = amp;
end
if nargout > 1
   phs = [];
   for j = 1:ds
     rowwd = [j j+ds];
     % phs: phase estimate
     phs(j,:) = atan(m(j, :)./m(j+ds,:));
   end
   phs = ([zeros(1,n); phs]);
   varargout{2} = phs;
end
if nargout > 2
   varargout{3} = vtar;
end

if nargout == 0 | disp_flg == 1
  disp_spec(amp, t, Fs, varargin{:});
end

function [t, Fs, Ns, smoothing, n_block, n_overlap, ...
    sigma, lambda, m0,v0, disp_flg, sparse] ...
  = parse_par(n, varargin)
% PARSE_PAR parses the parameters

% Set default values for the parameters: t, Fs, Ns
t = [];
Ns = 64; 
Fs = 6; 
spc = Fs/2 /Ns; 
Fs = (0+spc):spc:(Fs/2-spc);
% Parse the list of the input arguments: t, Fs, Ns
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
    case 't', t = args{i+1};
    case 'Ns', Ns = args{i+1};  
    case 'Fs',
      if  length(args{i+1})==1 
        Ftp = args{i+1}/2;
        spc = Ftp/Ns;
        Fs = (0+spc):spc:(Ftp-spc);
      else
        Ftp = [];
        Fs = args{i+1}; % In this case, Ns is not used.
        Fs = Fs(:)';
      end
  end
end
% set value for t if t is not specified by input arguments
if isempty(t)
  t = (1:n)/(Ftp*2); 
end

% Set default values for the parameters:sparse, smoothing, sigma, lambda, m0,v0, and disp_flg
sparse = 0;
smoothing = 0;
n_block = 1000;
n_overlap = 100;
n_bs = 2*length(Fs)+1;
sigma = 1*1e1*eye(n_bs);
lambda = 1e-2;
m0 = ones(n_bs,1); 
v0 = 1e2*eye(n_bs);
disp_flg = 0;
% Parse the list of the input arguments:smoothing, sigma, lambda, m0,v0
for i=1:2:nargs
  switch args{i},
  case 'sparse',  sparse = args{i+1};
    case 'smoothing',  smoothing = args{i+1};
    case 'n_block',  n_block = args{i+1};
    case 'n_overlap',  n_overlap = args{i+1};
    case 'var_pro_noise', sigma = args{i+1};  
    case 'var_obs_noise', lambda = args{i+1};  
    case 'mean_prior', m0 = args{i+1};
    case 'var_prior', V0 = args{i+1}; 
    case 'disp', disp_flg = args{i+1};
    case 't', ;
    case 'Ns',;  
    case 'Fs',;    
    case 'fig_sz',;
    case 'tp', ;  
    case 'te', ;
    case 'ts', ;
    case 'color',;
    case 'axis', ;
    case 'XTICK',;
    case 'YTICK', ;
    case 'FontSize',;
    case 'LogScale', ;
    otherwise,
       error(['invalid argument name: ' args{i}]);
  end
end

function disp_spec(amp, t, Fs, varargin)

% default value
tp =1; 
te = length(t);
ts = []; % the index of one time slice 
fig_sz = []; % fig_sz = [2.1 2.1]; 
color = 1;
log_scale = 1;
ax = []; xtick = []; ytick = []; font_sz = [];

% Parse the list of the input arguments
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
    case 'fig_sz', fig_sz = args{i+1};
    case 'tp', tp = args{i+1};  
    case 'te', te = args{i+1};  
    case 'ts', ts = args{i+1};  
    case 'color', color = args{i+1};  
    case 'axis', ax = args{i+1};  
    case 'XTICK', xtick = args{i+1};  
    case 'YTICK', ytick = args{i+1};  
    case 'FontSize', font_sz = args{i+1};
    case 'LogScale', log_scale = args{i+1};
  end
end

xa=[0 Fs];
figure, % 
if ~isempty(fig_sz)
  figsize(fig_sz);
end
if isempty(ts)
  if log_scale == 1
    imagesc(t(tp:te),xa,(20*log10(amp(:,tp:te)))),axis xy;
  else
    imagesc(t(tp:te),xa, amp(:,tp:te) ),axis xy;  
  end
  hy = ylabel('Frequency (Hz)');
  if ~isempty(font_sz), set(hy,'FontSize', font_sz), end
  hx = xlabel('Time (Second)');
  if ~isempty(font_sz), set(hx,'FontSize', font_sz), end
else
  if log_scale == 1
    plot(xa,  20*log10(abs(amp(:,ts))));
    hy = ylabel('Amplitude (dB)');  
  else
    plot(xa, abs(amp(:,ts)));
    hy = ylabel('Amplitude');
  end
  if ~isempty(font_sz), set(hy,'FontSize', font_sz), end
  hx = xlabel('Frequency (Hz)');
  if ~isempty(font_sz), set(hx,'FontSize', font_sz), end
end
colormap(jet);
grid on
if color == 0
  load grayflp;
  colormap(grayflp);
end

if ~isempty(ax), axis(ax), end
if ~isempty(xtick), set(gca,'XTick', xtick), end
if ~isempty(ytick), set(gca,'YTick', ytick), end
if ~isempty(font_sz), set(gca,'FontSize', font_sz), end
