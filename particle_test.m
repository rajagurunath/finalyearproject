function [m, vtar] =  kalman_est(x, t, Fs, Ns, smoothing, n_block, n_overlap, sigma, lambda, m0,v0, sparse)
% KALMAN_EST is called by PKALMAN 
% Yuan Qi
% April 2003

x = x(:);
n = length(x);
twopi = 2*pi;
m = zeros(length(m0), n);
C = [sin(twopi*Fs*t(1)) cos(twopi*Fs*t(1)) 1];

v0C = v0 *  C';
k1 = v0C * inv( C* v0C + lambda);
vt = v0 - k1*v0C';
m(:,1) = m0 + k1*(x(1) - C*m0);      

vtar = [];
if smoothing ==1 
  % preparing for smoothing
  vtar =[];
  n_block = min(n,n_block); 
  vtar = zeros([size(v0), n_block]);    
  ptar = zeros([size(v0), n_block - 1]);    
  vtar(:,:,1) = vt;
end

mdi = 1;
i_old = 0;
for i = 2:n
  inttp = t(i)-t(i-1);
  % update the state noise variance
  sigma_t = sigma * (inttp);  
  % if mod(i,150)==0,  vt = v0;  end
  pt = vt + sigma_t;
  
    C = [sin(twopi*Fs*t(i)) cos(twopi*Fs*t(i)) 1]; 
  ptC = pt *  C';
  kt = ptC * (1./( C * ptC + lambda));
  vtp = pt-kt*(C*pt);        
  mtp = m(:,i-1) + kt*(x(i)- C*m(:,i-1));

  if sparse == 1 & mod(i, 100) == 0  & i > 50 % used for obs4 % if mod(i, 40) == 0  & i > 50 % used for obssp
      % OBS pruning
      [m(:,i) vt] = obs4(mtp,vtp,v0);
      % vt = v0;   
  else
    m(:,i) = mtp; 
    vt = vtp;
  end
      
  if smoothing == 1 
    mdi = mdi + 1;
    if mdi == 0
      nn = i;
      ptar(:,:,n_block-1) = pt;
      vtar(:,:,n_block) = vt;
    elseif mdi == 1
      vtar(:,:,1) = vt;
    elseif mdi >= 2
      vtar(:,:,mdi) = vt;
      ptar(:,:,mdi-1) = pt;
    end
    if  (mdi == n_block) | (i == n)
      ind = i+1 - (n_block:-1:1);
      i_old = i;
      for j = mdi-1:-1:1 
        jt =  vtar(:,:,j)*inv( ptar(:,:,j)); % ss+1e-4*eye(size(ptar(:,:,j))) );
        m(:,ind(j)) = m(:,ind(j)) + jt * (m(:,ind(j)+1) - m(:,ind(j)));
        vtar(:,:,j) = vtar(:,:,j) +  jt * (vtar(:,:,j+1) - ptar(:,:,j))*jt';
      end
          
      vtar(:,:,1) = vtar(:,:,n_block - n_overlap + 1);
      for g = 2:n_overlap
        vtar(:,:,g) = vtar(:,:,n_block - n_overlap + g);
        ptar(:,:,g-1) = ptar(:,:,n_block - n_overlap + g - 1);
      end
      mdi = n_overlap;
    end
  end
end
 