function divfun=divfun(egVITout,lengthT)

global n

egVITouttimesno=egVITout.*repmat(0:n,lengthT,1);
divfun=sum(egVITouttimesno,2)./sum(egVITout,2)/n;