function bdb = genbnd11(w,wbd,W,A,B,C,D,Pnom,ph_r,info)
% Compute QFT bounds for the following closed-loop configuration
%
%                          ||A| + |B||G||
%                          |------------| <= WS
%                          |   C + DG   |
%
%     GENBNDS(11,W,WBD,WS,A,B,C,D,Pnom,PHS) computes bounds for the
%     above configuration.  W is the entire set of
%     frequencies and WBD is a subset of W designating bounds to compute.
%     WS is the performance specification. A, B, C and D are
%     the frequency response data sets, Pnom designates the nominal plant
%     for the configuration. PHS specifies at which phases (degrees)
%     to compute bounds.
%

% Author: Craig Borghesani
% 5/24/94
% Copyright (c) 2003, Terasoft, Inc.


m=length(ph_r); nbd=length(w);
pbds=qsubset(wbd,w);
final_state=ones(1,length(wbd));
state2=[];
myeps=1e-16;

[rma,cma]=size(A);[rmc,cmc]=size(C);
[rmb,cmb]=size(B);[rmd,cmd]=size(D);
[rmw,cmw]=size(W);
row = [rma,rmb,rmc,rmd,rmw];
maxr = max(row);

% preallocating matrix for more efficient programming
bmag=[ones(m,length(pbds))*myeps;ones(m,length(pbds))*1/myeps];

uPnom = abs(Pnom); vPnom = qatan4(Pnom);

%***
uA = abs(A); vA = 0;
uB = abs(B); vB = 0;
%***

uC = abs(C); vC = qatan4(C);
uD = abs(D); vD = qatan4(D);

% declaring sizes of replicating matricies
if repltest,
 u=ones(maxr,1); v=ones(1,m);
else
 u=ones(1,maxr); v=ones(1,m);
end

j=1; act=1; bct=1; cct=1; dct=1;
while j<=length(pbds),

 phi = ph_r(u,:) - (vPnom(pbds(j)));

 if cma>1, act=pbds(j); end
 if cmb>1, bct=pbds(j); end
 if cmc>1, cct=pbds(j); end
 if cmd>1, dct=pbds(j); end

 mA=uA(:,act); mB=uB(:,bct); pA=vA(:,act); pB=vB(:,bct);
 mC=uC(:,cct); mD=uD(:,dct); pC=vC(:,cct); pD=vD(:,dct);

%%%%%% V5 code
 if length(mA) == 1, mA = mA(u); end
 if length(mB) == 1, mB = mB(u); end
 if length(mC) == 1, mC = mC(u); end
 if length(mD) == 1, mD = mD(u); end
 if length(pA) == 1, pA = pA(u); end
 if length(pB) == 1, pB = pB(u); end
 if length(pC) == 1, pC = pC(u); end
 if length(pD) == 1, pD = pD(u); end

 psi1=pC-pD; psi1=psi1(:,v)-phi;
 psi2=pA-pB; psi2=psi2(:,v)-phi;
 Ws=W(:,pbds(j));

 A=Ws.^2 .* mD.^2 - mB.^2;

%***
 B=2*(Ws(:,v).^2 .*mC(:,v).*mD(:,v).*cos(psi1) - mA(:,v).*mB(:,v));
%***

 C=Ws(:,v).^2 .*mC(:,v).^2 - mA(:,v).^2;

 A=A(:,v);

 [g1,g2]=quadrtic(A,B,C); size_g1 = size(g1);
 cbdb=[[[g1';g2'],bmag(:,j)];ones(2,size_g1(1)+1)];
 [abvblw,state(j)]=sectbnds(cbdb,0);
 bmag(:,j)=abvblw(1:2*m);

 z=find(bmag(:,j)~=myeps & bmag(:,j)~=1/myeps & ...
        bmag(:,j)~=-248 & bmag(:,j)~=248 & ...
        bmag(:,j)~=-302 & bmag(:,j)~=302);
 if length(z), bmag(z,j)=bmag(z,j)*uPnom(pbds(j)); end

 frac=j/length(pbds);
 set(info(2),'xdata',[0,frac,frac,0]);
 set(info(3),'string',[num2str(100*frac),'%']);
 drawnow;

 j=j+1;
end
close(info(1));

bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302) = ...
   20*log10(bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302));

bnd=[bmag; w(pbds); 11*ones(1,length(pbds))];

mesgbnds(w,pbds,state,state2,ph_r,bnd);

[jk,t]=sort(w(pbds));
if nargout==0,
 plotbnds(bnd(:,t),[],180/pi*ph_r);
 title('GENBND11 Bounds'),xlabel('X: Phase (degrees)  Y: Magnitude (dB)')
else
 bdb=bnd(:,t);
end
