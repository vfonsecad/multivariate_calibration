function [tmcd,smcd,weights,initmean,initcov,rawrd,Hopt] = unimcd(y,h,varargin)

%UNIMCD computes the MCD estimator of a univariate data set.  This 
% estimator is given by the subset of h observations with smallest variance. 
% The MCD location estimate is then the mean of those h points,
% and the MCD scale estimate is their standard deviation (multiplied with a 
% consistency factor). Based on these raw estimates,
% a reweigthing step is applied as in the FASTMCD algorithm (see mcdcov.m).
% We recommend to use the function mcdcov.m, which calls unimcd.m.
%
% The MCD method was introduced in:
%
%   Rousseeuw, P.J. (1984), "Least Median of Squares Regression," 
%   Journal of the American Statistical Association, Vol. 79, pp. 871-881.
%
% The algorithm to compute the univariate MCD is described in
%
%   Rousseeuw, P.J., Leroy, A., (1988), "Robust Regression and Outlier Detection,"
%   John Wiley, New York.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust
%
% Written by: Katrien Van Driessen
% Revision by S. Verboven, M. Hubert, Wannes Van den Bossche
% Last update: 12/10/2014 

if nargin > 2
    centered=varargin{1};
else
    centered=0;
end

ncas=length(y);
len=ncas-h+1;
xorig=y;

if len==1
    if centered==1
        tmcd=0;
    else
        tmcd=mean(y);
    end
    smcd=sqrt(var(y));
    weights=ones(length(y),1);
else
    if centered==1
        [y,I]=sort(abs(y));
        Hopt=I(1:h);
        initmean=0;
        initcov=sum(y(1:h).^2)/(h-1);
    else
        [y,I]=sort(y);
        ay=zeros(1,len);
        ay(1)=sum(y(1:h));
        for samp=2:len
            ay(samp)=ay(samp-1)-y(samp-1)+y(samp+h-1);
        end
        ay2=ay.^2/h;
        sq=zeros(1,len);
        sq(1)=sum(y(1:h).^2)-ay2(1);
        for samp=2:len
            sq(samp)=sq(samp-1)-y(samp-1)^2+y(samp+h-1)^2-ay2(samp)+ay2(samp-1);
        end
        sqmin=min(sq);
        ii=find(sq==sqmin);
        Hopt = I(ii:ii+h-1);
        ndup=length(ii);
        slutn(1:ndup)=ay(ii);
        initmean=slutn(floor((ndup+1)/2))/h; %initial mean
        initcov=sqmin/(h-1); %initial variance
    end
    % calculating consistency factor
    res=(xorig-initmean).^2/initcov;
    sortres=sort(res);
    factor=sortres(h)/chi2inv(h/ncas,1);
    initcov=factor*initcov;
    res=(xorig-initmean).^2/initcov; %raw_robdist^2
    quantile=chi2inv(0.975,1);
    weights=res<=quantile; %raw-weights
    rawrd=sqrt(res);
    %reweighting procedure
    if size(weights,1)~=size(y,1)
        weights=weights';
    end
    if centered==1
        tmcd=0;
    else
        tmcd=sum(xorig.*weights)/sum(weights);
    end
    smcd=sqrt(sum((xorig-tmcd).^2.*weights)/(sum(weights)-1));
end


