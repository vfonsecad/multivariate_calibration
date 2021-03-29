% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% ------- Kennard-Stone algortihm as depicted in 
%         https://doi.org/10.1016/j.geoderma.2014.02.002 ------

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function Output = kennard_stone(X,Nout, varargin)

if Nout < size(X,1) 
    
    optN = length(varargin);
    args = cell(1,1);
    args{1,1} ='mahalanobis';


    aa = 0;
    while aa < optN
        aa = aa + 1;
        args{aa,1} = varargin{aa};
    end


    if args{1,1} == 'mahalanobis'
       X =  (X-mean(X))./std(X);
    end


    Nin = size(X,1);
    K = size(X,2);
    sample_selected = zeros(Nin,1);


    xcal_in = X(:,:);
    xcal_out = zeros(Nout,K);


    % Initialize
    in = 1;max_DD = 1000000;
    DD = pdist2(xcal_in,mean(xcal_in));
    [~,id] = min(DD);
    xcal_out(in,:) = xcal_in(id,:);
    sample_selected(id,1) = 1;

    %
    while and(in < Nout, max_DD > 0.00001)

        in = in + 1 ;
        DD = pdist2(xcal_in,xcal_out(sum(abs(xcal_out),2)>0,:));
        DD_row = min(DD,[],2);
        [max_DD,id] = max(DD_row);
        xcal_out(in,:) = xcal_in(id,:);
        sample_selected(id,1) = 1;   
        %disp(in)


    end


    Output.sample_id = sample_selected;
    Output.xout = xcal_out;

else 
    
     Output.sample_id = ones(size(X,1),1);
     Output.xout = X(:,:);
     
end






