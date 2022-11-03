function [Vhat, Vhatrobust] = Vhat_MLE(likelihood_f,params, varargin)
    
try
    [LLK, llk] = feval(likelihood_f,params,varargin{:});
catch
    error('There was an error evaluating the function.  Please check the arguements.');
end


hess = hessian_2sided(likelihood_f, params, varargin{:});

if rcond(hess) < 1e-15
    Vhat = NaN;
    Vhatrobust = NaN;
else


Vhat = hess^(-1);   % no need to multiply by (-1) because likelihood_f is usually negative log likelihood function

T = length(llk);
a = length(params);

if nargout == 2
   %h = eps.^(1/3)*max(abs(params),1e-2);
   h=min(abs(params/2),max(params,1e-2))*eps^(1/4);
   %h = 1e-5*ones(size(params));
   hplus=params+h;
   hminus=params-h;
   likelihoodsplus=zeros(T,a);
   likelihoodsminus=zeros(T,a);
   for i=1:a
      hparams=params;
      hparams(i)=hplus(i);      % only vary the i'th parameter
      [DUMMY, indivlike] = feval(likelihood_f,hparams,varargin{:});
      likelihoodsplus(:,i)=indivlike;
   end
   for i=1:length(params)
      hparams=params;
      hparams(i)=hminus(i);
      [DUMMY, indivlike] = feval(likelihood_f,hparams,varargin{:});         
      likelihoodsminus(:,i)=indivlike;
   end
   scores=(likelihoodsplus-likelihoodsminus)./(2*repmat(h',T,1));   % n rows with k columns, each row is the score for one observation
   scores=scores-repmat(mean(scores),T,1);      % make sure scores have exactly mean of zero
   B=scores'*scores;        % amounts to summing up the outer product of the scores
   Vhatrobust=Vhat*B*Vhat;      % Hamilton eq. 5.8.7
end

end