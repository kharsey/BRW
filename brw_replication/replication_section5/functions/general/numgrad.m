function [g, badg] = numgrad(fcn,x,varargin)
% function [g badg] = numgrad(fcn,x,varargin)
%
delta = 1e-5;
%delta=1e-2;
n=length(x);
tvec=delta*eye(n);
m = length(feval(fcn,x, varargin{:}));
g=zeros(m,n);

f0 = feval(fcn,x,varargin{:});
badg=0;
for i=1:n
   scale=1; 
   if size(x,1)>size(x,2)
      tvecv=tvec(i,:);
   else
      tvecv=tvec(:,i);
   end   
   g0 = (feval(fcn,x+scale*tvecv', varargin{:}) - f0) ...
         /(scale*delta);
   if abs(g0)< 1e15
       g(:,i)=g0;
      % disp('good gradient') % Jinill Kim
   else
      disp('bad gradient ------------------------') % Jinill Kim
      % fprintf('Gradient w.r.t. %3d: %10g\n',i,g0) %see above
      g(i)=0;
      badg=1;
      % return
      % can return here to save time if the gradient will never be
      % used when badg returns as true.
   end
end
