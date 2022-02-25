function [x,diagnostic] = cgls(A,x0,y,max_iter,tol)
%
% mimimize ||Ax-y||^2 with max_iter iterations, x0 is the intial sol
% and x is the final solution obtained with CG. tol is the tolerance
% M.D.Sacchi
% A: operator
% y: RHS
% x0: initial solution
% max_iter: maximum number of iteratiions
% tol : tol
% outputs are x: solution diagnostic: a varible that describes
% stopping criteria
% MD Sacchi for g431_531

 x  = x0;
 s=y-A*x;
 p=A'*s;
 r=p;
 q=A*p;
 old = r'*r;

 diagnostic.iter = max_iter;
 diagnostic.error = 0;

 error_old = 10.e10;
% CG loop

 for k=1:max_iter
     %sprintf('cgls iteration %d of %d',k,max_iter)
  alpha = (r'*r)/(q'*q);
  x= x +alpha*p;
  s = s-alpha*q;
  r= A'*s;
  new = r'*r;
  beta = new/old; 
  old = new;
  p = r + beta*p;
  q = A*p;
  error = s'*s; 
  if (abs(error-error_old)<tol);break;end
  error_old = error;
  diagnostic.iter = k;
  diagnostic.error = error;
 end

