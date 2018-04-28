% Copyright (C) 2001 Michel Juillard
%
% computes second order partial derivatives
% uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27 p. 884
% note that we have to reshape the result of this function. Also, if we use
% a minimizer, we have to multiply the result by -1.
%e.g.
%Hess = hessian('postfourier',xestimate,Tfinal,DATA, cutoff, periodigram);
%Hess = reshape(Hess,13,13);
%H= -eye(size(Hess))/Hess; 

function hessian_mat = hessian3(func,x,setup,data,block, params)
gstep_=1e-2;
func = str2func(func);
n=size(x,1);
%h1=max(abs(x),gstep_*ones(n,1))*eps^(1/3);
h1=max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/6);
h_1=h1;
xh1=x+h1;
h1=xh1-x;
xh1=x-h_1;
h_1=x-xh1;
xh1=x;
f0=feval(func,x,setup,data,block, params);
f1=zeros(size(f0,1),n);
f_1=f1;
for i=1:n
    xh1(i)=x(i)+h1(i);
    f1(:,i)=feval(func,xh1,setup,data,block, params);
    xh1(i)=x(i)-h_1(i);
    f_1(:,i)=feval(func,xh1,setup,data,block, params);
    xh1(i)=x(i);
   i=i+1;
end
xh_1=xh1;
hessian_mat = zeros(size(f0,1),n*n);
for i=1:n
   
    if i > 1
        k=[i:n:n*(i-1)];
        hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
    end 
    hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
    temp=f1+f_1-f0*ones(1,n);
        for j=i+1:n
        xh1(i)=x(i)+h1(i);
        xh1(j)=x(j)+h_1(j);
        xh_1(i)=x(i)-h1(i);
        xh_1(j)=x(j)-h_1(j);
        hessian_mat(:,(i-1)*n+j)=-(-feval(func,xh1,setup,data,block, params)-feval(func,xh_1,setup,data,block, params)+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
        xh1(i)=x(i);
        xh1(j)=x(j);
        xh_1(i)=x(i);
        xh_1(j)=x(j);
        j=j+1;
    end
    i=i+1;
end


% 11/25/03 SA Created from Hessian_sparse (removed sparse)
