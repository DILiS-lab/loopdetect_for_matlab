function[J]=numerical_jacobian(func,x,epsilon,method)
% NUMERICAL_JACOBIAN Determine Jacobian matrix via finite difference at
% certain variable values.
%
% J = NUMERICAL_JACOBIAN(FUNC,X) computes the Jacobian matrix of the ODE 
% y'=FUNC(y) (i.e. the partial derivatives of FUNC) at variable value 
% column vector X. FUNC is a function handle and depends only on the 
% variable values. Implemented is the finite difference approach in that
% df/dx_k is approximated by 1/e_k * (FUNC(X+h*e_k/2)-FUNC(X-h*e_k/2)) 
% (centered differences) with sufficiently small stepsize h. The size of 
% the matrix will be length(FUNC(X)) times length(X), and is usually 
% quadratic for ODE systems.
% 
% J = NUMERICAL_JACOBIAN(FUNC,X,EPSILON,METHOD) allows to tune the stepsize
% EPSILON (default: 1e-6) and to deviate from the default of centered 
% differences by setting METHOD to 'forward' for forward differences 
% (FUNC(X+EPSILON*e_k)-FUNC(X)) or to 'backward' for backward differences 
% in the derivative computation (FUNC(X)-FUNC(X-EPSILON*e_k)). Note that 
% centered differences are considered to be more exact in general than the 
% other two. However, if variable values are close to or at zero, forward 
% differences could be more appropriate to capture the function behavior. 
% 
% Example
%   J = numerical_jacobian(@(x)func_POSm4(1,x,klin,knonlin),[1 2 1 1])
%
% See also: numerical_jacobian_complex()

%set default value for the perturbation
if nargin<3
    epsilon= 1e-6;
end
if nargin<4
    method='centered';
end


    %for df/dxi:
for i=1:length(x)
    xplus=x;
    xminus=x;
    if strcmp(method, 'centered')
        xplus(i)=xplus(i)+epsilon/2;
        xminus(i)=xminus(i)-epsilon/2;
    end
    if strcmp(method, 'backward')
        xminus(i)=xminus(i)-epsilon;
    end
    if strcmp(method,'forward')
        xplus(i)=xplus(i)+epsilon;
    end
    J_temp=(feval(func,xplus)-feval(func,xminus))/epsilon;
    if i==1 %during the first run - initialize J 
        J=zeros(length(J_temp),length(x));
    end
    J(:,i)=J_temp; %fill the Jacobian for df/dxi
end

