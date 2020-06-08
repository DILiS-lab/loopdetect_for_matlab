function[J]=numerical_jacobian_complex(func,x,epsilon)
% NUMERICAL_JACOBIAN_COMPLEX Determine Jacobian matrix via finite 
% difference with complex step at certain variable values.
%
% Derivative approximation using a complex step approach (a small imaginary
% number is added) delivers more exact values than normal finite 
% differences that use real steps (Martins et al., 2003). It is exploited 
% that the partial deriavtive df/dx_k can be approximated by 
% Imag(f(x+h*e_k*i))/h for h close to zero. Thereby, no difference has to 
% be computed avoiding subtractive cancellation errors and thus improving 
% computation accuracy.
%
% J = NUMERICAL_JACOBIAN_COMPLEX(FUNC,X) computes the Jacobian matrix of 
% the ODE y'=FUNC(y) (i.e. the partial derivatives of FUNC) at variable 
% value column vector X using the complex step approach. FUNC is a function
% handle and depends only on the variable values. The size of the matrix 
% is length(FUNC(X)) times length(X), and is quadratic for ODE systems.
% 
% J = NUMERICAL_JACOBIAN(FUNC,X,EPSILON) allows to tune the stepsize
% EPSILON (default: 1e-6).
% 
% Example
%   klin=[165,0.044,0.27,550,5000,78,4.4,5.1];
%   knonlin=[0.3,2];
%   J = numerical_jacobian_complex(@(x)func_POSm4(1,x,klin,knonlin),...
%       [1 2 1 1])
%
% See also: numerical_jacobian()
% 
% References: Martins JRRA, Sturdza P, Alonso JJ. The complex-step 
% derivative approximation. ACM Trans Math Softw. 2003;29(3):245–62.


%set default value for the perturbation
if nargin<3
    epsilon= 1e-6;
end

    %for df/dxi:
for k=1:length(x)
    xplus=x;
    xplus(k)=xplus(k)+1i*epsilon;
    J_temp=imag(feval(func,xplus))/epsilon;
    if k==1 %during the first run - initialize J 
        J=zeros(length(J_temp),length(x));
    end
    J(:,k)=J_temp; %fill the Jacobian for df/dxi
end

