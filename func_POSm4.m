function[dx]=func_POSm4(t,x,klin,knonlin)
    dx=zeros(4,1);
    %J=jacob_POSm4(t,x,klin,knonlin);
    dx(1)= klin(1)-(klin(2)*(1 + (x(4)/knonlin(1))^knonlin(2))+klin(3))*x(1);
    dx(2)= klin(2)*(1 + (x(4)/knonlin(1))^knonlin(2))*x(1) - (klin(4) + klin(5))*x(2);
    dx(3)= klin(4)*x(2) - (klin(6) + klin(7))*x(3);
    dx(4)= klin(6)*x(3) - klin(8)*x(4);
end