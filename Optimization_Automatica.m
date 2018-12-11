%run SDP (dimension 2*2 here)
clear all
cvx_solver sedumi
n = 2;
tol = 10^(-3);
alpha = 0:0.1:3;
solution = [];
epsilon = 1;
vmax = 25;
Kpimin = 0;


for i = 1:size(alpha,2)
    
    alphai = alpha(i);

    cvx_begin sdp quiet
        variables gi Kvi Kvi0 Kpi0 c d 
        minimize -gi
        subject to:
            Jii0a = Kpimin*[-alphai alphai^2;-1 alphai]+[-alphai*(Kpi0) 1+alphai^2*(Kpi0)-alphai*(Kvi*(1+epsilon)+Kvi0*(1+epsilon));-Kpi0 alphai*(Kpi0)-Kvi*(1+epsilon)-Kvi0*(1+epsilon)];
            Jii0b = gi*[-alphai alphai^2;-1 alphai]+[-alphai*(Kpi0) 1+alphai^2*(Kpi0)-alphai*(Kvi*(1+epsilon)+Kvi0*(1+epsilon));-Kpi0 alphai*(Kpi0)-Kvi*(1+epsilon)-Kvi0*(1+epsilon)];
            Jii1a = Kpimin*[alphai -alphai^2; 1 -alphai]+[0 alphai*Kvi;0 Kvi];
            Jii1b = gi*[alphai -alphai^2; 1 -alphai]+[0 alphai*Kvi;0 Kvi];
            (Jii0b+Jii0b.')/2-c*eye(n) <= -eye(n)*tol;
            (Jii0a+Jii0a.')/2-c*eye(n) <= -eye(n)*tol;
            [d*eye(n) Jii1a;Jii1a.' d*eye(n)] >= 0;
            [d*eye(n) Jii1b;Jii1b.' d*eye(n)] >= 0;
            gi >= 0.1;
            Kvi >= 0.1;
            Kvi0 >= 0.1;
            Kpi0 >= 0.1;
            Kpi0 <= 0.6;
            c <= -tol;
            c+(1+epsilon)*d <= -tol;
        
    cvx_end 
    if cvx_optval ~= Inf 
        xa=max(eig((Jii0b+Jii0b.')/2));
        ya=max(svd(Jii1a));
        solution = [solution; alphai c d xa+(1+epsilon)*ya gi Kvi Kpi0 Kvi0];
    end

end
figure;
plot(solution(:,1),solution(:,4:end),'linewidth',1);
xlabel('alpha')
legend('Result','g_i','K_{vi}','K_{pi0}','K_{vi0}');