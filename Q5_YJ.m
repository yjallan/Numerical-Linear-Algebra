clc
clear all
close all

%% GENERATING PLOTS
for trial=1:3
    T1=[];T2=[];T3=[];
    for n=100:100:500        
        % creating random A based on d=single value information
        [U, ~] = qr(randn(n));   
        [V, ~] = qr(randn(n));    
        S = diag(3.^(-1:-1:-n));  
        A = U*S*V;    

        % Classical Gram-Schmidt Algorithm
        t1=cputime;
        [Q_CGS,R_CGS]=CGS_YJ(A);        
        T1=[T1,cputime-t1];               

        % Modified Gram-Schmidt Algorithm
        t2=cputime;
        [Q_MGS,R_MGS]=MGS_YJ(A);  
        T2=[T2,cputime-t2];
        
        % Householder Algorithm        
        t3=cputime;
        [Q_HH,R_HH]=Householder_YJ(A);
        T3=[T3,cputime-t3];        
        
    end

    hold on
    f=figure('units','normalized','outerposition',[0 0 1 1]);
    plot(100:100:500,T1,'-',100:100:500,T2,':',100:100:500,T3,'--');
    set(findall(gca, 'Type', 'Line'),'LineWidth',6);
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    xlabel('Dimension n-->');
    ylabel('CPU Time-->');
    legend({'Classical Gram Schmidt','Modified Gram Schmidt','Householder QR'},'FontSize',12);
    title(['CPU time Comparison for Different QR: Trial ',num2str(trial)]);
    saveas(f,['CPUTime_QR',num2str(trial),'.jpg']);
    hold off
    close all
end

%% TEST OF ORTHOGONALITY

clc
clear all
close all

n=500;
[U, ~] = qr(randn(n));   
[V, ~] = qr(randn(n));    
S = diag(3.^(-1:-1:-n));  
A = U*S*V;

[Q_MGS,R_MGS]=MGS_YJ(A); % modified gram Scmidt
[Q_HH,R_HH]=Householder_YJ(A); % Householder


T = Q_MGS; % assign whichever matrix you want to test.
T = Q_HH; % assign whichever matrix you want to test.
T_inv = inv(T);
T_tran = T';

tolerance=10^(-9);

error=abs(T_inv-T_tran);
if sum(error(:))<tolerance
   disp('Unitary')
else
   disp('NOT Unitary')
end