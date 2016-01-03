addpath('Codes');
addpath('Codes\BCT\');

N = 70;
m = 5;
weight = N/40;   %1
P = 0.02:0.04:0.20;
for j = 1:10
    j
    for pp = 1:length(P)
        pp
        A = Erdos_Renyi(N,P(pp));
        
        [BestNodes,~] = Optimize(A,m,weight,0);
        L = diag(sum(A))-A; for i = 1:m; L(BestNodes(i),BestNodes(i)) = L(BestNodes(i),BestNodes(i)) + weight; end;
        B = eig(L); R_opt(j,pp) = B(end)/B(1);
        D = sum(A); [D Di] = sort(D,'descend'); deg_ave(j,pp) = mean(D); deg_opt(j,pp) = mean(D(BestNodes)); deg_max(j,pp) = mean(D(Di(1:m)));
        L = diag(sum(A))-A; for i = 1:m; L(Di(i),Di(i)) = L(Di(i),Di(i)) + weight; end; B = eig(L); R_deg(j,pp) = B(end)/B(1);
        [EBC D] = edge_betweenness_bin(A); [D Di] = sort(D,'descend'); bet_ave(j,pp) = mean(D); bet_opt(j,pp) = mean(D(BestNodes)); bet_max(j,pp) = mean(D(Di(1:m)));
        L = diag(sum(A))-A; for i = 1:m; L(Di(i),Di(i)) = L(Di(i),Di(i)) + weight; end; B = eig(L); R_bet(j,pp) = B(end)/B(1);
        D = Closeness(A); [D Di] = sort(D,'descend'); clos_ave(j,pp) = mean(D); clos_opt(j,pp) = mean(D(BestNodes)); clos_max(j,pp) = mean(D(Di(1:m)));
        L = diag(sum(A))-A; for i = 1:m; L(Di(i),Di(i)) = L(Di(i),Di(i)) + weight; end; B = eig(L); R_clos(j,pp) = B(end)/B(1);
        D = eccentricity(A); [D Di] = sort(D,'descend'); ecc_ave(j,pp) = mean(D); ecc_opt(j,pp) = mean(D(BestNodes)); ecc_max(j,pp) = mean(D(Di(1:m)));
        L = diag(sum(A))-A; for i = 1:m; L(Di(i),Di(i)) = L(Di(i),Di(i)) + weight; end; B = eig(L); R_ecc(j,pp) = B(end)/B(1);
        D = clustering_coef_bu(A); [D Di] = sort(D,'descend'); clus_ave(j,pp) = mean(D); clus_opt(j,pp) = mean(D(BestNodes)); clus_max(j,pp) = mean(D(Di(1:m)));
        L = diag(sum(A))-A; for i = 1:m; L(Di(i),Di(i)) = L(Di(i),Di(i)) + weight; end; B = eig(L); R_clus(j,pp) = B(end)/B(1);
        
%         if R_deg(j,pp) < R_opt(j,pp)
%             disp('something wrong, Deg');
%         end
%         if R_bet(j,pp) < R_opt(j,pp)
%             disp('something wrong, Bet');
%         end
%         if R_clos(j,pp) < R_opt(j,pp)
%             disp('something wrong, Clos');
%         end
%         if R_ecc(j,pp) < R_opt(j,pp)
%             disp('something wrong, Ecc');
%         end
%         if R_clus(j,pp) < R_opt(j,pp)
%             disp('something wrong, Clus');
%         end
        
    end;
end;

Results.ER.N400.m10.Ropt = R_opt;
Results.ER.N400.m10.Rdeg = R_deg;
Results.ER.N400.m10.Rbet = R_bet;
Results.ER.N400.m10.Rclos = R_clos;
Results.ER.N400.m10.Recc = R_ecc;
Results.ER.N400.m10.Rclus = R_clus;
Results.ER.N400.m10.degree.opt = deg_opt;
Results.ER.N400.m10.degree.mean = deg_ave;
Results.ER.N400.m10.degree.max = deg_max;
Results.ER.N400.m10.bet.opt = bet_opt;
Results.ER.N400.m10.bet.mean = bet_ave;
Results.ER.N400.m10.bet.max = bet_max;
Results.ER.N400.m10.clos.opt = clos_opt;
Results.ER.N400.m10.clos.mean = clos_ave;
Results.ER.N400.m10.clos.max = clos_max;
Results.ER.N400.m10.ecc.opt = ecc_opt;
Results.ER.N400.m10.ecc.mean = ecc_ave;
Results.ER.N400.m10.ecc.max = ecc_max;
Results.ER.N400.m10.clus.opt = clus_opt;
Results.ER.N400.m10.clus.mean = clus_ave;
Results.ER.N400.m10.clus.max = clus_max;


