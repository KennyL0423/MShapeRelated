function [RI_test_combined,RI_test_alone,NMI2_test_combined,NMI2_test_alone] = compute_cluster_index(Y_true,Parameter,A_tp1,N_train)
Y_true_matrix=reshape_y_ture(Y_true,Parameter.C);
[~, Y_true_matrix2]=max(Y_true_matrix,[],1);

G = graph(sparse(A_tp1));  % 将稀疏矩阵 A_tpl 转换为无向图
A_y = conncomp(G);         % 计算连通分量
clusternum = max(A_y);  % 计算连通分量的数量
[~,~,NMI2_train_combined] = compute_nmi(A_y(:,(1:N_train))',Y_true_matrix2(:,(1:N_train))');
[~,~,NMI2_test_combined] = compute_nmi(A_y(:,(1+N_train:end))',Y_true_matrix2(:,(1+N_train:end))');
A_y2=reshape_y_ture(A_y',clusternum);
RI_train_combined = RandIndex(A_y2(:,(1:N_train)),Y_true_matrix(:,(1:N_train))) ;
RI_test_combined = RandIndex(A_y2(:,(1+N_train:end)),Y_true_matrix(:,(1+N_train:end))) ;

G = graph(sparse(A_tp1((1:N_train), (1:N_train))));  % 将稀疏矩阵转换为无向图
A_y = conncomp(G);  % 计算连通分量
clusternum = max(A_y);  % 获取连通分量的数量

% [clusternum, A_y]=graphconncomp(sparse(A_tp1((1:N_train),(1:N_train))));
[~,~,NMI2_train_alone] = compute_nmi(A_y',Y_true_matrix2(:,(1:N_train))');
A_y2=reshape_y_ture(A_y',clusternum);
RI_train_alone = RandIndex(A_y2,Y_true_matrix(:,(1:N_train))) ;

G = graph(sparse(A_tp1((1+N_train:end), (1+N_train:end))));  % 将稀疏矩阵转换为无向图
A_y = conncomp(G);  % 计算连通分量
clusternum = max(A_y);  % 获取连通分量的数量

% [clusternum, A_y]=graphconncomp(sparse(A_tp1((1+N_train:end),(1+N_train:end))));
A_y2=reshape_y_ture(A_y',clusternum);
[~,~,NMI2_test_alone] = compute_nmi(A_y',Y_true_matrix2(:,(1+N_train:end))');
RI_test_alone = RandIndex(A_y2,Y_true_matrix(:,(1+N_train:end))) ;
end