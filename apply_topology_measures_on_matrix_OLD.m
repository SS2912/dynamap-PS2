function results = apply_topology_measures_on_matrix(conn_mat,cfg)

% function to calculate and return several graph measure 
% input : a connectivity matrix as h2, or other connectivity methods
%       - conn_mat : matrix of connectivity
%       - cfg : structure contianing : 
%           - deg_thr : threshold for degree calculation
%
% results = apply_topology_measures_on_matrix(conn_mat,cfg)
%
% SMV 21/04/2023
% dependencies : Brain connectivity toolboxe
% https://sites.google.com/site/bctnet/network-construction

if isempty(cfg) || ~isfield(cfg,'deg_thr')
    cfg.deg_thr = 0.2;
end

% h2 = load('C:\Users\samuel\Documents\Data\hje360\EEG\64V_EMG_0004_resample_to_run_4.ades_algo-R2_hp-1Hz_lp-70Hz.mat');
% conn_mat = mean(h2.aw_h2,3);

% conn_mat = randn(20,20);
% conn_mat = conn_mat/max(abs(conn_mat(:)));

figure; 
imagesc(conn_mat);

nb_chan = size(conn_mat,1);

%%% faire attention matrice symetrique ou pas???
if issymmetric(conn_mat)
end
% normalized strength
new_conn_mat = weight_conversion(conn_mat, 'autofix');
norm_strength = sum(new_conn_mat,1)/(size(new_conn_mat,2)-1);

% degree
new_deg_mat = new_conn_mat>= cfg.deg_thr;
degrees = sum(new_deg_mat,1);

% betweeness
L = weight_conversion(new_conn_mat, 'lengths');
betweeness_centrality = betweenness_wei(L);
edge_betweeness = edge_betweenness_wei(L);

% clustering coefficient
clust_coeff = clustering_coef_wu(new_conn_mat);

% efficiency 
loc_efficiency = efficiency_wei(new_conn_mat);

% community 
com_louv = community_louvain(new_conn_mat,1,[],'negative_sym');
[B,I] = sort(com_louv);
reord_com_louv = new_conn_mat(I,I);
figure; 
subplot(1,2,1)
imagesc(new_conn_mat)
subplot(1,2,2)
imagesc(reord_com_louv)

% rich club
rich_club = rich_club_wu(new_conn_mat);

% shortest path length
[shortest_path_length,nb_shortest_path_length] = distance_wei(L);

% eigen vector centrality 
eig_vec_centrality = eigenvector_centrality_und(new_conn_mat);

%% output results 
results = [];
results.norm_strength = norm_strength';
results.degrees = degrees';
results.betweeness_centrality = betweeness_centrality;
results.edge_betweeness = edge_betweeness;
results.clust_coeff = clust_coeff;
results.loc_efficiency = loc_efficiency;
results.com_louv = com_louv;
results.reord_com_louv = reord_com_louv;
results.rich_club = rich_club;
results.eig_vec_centrality = eig_vec_centrality;
results.shortest_path_length = shortest_path_length;
results.nb_shortest_path_length = nb_shortest_path_length;
