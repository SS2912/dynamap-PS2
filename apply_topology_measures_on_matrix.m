function results = apply_topology_measures_on_matrix(conn_mat, cfg, symmetric, graphs)

tic
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
results = [];

for matrix_iter = 1 : size(conn_mat,3)

    clearvars -except matrix_iter symmetric results conn_mat cfg graphs

    Wmatrix = conn_mat(:,:,matrix_iter);
    nb_chan = size(Wmatrix,1);
    
    %%% faire attention matrice symetrique ou pas???
    if symmetric
        if ~issymmetric(Wmatrix)
            for i = 1 : length(Wmatrix)
                for j = 1 : length(Wmatrix)
                    if abs(Wmatrix(j,i))>=abs(Wmatrix(i,j))
                        Wmatrix(i,j) = Wmatrix(j,i);
                    end
                end
            end
        end
    end
    Wmatrix = abs(Wmatrix); 
    if graphs 
        figure; 
        imagesc(Wmatrix);
    end
    
    %% measure calculation
    
    % cleaning of the matrix 
    WU = weight_conversion(Wmatrix, 'autofix'); % deletes Inf and Nan in matrix, plus corrects rounding errors in symmetrization (if maatrix is symm)
    
    % normalized strength
    norm_strength = sum(WU,1)/(size(WU,2)-1)';  
    strenght = strengths_und(WU); % whats the difference between normalised and non-norm strength???
    
    % degree
    BU = WU>= cfg.deg_thr; % binarize matrix (BU= Binary Undirected matrix)
    degrees = sum(BU,1)';
    deg = degrees_und(BU);
    
    % centrality
        %   Conversion of connection weights to connection lengths is needed prior to computation of weighted distance-based measures, such as
        %   distance and betweenness centrality. In a weighted connection network, higher weights are naturally interpreted as shorter lengths. The
        %   connection-lengths matrix here is defined as the inverse of theconnection-weights matrix. E=find(W); W(E)=1./W(E); %E=index of non-zero elements in W
   
    L = weight_conversion(WU, 'lengths'); % converts weights to length 
    betweeness_centrality = betweenness_wei(L);
    edge_betweeness = edge_betweenness_wei(L);
    eig_vec_centrality = eigenvector_centrality_und(WU); % Eigenector centrality is a self-referential measure of centrality -- nodes have high eigenvector centrality if they connect to other nodes that have high eigenvector centrality.
    
    
    % clustering coefficient
    clust_coeff = clustering_coef_wu(WU); % The clustering coefficient is the fraction of triangles around a node and is equivalent to the fraction of node s neighbors that are neighbors of each other.
    
    % distance 
    [D, nb_edges] = distance_wei(L); % shortest weighted path matrix. The distance matrix contains lengths of shortest paths between all pairs of nodes. An entry (u,v) represents the length of shortest path from node u to node v. The average shortest path length is the characteristic path length of the network.
   
    % reachability and binary distance matrices
    [reachM, DB] = reachdist(BU); 
    
    % efficiency
    [char_path_len, glob_efficiency] = charpath(D); % The global efficiency is the average inverse shortest path length in the network, and is inversely related to the characteristic path length. Infinitely long paths (i.e. paths between disconnected nodes) are included in computations by default. This behavior may be modified with via the infinite_dist argument.
    loc_efficiency = efficiency_wei(WU,2);  % The local efficiency is the global efficiency computed on the neighborhood of the node, and is related to the clustering coefficient.
    
    % community 
    Q0 = -1; Q1 = 0;            % initialize modularity values
      while Q1-Q0>1e-5           % while modularity increases
          Q0 = Q1;                % perform community detection
          [com_louv, Q1] = community_louvain(WU,1,[],'modularity');
      end

    modularity = Q1;
%     [com_louv, modularity] = community_louvain(WU,1,[],'modularity');  % The optimal community structure is a subdivision of the network into nonoverlapping groups of nodes in a way that maximizes the number of within-group edges, and minimizes the number of between-group edges. The modularity is a statistic that quantifies the degree to which the network may be subdivided into such clearly delineated groups.
    [~,I] = sort(com_louv);
    reord_com_louv = WU(I,I);
%     figure; 
%     subplot(1,2,1)
%     imagesc(WU)
%     subplot(1,2,2)
%     imagesc(reord_com_louv)
%     xticks(1:length(I))
%     yticks(1:length(I))
%     xticklabels(num2str(I))
%     yticklabels(num2str(I))
%     title('louvain community')
    
    
    % assortativity and core structure
    assortativity = assortativity_wei(WU,1); % The assortativity coefficient is a correlation coefficient between the degrees of all nodes on two opposite ends of a link. A positive assortativity coefficient indicates that nodes tend to link to other nodes with the same or similar degree.
    rich_club = rich_club_wu(WU); % The rich club coefficient at level k is the fraction of edges that connect nodes of degree k or higher out of the maximum number of edges that such nodes might share
   
    
   
    %% output results 
    
    results.norm_strength(:, matrix_iter) = norm_strength';
    results.degrees(:, matrix_iter) = degrees;

    results.betweeness_centrality(:, matrix_iter) = betweeness_centrality;
    results.edge_betweeness(:, :, matrix_iter) = edge_betweeness;
    results.eig_vec_centrality(:,matrix_iter) = eig_vec_centrality;

    results.characteristic_path_length(matrix_iter) = char_path_len;
    results.distance_mat(:,:, matrix_iter) = D;
    results.nb_edges_distance_m(:,:, matrix_iter) = nb_edges;
    results.reachability_mat(:,:, matrix_iter) = reachM;
    results.distance_binary(:,:, matrix_iter) = DB;

    results.clust_coef(:, matrix_iter) = clust_coeff;
    results.loc_efficiency(:, matrix_iter) = loc_efficiency;
    results.global_efficiency(matrix_iter) = glob_efficiency;

    results.com_louv(:, matrix_iter) = com_louv;
    results.reord_com_louv(:,:, matrix_iter) = reord_com_louv;
    results.modularity(matrix_iter) = modularity;
    results.rich_club(:, matrix_iter) = rich_club';
    results.assortativity(matrix_iter) = assortativity;
    
end

toc
end

