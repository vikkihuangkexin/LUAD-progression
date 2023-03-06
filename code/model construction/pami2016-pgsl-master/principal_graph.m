function [C, W, P, objs] = principal_graph(X, C0, G, params)
%INPUT:
% X: the input data DxN
% G: graph matrix with side information where cannot-link pair is 0
% C0: the initialization of centroids

% params = struct('maxiter', 10, ...
%                 'eps', 1e-5, ...
%                 'gstruct', 'l1-graph',...
%                 'lambda', 1.0, ...
%                 'gamma', 0.5, ...
%                 'sigma', 0.01, ...
%                 'nn', 5, ...
%                 'verbose',true);

%OUTPUT:
% C: centers for principal graph
% W: principal graph matrix
% P: cluster assignment matrix

C = C0;
[D,K] = size(C);

% construct only once for all iterations
if strcmp(params.gstruct,'l1-graph')
    
    % low triangular sum_i sum_{j < i}
    [row, col] = find(tril(G));
    nw = length(row);
    nvar = nw + K*D;

    rc = java.util.HashMap;
    for i=1:nw
        key_ij = row(i) + col(i)*K;
        rc.put(key_ij, i);
    end
    
    % construct A and b
    A=sparse([]);
    b=sparse([]);
    for i=1:K
        nn_i = find(G(:,i)==1);
        a = sparse(zeros(2*D,nvar));
        for jj=1:length(nn_i)
            j = nn_i(jj);
            key_ij = i+j*K;
            if i<j
                key_ij = j + i*K;
            end
            pos_ij = rc.get(key_ij);
            a(:,pos_ij) = [-X(:,j);X(:,j)];
        end
        start_i = nw + (i-1)*D + 1;
        end_i = start_i + D-1;
        a(:, start_i:end_i ) = -[ eye(D,D); eye(D,D) ];
        A = [A;a];
        b = [b;-X(:,i);X(:,i)];
    end
end

objs = [];
lp_vars = []; 
for iter=1:params.maxiter
    
    norm_sq = repmat(sum(C.^2, 1), K, 1);
    Phi = norm_sq + norm_sq' - 2 .* C' * C;
    switch(params.gstruct)
        case 'l1-graph'
            val = zeros(nw,1);
            for i=1:nw
                val(i) = Phi(row(i), col(i));
            end
            f = [2.*val; params.lambda .* ones(K*D,1)];
            
            % MATLAB solver
            options = optimset( 'Display', 'off','Algorithm','interior-point');
            [w_eta, obj_W] = linprog(f, A, b, [], [], [zeros(nw, 1); -Inf.*ones(K*D,1)], [], lp_vars, options);
            
%             % Mosek solver
%             prob.c = f; prob.a = A;
%             prob.buc = b;
%             prob.blx = sparse( [zeros(nw, 1); -Inf.*ones(K*D,1)] );
%             [r,res] = mosekopt('minimize echo(0)',prob);  
%             w_eta = res.sol.bas.xx;
%             obj_W = f'*w_eta;
            
            % recover results
            w = w_eta(1:nw);
            W_tril = sparse(row, col, w, K, K);
            W = W_tril + W_tril';
            
            % warm start
            lp_vars = w_eta;
            
        case 'span-tree' 
            stree = graphminspantree(sparse(tril(Phi)),'Method','Kruskal');
            stree = stree + stree';
            W = stree ~= 0;
            obj_W = sum(sum(stree));
            
        otherwise
            warning('graph structure %s is not supported yet.', params.gstruct);
            return            
    end

% for intermediate results output
%     eps_filename = sprintf('results/converge/tree_300_%d.eps', iter);
%     plot_graph(X,C,W,eps_filename);
    
    [P, obj_P] = soft_assignment(X, C, params.sigma);
    
    obj = obj_W + params.gamma * obj_P;
    objs = [objs; obj];
    if params.verbose
        fprintf('iter=%d, obj=%f\n',iter, obj);
    end
    
    if iter > 1
        relative_diff = abs( objs(iter-1) - obj) / abs(objs(iter-1));
        if relative_diff < params.eps
            if params.verbose
                fprintf('eps=%f, converge.\n', relative_diff);
            end
            break;
        end
        if iter >= params.maxiter
            if params.verbose
                fprintf('eps=%f, reach maxiter.\n', relative_diff);
            end            
        end
    end
    
    C = generate_centers(X, W, P, params.gamma);
end

function C = generate_centers(X, W, P, gamma)

[D, N] = size(X);
K = size(W,1);
% prevent singular
Q = 2.*( diag(sum(W)) - W ) + gamma .* diag(sum(P));% + 1e-10.*eye(K,K);
B = gamma .* X * P;
C = B/Q;

function [P, obj] = soft_assignment(X, C, sigma)
% sigma: bandwidth parameter
% C: centers DxK
% X: input data DxN

[D, N] = size(X);
K = size(C,2);
norm_X_sq = repmat(sum(X.^2,1)', 1, K);
norm_C_sq = repmat(sum(C.^2,1), N, 1);
dist_XC = norm_X_sq + norm_C_sq - 2.* X' * C;

% %% direct compute P encountering 0/0 numerical issue
% Phi_XC = exp(- dist_XC ./ sigma);
% P = Phi_XC ./ repmat(sum(Phi_XC, 2), 1, K);
% tmp_obj = log( sum( exp(- dist_XC./sigma) ,2) );
% obj = - sigma * sum( tmp_obj );

%% handle numerical problems 0/0 for P
min_dist = min(dist_XC,[], 2);
dist_XC = dist_XC - repmat(min_dist, 1, K );
Phi_XC = exp(- dist_XC ./ sigma);
P = Phi_XC ./ repmat(sum(Phi_XC, 2), 1, K);

obj = - sigma * sum( log( sum( exp(- dist_XC./sigma) ,2) ) ...
        - min_dist(:,:) ./ sigma );

% fprintf('obj=%f, obj1=%f\n',obj,obj1);
