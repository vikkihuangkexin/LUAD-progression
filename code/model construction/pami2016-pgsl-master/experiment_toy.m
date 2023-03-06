function experiment_toy
addpath('toy')

filenames = {'Circle','two_moon','tree_300','Spiral','three_clusters','DistortedSShape'};
for i=1:length(filenames)
    filename = filenames{i};

    % l1-graph
    eps_filename = sprintf('results/l1_%s.eps',filename);
    params = struct('maxiter', 20, ...
                    'eps', 1e-5, ...
                    'gstruct', 'l1-graph',...
                    'lambda', 1.0, ...
                    'gamma', 0.5, ...
                    'sigma', 0.01, ...
                    'nn', 5, ...
                    'verbose',true);
    
    run_our_method(filename, eps_filename, params);
    
    % span-tree
    eps_filename = sprintf('results/sp_%s.eps',filename);
    params = struct('maxiter', 20, ...
                    'eps', 1e-5, ...
                    'gstruct', 'span-tree',...
                    'lambda', 1.0, ...
                    'gamma', 0.5, ...
                    'sigma', 0.01, ...
                    'nn', 5, ...
                    'verbose',true);

    run_our_method(filename, eps_filename, params);
    
end

function run_our_method(filename, eps_filename, params)

switch (filename)
    case 'Circle'
        % circle
        load('Circle'); X = X';
        params.sigma = 0.1; params.nn=10;
        
    case 'two_moon'
        % two_moon
        load('two_moon'); X = X';
        params.lambda = 0.1; params.gamma = 3;
        
    case 'tree_300'
        % tree
        load('tree_300'); X = X';
        params.gamma = 10; 
            
    case 'Spiral'
        % spiral
        load('Spiral_SUN'); X = X';
        params.nn=10;
        
    case 'three_clusters'
        % three_clusters
        load('three_clusters'); X = X';
        params.lambda = 0.1;
        
    case 'DistortedSShape'
        % distrotedSShape
        load('DistortedSShape'); X = X';
        
    otherwise
        warning('unexpected data settings for %s', name);
        return
end
[D,N] = size(X);

Z = X;
C0 = Z;
Nz = size(C0,2);

if params.nn<N
    G = get_knn(C0, params.nn);
else
    G = ones(Nz,Nz) - eye(Nz,Nz);
end

% run the proposed method
start_time = cputime;
[C, W, P, objs] = principal_graph(X, C0, G, params);
elapse_time = cputime - start_time;
fprintf('time cost is: %f sec\n',elapse_time);

% plot results
plot_graph_toy(X, C, W, eps_filename)

%% print convergence
% h = figure; 
% plot(objs,'LineWidth',2);
% set(gca, 'FontSize',16);
% xlabel('iteration')
% ylabel('objective value');
% 
% converge_filename = sprintf('results/converge/converge_%s.eps',filename);
% print(h, '-depsc',  converge_filename);
% close(h);

function plot_graph_toy(X, C, W, eps_filename)

W(W <1e-5) = 0;

[iidx, jidx, val] = find(sparse(W));
% fprintf ('%d, %d, %f\n', [iidx, jidx, val]');
h=figure;
box on;
hold on;
for i=1:length(iidx)
    if ~isnan(val(i))
        plot( [C(1, iidx(i)), C(1, jidx(i))], [C(2, iidx(i)), C(2, jidx(i))],...
            'ko-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k' );
    end
end
plot(C(1,:), C(2,:),'ko');
plot(X(1,:), X(2,:), '.b','MarkerSize',8);
xlim([ min(X(1,:)), max(X(1,:)) ]);
ylim([ min(X(2,:)), max(X(2,:)) ]);
set(gca, 'FontSize',16);

print(h, '-depsc',  eps_filename);

