function plot_3Dgraph(X, C, W, y)

W(W <1e-5) = 0;

[iidx, jidx, val] = find(sparse(W));
% fprintf ('%d, %d, %f\n', [iidx, jidx, val]');

figure;
hold on;
plot3(C(1,:), C(2,:), C(3,:),'or');
for i=1:length(iidx)
    
    plot3( [C(1, iidx(i)), C(1, jidx(i))], [C(2, iidx(i)), C(2, jidx(i))],...
        [C(3, iidx(i)), C(3, jidx(i))], '-ok','LineWidth',0.5);%, 'LineWidth', W(iidx(i),jidx(i))*10 );
    hold on;

end


cls = unique(y);
ncls = length(cls);

colors = distinguishable_colors(ncls);
names={};
for c = 1:ncls
    idx = find(y==cls(c));
    names{c}=int2str(cls(c));
    h(c) = plot3(C(1,idx), C(2,idx), C(3,idx),'o','Color',...
        colors(c,:), 'MarkerSize', 10,'MarkerFaceColor',colors(c,:));
end
legend(h, names,'FontSize', 26,'Location','bestoutside','Orientation','vertical');