function plot_graph(tarray,arry, pfunc,str_title)

plot(tarray{4},pfunc{4},tarray{2},arry{2},tarray{3},arry{3},tarray{4},arry{4},tarray{5},arry{5},tarray{6},arry{6});
 legend('Exact Solution','t=0.5', 't=0.25','t=0.125','t=0.0625','t=0.0313');
 title(str_title);
 axis([0 5 0 20]);
end