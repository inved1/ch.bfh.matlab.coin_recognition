%Let's make some fake data with two groups
n=150; %sample size
x=[randn(n,1)+3;randn(n,1)+3];
y=[randn(n,1)+3;randn(n,1)+3];

%plot the data
subplot(1,2,1)
plot(x,y,'ok','MarkerFaceColor','k')

%now divide into two using k-means and plot the results
data2=[x,y];
IDX=kmeans(data2,2);

%plot the k-means results
subplot(1,2,2)
scatter(x,y,50,IDX,'filled')