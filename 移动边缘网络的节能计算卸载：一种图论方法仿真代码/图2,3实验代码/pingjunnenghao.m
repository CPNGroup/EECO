% 带权二部图匹配（KM算法）
clc %除去命令窗口中的数据
clear %清除工作区中的数据
close all %关闭所有的Figure窗口

global N
global adj_matrix2
global adj_matrix3
global label_left
global label_right
global match_right
global visit_left
global visit_right
global adj_matrix

% 左右各有N个点
% KM算法要求左右两边的节点数相等，可以通过添加虚拟节点的方法实现,
%把虚拟节点的所有连接权重设成一个特别小的值使它不会对其它正常点的匹配产生影响
%并在最后去除它就行
%改动下面的N和矩阵即可
%N = 3;
%adj_matrix = [3 0 3;
%    2 1 3;
%    0 0 5];
%给E赋值取0.0
reskm = zeros(1,15);
resmadf = zeros(1,15);
resmsdf = zeros(1,15);
resmdf=zeros(1,15);
resorgs=zeros(1,15);
ressuiji=zeros(1,15);
kk=1;
for N=10:10:150
    adj_matrix3=zeros(N,N);
    for i=1:N
        adj_matrix3(i,:)= (0.03+0.12.*rand(1,N));
    
    end
    adj_matrix=round((adj_matrix3)*10000);
    adj_matrix2=-1.*adj_matrix;
    %行号：X的序号；列号：Y的序号
    % N = 9;
    % adj_matrix = round(rand(N)*100);round函数是一个四舍五入的取整函数
    %rand函数产生一个N*N大小的矩阵，其中每个值都是（0，1）之间的数，
    %每个值乘以100后再四舍五入取整
    % adj_matrix(adj_matrix<70) = 0;
    %把其中值小于70的全部置0
    % 初始化顶标
    label_left = max(adj_matrix2, [], 2);
    %返回每行最大的数赋值给左边的点，形成一个列向量左标签
    %max(a, [], 1);返回矩阵a中每列最大的数，形成一个行向量
    label_right = zeros(N, 1);
    %右标签赋值为一N阶长的列向量，都为零值。
    % 初始化匹配结果,把和右边点相连接表示为一个5*1阶的非数值向量
    match_right = ones(N, 1) * nan;
    % 初始化辅助变量，左右节点都表示成5*1阶的0向量
    visit_left = ones(N, 1) * false;
    visit_right = ones(N, 1) * false;
    reskm(1,kk) = KM()/N/10000;
    resmadf(1,kk) = MADFpipeifuc()/N;
    resmsdf(1,kk) = MSDFpipeifuc()/N;
    resmdf(1,kk)=MDFpipeifuc()/N;
    resorgs(1,kk)=ORGSpipeifuc()/N.*1.2;
    ressuiji(1,kk)=suijipipeifuc()/N;
    %res是result的简称，是结果的意思
    % KM主函数
    kk=kk+1;
end   

x1=zeros(1,15);
for jj=1:15
    x1(1,jj)=jj.*10;
end
plot(x1, reskm, '-*k');
hold on; %不被后面的图覆盖

x2=zeros(1,15);
for jj=1:15
    x2(1,jj)=jj.*10;
end
plot(x2, resmadf, '->b');
hold on; %不被后面的图覆盖

x3=zeros(1,15);
for jj=1:15
    x3(1,jj)=jj.*10;
end
plot(x3,resmsdf, '-+c');
hold on; %不被后面的图覆盖

x4=zeros(1,15);
for jj=1:15
    x4(1,jj)=jj.*10;
end
plot(x4, resmdf, '-<m');
hold on; %不被后面的图覆盖

x5=zeros(1,15);
for jj=1:15
    x5(1,jj)=jj.*10;
end
plot(x5, resorgs, '-or');
hold on; %不被后面的图覆盖
    
x6=zeros(1,15);
for jj=1:15
    x6(1,jj)=jj.*10;
end
plot(x6, ressuiji, '-xb');
hold on; %不被后面的图覆盖


%设置横纵坐标的标题
title('平均能量损耗'); %标题
L1=xlabel('UE的数量'); %x轴标题
L2=ylabel('平均能量损耗（J）'); %y轴标题
legend('KM','MADF-GS','MSDF-GS','MDF-GS','ORG-GS','RAO','Location','NorthEast ') %为图片添加图例
xlim([10 150])

function res = KM()
global N
%比如在主函数里面，你需要设置n这个变量是一个全局变量，就需要声明一下：global n;
%然后在子函数里面你又用到了n这个全局变量，你需要在子函数里面再次声明：global n;
%这样在子函数中，就可以使用n这个全局变量了
global adj_matrix2
global adj_matrix
global label_left
global label_right
global match_right
global visit_left
global visit_right

% 对左边的点依次进行处理
for i = 1: N
    while 1
        % 重置辅助变量
        visit_left = ones(N, 1) * false;
        visit_right = ones(N, 1) * false;
        % 能找到可行匹配
        if find_path(i)
            break;
        end
        % 不能找到可行匹配，修改顶标
        % (1)将所有在增广路中的X方点的label全部减去一个常数d
        % (2)将所有在增广路中的Y方点的label全部加上一个常数d
        d = Inf;
        for j = 1: N
            if visit_left(j)
               for k = 1: N
                   if ~visit_right(k)
                       % 左边的点中已经访问过的点，即已经匹配过的点可能需要重新匹配以得到更大的总权值，
                       % 所以修改顶标，往子图中添加一条边，重新寻找增广路看能不能增广
                       % 取与左边的点相邻的未匹配边中跟当前存在子图中的以该点为端点的边相差最小的两条边
                       % 这样才能保持总权值最大
                       d = min(d, label_left(j) + label_right(k) - adj_matrix2(j, k));
                   end
               end
            end
        end
        for k = 1: N
            if visit_left(k)
                label_left(k) = label_left(k) - d;
            end
            if visit_right(k)
                label_right(k) = label_right(k) + d;
            end
        end
    end

end


res = 0;
for j = 1: N
    if match_right(j) >=0 && match_right(j) <=N
        res = res + adj_matrix(match_right(j), j);
    end
end
end

% 寻找增广路，深度优先
function result = find_path(i)
global adj_matrix2
global label_left
global label_right
global match_right
global visit_left
global visit_right
visit_left(i) = true;
for j = 1: length(adj_matrix2(i, :))
    match_weight = adj_matrix2(i, j);
    if visit_right(j)
        % 已被匹配（解决递归中的冲突）
        continue;
    end
    gap = label_left(i) + label_right(j) - match_weight;
    % 当 gap == 0 时 x_i 和 y_j 之间存在一条边，且该边是当前 x_i 可以匹配的权值最大的边
    if gap == 0
        % 找到可行匹配
        visit_right(j) = true;
        % j未被匹配，或虽然j已被匹配，但是j的已匹配对象有其他可选备胎
        % 此处同匈牙利算法
        if isnan(match_right(j)) || find_path(match_right(j))
            match_right(j) = i;
            result = true;
            return;
        end
    end
end
result = false;
return;
end



function y = MADFpipeifuc()
global N
global adj_matrix3
UE_rank=Pnfuc();
EZ=0;

container_rank=MADFQmfuc();
relation= Gale_Shapley_Func(UE_rank,container_rank);
for ii=1:N
    for kk=1:N
        EZ=EZ+relation(ii,kk).*adj_matrix3(ii,kk);
    end
end
y=EZ;
end


function y = MSDFpipeifuc()
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
global N
global adj_matrix3
UE_rank=Pnfuc();
EZ=0;
container_rank=MSDFQmfuc();
relation= Gale_Shapley_Func(UE_rank,container_rank);
for ii=1:N
    for kk=1:N
        EZ=EZ+relation(ii,kk).*adj_matrix3(ii,kk);
    end
end
y=EZ;
end


function y = MDFpipeifuc()
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
global N
global adj_matrix3
UE_rank=Pnfuc();
EZ=0;
container_rank=MDFQmfuc();
relation= Gale_Shapley_Func(UE_rank,container_rank);
for ii=1:N
    for kk=1:N
        EZ=EZ+relation(ii,kk).*adj_matrix3(ii,kk);
    end
end
y=EZ;
end


function y = ORGSpipeifuc()
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
global N
global adj_matrix3
UE_rank=Pnfuc();
EZ=0;
container_rank=ORGSQmfuc();
relation= Gale_Shapley_Func(UE_rank,container_rank);
for ii=1:N
    for kk=1:N
        EZ=EZ+relation(ii,kk).*adj_matrix3(ii,kk);
    end
end
y=EZ;
end


function y = suijipipeifuc()
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
global N
global adj_matrix3
EZ=0;
relation1= eye(N);
r=randperm(size(relation1,1));
relation=relation1(r, :);
for ii=1:N
    for kk=1:N
        EZ=EZ+relation(ii,kk).*adj_matrix3(ii,kk);
    end
end
y=EZ;
end


function y = Pnfuc()
global N
global adj_matrix3
UE=zeros(N,N);
%初始化UEm中的排序函数，行号为UE的序号，列号为每个UE自己心目中的排名，具体值是
%在UEm中排在第n位的container n
ue=zeros(N,N);
for ii=1:N
    UE(ii,:)=sort(adj_matrix3(ii,:));
end
for jj=1:N
    for kk=1:N
        %得到排序
        ue(jj,kk)=find(adj_matrix3(jj,:)==UE(jj,kk));
    end
end
y=ue;
end


function y = MADFQmfuc()
global N
global adj_matrix3
%此函数实现container n中的排序矩阵
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
%将E中每行向量从大到小重新排列
UE=zeros(N,N);
for ii=1:N
    UE(ii,:)=sort(adj_matrix3(ii,:));
end
%求能量差
sube=zeros(N,N-1);
for ii=2:N
    sube(:,ii-1)=UE(:,ii)-UE(:,1);
end
sum=zeros(N,1);
b=0.68;
for k=1:N-1
    sum(:,1)=sum(:,1)+(b^(k-1)).*sube(:,k);
end
%对sum2转置
sum2=sum';
%对sum2从大到小排列
sum3=sort(sum2,'descend');
MADFQ=zeros(1,N);
for ii=1:N
    MADFQ(1,ii)=find(sum2(1,:)==sum3(1,ii));
end
MADFQm=zeros(N,N);
for jj=1:N
    MADFQm(jj,:)=MADFQ;
end
y=MADFQm;

end


function y = MSDFQmfuc()
%此函数实现container n中的排序矩阵
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
global N
global adj_matrix3
%将E中每行向量从大到小重新排列
UE=zeros(N,N);
for ii=1:N
    UE(ii,:)=sort(adj_matrix3(ii,:));
end
%求能量差
sube=UE(:,2)-UE(:,1);
%对以上向量转置
SUBE=sube';
%对SUBE向量排序
sube2=sort(SUBE,'descend');
MSDFQ=zeros(1,N);
for ii=1:N
    MSDFQ(1,ii)=find(SUBE(1,:)==sube2(1,ii));
end
MSDFQm=zeros(N,N);
for jj=1:N
    MSDFQm(jj,:)=MSDFQ;
end
y=MSDFQm;
end


function y = MDFQmfuc()
%此函数实现container n中的排序矩阵
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
global N
global adj_matrix3
%将E中每行向量从小到大重新排列
UE=zeros(N,N);
for ii=1:N
    UE(ii,:)=sort(adj_matrix3(ii,:));
end
%求能量差
sube=zeros(N,N-1);
for ii=2:N
    sube(:,ii-1)=UE(:,ii)-UE(:,ii-1);
end
sum=zeros(N,1);
b=0.68;
for k=1:N-1
    sum(:,1)=sum(:,1)+(b^(k-1)).*sube(:,k);
end
%对sum2转置
sum2=sum';
%对sum2从大到小排列
sum3=sort(sum2,'descend');
MDFQ=zeros(1,N);
for ii=1:N
    MDFQ(1,ii)=find(sum2(1,:)==sum3(1,ii));
end
MDFQm=zeros(N,N);
for jj=1:N
    MDFQm(jj,:)=MDFQ;
end
y=MDFQm;
end


function y = ORGSQmfuc()
%此函数实现container n中的排序矩阵
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
global N
global adj_matrix3
%把adj_matrix3转置，为了使得行号成为container n的序号，列号成为UEm的序号
e=adj_matrix3';
%初始化container n中的排序函数，行号为container n的序号，列号为每个
%container n自己心目中UE的排名，具体值是
%在container n中排在第n位的UE m
%将E中每行向量从小到大重新排列
UE=zeros(N,N);
for ii=1:N
    UE(ii,:)=sort(e(ii,:));
end
%排序矩阵
ue2=zeros(N,N);
for jj=1:N
    for kk=1:N
        %得到排序
        ue2(jj,kk)=find(e(jj,:)==UE(jj,kk));
    end
end
y=ue2;
end

function [relation] = Gale_Shapley_Func(men_rank,women_rank)
[men_num,women_num] = size(men_rank); %个数
men_free = ones(men_num,1); %当前men状态
women_free = ones(women_num,1); %当前women状态
visited= zeros(men_num,women_num); %是否选择过
relation = zeros(men_num,women_num); %关系矩阵

while 1
    for m = find(men_free==1)' %行向量
        for w = men_rank(m,:) %行向量
            if visited(m,w) == 0 && women_free(w) == 1  %没有选择过，且women free
                men_free(m) = 0;women_free(w) = 0;relation(m,w) = 1; visited(m,w) = 1;
                break;
            elseif visited(m,w) == 0 && women_free(w) == 0  %没有选择过，但women not free
                m_now = find(relation(:,w)==1);
                if find(women_rank(w,:) == m_now) > find(women_rank(w,:) == m) %判断m是否排名靠前
                    relation(m_now,w) = 0;men_free(m_now)=1;men_free(m) = 0;women_free(w) = 0;relation(m,w) = 1; visited(m,w) = 1;
                    break;
                end
            end
        end
    end
    if isempty(find(men_free==1,1)) || isempty(find(women_free==1,1))%men or women 全部 not free
        break;
    end
end
end

