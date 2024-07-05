clear;
%绘制MSDF算法排序函数下的能耗图
x1=zeros(1,26);
B=1;
for kk=0:0.04:1
    x1(1,B)=kk;
    B=B+1;
end
EZMSDF=MSDFpipeifuc();
plot(x1, EZMSDF, '--b');
hold on; %不被后面的图覆盖

%绘制MADF算法排序函数下的能耗图
x4=zeros(1,26);
B=1;
for kk=0:0.04:1
    x4(1,B)=kk;
    B=B+1;
end
EZMDF=MADFpipeifuc();
plot(x4, EZMDF, '-vg');
hold on; %不被后面的图覆盖

%求均值
sumEZMDF=0;
for ii=1:26
    sumEZMDF=EZMDF(1,ii)+sumEZMDF;
end
averageEZMDF=sumEZMDF/26;
avgMDF=zeros(1,26);
B=1;
for kk=0:0.04:1
    avgMDF(1,B)=averageEZMDF;
    B=B+1;
end


%绘制MADF的均值
x3=zeros(1,26);
B=1;
for kk=0:0.04:1
    x3(1,B)=kk;
    B=B+1;
end
plot(x3, avgMDF, '--g');
hold on; %不被后面的图覆盖


%绘制MDF算法排序函数下的能耗图
x4=zeros(1,26);
B=1;
for kk=0:0.04:1
    x4(1,B)=kk;
    B=B+1;
end
EZMDF=MDFpipeifuc();
plot(x4, EZMDF, '-<r');
hold on; %不被后面的图覆盖

%求均值
sumEZMDF=0;
for ii=1:26
    sumEZMDF=EZMDF(1,ii)+sumEZMDF;
end
averageEZMDF=sumEZMDF/26;
avgMDF=zeros(1,26);
B=1;
for kk=0:0.04:1
    avgMDF(1,B)=averageEZMDF;
    B=B+1;
end


%绘制MDF的均值
x5=zeros(1,26);
B=1;
for kk=0:0.04:1
    x5(1,B)=kk;
    B=B+1;
end
plot(x5, avgMDF, '--r');
hold on; %不被后面的图覆盖

%设置横纵坐标的标题
title('能量损耗与γ'); %标题
L1=xlabel('γ'); %x轴标题
L2=ylabel('能量损耗加和(J)'); %y轴标题
legend('EZMSDF','EZMADF','(AVG.)MADF-GS','EZMDF','(AVG.)MDF-GS','Location','SouthWest ') %为图片添加图例
