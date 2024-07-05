function y = MDFpipeifuc()
%先给定Em,n,行号代表UEm 的序号，列号代表contianer n的序号
E=[8.9173   12.0292    7.3664    8.2316    8.7052   11.3980   10.4222   10.7436    8.5403    8.0174
    8.4361    7.5462   12.1291   12.2768   11.1483    8.4183   10.6947    8.8682    9.7114    7.9806
    8.3631    7.2585    8.7892    9.8404    9.0622   11.9727    7.8564    7.5247   12.4383   12.4225
    8.4074   10.0468   10.3038   10.8300   10.4367   10.1888   12.5146    7.0956    8.2440   10.4965
   11.4474   10.1586   12.7466   12.9157   11.7905   11.3635    7.9383   12.0940    7.3093   10.8656
    9.6513    7.0799    7.8882   10.7076   12.6019   10.8315   12.5196   11.4447    9.2719    9.1685
   12.5600    9.4477    7.0966    7.5281   11.2354   12.3965   11.9796    9.9102    8.6493    9.7876
   11.2181   11.5781   12.9494   11.3188    9.7418   10.0832    9.1313    9.9075   11.2673   11.9089
   11.4248   12.6884    9.1063   11.0073   12.3290   11.1360   11.5766   10.3139    7.4051   10.4496
    8.1963   10.1092   12.7087    7.7096   12.8117   12.1638    7.3382   11.8611    8.6457    7.1091];
UE_rank=Pnfuc();
EZ=zeros(1,26);
K=1;
for b=0:0.04:1
    container_rank=MDFQmfuc(b);
    relation= Gale_Shapley_Func(UE_rank,container_rank);
    for ii=1:10
        for kk=1:10
            EZ(1,K)=EZ(1,K)+relation(ii,kk).*E(ii,kk);
        end
    end
    K=K+1;
end

y=EZ;



