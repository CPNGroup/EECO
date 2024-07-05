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