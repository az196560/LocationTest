function Z=ICA3(X)
%-----------去均值---------
[M,T] = size(X); %获取输入矩阵的行/列数，行数为观测数据的数目，列数为采样点数      
average= mean(X')';  %均值
for i=1:M
    X(i,:)=X(i,:)-average(i)*ones(1,T); 
end

%---------白化/球化------
Cx = cov(X',1);    %计算协方差矩阵Cx
[eigvector,eigvalue] = eig(Cx); %计算Cx的特征值和特征向量
W=eigvalue^(-1/2)*eigvector';   %白化矩阵
Z=W*X;   %正交矩阵

%----------迭代-------
m=M;                %需要估计的分量的个数
W=rand(m);
W1=zeros(4,4);

WP1=W(:,1);
WP1=WP1/norm(WP1);
WP2=W(:,2);
WP2=WP2/norm(WP2);
WP3=W(:,3);
WP3=WP3/norm(WP3);
WP4=W(:,4);
WP4=WP4/norm(WP4);

count=0;


W5=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
while count<70
        count=count+1;
        count1=0;
        WP1=1/T*Z*(((WP1)'*Z).^3)'-3*WP1;
        WP2=1/T*Z*(((WP2)'*Z).^3)'-3*WP2;
        WP3=1/T*Z*(((WP3)'*Z).^3)'-3*WP3;
        WP4=1/T*Z*(((WP4)'*Z).^3)'-3*WP4;
        W=[WP1,WP2,WP3,WP4];
        while norm(W*W'-W5)>0.0001
           W=1.5*W-0.5*W*W'*W;
           count1=count1+1;
        end
%        W=(W*W')^(-0.5)*W;
        WP1=W(:,1);
        WP2=W(:,2);
        WP3=W(:,3);
        WP4=W(:,4);
        if( abs(W1-W)&abs(W1+W)<0.001)
            disp('提前调出循环了');
            disp(count);
            return;
        else
            W1=W;
        end
end
Z=W'*Z;
disp('ok?');
