function Z=ICA3(X)
%-----------ȥ��ֵ---------
[M,T] = size(X); %��ȡ����������/����������Ϊ�۲����ݵ���Ŀ������Ϊ��������      
average= mean(X')';  %��ֵ
for i=1:M
    X(i,:)=X(i,:)-average(i)*ones(1,T); 
end

%---------�׻�/��------
Cx = cov(X',1);    %����Э�������Cx
[eigvector,eigvalue] = eig(Cx); %����Cx������ֵ����������
W=eigvalue^(-1/2)*eigvector';   %�׻�����
Z=W*X;   %��������

%----------����-------
m=M;                %��Ҫ���Ƶķ����ĸ���
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
            disp('��ǰ����ѭ����');
            disp(count);
            return;
        else
            W1=W;
        end
end
Z=W'*Z;
disp('ok?');
