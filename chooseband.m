function [Fj,J] = chooseband(x)
%�ú�������ںϵ���������Fj����ѡ����Ӵ�J
    xij_norm=zeros();
    Wi= weight(x);
%��һ����������    
for i = 1:size(x,1)   
    for j =1:size(x,2)
        xij_norm(i,j)=(x(i,j)-min(x(i,:)))/(max(x(i,:))-min(x(i,:)));
    end
end
    Fj=zeros();
%����ں���������
    for j=1:size(x,2)
        Fj(j)=Wi*xij_norm(:,j);
    end
%��ѡ���ֵ��Ӧ���Ӵ� ��J�Ĵ�С
    [~,J]=max(Fj);  %~Ϊ�������ֵ��JΪ���ֵ���±� ����main�з�ΧΪ1:1:10 �򷵻ص�λ�ü�ΪQֵ
end

