function Wi= weight(x)
     %xΪ�������� yΪ���󳤶�
     %���Ȩ��Wi WiΪ����Ϊ������Ŀ��������
     [~, Ei] = info_entropy(x);
     Wi=zeros();
     for i=1:size(x,1)
        Wi(i)=1-Ei(i)/log(size(x,2));
     end
end

