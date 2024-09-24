clc;clear

% matrix=cell(2,2,2);
% matrix{1,1,1}=[1;2;3];
% matrix{1,1,2}=[1;2;3];
% matrix{1,2,1}=[1;2;3];
% matrix{1,2,2}=[1;2;3];
% matrix{2,1,1}=[1;2;3];
% matrix{2,1,2}=[1;2;3];
% matrix{2,2,1}=[1;2;3];
% matrix{2,2,2}=[1;2;3];
% a=zeros(0);
% b=zeros(0);
% c=zeros(0);
%  for i=1:size(matrix,1)
%     for j=1:size(matrix,2)
%         for k=1:size(matrix,3)
%             a(end+1)=matrix{i,j,k}(1);
%             b(end+1)=matrix{i,j,k}(2);
%             c(end+1)=matrix{i,j,k}(3);
%         end
%     end
%  end
% matrix_index=[a;b;c];

b=zeros(0);
Jmax=5;J_range=1:4;
for j=1:Jmax
    if j<=length(J_range)
        b(j)=1;
    else 
        b(j)=0;
    end
end