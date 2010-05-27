clear m
for i = 0:9
  m(i+1,:)=sprintf('0.%id0',i);
end
m(11,:)='1.0d0'
ca=['1.d-3'; '2.d-3'; '5.d-3'; '1.d-2'; '2.d-2']

mlamda1 = zeros(size(ca, 1),size(m, 1));
mlamda1 = zeros(size(ca, 1),size(m, 1));

for i = 1:size(ca, 1)
  for j = 1:size(m, 1)
    cd(strcat('~/Desktop/Analysis/cell/batch/planar-axi-nd5m',m(j,:),'g',ca(i,:)))
    load stretches_02000.txt
    load shpa0001.txt
    mlambda1(i,j) = shpa0001(:,4)'*stretches_02000(:,2)/sum(shpa0001(:,4));
    mlambda2(i,j) = shpa0001(:,4)'*stretches_02000(:,3)/sum(shpa0001(:,4));
  end
end
