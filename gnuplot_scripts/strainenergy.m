clear r
for i = 0:9
  r(i+1,:)=sprintf('0.%id0',i);
end
r(11,:)='1.0d0'
ca=['1.d-2'];
w = zeros(size(ca, 1),size(r, 1));
for i = 1:size(ca, 1)
  for j = 1:size(r, 1)
%    cd(strcat('~/Documents/Linear flow/mix20100507-psc-shadow/batch/bp-ob-l5r',r(j,:),'g',ca(i,:), 't0.00d0p0.50d0'))
    cd(strcat('~/Documents/Linear flow/mix20100507-psc-shadow/batch/bp-sp-l5r',r(j,:),'g',ca(i,:)))
    load stretches_02000.txt
    load shpa0001.txt
    for k=1:size(stretches_02000,1)
      w(i,j) = w(i,j) + (stretches_02000(k,2)^2+stretches_02000(k,3)^2+stretches_02000(k,2)^-2*stretches_02000(k,3)^-2-3)*shpa0001(k,4)/3;
    end
    w(i,j) = w(i,j)/sum(shpa0001(:,4));
  end
end
cd ../../gnuplot_scripts
w'
