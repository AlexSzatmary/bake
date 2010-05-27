clear m
paths={'~/Desktop/Analysis/mix20091209/batch/planar-axi6m0.0d0g2.d-2';... 
'~/Desktop/Analysis/mix20091209/batch/planar-axi6m0.5d0g2.d-2';... 
'~/Desktop/Analysis/mix20091209/batch/planar-axi6m1.0d0g2.d-2';... 
'~/Desktop/Analysis/mix20091209/batch/planar-biaxi6m0.5d0g2.d-2';...
'~/Desktop/Analysis/mix20091209/batch/planar-biaxi6m1.0d0g2.d-2'}
titles={'Planar'; 'Uniaxial-Planar r=0.5'; 'Uniaxial';...
        'Biaxial-Planar r=0.5'; 'Biaxial'}
isubplot=[233 232 231 235 234];

for i = 1:size(paths, 1)
  cd(paths{i});
  load cappro0001_002000.txt;
  load uvwpdump__02000.txt;
  center = mean(cappro0001_002000,1);
  theta = atan((cappro0001_002000(:,2)-center(2))./(cappro0001_002000(:,1)-center(1)));
  for j = 1:size(cappro0001_002000(:,1))
    if cappro0001_002000(j,1)< center(1)
      theta(j) = theta(j) + pi;
    end
  end
  yz = sortrows([cappro0001_002000 theta], 3);
  range=2:2:64;
  range2=1:64;
  [xr, yr] = meshgrid(range,range);
  [xr2, yr2] = meshgrid(range2+0.5,range2+0.5);
  p=reshape((uvwpdump__02000(32*64*64+1:33*64*64,7)+uvwpdump__02000(32*64*64+1:33*64*64,7))/2,64,64);
  v=reshape((uvwpdump__02000(32*64*64+1:33*64*64,5)+uvwpdump__02000(32*64*64+1:33*64*64,5))/2,64,64);
  w=reshape((uvwpdump__02000(32*64*64+1:33*64*64,6)+uvwpdump__02000(32*64*64+1:33*64*64,6))/2,64,64);
  k = convhull(yz(:,1), yz(:,2));
  subplot(isubplot(i));
  plot(yz(k,1), yz(k,2)+1, 'k')
  title(titles{i})
  axis([1,64,1,64])
  axis off
  hold on
  for j = 1:63
    for k = 1:63
      p2(j, k) = mean([p(j+1, k+1), p(j, k+1), ...
        p(j+1, k), p(j, k)]);
    end
  end
  for j = 1:63
    p2(j, 64) = mean([p(j+1, 1), p(j, 1), ...
      p(j+1, end), p(j, end)]);
  end
  for k = 1:63
    p2(64, k) = mean([p(1, k+1), p(end, k+1), ...
      p(1, k), p(end, k)]);
  end
  p2(64, 64) = mean([p(1,1), p(1,end), p(end,1), p(end, end)]);
%[c, h] = contour(sqrt(v.*v+w.*w));
[c, h] = contour(p2);
%  clabel(c, h, 'Rotation', 0);
  quiver(xr, yr, v(range, range), w(range, range), sqrt(max(max(v.*v+w.*w))/8.318074524560400e-03), 'k')
  hold off
  get(gca,'PlotBoxAspectRatio');
  set(gca,'PlotBoxAspectRatio', ans)
  img = getframe;
  paths{i}
end
