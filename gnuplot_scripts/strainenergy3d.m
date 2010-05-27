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
  load elm_0001.txt;
  load solidnodes02000.txt;
  stretches = load('stretches_02000.txt');
  w = (stretches(:,2).^2+stretches(:,3).^2+stretches(:,2).^-2.*stretches(:,3).^-2-3)/3;
  trisurf(elm_0001, solidnodes02000(:,1), solidnodes02000(:,2), solidnodes02000(:,3), w, 'EdgeColor', 'none');
  axis vis3d;
  subplot(isubplot(i));
end
