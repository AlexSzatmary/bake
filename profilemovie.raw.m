t_step=$smalldumpint$;
t_end=$nstep$;
ncap=$ncap$;
for k = 0:t_step:t_end
	hold off;
	for l = 1:ncap
		clear pro;
		pro=load(sprintf('cappro%04d_%06d.txt',l,k), '-ascii');
		plot(pro(:,1),pro(:,2),'.');
		axis([0 32 0 32])
		hold on;
	end
	M(k/t_step+1)=getframe;
end

