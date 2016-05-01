load movie.mat A;
axes('pos',[0 0 1 1],'visible','off'); % hide annoying  additional axis
movie(A,5,8);   % play movie A 5 times with 8 fps
