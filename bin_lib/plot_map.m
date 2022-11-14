clear all
intersections;
streets;

plot(intersections(:,1),intersections(:,2),'*')
xlim([-1000,1000])
ylim([-1000,1000])

len_inters = length(intersections(:,1));
for i = 1:len_inters
  text(intersections(i,1)-40, intersections(i,2)-40,int2str(i-1))
end
  
len_streets = length(streets(:,1));
for i = 1:len_streets
  A = intersections(streets(i,1)+1,:);
  B = intersections(streets(i,2)+1,:);
  line([A(1) B(1)], [A(2) B(2)]);
end

grid on