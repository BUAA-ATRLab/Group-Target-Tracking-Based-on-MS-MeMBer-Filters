load truth.mat;
for k=1:truth.K
   if ~isempty(truth.X{k}) 
       truth.X{k}=truth.X{k}([1 3 2 4],:);
   end
end