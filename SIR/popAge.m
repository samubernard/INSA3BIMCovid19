function pop = popAge(cbound)
%% Population divided in age classes defined in cbound

popfr = readtable('./popfr.txt');
popfr.Properties.VariableNames = {'birth', 'age', 'M', 'F', 'total'};

n = length(cbound);
M = -1;
pop = zeros(n+1,1);
for i = 1:n
    m = M;
    M = cbound(i);
    cat = popfr((popfr.age > m) & (popfr.age <= M) ,:);
    pop(i) = sum(cat.total);
end
cat = popfr(popfr.age > cbound(end) ,:);
pop(n+1) = sum(cat.total);

end