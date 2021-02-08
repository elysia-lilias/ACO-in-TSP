function Chaosp = chaosinit(n) %D:distance matrix

a = rand();
Chaos_element = zeros(1,100);
Chaos_element(1) = 3.8*a*(1-a);
for i=2:100
Chaos_element(i) = 3.8* Chaos_element(i-1) * (1-Chaos_element(i-1));
end
Chaos_element = Chaos_element;
Phem = zeros(n);

%DVC idea based on the essay "Solving Traveling Salesman Problem by Chaos Ant
% Colony Optimization Algorithm" by GAO Shang

%randomly choose 100 Full Permutation based on chaos series

D = zeros(1,100);    %index of the Permutations in lexicographical order
d1 = zeros(100,n-1);
V = zeros(100,n-1);  %Used in finding the Permutation
C = zeros(100,n);  %Permutation

%From D to V
D = ceil(Chaos_element.*factorial(n));
d1(:,1) = n*Chaos_element;
V(:,1) = ceil(d1(:,1));
for i = 2:n-1
   d1(:,i) = (n-i+1).*(d1(:,i-1)-V(:,i-1)+1);
   V(:,i) = ceil(d1(:,i));
end


%From V to C
order = zeros(1,n);
temporder = zeros(1,n);
for i = 1:n
    order(i) = i;
end
for num = 1:100
    temporder = order;
    for i=1:n-1
        C(num,i) = temporder(V(num,i));
        temporder(V(num,i)) = [];
    end
      C(num,n) = temporder(1);
end
Chaosp = C;

end