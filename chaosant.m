function chaosant()
clear all;
close all; 
clc;

%variables
number_ants = 50;
Alpha = 1; %power of Phem 
Beta = 0.5; %power of 1/Distance
Phemloss = 0.2; %lose percentage of Phem every trip
Stepmax = 100; %number of steps
Chaos_element_weight = 0.1;
Cities = [
    4214 5325;
    3639 1315;
    4177 2244;
    3712 1399;
    3488 1535;
    3326 1556;
    3238 1229;
    4196 1004;
    5252 3525;
    2251 3225;
    2623 1222;
    2552 2533;
    2525 5435;
    1244 5235;
    2364 7567;
    245 1425;
    636 346;
    14 366;
    1523 6346;
    1242 6363;
    1523 5254;
    5234 12;
    1525 3654;
    ];

n = size(Cities,1); %number of cities
D = zeros(n,n); %Complete graph
Memory = n; %how many place can each ant remember, no use here
            %just let it easier to write ant alogrithm for finding image 
            %   edge based on this program

%distance between cities
for i=1:n
    for j = 1:n
        if(i==j)
         D(i,j) = eps;
        else
            D(i,j) = ((Cities(i,1)-Cities(j,1))^2+(Cities(i,2)-Cities(j,2))^2)^0.5;
        end
         D(j,i)= D(i,j);
    end
end
%longer path,less phem
Eta = 1./D;
%Q = min(min(D))
Q = 100;
Chaosp = chaosinit(n); % Permutation based on chaos series
Phem = ones(n)*eps;
%chaos initialization of Phem
Chaossum = zeros(100,1);
for i=1:100
    Chaossum(i) = sum(Chaosp(i,:));
end
for i=1:100
    for j=1:n-1    
    Phem(Chaosp(i,j),Chaosp(i,j+1)) = Phem(Chaosp(i,j),Chaosp(i,j+1)) + Q/Chaossum(i);
    end
end
Phem;

currentstep = 1;
Path = zeros(number_ants,Memory);
%ants start traveling
AllBestDist = zeros(Stepmax,1);
AllBestPath = zeros(Stepmax,n);
AveDist = zeros(Stepmax,1);
wa=waitbar(0,'wait');
while currentstep <= Stepmax
    initpos = [];
     c = 0;
 while c<number_ants;
 initpos = [initpos,randperm(n)];
 c = c+n;
 end
 %at the beginning of each step, place ants randomly
 waitbar((n*(currentstep-1)+1)/(n*Stepmax));
 Path(:,1) = (initpos(1:number_ants))';
 for currentpos = 2:n
     waitbar((n*(currentstep-1)+currentpos)/(n*Stepmax))
     for current_ant = 1:number_ants
         visited = Path(current_ant,1:(currentpos-1)); %cities already visited
         nonvisited = zeros(1,(n-currentpos+1));      %cities non visited (init)
    pnv = nonvisited; %possibility of go to each nonvisited cities (init)
    temp_0 = 1;
    for temp_1 = 1:n;
        if(ismember(temp_1,visited)==0)
            nonvisited(temp_0) = temp_1;  %record each nonvisited cities
            temp_0 = temp_0 + 1;
        end
    end
    
    for temp_2 = 1:n-currentpos+1;
        pnv(temp_2) = (Phem(visited(end),nonvisited(temp_2))^Alpha)*(Eta(visited(end),nonvisited(temp_2))^Beta)/D(visited(end),nonvisited(temp_2));
    end
    %calculate possibility
    pnv = pnv/(sum(pnv));
    pnc = cumsum(pnv);
    temp_3 = find(pnc>=rand);
    next_city = nonvisited(temp_3(1));
    Path(current_ant,currentpos) = next_city;  %next city
     end
    Path
 end
  %find best path in this step
     Dist = zeros(number_ants,n-1); %distance
     for current_ant = 1:number_ants;
     for currentpos = 1:n-1
         Dist(current_ant,currentpos) = D(Path(current_ant,currentpos),Path(current_ant,currentpos+1));
     end
     end
     Totaldist = zeros(number_ants,1);
     for current_ant = 1:number_ants;
         Totaldist(current_ant) = sum(Dist(current_ant,:));
     end
     CurrentBest = min(Totaldist);
     AllBestDist(currentstep) = CurrentBest;
     CurrenBestPath = find(Totaldist== CurrentBest);
     AllBestPath(currentstep,:) = Path(CurrenBestPath(1),:);
     AveDist(currentstep) = mean(Totaldist);
%choose those path smaller than average
     ChosenPathNo = find(Totaldist <= AveDist(currentstep));
     DeltaPhem = zeros(n,n);
%how much phem each chosen ant path will be added
     for tempp = 1:length(ChosenPathNo)
     DeltaPhem(tempp) = Q/Totaldist(ChosenPathNo(tempp));
     end
%chaos variable
a = rand();
chaosa = 3.8*a*(1-a);
chaosmat = zeros(n,n);
%phem update
  tempdelt = zeros(n,n);
    for PUP = 1:length(ChosenPathNo) 
        ChosenAnt = ChosenPathNo(PUP);     
        tempPath = Path(ChosenAnt,:);
        for currentpos = 1:n-1
        tempdelt(tempPath(currentpos),tempPath(currentpos+1)) = tempdelt(tempPath(currentpos),tempPath(currentpos+1)) + DeltaPhem(PUP);
        if(chaosmat(tempPath(currentpos),tempPath(currentpos+1)==0))
            chaosmat(tempPath(currentpos),tempPath(currentpos+1)) = chaosa;
            chaosa = 3.8*chaosa*(1-chaosa);
        end
        end
    end
    Phem = (1-Phemloss)*Phem + tempdelt + chaosmat*Chaos_element_weight;
    currentstep = currentstep + 1;
    Path = zeros(number_ants,Memory);
end
close(wa);
Fin = find(AllBestDist==min(AllBestDist));
Shortest_Path = AllBestPath(Fin(1),:);
Shortest_Length = AllBestDist(Fin(1))
%plot the result
figure
%best distance as time goes by
subplot(1,2,1)
plot(AllBestDist)
%ave distance as time goes by
hold on
plot(AveDist,'r')
%best path
subplot(1,2,2)
scatter(Cities(:,1),Cities(:,2));
hold on
plot([Cities(Shortest_Path(1),1),Cities(Shortest_Path(n),1)],[Cities(Shortest_Path(1),2),Cities(Shortest_Path(n),2)],'r')
hold on
for iz = 2:n
   plot([Cities(Shortest_Path(iz-1),1),Cities(Shortest_Path(iz),1)],[Cities(Shortest_Path(iz-1),2),Cities(Shortest_Path(iz),2)],'r')
   hold on
end




end