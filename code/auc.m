function overallauc=positiontooverallaucfold(position)
position1=position;
for i=1:100
    position=position1(i,:);
interaction = xlsread('Protein_Disease_Associations.xlsx');
[n,m]=size(interaction);
sID=textread('Protein_Disease_adj.txt');
[pp,qq]=size(sID);
end
for k=1:m*n-floor(pp/5)*4
    tp=0;
    for t=1:pp
        if position(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
fp=k*pp-tp;
  
     fpr(1,k)=fp/(floor(pp/5)*4*(m*n-pp+floor(pp/5)-1)+(pp-floor(pp/5)*4)*(m*n-floor(pp/5)*4-1));
end
plot(fpr,tpr)