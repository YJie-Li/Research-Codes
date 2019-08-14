field1='set';field2='m';
value1=[];value2=[];
E=struct(field1,value1,field2,value2);
%%% Input  variables. 
%%% "evidence_num". The number of evidence. Data type: double.
%%% "set". Data type: double. format:[1,0,0;0,1,1]. Num is the number of elements in framework.
%%% "m". The BPAs of each evidence.
singleelement_num=3;
Keys={'1','2','3'};%The number of singletons. Here we use three single elements.
d_threshold=0.5;%Determine the threshold value according to the actual situation.


conflict=zeros(length(E));%The initialization of the measure of conflicts.
for i=1:length(E)
    for j=1:length(E)
        conflict(i,j)=confl_alter(E(i).m,E(j).m,E(i),E(j),singleelement_num,Keys);
    end
end
similarity=1-conflict;
support(1,[1:length(E)])=0;
for i=1:length(E)
    for j=1:length(E)
        if i~=j
            support(i)=support(i)+similarity(i,j);
        end
    end
end
crd=normalization(E,support);

 
 mnew=zeros(length(E),length(E(1).set));
 for i=1:length(mnew)
     mnew(i,:)=E(i).m;
 end
 [row,column]=size(mnew);
 for i=1:row
     for j=1:(column-1)
         mnew(i,j)=mnew(i,j)*crd(i);
     end
 end
 for i=1:row
     mnew(i,column)=crd(i)*mnew(i,column)+(1-crd(i));
 end
 newstruct=E;
 for i=1:length(newstruct)
     newstruct(i).m=mnew(i,:);
 end
 
for i=1:length(E)
     for j=1:length(E)
Dbpa(i,j) = EviDistance_TwoEvi(newstruct(i),newstruct(j),singleelement_num);
end
end
%Compare Dbpa
confl2vec=squareform(Dbpa);
preclass=linkage(confl2vec,'single');
 category=cluster(preclass,'cutoff',d_threshold,'criterion','distance');%category
 numofclass=numel(unique(category));
 
for i=1:length(newstruct)
    newstruct(i).class=category(i);
end

count=[newstruct.class];
field1='set';field2='m';field3='class';
value1=[];value2=[];value3=[];
classofevi=struct(field1,value1,field2,value2,field3,value3);%The initialization of clusters
sumcrd=zeros(1,numofclass);
for i=1:numofclass
    if (sum(count==i)>1)
        classofevi(i).m=fuscond(newstruct(find(category==i)));
        classofevi(i).set=newstruct(i).set;
        classofevi(i).class=i;
        sumcrd(1,i)=sum(crd(find(category==i)));
    end
    if (sum(count==i)==1)
        classofevi(i)=newstruct(find(category==i));
        sumcrd(1,i)=crd(i);
    end
end
 
if (numofclass>1)
    massmul = fusconfkconf(classofevi,singleelement_num,Keys);
else
    massmul = fuscond(classofevi);
end
