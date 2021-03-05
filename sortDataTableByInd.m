function Snew=sortDataTableByInd(S,indFld)
% sorting pass data
nStd=size(S{1}.(indFld),2);
nSubj=size(S{1}.(indFld),3);
vals=cell(nStd,1);
for i = 1:numel(S)
for s = 1:nStd
    vals{s}=[vals{s}; S{i}.(indFld)(:,s,:)];
end
end
n=zeros(size(S{i}.(indFld),2),1);
for s = 1:size(S{i}.(indFld),2)
    vals{s}=sort(unique(vals{s}));
    n(s)=size(vals{s},1);
end
n=max(n);

Snew=cell(size(S));
for i =1:numel(S)

    %ind=num2cell(INDS(j,:));
    for s = 1:nStd
    for j = 1:nSubj
        A=S{i}.(indFld)(:,s,j);
        indsA{s,j}=cell2mat(arrayfun(@(y) find(A==y),vals{s},'UniformOutput',false));
        indsB{s,j}=cell2mat(arrayfun(@(x) find(vals{s}==x),A,'UniformOutput',false));
    end
    end

    Snew{i}=struct();
    flds=fieldnames(S{i});
    nFld=length(flds);
    for f = 1:nFld
        fld=flds{f};
        if size(S{i}.(fld),1) ~= size(S{i}.(indFld),1)
            continue
        end
        sz=size(S{i}.(fld));
        sz=sz(2:end);
        Snew{i}.(fld)=nan([n,sz]);
        for s = 1:nStd
        for j = 1:nSubj
            Snew{i}.(fld)(indsB{s,j},s,j) =S{i}.(fld)(indsA{s,j},s,j);
        end
        end
    end
end
%14880
%14865
