function Snew=combineStdsDataTable(S,indFld)
    nStd=numel(S{1}.(indFld));
    nSubj=size(S{1}.(indFld){1},2);
    nExp=numel(S);

    flds=cellStructGetSameFlds(S);
    Snew=cell(nStd,1);
    for s=1:nStd
        Snew{s}=struct();
        for f=1:length(flds)
            fld=flds{f};
            n=size(S{1}.(fld){s},1);
            Snew{s}.(fld)=zeros(n,nSubj,nExp);
        end
    end

    for s= 1:nStd
    for i=1:nExp
    for f = 1:length(flds)
        fld=flds{f};
        n=size(S{i}.(fld){s},1);

        size(S{i}.(fld){s})
        Snew{s}.(fld)(:,:,i)=S{i}.(fld){s};
    end
    end
    end
end
