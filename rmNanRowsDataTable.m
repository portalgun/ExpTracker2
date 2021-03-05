function Snew=rmNanRowsDataTable(S,indFld)
% rms nans and places by making new cell in std
    nStd=size(S{1}.(indFld),2);
    nSubj=size(S{1}.(indFld),3);
    n=size(S{1}.(indFld),1);
    ind=false(n,nStd);
    for s=1:nStd
    for i=1:numel(S)
        k=any(isnan(S{i}.(indFld)(:,s,:)),3);
        ind(:,s)=ind(:,s) | k;
    end
    end

    Snew=cell(numel(s),1);
    for i = 1:numel(S)
        flds=fieldnames(S{i});
        Snew{i}=struct();
        for f=  1:length(flds)
            fld=flds{f};
            Snew{i}.(fld)=cell(nStd,1);
            for s = 1:nStd
                Snew{i}.(fld){s}=squeeze(S{i}.(fld)(~ind(:,s),s,:));
            end
        end
    end
end
