function obj=runSubj(subj,Aliases,rule,mode,std,blk);
%function obj=runSubj(subj,optional Aliases,rule,mode,std,blk);
% Aliases,rule,mode,std,blk and optional
    if ~exist('Aliases','var')
        Aliases=[];
    end
    if ~exist('rule','var') || isempty(rule)
        rule='rand_all';
    end
    if ~exist('mode','var')
        mode=[];
    end
    if ~exist('std','var')
        std=[];
    end
    if ~exist('std','var')
        blk=[];
    end
    obj=EAobj.run_loop_all(subj,Aliases,rule,mode,std,blk);
end
