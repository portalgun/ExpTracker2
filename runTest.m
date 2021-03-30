function obj=runTest(alias,mode,std,blk,S);
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
    if ~exist('blk','var')
        blk=[];
    end
    if ~exist('S','Var')
        S=[];
    end
    if evalin('base','exist(''Stest'',''var'');') 
        S=evalin('base','Stest;');
    end
    obj=Eobj.load(alias);
    obj.run('DNW',mode,std,blk,0,1,S);
end
