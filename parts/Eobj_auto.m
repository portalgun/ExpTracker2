classdef Eobj_auto < handle
methods
    function stds=get_std_nums(obj)
        if isfield(obj.methodVars,'stdXunqAll')
            stds=obj.methodVars.stdXunqAll;
            if iscell(stds)
                stds=distribute(stds{:});
            end
        else
            stds=0;
        end
    end
    function stds=get_std_flds(obj)
        if isfield(obj.methodVars,'stdXunqAll')
            stds=obj.methodVars.stdXunqAll;
            if iscell(stds)
                stds=distribute(stds{:});
                s=cell(size(stds));
                for i = 1:size(stds,2)
                    s(:,i)=num2fldstr(stds(:,i));
                end
                stds=s;
                if size(stds,2) > 1
                    stds=strrep(join(stds,2),' ','_');
                end
            else
                stds=num2fldstr(stds);
            end
        else
            stds='p0';
        end
    end
    function std=get_std_fld(obj,num)
        stds=obj.get_std_flds;
        std=stds{num};
    end
    function out=get_all_data_combs(obj)
        mthds={'test','train','pilot'};
        stds=obj.get_std_flds();
        blks=1:obj.nBlk;
        out=distribute(obj.subjs,mthds,stds,blks);
    end
%end
%methods(Access=private)
    function [subj,mode,std,blk]=auto_fld(obj,subj,mode,std,blk)
        subj=obj.auto_subj(subj);
        if nargin > 2 && nargout > 1
            mode=obj.auto_mode(mode);
        end
        if nargin > 3 && nargout > 2
            std=obj.auto_std_fld(std);
        end
        if nargin > 4 && nargout > 3
            blk=obj.auto_blk_num(blk);
        end
    end
    function subj=auto_subj(obj,subj)
        if ischar(subj) && ismember(subj,obj.subjs)
            return
        elseif ischar(subj) && ~ismember(subj,obj.subjs)
            error(['Invalid subject name: ' subj]);
        elseif all(isint(subj)) && all(subj > 0) && all(subj <= length(obj.subjs))
            subj=obj.subjs{subj};
        elseif subj <= 0 || subj > length(obj.subjs)
            error(['Subject index out of bounds: ' num2str(subj)]);
        elseif ~isint(subj)
            error(['Subject index must be positive: ' num2str(subj)]);
        end
    end
    function subj=auto_subj_num(obj,subj)
        if all(isint(subj())) && all(ismember(subj,1:numel(obj.subjs)))
            return
        elseif all(isint(subj))
            error('Subj ind out of bounds')
        else
            subj=find(ismember(obj.subjs,subj));
        end

    end
    function mode=auto_mode_num(obj,mode)
        if isempty(mode)
            mode=1;
        elseif all(isint(mode)) && all(ismember(mode,1:3))
            return
        elseif all(isint(mode))
            error('Mode ind out of bounds')
        else
            mode=find(ismember({'test','trian','pilot'},mode));
        end

    end
    function mode=auto_mode(obj,mode)
        if ~exist('mode','var') || isempty(mode)
            mode='test';
        elseif ischar(mode) && ismember(mode,{'test','train','pilot'})
            return
        elseif ischar(mode) && ~ismember(mode,{'test','train','pilot'})
            error(['Invalid mode name: ' mode]);
        elseif (isint(mode) && (mode==1 || mode==0)) || isempty(mode)
            mode='test';
        elseif  isint(mode) && mode==2
            mode='train';
        elseif  isint(mode) && mode==3
            mode='pilot';
        elseif mode < 0 || mode > 3
            error(['Mode index out of bounds: ' num2str(mode)]);
        elseif ~isint(mode)
            error(['Mode index must be positive: ' num2str(mode)]);
        end
    end
    function std=auto_std_fld(obj,std)
        % XXX make more like auto_std_num
        stdsf=obj.get_std_flds;
        stds=obj.get_std_nums;
        if iscell(stds)
            stds=distribute(stds{:});
        end
        stdsi=1:length(stds);
        bStds=isfield(obj.methodVars,'stdXunqAll');
        if ~bStds && ((~exist('std','var')  || isempty(std)) || std==1)
            std='p0';
        elseif ischar(std) && ismember(std,stdsf)
            return
        elseif isnumeric(std) && size(stds,2)> 1 && length(std)==size(stds,1) && ismember(std,stds,'rows')
            std=num2fldstr(std);
        elseif isnumeric(std) && size(stds,2)==1 && numel(std)== 1 && ismember(std,stds)
            std=num2fldstr(std);
        elseif isint(std) && ismember(std,stdsi)
            stdX=obj.methodVars.stdXunqAll;
            if iscell(stdX)
                std=join(stdsf(std,:),'_');
                std=std{1};
            else
                std=num2fldstr(stdX(std));
            end
        elseif isint(std)
            error(['Std index/value out of bounds: ' num2str(std)]);
        elseif isnumeric(std)
            error(['Std number out of bounds: ' num2str(std)]);
        elseif ischar(std)
            error(['Invalid std field-string' std]);
        end
    end
    function std=auto_std_ind(obj,std)
        stds=obj.get_std_nums();
        flds=obj.get_std_flds();
        if iscell(std) && all(ismember(str2double(std),stds))
            std=str2double(std);
            std=find(ismember(stds,std));
            
        elseif all(isint(std)) && all(ismember(std,1:length(stds)))
           return
        elseif isnumeric(std) && ismember(std,stds)
            std=find(stds==std);
        elseif ischar(std) &&  ismember(std,flds)
            std=find(ismember(flds,std));
        elseif ischar(std) && ismember(str2double(std),stds)
            std=find(ismember(stds,str2double(std)));
        else
            error(['Invalid std value' std]);
        end
    end
    function std=auto_std_str(obj,std)
        std=strrep(num2strSane(obj.auto_std_num(std)),',',' ');

    end
    function std=auto_std_num(obj,std)
        stds=obj.get_std_nums();
        flds=obj.get_std_flds();

        if iscell(std) && all(ismember(str2double(std),stds))
            dk
        elseif all(isint(std)) && all(ismember(std,1:length(stds)))
            std=stds(std,:);
        elseif all(isnumeric(std)) && all(ismember(std,stds))
            return
        elseif all(isint(std)) && all(ismember(std,1:length(stds)))
            std=stds(std);
        elseif ischar(std) &&  ismember(std,flds)
            std=fldstr2num(std);
        else
            error(['Invalid std value' std]);
        end
    end
    function blk=auto_blk_num(obj,blk)
        if ischar(blk) && ismember(num2double(blk),1:obj.nBlk)
            blk=num2str(blk);
        elseif ischar(blk) && ~ismember(num2double(blk),1:obj.nBlk)
            error(['Invalid block string ' blk])
        elseif isint(blk) && ismember(blk,1:obj.nBlk)
            return
        elseif isint(blk) && ~ismember(blk,1:obj.nBlk)
            error(['Invalid block number ' num2str(blk)])
        end
    end
    function blk=auto_blk_str(obj,blk)
        if ischar(blk) && ismember(num2double(blk),1:obj.nBlk)
            blk=sprintf('%03i',num2double(blk));
        elseif ischar(blk) && ~ismember(num2double(blk),1:obj.nBlk)
            error(['Invalid block string ' blk]);
        elseif isint(blk) && ismember(blk,1:obj.nBlk)
            blk=sprintf('%03i',blk);
        elseif isint(blk) && ~ismember(blk,1:obj.nBlk)
            error(['Invalid block number ' num2str(blk)]);
        end
    end
    function r1=auto_int_str(obj,r1)
        if isint(r1)
            r1=num2str(r1);
        elseif ischar(r1) && isint(num2str(r1))
            return
        end
    end
    function str=auto_pass_str(obj)
        if isempty(obj.pass)
            pass=1;
        else
            pass=obj.pass;
        end
        str=['pass' num2str(pass)];
    end
    function str=auto_ind_str(obj)
        str=strrep(num2strSane(obj.prjInd),',','-');
    end
    function str=auto_indd_str(obj)
        str=strrep(num2strSane(obj.prjInd),',','d');
    end
    function fld=auto_out_fld(obj,fld)
        fld=makeLowerCase(fld);
        fld=strrep(fld,'prunedp','prunedP');
    end
end
methods(Static)
    function fld=stdXunq2string(val)
        if iscell(val)
            fld=strrep(num2strSane(val),',','_');
        else
            fld=num2strSane(val);
        end
    end
end
end
