function testDataAlignment(fname1,fname2)
f1=load_fun(fname1);
f2=load_fun(fname2);

mi=min([f1.cmpIind(:); f2.cmpIind(:)]);
ma=max([f1.cmpIind(:); f2.cmpIind(:)])-mi+1;

disp(newline)
disp('NO. MATCHING')

% R
t1=table_fun_chs(f1,ma);
t2=table_fun_chs(f2,ma);
[t12,t21]=Rcomp_fun(t1,t2);
disp(['RcmpChs'])
disp_fun(t12,t21,t1,t2);

% ALL
t1=table_fun_All(f1,ma);
t2=table_fun_All(f2,ma);
[t12,t21]=comp_fun(t1,t2);
disp(['CmpIind & stdIind & CmpX & stdX'])
disp_fun(t12,t21,t1,t2);

% INDS
t1=table_fun_ind(f1,ma);
t2=table_fun_ind(f2,ma);
[t12,t21]=comp_fun(t1,t2);
disp(['CmpIind & stdIind'])
disp_fun(t12,t21,t1,t2);

% X
t1=table_fun_X(f1,ma);
t2=table_fun_X(f2,ma);
[t12,t21]=comp_fun(t1,t2);
disp(['CmpX & stdX'])
disp_fun(t12,t21,t1,t2);
disp(newline)

function [t12,t21]=Rcomp_fun(t1,t2)
    t12=sum(~any(isnan(t1),2) & t1==t2);
    t21=sum(~any(isnan(t1),2) & t2==t1);
end
function [t12,t21]=comp_fun(t1,t2)
    t12=sum(~any(isnan(t1),2) & ismember(t1,t2,'rows'));
    t21=sum(~any(isnan(t2),2) & ismember(t2,t1,'rows'));
end

function disp_fun(t12,t21,t1,t2)
    n1=sum(~any(isnan(t1),2));
    n2=sum(~any(isnan(t1),2));
    disp(['   1-2: ' num2str(t12) ' / ' num2str(n1) ]);
    disp(['   2-1: ' num2str(t21) ' / ' num2str(n2) ]);
end

function table=table_fun_ind(f,ma)
    table=nan(ma,2);
    ind=f.cmpIind(:)-mi+1;
    table(ind,1)=f.cmpIind(:);
    table(ind,2)=f.stdIind(:);
end
function table=table_fun_X(f,ma)
    table=nan(ma,2);
    ind=f.cmpIind(:)-mi+1;
    table(ind,1)=f.cmpX(:);
    table(ind,2)=f.stdX(:);
end
function table=table_fun_All(f,ma)
    table=nan(ma,4);
    ind=f.cmpIind(:)-mi+1;

    table(ind,1)=f.cmpIind(:);
    table(ind,2)=f.stdIind(:);
    table(ind,3)=f.cmpX(:);
    table(ind,4)=f.stdX(:);
end
function table=table_fun_chs(f,ma)
    table=nan(ma,1);
    ind=f.cmpIind(:)-mi+1;

    table(ind,1)=f.RcmpChosen;
end

function out=load_fun(fname)
    exp=load(fname);
    out=struct();

    if isfield(exp,'S')
        exp=exp.S;
    end
    if isfield(exp,'raw')
        exp=exp.raw;
    end

    % R
    if isfield(exp,'RSP')
        out.R=exp.RSP.responses;
        if size(R,2) > 1
            out.R(:, ~any(out.R,1))=[];
        end
    elseif isfield(exp,'R')
        out.R=exp.R;
    end

    % RcmpChoesn
    if ~isfield(exp,'RcmpChosen')
        out.RcmpChosen=out.R==exp.cmpIntrvl;
    elseif isfield(exp,'cmpIntrvl')
        out.RcmpChosen=exp.RcmpChosen;
    end

    %cmpIind
    if isfield(exp,'cmpIind')
        out.cmpIind=exp.cmpIind;
    end

    %stdIind
    if isfield(exp,'stdIind')
        out.stdIind=exp.stdIind;
    end

    %stdX
    if isfield(exp,'stdX')
        out.stdX=exp.stdX;
    end

    %cmpX
    if isfield(exp,'cmpX')
        out.cmpX=exp.cmpX;
    end
end
end
