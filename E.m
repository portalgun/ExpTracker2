function out=E(name);
if ~exist('name','var')
    name=[];
end
out=Eobj.load(name);
if nargout < 1
    assignin('base','obj',out);
end
