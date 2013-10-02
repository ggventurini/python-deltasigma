function [sp, op] = dsexpand(s,c,k,n,o)
%function [sp, op] = dsexpand(s,c,k,n,o). Expand points s outward from c
%by a factor 1+k and update the offsets associated with normals n.
shift = c(:,ones(1,size(s,2)));
sp = shift + (1+k)*( s-shift );
if(nargin>3 & nargout==2)
    op = o + k*(o+c'*n);
end
