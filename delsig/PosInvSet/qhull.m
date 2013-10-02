function [vertices, edges, normals, offsets] = qhull(points,str)
% function [V E N O] = qhull(P,str): Convex hull finder based on qhull
% P is an nxm list of m n-dimensional points.
% V is the vertices of the hull.
% E is the edges of the hull: pairs of indices into the V array.
% N is the normals for the facets and O is their offsets.
% Points inside the hull are characterized by x'*N + O < 0
%
% str is the argument string to qhull, the default is 
%  'qhull(mex) Qcx C0.001 A0.999'.
%
% See the documentation on qhull for more information.
%
% NOTE: In 2D, qhull() modifies those points in P that are vertices of 
%  merged facets.
