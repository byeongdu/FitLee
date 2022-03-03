function [y, v2] = schultzpolyhedralFun(q, r0, sig, particleshape, truncationlevel)
% function [y, v2] = schultzpolyhedralFun(q, r0, sig, particleshape, truncationlevel)
% 
    switch particleshape
        case {'rhombicdodecahedron', 'rhombic dodecahedron'}
            anisofile = 'rhombicdodecahedron.mat';
        case 'octahedron'
            anisofile = 'octahedron.mat';
        case 'truncatedoctahedron'
            switch truncationlevel 
                case 0.1
                    anisofile = 'octhedron_tr10p.mat';
                case 0.2
                    anisofile = 'octhedron_tr20p.mat';
                case 0.3
                    anisofile = 'octhedron_tr30p.mat';
                case 0.4
                    anisofile = 'octhedron_tr40p.mat';
                case 0.5
                    anisofile = 'octhedron_tr40p.mat';
                case 0.6
                    anisofile = 'octhedron_tr40p.mat';
            end
        case 'truncatedcube'
            switch truncationlevel 
                case 0.1
                    anisofile = 'cube_tr10p.mat';
                case 0.2
                    anisofile = 'cube_tr20p.mat';
                case 0.3
                    anisofile = 'cube_tr30p.mat';
                case 0.4
                    anisofile = 'cube_tr40p.mat';
                case 0.5
                    anisofile = 'cubooctahedron.mat';
            end
        case 'cube'
            anisofile = 'cube.mat';
        case 'concavecube'
            anisofile = 'concavecube.mat';
        case 'convexcube'
            anisofile = 'convexcube.mat';
        case 'cylinder'
            anisofile = [];
            %anisofile = 'disk2.mat';
        case 'THH'
            anisofile = 'THH.mat';
        otherwise
            anisofile = [];
    end
if isempty(anisofile)
%    y = []; v2 = [];
    p{1}.shape = particleshape;
    p{1}.edgelength = r0;
    p{1}.edgelength_sig = sig;
    [~, y] = AnisoFormFactor(p, q);
    return
end
    a = load(anisofile);
    dt = a.formfactortoscale.data;
    sz = a.formfactortoscale.size;
    [y, v2] = formfactorscale(dt(:,1), dt(:,2), sz, q, r0, sig, 1);
    