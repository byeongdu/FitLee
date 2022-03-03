% Usage:
%   edensity('chemicalformula', density)
%       density in g/cm3
%       example: edensity('Al2O3', 3.95)
%   edensity(n, N, A, density)
%       n : number of atoms of each indivisulal element
%       N : number of electrons of a single atom
%       A : atomic weight 
%       density : density (g/cm3)
%       example : MIBK
%           chemical formula :C6 H12 O
%           density : 0.7965g/cm3 at 25oC
%           n = [6, 12, 1]
%           N = [6, 1, 8]
%           A = [12.011, 1.008, 16]
%
%           edensity = 2.6808E-1 e/Angstrom3
%


function eden = edensity(varargin)
switch numel(varargin)
    case 4
        n = varargin{1};
        N = varargin{2};
        A = varargin{3};
        density = varargin{4};
        
        z = (n*N')/(n*A');
        na = (n./sum(n.*A))*density*6.02E23*1E-24;

    case 2
        CFM = varargin{1};
        density = varargin{2};
        mol = molecule(CFM);
        z = sum(mol.Atom)/mol.Mw;
    case 1
        CFM = varargin{1};
        mol = atom(CFM);
        z = sum(mol.AtomicNumber)/mol.Mass;
        density = mol.density_solid{1}(1);
    otherwise
end
eden = z*density*6.02E23*1E-24;