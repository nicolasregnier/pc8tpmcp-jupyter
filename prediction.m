function yp=prediction(y,u,A,B,r,onestep)
% PREDICTION des sorties d'un modèle ARX 
%
%      -1            -1   -r
%   A(q  ).y(k) = B(q  ).q  .u(k)
%             -1            -1      -2           -na
%          A(q  ) = [1 -a1.q   -a2.q   ... -ana.q   ]
%             -1            -1      -2           -nb
%          B(q  ) = [0  b1.q    b2.q   ...  bnb.q   ]
%
%   yp=prediction(y,u,A,B,r)
% 
%   y : vecteur des sorties (prétraitées)
%   u : vecteur des entrées (prétraitées)
%   A,B : vecteurs des coefficients des polynômes définissant le modèle (estimations)
%          A = [1 -a1 -a2 ...] 
%          B = [0  b1  b2 ...]
%   r : retard pur de l'entrée
%   yp : vecteur des sortie prédites par le modèle
%
%   Exemple : >> yp = prediction(y, u, A, B, 1) calcule les sorties prédite par un modèle ARX définit par A, B, et r=1.

Le résidu d’estimation $\varepsilon = {y - y}_{P}$ peut être calculé à partir de ces données, de la façon suivante : **eps = y-yp**.


% Vérification des arguments

if nargin<6 || isempty(onestep)
    onestep = true; % recalage
end
if nargin<5
    error('PREDICTION --> nombre incorrect d''arguments d''entree (tapez ''help prediction'')');
end

% Calcul de la sortie prédite

A = A(:);
B = B(:);
theta = [-A(2:end); B(2:end)];
na = length(A) - 1;
nb = length(B) - 1;
nt = na + nb;
k0 = max([na nb+r]) + 1;
N = length(u);
yp = y(1:k0-1);

if onestep==true % calcul avec recalage
    for i=k0:N
        x = [y(i-1:-1:i-na); u(i-1-r:-1:i-nb-r)];
        yh = x'*theta;
        yp = [yp; yh];
    end
end

if onestep==false % calcul sans recalage
    x = [y(k0-1:-1:k0-na); u(k0-1-r:-1:k0-nb-r)];
    for i=k0:N
        yh = x'*theta;
        yp = [yp;yh];
        x = [yp(i:-1:i-na+1); u(i-r:-1:i-nb+1-r)];
    end
end
