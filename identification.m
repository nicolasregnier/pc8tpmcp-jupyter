function varargout=identification(y,u,varargin)
% IDENTIFICATION d'un modèle ARX par la méthode des moindres carrés récurrents
% Renvoie les valeurs des coefficients du modèle et leurs intervalles de confiance,
% ou la variance des résidus
%
%   [A,icA,B,icB]=identification(y,u,na,nb,r)
%   [A,icA,B,icB]=identification(y,u,n,r)
%   [A,icA,B,icB]=identification(y,u,[na nb r])
%   [A,icA,B,icB]=identification(y,u,[n r])
%   Renvoie les coefficients estimés et leurs intervalles de confiance
%   
%   y : vecteur des sorties (prétraitées)
%   u : vecteur des entrées (prétraitées)
%   na, nb : ordres des polynômes A et B. 
%   Si une seule valeur n est fournie, alors na=nb=n
%   r : retard pur de l'entrée
%   A,B : vecteurs des coefficients des polynômes définissant le modèle (estimations)
%          A = [1 -a1 -a2 ...] 
%          B = [0  b1  b2 ...]
%   icA, icB : intervalles de confiance à 97.5% associés aux estimations
%          icA = [0 ic(a1) ic(a2) ...] 
%          icB = [0 ic(b1) ic(b2) ...] 
%          (ic(p) : intervalle de confiance à 97.5% associé au paramètre p)
%
%   vareps=identification(y,u,...)
%   vareps=identification(y,u,n1:n2,r1:r2)
%   vareps=identification(y,u,[na1 nb1 r1 ; na2 nb2 r2 ; ...])
%   Renvoie la variance des résidus
%
%   vareps : variance des résidus d'estimation, écarts y-yp entre les valeurs
%            de sortie y mesurées et yp prédites par le modèle
%   Dans ce cas, les structures de modèles à tester peuvent être spécifiées de plusieurs façons :
%          - soit à l'aide d'un vecteur d'ordres (na=nb) n1:n2 et d'un vecteur de retards r1:r2
%            vareps est alors une matrice dans laquelle chaque ligne correspond à un ordre et 
%            chaque colonne à un retard. Pour tracer la courbe, tapez plot(vareps).
%          - soit à l'aide d'une matrice dont chaque ligne est un triplet [na nb r]
%            vareps est alors un vecteur comportant autant de lignes. Pour tracer les 3 courbes
%            correspondant chacune à une valeur du retard en fonction de l'ordre, tapez plot(0:5, vareps).
%
%   Note : un ordre nul correspond à l'absence de modèle, et la variance retournée est alors celle de y.
%
%   Exemples : >> [A,icA,B,icB]=identification(y,u,[2 2 1])
%              Estimation d'un modèle ARX avec na=nb=2 et r=1.
%              >> vareps=identification(y,u,[1 1 1; 2 2 1; 3 3 1])
%              Estimation de 3 modèles ARX de structures différentes
%              avec na=nb=1 et r=1, na=nb=2 et r=1, na=nb=3 et r=1.
%              >> vareps=identification(y,u,1:5,0:2)
%              Estimation de 15 modèles ARX de structures différentes
%              avec na=nb=1, 2, 3, 4, ou 5 et r=0, 1, ou 2.

    % Vérification des arguments
    
    if length(varargin)==1
        if size(varargin{1},2)==2
            structure = [varargin{1}(1) varargin{1}(1) varargin{1}(2)];
        elseif size(varargin{1},2)==3
            structure = varargin{1};
        else
            error('IDENTIFICATION --> taille incorrecte d''arguments d''entrée (tapez ''help identification'')');
        end
        delay = [];
    elseif length(varargin)==2
        structure = repmat(varargin{1}(:),1,2);
        delay = (varargin{2}(:))';
    elseif length(varargin)==3
        structure = [varargin{1} varargin{2} varargin{3}];
        delay = [];
    end
    
    Norders = size(structure,1);
    Ndelays = size(delay,2);
    
    if Norders>1 && nargout>1
        error('IDENTIFICATION --> nombre incorrect d''arguments de sortie (tapez ''help identification'')');
    end
    
    nc = 0;

    uu = u;
    yy = y;
    N = length(u);

    % Normalisation des données

    ectu = std(u);
    ecty = std(y);

    vareps = zeros(Norders,max(1,Ndelays));
    
    for idelay=1:max(1,Ndelays)
        for iorder=1:Norders

            na = structure(iorder,1);
            nb = structure(iorder,2);
            if Ndelays>0
                r = delay(idelay);
            else
                r = structure(iorder,3);
            end

            % Modèle fictif d'ordre zéro
            if na*nb==0
                vareps(iorder,idelay) = var(yy);
                theta = [];
                precisions = [];
                continue;
            end
            
            %%%%%%%%%%%%%%%%%%%%%
            % Calcul paramètres
            %%%%%%%%%%%%%%%%%%%%%

            % Initialisations

            u = uu/ectu;
            y = yy/ecty;

            nt = na + nb + nc;
            theta = zeros(nt,1);
            p = eye(nt)*1e+06;
            k0 = max([na nb+r nc]) + 1;
            net = zeros(N,1);

            % Calcul récurrent de theta
            for i=k0:N
                x = [y(i-1:-1:i-na); u(i-1-r:-1:i-nb-r); net(i-1:-1:i-nc)];
                K = p*x/(1+x'*p*x);
                p = p - K*x'*p;
                theta = theta + K*(y(i)-x'*theta);
                net(i) = y(i) - x'*theta;
            end

            % Dénormalisation des paramètres
            theta=[theta(1:na); theta(na+1:na+nb)*ecty/ectu]';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculs variances & précisions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            u = uu;
            y = yy; 

            % Calcul de la sortie prédite
            yp = zeros(N-k0+1,1);
            for i=k0:N
                x = [y(i-1:-1:i-na); u(i-1-r:-1:i-nb-r)];
                yp(i-k0+1) = x'*theta';
            end

            % Calcul du résidu d'équation, de sa variance et des intervalles de confiance
            eps = y(k0:N) - yp;

            d = sqrt(diag(p));
            d = [d(1:na)/ecty; d(na+1:na+nb)/ectu];
            precisions = 2*std(eps)*d'; 

            vareps(iorder,idelay) = var(eps);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mise en forme et retour des résultats
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargout<=1
        varargout{1} = vareps; % vareps
    end

    if nargout>=4
        varargout{1} = [1 -theta(1:na)]; % A
        varargout{2} = [0 precisions(1:na)]; % icA
        varargout{3} = [0 theta(na+1:na+nb)]; % B
        varargout{4} = [0 precisions(na+1:na+nb)]; % icB
    end
    
    if nargout==5
        varargout{5} = vareps; % vareps
    end

end