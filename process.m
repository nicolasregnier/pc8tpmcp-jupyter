function varargout = process(varargin)
% SIMULATION d'un procédé
%
%   [y,t] = process(u)
%   Simule la réponse y du procédé à un signal d'entrée u quelconque
%
%      u : vecteur contenant les valeurs échantillonnées de la variable de commande à appliquer à l'entrée du procédé
%          Valeurs appliquées selon la période d'échantillonnage Te du procédé
%      y : vecteur contenant les valeurs échantillonnées de la réponse du procédé (variable de sortie mesurée)
%          Dimension égale à celle du vecteur u
%      t : vecteur contenant les valeurs du temps en secondes 
%          Valeur initiale : 0, incrément : période d'échantillonnage Te du procédé
%          Dimension égale à celle du vecteur u
%
%  [Te,T98,T95] = process()
%  Fournit les caractéristiques temporelles du système
%
%      Te : période d'échantillonnage choisie pour l'échantillonnage des signaux, en secondes 
%      T98 : temps de réponse à 98% estimé d'après la réponse indicielle du procédé, en secondes
%      T95 : temps de réponse à 95% estimé d'après la réponse indicielle du procédé, en secondes

    str = getenv("JUPYTERHUB_USER");
    code = hash('md5', str); % 32 caractères hexadécimaux
    code = reshape(code,8,4); % 8x4 caractères hexadécimaux
    code = hex2dec(code)./65535.; % 8 nombre réels [0,1]
    proc.n = round(code(1)) + 1;
    proc.tau = 50*code(2)+30;
    proc.K = 0.2*code(3)+0.2;
    proc.b = 0.10*code(5)+0.04;
    proc.y0 = 7*code(6)+18;
    % Période d'échantillonnage, temps de réponse, et fonction de transfert
    if proc.n==1
        proc.Te = round(-proc.tau*log(0.05)/15);
        proc.T95 = round(-proc.tau*log(0.05));        
        proc.T98 = round(-proc.tau*log(0.02));
        proc.H = tf(proc.K,[proc.tau 1]);
    else
        w = lambertw(-1,-0.02*exp(-1));
        proc.Te = round(-proc.tau*(lambertw(-1,-0.05*exp(-1))+1)/15);
        proc.T95 = round(-proc.tau*(lambertw(-1,-0.05*exp(-1))+1));        
        proc.T98 = round(-proc.tau*(lambertw(-1,-0.02*exp(-1))+1));
        proc.H = tf(proc.K,[proc.tau^2 2*proc.tau 1]);
    end
    proc.taur = round(code(4)*proc.Te); 
    
    % varargin = signal d'entrée --> vargout = signal de sortie et temps
    if length(varargin)==1 && isfloat(varargin{1}) && nargout<=2
        u = varargin{1};
        randn("seed", "reset");
        % Calcul de la réponse
        t = 0:proc.Te:(length(u)-1)*proc.Te;
        y = lsim(proc.H,u,t);
        % Prise en compte du retard
        if proc.taur>0
            delay = thiran(proc.taur, proc.Te);
            y = lsim(delay,y);
        end
        % Ajout de la valeur initiale et du bruit de mesure
        noise = randn(length(u),1)*proc.b;
        y = y + proc.y0 + noise;
        % Arguments de sortie
        varargout{1} = y;
        if nargout>1
            varargout{2} = t;
        end
        return;
    end
    
    % varargin = vide --> vargout = période d'échantillonnage et temps de réponse à 98%
    if length(varargin)==0 && nargout<=3
        varargout{1} = proc.Te;
        if nargout>1
            varargout{2} = proc.T98;
        end
        if nargout>2
            varargout{3} = proc.T95;
        end        
        return;
    end
    
    error("PROCESS -> nombre d'arguments d'entrée ou de sortie incorrect. Consultez l'aide (help process)");
    
end