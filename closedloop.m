function [y,u,t] = closedloop(w, R, S, T)
% SIMULATION d'un système en boucle fermée
%
%   [y,u,t] = closedloop(w, R, S, T)
%
%      w : vecteur contenant les valeurs échantillonnées de la consigne
%          Valeurs appliquées selon la période d'échantillonnage Te du procédé
%      R, S, T : vecteurs contenant les coefficients des polynômes qui définissent le régulateur
%      y : vecteur contenant les valeurs échantillonnées de la réponse du procédé (variable de sortie mesurée)
%          Dimension égale à celle du vecteur w
%      u : vecteur contenant les valeurs échantillonnées de la variable de commande appliquée à l'entrée du procédé
%          imension égale à celle du vecteur w
%      t : vecteur contenant les valeurs du temps en secondes 
%          Valeur initiale : 0, incrément : période d'échantillonnage Te du procédé
%          Dimension égale à celle du vecteur w

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
    
    % Construction de la boucle fermée
    Hd = c2d(proc.H,proc.Te);
    Hbf = T*feedback(tf(1,R,proc.Te)*Hd, tf(S,1,proc.Te));

    % Calcul de la sortie
    t = 0:proc.Te:(length(w)-1)*proc.Te;
    y = lsim(Hbf,w-proc.y0,t);

    % Prise en compte du retard
    if proc.taur>0
        delay = thiran(proc.taur, proc.Te);
        y = lsim(delay,y);
    end
        
    % Ajout de la valeur initiale et du bruit de mesure
    randn("seed", "reset");
    noise = randn(length(w),1)*proc.b;
    y = y + proc.y0 + noise;   
    
    % Calcul de l'entrée
    Hu = T*feedback(tf(1,R,proc.Te), Hd*tf(S,1,proc.Te));
    u = lsim(Hu,w-proc.y0,t);
    
    % Ajout du bruit de mesure
end