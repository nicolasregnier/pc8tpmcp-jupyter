function x=sbpa(n, k, amplitude, N)
% Génération d'une séquence binaire pseudo aléatoire (SBPA)
%
%   u = sbpa(n, k, LEVELS)
%
%   n, k : paramètres de la séquence, respectivement nombre de bistables et periode. 
%   LEVEL : vecteur dont les deux éléments représentent respectivement la valeur inférieure 
%   et la valeur supérieure de l’amplitude du signal.
% 
%   Exemples : >> u = sbpa(7, 3, [50-20 50+20])
%              Génération d'une SBPA définie par n=7 et k=3, et dont l’amplitude varie de 30 à 70.
%              (plus ou moins 20 autour de 50, qui est alors la valeur moyenne de la séquence)


    if nargin<3 || isempty(amplitude)
        amplitude = [-1 1];
    end
    if nargin<2 || isempty(k)
        k = 1;
    end    
    if nargin<1
        error('PREDICTION --> nombre incorrect d''arguments d''entree (tapez ''help sbpa'')');
    end

    if nargin<4 || isempty(N)
        N = k*(2^n-1)+n*k;
    end
    
    freqmax = 1/k; 

    % Positions
    feedback = hex2dec([...
        '0000000000',
        '0000000003',
        '0000000005',
        '0000000009',
        '0000000012',
        '0000000021',
        '0000000041',
        '000000008e',
       %'00000000b1',
        '0000000108',
        '0000000204',
        '0000000402',
        '0000000829',
        '000000100d',
        '0000002015',
        '0000004001',
        '0000008016',
        '0000010004',
        '0000020013',
        '0000040013',
        '0000080004',
        '0000100002',
        '0000200001',
        '0000400010',
        '000080000d',
        '0001000004',
        '0002000023',
        '0004000013',
        '0008000004',
        '0010000002',
        '0020000029',
        '0040000004',
        '0080000057',
        '0100000029',
        '0200000073',
        '0400000002',
        '080000003b',
        '100000001f',
        '2000000031',
        '4000000008'
     ]);

    T = floor(1/freqmax);
    n = round(log2(N/T+1)); % par arrondi
    n = min([n, length(feedback)]);
    feedback = str2num(char(cellstr(dec2bin(feedback(n),n)')))';

    % Génération
    len = ceil(N/T);
    x = ones(1,len+n);
    xx = zeros(1,N);
    for i=1+n:1:len+n
        x(i) = mod(sum(and(feedback,x(i-n:i-1))),2);
        xx((i-n-1)*T+1:(i-n)*T) = ones(1,T)*x(i);
    end

    x = xx;

    % Mise à l'échelle
    x = amplitude(1) + (amplitude(2)-amplitude(1))*(x-min(x))/(max(x)-min(x));

    x = x(:);
