clear;
%ALEJO GARCIA 31263 

%%%%%%%%%%%%%%%%%%%%%
%%% PHASE VOCODER %%%
%%%%%%%%%%%%%%%%%%%%%

% PARA VARIACIONES DE PROCESAMIENTO:
% linea 107 max_peak
% wL
% p
% PARA MODIFICAR FRECUENCIA:
% >>> SHIFT 
% >>> HOPFACTOR = SHIFT
% PARA MODIFICAR TEMPO:
% >>> HOP OUT 
% >>> SHIFT = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%shift=0.25;
shift=0.5;
%shift=0.75;
%shift=1;                                  %ajuste freq shift
%shift=1.5;
%shift=2;
%shift=4;

%hopfactor=4;
%hopfactor=2;
%hopfactor=1.5;
%hopfactor=1;                              %ajuste factor de hop
%hopfactor=0.75;
hopfactor=0.5;
%hopfactor=0.25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
in=['VOICE_MALE.wav'];
out=['VOICE_MALE_F_x05.wav'];
  
[x,fS] = audioread(in);                     % wav mono
                                            %Fs freq de sampleo
                                              

[mx , nx ]=size(x);
if nx==2
        x=sum(x,2)/size(x,2);               %pasaje a mono si era estereo
end

xL=length(x);                             
                             
tIn=(0:1/fS:(xL-1)/fS);                       % vector tiempo

piArray=2*pi*(0:100);                       %usado para el unwrapping en el ajuste de fase




%wL=2^8; 
%wL=2^10; 
wL=2^11;                                   %longitud ventana en samples (Window Length)
%wL=2^12;
%wL=2^13;

%w=hanning(wL,'periodic');                   %tipo de ventana PROBAR VARIOS
w=hamming(wL,'periodic');
%w=blackman(wL, 'periodic');
%w=hann(wL, 'periodic');

p=2;                                        %parámetro de solapamiento de salto
%p=4;                                                    %si p<4 se empiezan a notar los saltos en tiempo
%p=8;                                                     % mientras mas grande es p mejor funciona bajando el tempo
p=16;

hopLIn=wL/p;                                %longitud salto carga (Hop Length In)

hopLOut=hopLIn*hopfactor;
%longitud salto sintesis 
%la relacion HopLIn HopLOut determina
%el estiramiento o achicamiento temporal de la señal

tOut=zeros(1,ceil(hopLOut/hopLIn)*xL);                         %vector tiempo de salida, varia con la relacion de hop.

stftF=(0:fS/wL:fS/2);                           %vector de saltos de frecuencia de la stft       

dTIn=hopLIn/fS;                             %salto temporal entre ventanas 

dTOut=hopLOut/fS;                           %salto temporal entre ventanas salida

Col=1+fix((xL-wL)/hopLIn);                      % las columnas es xL-wL por que el ultimo ventaneo no que queda por fuera de xL no tiene mas saltos.

stftL=ceil((1+wL)/2);                       % largo del array stft definido por el ancho de ventana  

phAux=zeros(wL/2+1);
phadvance=zeros(wL/2+1);                    % arrays auxiliares para comparar fases en el ajuste

phAdv=zeros(stftL+1);


for aux1=1:Col-1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ANALISIS STFT 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    wS=(((aux1-1)*hopLIn+1):((aux1-1)*hopLIn+wL)); %porcion de audio que se va a procesar
    
    xW=transpose(w.*x(wS));
   

    stft=fft(xW);
    
    mag=abs(stft(1:stftL+1));                       %magnitud
    
    ph=angle(stft(1:stftL+1));                      %fase
    
    %%%%%%%%%%%%%%%%%%%%%
    
    % AJUSTE DE FASE
    
    %%%%%%%%%%%%%%%%%%%%%
  
    
    max_peak=100;                              % maxima cantidad de picos para hacer el seguimiento de fase
                                               % probar como varia la
                                               % resolucion de salida
                                               % aumentando la cantidad de
                                               % picos
    eps_peak=0.005;                            % alto de picos
    
    p=FindPeaks4(mag, max_peak, eps_peak, stftF); %funcion del estilo findpeaks pero con parametros mas utiles para el uso que le doy
                                                  %recomendada por http://sethares.engr.wisc.edu/vocoders/Transforms.pdf
                                                  %FindPeaks4 devuelve una matriz p, donde la columna central es la posicion de el pico y
                                                  %las columnas 1 y 3 determinan un Q del pico encontrado.

    [nouse , peakIndex ] = sort(mag(p(:,2)));     %ordeno los picos por nivel, y me quedo con un vector indice de su posicion con respecto a los otros.                                        
                                                  
    pOrden=p(peakIndex,:);                        %el mismo vector p, pero ahora los picos estan en orden.
    
    pUse=pOrden(:,2);                             %me quedo solo con la frecuencia central.   
    
    
    for aux2=1:length(pUse)                                               % 
        
        dph=(ph(pUse(aux2))-phAux(pUse(aux2)))+piArray;                   %diferencial de fase phAux guarda la fase del dT anterior.
        
        dF=dph./(dTIn*2*pi);                                              %calculo de diferencial de frec a partir de diferencia de fase
        
        [nouse2 , indexF ] = min(abs(stftF(pUse(aux2))-dF));              %encuentra el minimo diferencial de frecuencia para cada salto
        
        fAj(aux2)=dF(indexF);                                             %mejor frecuencia estimada de cada fila
        
    end  
    
    magOut=mag;
    phOut=ph;
    
    for aux2=1:length(pUse)
        
        fOut=fAj(aux2);
        
        indexFQ=(pOrden(aux2,1):pOrden(aux2,3));        %Q's de frecuencia
         
        phAdv(pUse(aux2))=phAdv(pUse(aux2))+2*pi*fOut*dTOut;  %fase del pico + fase de ajuste, para cada pico
        
        piZero=pi*ones(1,length(indexFQ));
        
        pUseEnt=pUse(aux2)-pOrden(aux2,1)+1;               %Q de frec /2
        
        indpUse=(2-mod(pUseEnt,2)):2:length(indexFQ);       %vector del q / 2 en pasos de 2
        
        piZero(indpUse)=zeros(1,length(indpUse));          %vector de pi y cero segun Q de F sea par o inpar
        
        phOut(indexFQ)=phAdv(pUse(aux2))+piZero;           %%%%% SACANDO ESTE TERMINO EL VOCODER FUNCIONA
                                                                %SIN CORRECCION DE FASE
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    
    % SINTESIS
    
    %%%%%%%%%%%%%%%%%%%%%
    
    
        phAux=ph;
        
        
        c=magOut.*exp(sqrt(-1)*phOut);            %armado de array complejo de magnitud y fase para hacer la ISTFT
        
        c((wL/2)+1)=stft((wL/2)+1);               %el ultimo es igual al de entrada por que no tengo fase con que compararlo
        
        c=[c,fliplr(conj(c(2:(wL/2))))];
        
        
        wav=real(ifft(c));                        %genero el array de amp. vs t.
            
        outTL=round(((aux1-1)*hopLOut+1):(((aux1-1)*hopLOut)+1+wL));
        
        tOut(:,outTL)=tOut(:,outTL)+ wav  ;
end   
    
tOut=0.8*tOut/max(max(abs(tOut)));   %normalizo

tOOO=(0:1/fS:((length(tOut)-1)/fS));  %vector tiempo para graficar

wavwrite(tOut,fS*shift,16,out);
