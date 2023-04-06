function [bits, table]=huffman_code(x, freqs, tabella)
% HUFFMAN_CODE           Codifica di Huffman
%
%      Synopsys:
%
%            [BITS, TABLE]=HUFFMAN_CODE(X, FREQS)
%
%            TABLE = HUFFMAN_CODE([], FREQS)
%
%            BITS  = HUFFMAN_CODE(X,  [], TABLE)
%
%      Parameters:
%
%           BITS  = Stringa di "0" e "1" che rappresenta il bitstream 
%                   in uscita dalla codifica di Huffman
%
%           TABLE = Tabella che descrive il codice di Huffman
%
%           X     = Sequenza di interi non negativi da codificare
% 
%           FREQS = Vettore di probabilita` su cui costruire la tabella di
%                   Huffman 
%
%      Description:
%
%           Codifica il vettore di interi non negativi X usando un codice
%           di Huffman.  In uscita BITS e` una stringa di caratteri "0" e
%           "1" che rappresenta X codificato e TABLE e` una tabella che
%           descrive il codice.  BITS e TABLE possono essere passati
%           direttamente a huffman_decode per ricostruire X
%
%              X = huffman_decode(BITS, TABLE)
%
%           Il parametro FREQS contiene le probabilita` su cui costruire
%           il codice (FREQS(n) e` la probabilita` del simbolo n-1).  Se
%           FREQS e` omesso viene calcolato sulla base delle frequenze
%           campionare di X.
%
%           Nella seconda forma
%
%                    TABLE = HUFFMAN_CODE([], FREQS)
%
%           crea e restituisce in TABLE la tabella di Huffman basata sulle
%           frequenze in FREQS.  Tale tabella puo` essere usata chiamando
%           HUFFMAN_CODE nella terza forma
%
%                    BITS  = HUFFMAN_CODE(X,  [], TABLE)
%
%
%      Defaults:
%
%           FREQS = frequenze campionare di X
%
%      See also: 
%
%           HUFFMAN_DECODE
%

% 
% 25/2/2004 Versione 1.1: 
%
%       * Reso compatibile con octave
%
%       * Aggiunta la possibilita` di separare la creazione della tabella
%         dalla codifica
%
% ---
% Autore: R. Bernardini
%


%%
%% Default handling
%%

%
% Call parsing
%


%
% Default values
%


%%
%% True code
%%

empty_list_elements_ok=1;  % Per Octave, no-op per Matlab

if (nargin < 2)
  %
  % FREQS e` stato omesso.  Calcolalo a partire da X.
  %
  if (isempty(x))
    error('X non puo` essere vuoto se manca FREQS');
  end
  
  freqs=hist(x, 0:(max(x)));
  freqs=freqs/length(x);
end

if (isempty(x))
  max_simb = length(freqs)-1;
else
  max_simb = max(x);
end

%
% Se FREQ e` NaN, vuol dire che l'utente ti ha passato una tabella
% precostruita.  In caso contrario, ricavati la tabella sulla base di
% freq. 
%

if (isempty(freqs))
  if (nargin < 3)
    error('FREQ e` vuota, ma manca il terzo parametro')
  end
  
  table = tabella;  % Usa la tabella passata dall'utente
else
  %
  % Oooh... bene... per prima cosa devo costruire il codice di Huffman..  
  % 
  % Conviene spendere due parole sull'algoritmo che intendo usare..
  % Il codice di Huffman puo` essere rappresentato da un albero in cui
  %
  %   1. le foglie sono i simboli dell'alfabeto
  %
  %   2. i nodi intermedi sono "simboli di comodo" introdotti durante la
  %      costruzione del codice 
  %
  %   3. ad ogni ramo e` associato il bit 0 o il bit 1
  %
  %   4. il codice associato ad ogni simbolo e` la stringa di bit associata
  %      al percorso che va dalla radice alla foglia.
  %
  % Per gestire tale albero in Matlab, costruiro` una tabella rappresentata
  % da una matrice con N righe (N=numero totale di nodi) e 4 colonne.  Ogni
  % riga rappresenta un nodo, mentre il significato delle colonne e` il
  % seguente: 
  %
  %   1. La prima colonna contiene la probabilita` associata al nodo (la
  %      probabilita` del  simbolo se il nodo e` una foglia, la somma delle
  %      probabilita` dei suoi figli, se il nodo e` un nodo intermedio)
  %
  %   2. La seconda colonna il simbolo (rappresentato da un intero non
  %      negativo) associato al nodo.  
  %
  %   3. La terza colonna contiene il bit associato al ramo che lega il nodo
  %      corrente al padre
  %
  %   4. La quarta colonna contiene il simbolo associato al padre.  Per il
  %      nodo radice (che non ha antenati) la quarta colonna contiene -1
  %
  % L'algoritmo, a grandi linee e` il seguente:
  %
  %   1. Inizializzo la tabella ed ordino le righe in ordine di probabilita`
  %      crescente 
  %
  %   2. Prendo le prime due righe della tabella (che corrispondono ai due
  %      simboli meno probabili), creo un nodo di comodo che uso come padre
  %      dei due nodi appena estratti ed assegno '0' ad un ramo ed '1'
  %      all'altro 
  %
  %   3. Elimino dalla tabella i due nodi appena usati e riprendo dal passo
  %      numero 2.
  %
  % Ad onor del vero, nel codice Matlab che segue non elimino i due nodi
  % dalla tabella, ma aggiungo 1 alla loro probabilita` e riordino la
  % tabella in ordine di probabilita` crescente.  In questo modo i due nodi
  % appena usati finiscono in fondo e non rompono piu`.  Mi accorgero` di
  % aver finito quando il primo nodo avra` probabilita` pari ad 1 (e` la
  % radice e ha come discendenti tutte le foglie) ed il secondo probabilita`
  % maggiore di uno. 
  %
  % Ok, si comincia!
  %
  
  prob  = 1;   % Alcuni nomi simbolici per le colonne della tabella,
  simb  = 2;   % in questo modo evito di confondere una colonna con
  bit   = 3;   % l'altra
  padre = 4;
  depth = 5;

  %
  % Crea un vettore contenente i simboli dell'alfabeto.  Rappresentiamo ogni
  % simbolo con un intero non negativo.
  %
  alfabeto=0:(length(freqs)-1);

  %
  % Inizializza la tabella con i cinque campi di cui sopra
  %
  tabella = [freqs(:), alfabeto(:), zeros(length(alfabeto), 2)-1];
  tabella = [tabella zeros(length(alfabeto), 1)];

  %
  % I simboli a frequenza nulla rompono le scatole quando si arriva a
  % costruire l'albero di Huffman.  Trova quali simboli hanno frequenza non
  % nulla e butta via gli altri
  %
  buoni = find(freqs ~= 0);
  tabella = tabella(buoni, :);

  %
  % Ok, ora ordina la tabella in ordine di probabilita` crescente
  %
  tabella = sortrows(tabella, prob);

  %
  % Siamo pronti.  Possiamo partire..
  %

  next_aux = max_simb+1;  % Intero associato al prossimo simbolo di comodo
  max_depth = 0;
  %
  % Vedi il commento piu` sopra per il test di terminazione
  %
  while (tabella(2,prob) < 1)
    %
    % Crea la riga relativa al padre dei due simboli meno probabili
    %
    %             Probabilita`              simbolo   bit  padre
    %                  |                       |       |     |
    %                  |                       |       |     |
    %                  V                       V       V     V
    new_entry = [sum(tabella(1:2,prob)),   next_aux,   -1,   -1];
    new_entry = [new_entry,  max(tabella(1:2, depth))+1];

    max_depth = max(max_depth, new_entry(depth));
    
    %
    % Aggiungi 1 alle probabilita` dei due nodi che stai elaborando, in modo
    % da spedirli in fondo alla tabella al prossimo giro.
    %
    tabella(1:2, prob)  = tabella(1:2,prob)+1;
    
    %
    % Assegna il bit 0 ad un ramo ed il bit 1 all'altro
    %
    tabella(1,   bit)   = 0;
    tabella(2,   bit)   = 1;
    
    %
    % Dichiara il nuovo nodo padre dei due nodi correnti
    %
    tabella(1:2, padre) = next_aux;
    next_aux=next_aux+1;
    
    %
    % Aggiungi il nuovo nodo alla tabella
    %
    tabella = [new_entry ; tabella];
    
    %
    % Riordina la tabella.  Ora i due nodi appena usati finiscono in fondo. 
    %
    tabella = sortrows(tabella, prob);
  end



  %
  % Molto, molto bene... Ora che abbiamo l'albero possiamo costruire il
  % codice vero e proprio.  Le parole di codice verranno tenute come
  % stringhe in un "cell array" con la convenzione che la parola associata
  % al simbolo n >= 0 sara` nella cella numero n+1.
  %
  % Per costruire il codice sara` sufficiente attraversare l'albero appena
  % costruito.  L'algoritmo, a spanne, e` il seguente
  %
  %   1. Crea un vettore di "padri" ancora da visitare ed inizializzalo
  %      mettendoci la radice dell'albero
  %   
  %   2. Fintanto che il vettore di padri non e` vuoto
  %
  %        2.1 Estrai il primo nodo della lista dei padri
  %
  %        2.2 Cerca i due nodi che sono figli del nodo appena estratto
  %
  %        2.3 Se il nodo appena estratto non ha figli, allora e` una
  %            foglia --> non fare niente 
  %
  %        2.4 Se il nodo ha figli (necessariamente due)
  %
  %             2.4.1 Assegna ad ogni figlio la parola associata al padre
  %                   seguita dal bit associato al figlio
  % 
  %             2.4.2 Metti i figli nella lista dei padri
  %
  % Ok, al lavoro!
  %

  %
  % Per prima cosa allochiamo la tabella delle parole di codice
  %
  n_simb = max(tabella(:,simb)+1);
  stringhe = reshape(blanks(n_simb*max_depth), n_simb, max_depth);
  vuota = blanks(max_depth);

  %
  % Trova il nodo radice.  E` il simbolo con padre==-1
  %
  idx = find(tabella(:,padre)==-1);  % Riga associata al nodo radice
  radice = tabella(idx, simb);       % Simbolo associato al nodo radice
  stringhe(radice,:) = vuota;        % Parola di codice (vuota) della radice

  %
  % Inizializza la lista dei padri
  %
  padri = [radice];       

  %
  % Itera fintanto che ci sono padri
  %
  while (~isempty(padri))
    %
    % Estrai il primo padre e cancellalo dalla lista
    %
    babbo = padri(1);  % Non sarebbe piu` carino in Ruby? ;-)
    padri(1) = [];     %    babbo = padri.shift
    
    %
    % Trova le righe dei nodi che hanno il nodo corrente come padre
    %
    idx = find(tabella(:,padre) == babbo);
    if (~isempty(idx))
      %
      % Il nodo corrente ha dei figli
      %
      figli = tabella(idx, simb);  % Simboli associati ai nodi figli

      %
      % Crea le stringhe associate ai figli concatenando quella del padre
      % con 0 o 1
      %
      for n=1:length(figli)
	s=stringhe(babbo+1, :);
	primo_libero=min(find(s==' '));
	s(primo_libero) = sprintf('%d', tabella(idx(n),bit));
	stringhe(figli(n)+1,:) = s;
      end
      
      %
      % Inserisci i figli nella lista dei padri da visitare
      %
      padri = [padri ; figli];
    end
  end


  %
  % Infine, elimina dalla tabella le parole di codice associate ai simboli
  % di comodo.   Nota che un simbolo e` di comodo se e solo se > max(X)
  %
  table=stringhe(1:(max_simb+1),:);
  
  %
  % Butta via le colonne identicamente vuote
  %
  da_buttare = find(all(table==' '));
  if (~isempty(da_buttare))
    table = table(:, 1:(min(da_buttare)-1));
  end
end

%
% Controlla che l'utente non sia interessato solo alla tabella
%
if (isempty(x))
  %
  % L'utente vuole solo la tabella: ritorna dopo aver copiato la tabella
  % nel primo parametro di uscita.
  %
  bits = table;
  return;
end

%
% Ora dobbiamo codificare X.  Questa e` la parte piu` facile: basta
% concatenare le parole di codice.
%

%
% Flag introdotto durante la fase di ottimizzazione per fare prove
% di velocita`.  Provate a codificare 30000 campioni (da un
% alfabeto di 256 simboli) ponendo prima fast=0 e poi fast=1 per
% apprezzare la differenza di velocita` tra i cicli for e le
% istruzioni vettorizzate in Matlab...
%
fast=1;

if (fast)
  bits = table(x+1,:);        % Estrai le parole di codice
  bits = bits';               %
  bits = bits(:)';            % Raddrizza la matrice in una riga
  
  %
  % Ok, infine buttiamo via gli spazi vuoti in bits.
  %
  da_cancellare = find(bits == ' ');
  
  if (exist('crash_dumps_octave_core'))
    bits(da_cancellare)='';    % Per octave
  else
    bits(da_cancellare)=[];    % Per matlab
  end
else
  % x=x(:)';   % Forza X a vettore riga (per usarlo nel for)
  % 
  % bits='';   % Inizializza la stringa di bit
  % for n=x    % Cicla sui simboli in X
  %   bits = [bits stringhe{n+1}];
  % end
end


%
% That's all, folks!
%


%%% May be obsolete
  %
  % Questa porzione di codice (nome pomposo per 1 (!) riga di
  % Matlab) ha un comportamento a dir poco bislacco per quanto
  % riguarda i tempi di esecuzione.  In particolare, il tempo
  % richiesto dipende dalla dimensione dell'alfabeto max(X) e,
  % soprattutto, l'esecuzione e` piu` lenta quando l'alfabeto e`
  % piu` piccolo?!? (a parita` di lunghezza di X)
  %
  % Si provi, per esempio, ad eseguire il codice che segue 
  %
  %       a=floor(rand(1,40000)*256);
  %       moduli=[4, 8, 16, 32, 50, 64, 100, 128, 200, 256];
  %       
  %       for m=moduli
  %         sprintf('Dim alfabeto=%d', m);
  %
  %         x=mod(a,m);   % I valori di x sono tra 0 e m-1
  %         
  %         tic
  %         b=huffman_code(x);
  %         toc
  %       end
  %
  % Il grafico loglog dei tempi ottenuti con questo codice mostra
  % che il tempo, per alfabeti piccoli, e` quasi inversamente
  % proporzionale alla dimensione dell'alfabeto?!?
  %
  % Non riesco ad immaginarmi nessuna spiegazione razionale...
  %
  % tic
