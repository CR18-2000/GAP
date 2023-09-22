#Hamming
Print("Richiamare la funzione hamming(code) passando al posto di code \n il vettore binario da inviare (usa ad esempio [1,0,0,1]) \n");


isWord := function(vec, F)
    local g;
    #la funzione UnivariatePolynomial mi costruisce un polinomio in una variabile sul campo theField e con coefficienti quelli passati nel vettore
    g := UnivariatePolynomial(F, Reversed(vec)*One(F));
    #valutiamo il polinomio in una radice del polinomio che definisce il campo; se è uno zero del campo ritorna true, altrimenti false
    return Value(g, RootOfDefiningPolynomial(F)) = Zero(F);
end;


errorLocator := function(vec, F)
    local g, error, disLog;
    #se è una parola del dizionario allora non ho un errore
    if (isWord(vec,F)) then return -1; fi;
    g := UnivariatePolynomial(F, Reversed(vec));
    error:= Value(g, RootOfDefiningPolynomial(F));
    #LogFFE ritorna il logaritmo discreto dell'elemento error rispetto alla radice del polinomio che definisce il campo
    disLog := LogFFE(error, RootOfDefiningPolynomial(F));
    while disLog < 0 do disLog := disLog + Size(F) - 1; od;
    return disLog;
end;


padToLength := function (vec, length)
    local res;
    #fa una copia modificabili (senza puntatore)
    res := ShallowCopy(vec);
    while Length(res) < length do
      Add(res, 0);
    od;
    return res;
end;


encode := function (vec, F) 
    local res, limit;
    limit := Size(F) - Length(FactorsInt(Size(F)));
    if Length(vec) > limit then
        Print("Il vettore è troppo lungo: la lunghezza deve essere al massimo ", limit, "\n");
    return; fi;
    res := padToLength(vec, Size(F) - 1);
    res := UnivariatePolynomial(F, Reversed(res));
    #aggiunge il resto della divisione tra il polinomio che sto considerando per il polinomio che definisce il campo
    res := res + EuclideanRemainder(res,DefiningPolynomial(F));
    res := Reversed(CoefficientsOfUnivariatePolynomial(res));
    res := List(res, IntFFE);
    return padToLength(res, Size(F) - 1);
end;


decode := function (vec, F)
    local errLoc, res, limit;
    limit := Size(F) - 1;
    if Length(vec)>limit then Print("Il vettore è troppo lungo, la lunghezza deve essere al massimo ", limit, "\n"); return; fi;
    res := padToLength(vec, limit);
    errLoc := errorLocator(res, F);
    if errLoc = -1 then Print("Messaggio inviato in modo giusto \n"); return vec; fi;
    errLoc := Size(F) - 1 - errLoc;
    Print("Bit #", errLoc, "corretto \n");
    res[errLoc] := (res[errLoc] + 1) mod 2;
    return res;
end;


hamming := function(code)
    local n, r, grado, bool, res, limit, vec, errore, F, theField, i;
    #n è il numero di bit del messaggio che vogliamo mandare
    n := Length(code);
    bool := false;
    #r è il numero di bit di ridondanza
    r := 1;
    while not bool do
        if (n+1)<=(2^r-r) then bool := true; 
        else r := r+1; fi;
    od;
    #r := r-1;
    #magari metti di controllare che n+r sia una potenza di 2 -1
    grado := LogInt(r+n+1, 2);
    #Print("n=", n, " r=", r, " grado=", grado, "\n");
    if r+n+1 <> 2^grado then Print("La somma del numero di bit trasmessi e di ridondanza +1 non è una potenza di 2 \n"); return; fi;
    Print("Il numero di bit di ridondanza è: ", r, "\n");
    #campo con due elementi {0,1}
    F := GF(2);
    #indeterminata del campo
    #x := Indeterminate(F);
    #costruisco il campo a partire dal campo con due elementi e da un polinomio di grado r in F[x]
    theField := GaloisField(F, ConwayPolynomial(2,r));
    Print("Costruisco un campo a partire da F[x], con F={0,1} e dal polinomio di Conway di grado i, con i=LogInt(r+n+1, 2) \n");
    Print("Il polinomio su cui costruisco il campo è: ", ConwayPolynomial(2,grado), "\n");
    Print("Codifico ora il vettore passato per parametro \n");
    res := encode(code, theField);
    Print("Proviamo che non trovi errori in un messaggio non perturbato \n");
    if errorLocator(res, theField) = -1 then Print("Codice trasmesso senza errori \n"); 
    else Print("Errore nella trasmissione \n"); return; fi;
    Print("Il vettore con i bit di controllo risulta: ", res, "\n");
    Print("Senza bit di ridondanza è: ");
    for i in [1..n] do 
        Print(res[i], " ");
    od;
    Print("\n");
    Print("Poniamo ora che ci sia stato un errore nella trasmissione \n");
    errore := Random([1..n+r]);
    res[errore] := (res[errore] + 1) mod 2;
    Print("Errore al bit ", errore, "; la nuova stringa è: ", res, "\n");
    Print("Vediamo che passando il vettore alla funzione decodifica mi trova l'errore \n");
    if errorLocator(res, theField) = -1 then Print("Codice trasmesso in modo esatto \n"); return; 
    else Print("Errore nella trasmissione della stringa \n"); fi; 
    vec := decode(res, theField);
    Print("La funzione decode mi ritorna ", vec, "\n");
    Print("Senza bit di ridondanza è: ");
    for i in [1..n] do 
        Print(vec[i], " ");
    od;
end;