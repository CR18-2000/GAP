Print("Richiamare la funzione bch passando come parametro il vettore che voglio trasmettere \n Il vettore deve avere lunghezza 7 bit \n");

#Funzione per fare padding di 0
padToLength := function (vec, length)
    local res;
    #fa una copia modificabili (senza puntatore)
    res := ShallowCopy(vec);
    while Length(res) < length do
      Add(res, 0);
    od;
    return res;
end;


encode := function (vec, F, pol) 
    local res, limit;
    limit := Size(F)-1;
    if Length(vec) > limit then
        Print("Il vettore è troppo lungo: la lunghezza deve essere al massimo ", limit, "\n");
    return; fi;
    res := padToLength(vec, Size(F) - 1);
    res := UnivariatePolynomial(F, Reversed(res));
    #aggiungo il resto per la divisione tra il polinomio che sto trasmettendo per il prodotto tra il polinomio che definisce il campo e il polinomio minimo di alpha^3
    res := res + EuclideanRemainder(res,DefiningPolynomial(F)*pol);
    res := Reversed(CoefficientsOfUnivariatePolynomial(res));
    res := List(res, IntFFE);
    return padToLength(res, Size(F) - 1);
end;


errorLocator := function(vec, F, alpha, mat)
    local g, a, b, i, j, k, vettore, appo, unacol, bool, sol;
    #Polinomio i cui coefficienti sono ottenuti a partire dal vettore che sto passando al contrario
    g := UnivariatePolynomial(F, Reversed(vec)*One(F));
    a := Value(g, alpha*One(F));
    b := Value(g, alpha^3*One(F));
    vettore := [];
    for i in [1..8] do
        appo := 0;
        for j in [1..15] do
            appo := appo + mat[i][j]*vec[j];
        od;
        vettore[i] := Int(ZmodnZObj(appo,2));
    od;
    bool := true;
    #si può ottimizzare togliendo il ciclo mettendo in quello sopra, ma per un ciclo di 8 non mi cambia
    for i in [1..8] do 
        if vettore[i] <> 0 then bool := false; fi;
    od;
    if bool then return -1; fi;
    unacol := false;
    #se vale una tra a=b^3 e b=a^3 allora ho un solo errore
    if a = b^3 then unacol := true;
    else if b = a^3 then unacol := true;
        fi;
    fi;
    
    if unacol then 
        #se ho un errore il vettore ottenuto è uguale a una colonna
        Print("Errore di una colonna \n");
        sol := [];
        for i in [1..15] do
            bool := true;
            for j in [1..8] do
                if vettore[j] <> mat[j][i] then bool := false;
                fi;
            od;
            if bool then Print ("Errore al bit ", i, "\n");
            sol[1] := i; return sol; fi;
        od;
    else 
        #se ho due errori il vettore ottenuto sarà la somma delle due solonne corrispondenti ai bit errati
        sol := [];
        for i in [1..15] do
            for j in [1..15] do
                bool := true;
                for k in [1..8] do
                    if i <> j then 
                        if (mat[k][i] + mat[k][j]) mod 2 <> vettore[k] then bool := false; fi;
                    else bool := false;
                    fi;
                od;
                if bool then
                    Print("Errore in due colonne \n");
                    sol[1] := i;
                    Print("Errore ai bit ", i, " e ");
                    sol[2] := j;
                    Print(j, "\n");
                    return sol;
                fi;
            od;
        od;
    fi;
    return -2;
end;


bch := function(code)
    local n, r, alpha, bool, pol, mat, res, limit, vec, err, errore, errore2, F, theField, i;
    #n è il numero di bit del messaggio che vogliamo mandare
    n := Length(code); #Deve essere 7 in questo caso
    if n <> 7 then Print("Lunghezza non 7 \n"); return fail; fi; 
    r := 4; #8 bit di ridondanza, 2^4=7+8+1
    F := GF(2);
    #indeterminata del campo
    #costruisco il campo a partire dal campo con due elementi e da un polinomio di grado r in F[x]
    theField := GaloisField(F, ConwayPolynomial(2,r));
    #if (isWord(code, theField)) = false then Print("Il vettore non è nel dizionario \n"); return; fi;
    Print("Costruisco un campo a partire da F[x], con F={0,1} e dal polinomio di Conway di grado 4 \n");
    Print("Il polinomio su cui costruisco il campo è: ", ConwayPolynomial(2,r), "\n");
    Print("Codifico ora il vettore passato per parametro \n");
    alpha := RootsOfPolynomial(theField, ConwayPolynomial(2,r))[1];
    pol := MinimalPolynomial(F, alpha^3);
    res := encode(code, theField, pol);
    mat := [[1,1,1,1,0,1,0,1,1,0,0,1,0,0,0], [0,1,1,1,1,0,1,0,1,1,0,0,1,0,0], [0,0,1,1,1,1,0,1,0,1,1,0,0,1,0], [1,1,1,0,1,0,1,1,0,0,1,0,0,0,1], [0,1,0,1,0,0,1,0,1,0,0,1,0,1,0], [1,1,0,0,0,1,1,0,0,0,1,1,0,0,0], [0,1,1,0,0,0,1,1,0,0,0,1,1,0,0], [0,1,0,0,1,0,1,0,0,1,0,1,0,0,1]];
    Print(res, "\n");
    Print("Proviamo che non trovi errori in un messaggio non perturbato \n");
    if errorLocator(res, theField, alpha, mat) = -1 then Print("Codice trasmesso senza errori \n"); 
    else Print("Errore nella trasmissione \n"); return; fi;
    Print("Il vettore con i bit di controllo risulta: ", res, "\n");
    Print("Senza bit di ridondanza è: ");
    for i in [1..n] do 
        Print(res[i], " ");
    od;
    Print("\n");
    Print("Poniamo ora che ci siano stati errori nella trasmissione \n");
    errore := Random([1..2^r-1]);
    res[errore] := (res[errore] + 1) mod 2;
    errore2 := Random([0..2^r-1]);
    while errore2 = errore do
        errore2 := Random([0..2^r-1]);
    od;
    if errore2 <> 0 then res[errore2] := (res[errore2] + 1) mod 2; fi;
    Print("Errore ai bit ", errore, " ", errore2, "; la nuova stringa è: ", res, "\n");
    Print("Vediamo che passando il vettore alla funzione decodifica mi trova l'errore \n");
    err := errorLocator(res, theField, alpha, mat);
    if err = -1 then Print("Codice trasmesso in modo esatto \n"); return; 
    else if err = -2 then Print("Codice con più di due errori \n"); return;
    else Print("Errore nella trasmissione della stringa \n"); fi;
    fi; 
    if Length(err) = 1 then 
        res[err[1]] := (res[err[1]] + 1) mod 2; 
    fi;
    if Length(err) = 2 then
        res[err[1]] := (res[err[1]] + 1) mod 2; 
        res[err[2]] := (res[err[2]] + 1) mod 2; 
    fi;
    #vec := decode(res, theField);
    Print("La funzione decode mi ritorna ", res, "\n");
    #Print("Senza bit di ridondanza è: ");
    #for i in [1..n] do 
    #    Print(vec[i], " ");
    #od;
end;