#Diffie-Helman
Print("Richiamare la funzione algoritmo() passando come parametro un naturale \n");

alice1 := function(p, g)
    #alice sceglie un int e manda a bob g alla int mod p
    local ret, list, a;
    Print("Nel primo passaggio Alice sceglie un a e calcola g^a(mod p)", "\n");
    Print ("Primo passaggio di Alice: ");
    a := Random([100..100000]);
    Print("a=",a);
    ret := Int(ZmodnZObj(g^a, p));
    Print(", g^a(mod p)=", ret, "\n");
    list := [a,ret];
    return list;
end;

bob1 := function(p, g)
    #bob sceglie un int e manda a bob g alla int mod p
    local ret, list, b;
    Print("Nel primo passaggio Bob sceglie un b e calcola g^b(mod p)", "\n");
    Print ("Primo passaggio di Bob: ");
    b := Random([100..100000]);
    Print("b=", b);
    ret := Int(ZmodnZObj(g^b, p));
    Print(", g^b(mod p)=", ret, "\n");
    list := [b,ret];
    return list;
end;

alice2 := function (p, a, B)
    local s;
    Print("Nel secondo passaggio Alice calcola B^a (mod p), con B=g^b(mod p)", "\n");
    s := Int(ZmodnZObj(B^a, p));
    return s;
end;

bob2 := function (p, b, A)
    local s;
    Print("Nel secondo passaggio Bob calcola A^b (mod p), con A=g^a(mod p)", "\n");
    s := Int(ZmodnZObj(A^b,p));
    return s;
end;

algoritmo := function(p)
    local g, list, appo, a2, b2;
    if (IsPrime(p)=false) then p := NextPrimeInt(p); fi; 
    Print("Il numero primo che viene usato è ", p, "\n");
    g := PrimitiveRootMod(p); #in caso posso mettere da dove cercarla se la voglio più alta
    Print("La radice primitiva modulo p è g=", g, "\n");
    list := alice1(p, g);
    appo := bob1(p, g);
    list[3] := appo[1];
    list[4] := appo[2];
    a2 := alice2(p, list[1], list[4]);
    b2 := bob2(p, list[3], list[2]);
    Print("Le chiavi segrete condivise sono ", a2, " ", b2);
end;