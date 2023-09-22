#funzione per RSA 
Print("Richiamare la funzione algoritmo passando per parametro due numeri naturali e il messaggio da trasmettere \n");

#funzione per lo scambio di chiavi
chiavi := function (p,q)
	local n, feul, i, r, s, t, list, bool;
	#se i numeri passati come parametro non sono primi trova i primi successivi (diversi)
	if(IsPrime(p)=false) then p := NextPrimeInt(p); fi;
	if(IsPrime(q)=false) then q := NextPrimeInt(q); fi;
	if(p=q) then q := NextPrimeInt(p); fi;
	Print("p=", p, ", q=", q, "\n");
	n := p*q;
	Print("n=", n, "\n");
	#feul è la funzione di eulero di n
	feul := (p-1)*(q-1);
	Print("Funzione di Eulero di n=", feul, "\n");
	#migliora se vado dal basso in alto con numeri solo primi? Oppure posso scendere sui numeri primi anziche salire?
	#in realtà non mi serve che i sia primo (quindi come lo scelgo? Randomicamente?)
	bool := false;
	#cerco un r coprimo con la funzione di eulero
	while (bool=false) do
		r := Random([feul/2..3*feul/4]);
		Print("Tentativo di r: ", r, "\n");
		if(Gcd(feul,r)=1) then bool := true; fi;
	od;
	#la chiave pubblica sarà (r,n)
	Print("La chiave pubblica sarà (r,n)= (", r, ",", n, ")", "\n");
	list := GcdRepresentation(r, feul);
	s := list[1];
	t := list[2];
	while s<0 do
		bool := false;
		#cerco un r coprimo con la funzione di eulero
		while (bool=false) do
			r := Random([feul/2..3*feul/4]);
			Print("Tentativo di r: ", r, "\n");
			if(Gcd(feul,r)=1) then bool := true; fi;
		od;
		#la chiave pubblica sarà (r,n)
		Print("La chiave pubblica sarà (r,n)= (", r, ",", n, ")", "\n");
		list := GcdRepresentation(r, feul);
		s := list[1];
		t := list[2];
	od;
	Print("La chiave privata sarà (s,n)=(", s, ",", n, ")", "\n");
	#la chiave privata sarà (s,n)
	list := [n,r,s];
	return list;
end;


#alice deve inviare un messaggio a bob
alice254 := function (r,n,parola)
	local l, bool, lBlocchi, list, k, nBlocchi, sum, letter, j, app, a, var, ret;
	l := Length(parola);
	bool := false;
	#lBlocchi è la lunghezza dei blocchi
	#divido in blocchi lunghi in modo tale avere n<254^lBlocchi
	lBlocchi := 1;
	while (bool=false) do
		if(n<254^lBlocchi) then bool := true; 
		else lBlocchi := lBlocchi+1;
		fi;
	od;
	lBlocchi := lBlocchi-1;
	Print("Lunghezza blocchi: ", lBlocchi, "\n");
	#divido la parola in blocchi lunghi lBlocchi
	if (Int(ZmodnZObj(l,lBlocchi))=0) then nBlocchi := l/lBlocchi;
	else nBlocchi := QuoInt(l,lBlocchi)+1; fi;
	Print("Numero blocchi: ", nBlocchi, "\n");
	list := [];
	k := 0;
	sum := 0;
	letter := 1;
	j := 1;
	#converto in base 254 ogni blocco usando per ogni carattere il corrispondente valore sulla tabella ascii
	#faccio padding con spazi vuoti
	while (letter<=lBlocchi*nBlocchi) do
		if (letter<=l) then
			sum := sum + IntChar(parola[letter]) * 254^k;
		else 
			sum := sum + IntChar(' ') * 254^k;
		fi;
		k := k+1;
		letter := letter + 1;
		if (k=lBlocchi) then k := 0; list[j] := sum; sum := 0; j := j+1; fi;
	od;
	if (Length(list) <> nBlocchi) then return fail; fi;
	#il comando Int(ZmodnZObj(a,n)) mi restituisce l'intero congruo ad a modulo n
	for a in [1..nBlocchi] do
		app := list[a]^r;
		list[a] := Int(ZmodnZObj(app, n));
	od;
	#ora nella mia lista ho il chiper text
	Print("Vettore di parola cifrata ", list, "\n");
	#devo ritornare la lista e la lunghezza dei blocchi
	ret := [];
	ret[1] := list;
	ret[2] := lBlocchi;
	return ret;
end;



#devo ora creare la funzione che rappresenta bob, che con la chiave privata decifra il messaggio
bob254 := function(s,n,stringa, lunBlocchi)
	local list, i, appo, j, inv;
	i := 1;
	list := [];
	inv := [];
	#faccio l'inverso di ogni intero inviato
	if s<0 then
		while i<=Length(stringa) do
			#if (Gcd(n,stringa[i]) <> 1) then return fail; fi;
			inv := GcdRepresentation(n,stringa[i]);
			stringa[i] := inv[2];
			#Print("Inverso ", i, " ", stringa[i], "\n");
			i := i+1;
		od;
		s := AbsInt(s);
	fi;
	Print("lunghezza ", Length(stringa), "\n");
	i := 1;
	#elevo gli inversi alla s mod n (ritornando da base 254 a base 10)
	while i<=Length(stringa) do
		appo := Int(ZmodnZObj(stringa[i]^s, n));
		for j in [1..lunBlocchi] do
			list[(i-1)*lunBlocchi+j] :=  RemInt(appo, 254);
			appo := QuoInt(appo, 254);
		od;
		i := i+1;
		od;
		Print("Prima decodifica ", list, "\n");
		for j in [1..Length(list)] do 
			list[j] := CharInt(list[j]);
	od;
	Print ("Il messaggio in chiaro è ", list);
end;


algoritmo := function(p,q,parola)
	local list, appo, i;
	list := [];
	appo := chiavi(p,q);
	for i in [1..Length(appo)] do
		list[i] := appo[i];
	od;
	#quindi list[1]=n list[2]=r, list[3]=s
	appo := alice254(list[2], list[1], parola);
	for i in [1..Length(appo)] do
		list [i+3] := appo[i];
	od;
	#quindi list[4] è la parola cifrata e list[5]=lBlocchi
	bob254(list[3], list[1], list[4], list[5]); 
end;