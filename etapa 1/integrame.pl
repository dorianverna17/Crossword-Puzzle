:- ensure_loaded('checker.pl').

%test_mode(detailed).

% Considerăm următoarele reprezentări:
%
% O integramă este reprezentată prin structura (compusul)
% integ(H, W, Lista, Vocab), unde:
% H este înălțimea integramei
% W este lățimea integramei
% Lista este o listă de tupluri (Poz, Valoare), unde
%   Poz este un tuplu (R, C) conținând rândul și coloana (0-based)
%   Valoare este una dintre:
%     x - dacă celula este neagră (nu poate fi completată cu litere)
%     o literă, dacă celula este completată cu o literă
%     o listă de întrebări, reprezentate ca tupluri (Text, Dir, ID), cu
%       Text - un srting, textul întrebării
%       Dir - una dintre valorile j sau d, indicând direcția întrebării
%       ID - un identificator numeric al întrebării
% Vocab este o listă de stringuri reprezentând cuvinte disponibile
% pentru a rezolva întrebarea.
%
% În ieșirea predicatului intrebări, o întrebare este reprezentată ca
% ((R, C), Text, Dir, ID), unde
% R este rândul căsuței cu întrebarea (0-based)
% C este coloana căsuței cu întrebarea (0-based)
% Text este textul întrebării (un string)
% Dir este j sau d, reprezentând direcția în care trebuie să fie plasat
% răspunsul (jos sau dreapta)
% ID este un identificator numeric al întrebării.

% Puteți vizualiza integramele cu:
% integrama(0, W), print_integrama(W).
% integrama(1, W), print_integrama(W).
% integrama(2, W), print_integrama(W).
% integrama(3, W), print_integrama(W).
%
% Testați cu
% vmtest.
% Testați teste individuale (vedeți predicatul tt din checker.pl) cu
% vmtest(Test).
% de exemplu cu vmtest(intrebari).


% intrebari/2
% intrebari(integ(+H, +W, +Lista, +Vocab), -Lista_intrebari)
% Este adevărat atunci când Lista_intrebari este o lista de tupluri
% ((R, C), Text, Dir, ID), fiecare tuplu corespunzând unei întrebări din
% integramă (rândul, coloana, textul întrebării, direcția (j/d),
% identificatorul).
% BONUS: intrebari are o singură soluție (o singură listă) pentru o
% anumită integramă.

getTuples(integ(_, _, L1, _), L1).

% diferenta dintre doua multimi:
setDiff([], _, []).
setDiff([H1|T], L2, [H1|T1]):- \+member(H1, L2), setDiff(T, L2, T1).
setDiff([H1|T], L2, L3):- member(H1, L2), setDiff(T, L2, L3).

getTuplesWithNoX(_, []).
getTuplesWithNoX(T, [(H, V)|L]):- member((H, V), T), is_list(V), !,%\+ V = x, !,
	setDiff(T, [(H, V)], T1), getTuplesWithNoX(T1, L).

tuples(integ(_, _, L1, _), V):- getTuplesWithNoX(L1, V).

numberAll([], 0).
numberAll([(_, V)|T], Nr):- is_list(V), numberAll(T, N1), Nr is N1 + 1.
numberAll([(_, V)|T], Nr):- \+ is_list(V), numberAll(T, Nr).

tuplesWithNoX(integ(H, W, L1, Vocab), L2) :- tuples(integ(H, W, L1, Vocab), L2),
	numberAll(L1, X), length(L2, X).

transformTuples1([], []).
transformTuples1([((R, C), L1)|T], [((R, C), L2)|L]) :- length(L1, 1), L1 = [L2], transformTuples1(T, L).
transformTuples1([((R, C), L1)|T], [((G1, F1), K1)|L]) :- length(L1, 2), L1 = [A1|[A2]], L = [((G2, F2), K2)|Laux],
	((G1, F1), K1) = ((R, C), A1), ((R, C), A2) = ((G2, F2), K2), transformTuples1(T, Laux).

transformTuples2([], _, []).
transformTuples2([((R, C), (T, (D1, _)))|L1], Acc, [H2|L2]):- H2 = ((R, C), T, D1, Acc),
	Acc1 is Acc + 1, transformTuples2(L1, Acc1, L2).

% urmatoarea functie intoarce rezultatul dorit dar in forma ((R, C), [(T, D)])
intrebari(integ(H, W, L1, Vocab), L) :- tuplesWithNoX(integ(H, W, L1, Vocab), X),
	transformTuples1(X, Y), transformTuples2(Y, 0, L).

% id_intrebare/2
% id_intrebare(+Integ, ?Intrebare, ?Q_ID)
% Este adevărat dacă în integrama reprezentată ca integ(...), Intrebare
% este un text iar Q_ID este un identificator care corespund aceleași
% întrebări.

equalLists([], []).
equalLists([H1|L1], [H1|L2]):- length(L1, X), equalLists(L1, L2), length(L2, X).

equalStrings(S1, S2):- atom_chars(S1, L1), equalLists(L1, L2), atom_chars(S2, L2).

findQinIntegrams([], _, _):- false.
findQinIntegrams([((_, _), _, _, Id)|L], Q, Q_ID):- findQinIntegrams(L, Q, Q_ID), \+ Id = Q_ID. 
findQinIntegrams([((_, _), T, _, Q_ID)|_], Q, Q_ID):- equalStrings(T, Q).

id_intrebare(_, _, _):- false.
id_intrebare(W, Q, Q_ID):- intrebari(W, L), findQinIntegrams(L, Q, Q_ID).

% completare/3
% completare(+Integ, +Sol, -Integrama)
% Predicatul produce Integrama, o structură de forma integ(...),
% pornind de la Integ, în care au fost completate celule conform cu
% soluția Sol.
% Soluția este reprezentată ca o listă de perechi (Întrebare, Răspuns),
% unde Întrebarea este textul unei întrebări, iar Răspuns este un cuvând
% de completat; ambele sunt stringuri.
% De exemplu, o soluție parțială pentru integrama 0 poate fi:
% [('Din care plouă', 'NOR'), ('Al doilea număr', 'DOI')]
% BONUS: lungime_spatiu are o singură soluție pentru o anumită
% întrebare.
% Puteți testa manual predicatul cu o interogare de forma:
% integrama(0, W), solutie(0, Sol), completare(W, Sol, W2),
%   print_integrama(W2).

%integrama(0, W), solutie(0, Sol), completare(W, Sol, W2), print_integrama(W2).

findDirAndPosInIntegram([], _, _):- false.
findDirAndPosInIntegram([((R, C), T, D, _)|_], S1, ((R, C), D)) :-
	atom_chars(T, X), atom_chars(S1, Y), equalLists(X, Y).
findDirAndPosInIntegram([((_, _), T, _, _)|L], S1, ((R1, C1), D1)) :-
	atom_chars(T, X), atom_chars(S1, Y), \+ equalLists(X, Y),
	findDirAndPosInIntegram(L, S1, ((R1, C1), D1)).

putSolutions(I, ((_, _), d), _, 0, I).
putSolutions(I, ((_, _), j), _, 0, I).
putSolutions(integ(H, W, L, Vocab), ((R, C), d), [H2|S2], N, I3):- N > 0,
	C1 is C + 1, setDiff(L, [((R, C1), _)], Laux),
	I2 = integ(H, W, [((R, C1), H2)|Laux], Vocab), Nr is N - 1,
	putSolutions(I2, ((R, C1), d), S2, Nr, I3).
putSolutions(integ(H, W, L, Vocab), ((R, C), j), [H2|S2], N, I3):- N > 0,
	R1 is R + 1, setDiff(L, [((R1, C), _)], Laux),
	I2 = integ(H, W, [((R1, C), H2)|Laux], Vocab), Nr is N - 1,
	putSolutions(I2, ((R1, C), j), S2, Nr, I3).

putSolStringInI(I, S1, S2, Isol):- intrebari(I, L),
	findDirAndPosInIntegram(L, S1, ((R, C), D)),
	atom_chars(S2, L2), length(L2, Nr),
	putSolutions(I, ((R, C), D), L2, Nr, Isol).

putSolInI(I, (S1, S2), Isol):- putSolStringInI(I, S1, S2, Isol).

completare(I, [], I).
completare(I, [H], I3):- putSolInI(I, H, I2), completare(I2, [], I3).
completare(I, [H|S], I3) :- putSolInI(I, H, I2), completare(I2, S, I3), !.

%integrama(0, W), W = integ(A, B, C, D), completare(W, [('Afirmativ', 'DA'), ('Al doilea număr', 'DOI'), ('Primii 3 din artă', 'ART'),('Din care plouă', 'NOR')],integ(A1, B1, C1, D1)).

% lungime_spatiu/3
% lungime_spatiu(integ(+H, +W, +Lista, +Vocab), +Intrebare, -Lungime)
% Returnează lungimea spațiului asociat întrebării date.
% Întrebarea este indicată prin textul ei. De exemplu:
% lungime_spatiu pentru integrama 0 și întrebarea 'Al doilea număr'
% trebuie să lege Lungime la 3.
% BONUS: lungime_spatiu are o singură soluție pentru o anumită
% întrebare.
% Puteți testa manual predicatul cu o interogare de forma:
% integrama(0, W), id_intrebare(W, Text, 3), lungime_spatiu(W, Text, X).
lungime_spatiu(_, _, _) :- false.

% intersectie/5
% intersectie(integ(+H, +W, +Lista, +Voc), +I1, -Poz1, +I2, -Poz2)
% Pentru o integramă și două întrebări date prin textul lor (I1 și I2),
% al căror răspunsuri se intersectează, întoarce în Poz1 indicele din
% răspunsul la I1 la care este intersecția, și în Poz2 indicele din
% răspunsul la I2 la care este intersecția. Indecșii incep de la 0.
%
% De exemplu, în integrama 0:
%  █       █       2↓      3↓      █
%  █       0↓,1→   -       -       █
%  4→      -       -       -       █
%  5→      -       -       -       █
%  █       █       █       █       █
%
%  Întrebările 'Primii 3 din artă' și 'Afirmativ' (3, respectiv 1) se
%  intersectează la pozițiile 0, respectiv 2 (va fi litera A, de la
%  ART, respectiv DA).
intersectie(_, _, _, _, _) :- false.

% solutii_posibile/2
% solutii_posibile(integ(+H, +W, +Lista, +Vocabular), -Solutii)
% Formează o listă Solutii, conținând perechi de forma
% (Întrebare, Cuvinte), unde
% Întrebare este textul unei întrebări din integramă, iar Cuvinte este o
% listă de cuvinte sunt din Vocabular și au lungimea corectă pentru a fi
% răspuns la întrebare. Solutii conține câte o pereche pentru fiecare
% întrebare din integramă.
% Cuvintele sunt reprezentate ca liste de stringuri, fiecare string
% având lungime 1 (o singură literă).
% De exemplu, pentru integrama 0, Solutii conține 6 perechi, două dintre
% ele fiind:
% ('Afirmativ', [['D', 'A'], ['N', 'U']])
% ('Din care plouă',
% [['N','O','R'],['A','R','T'],['U','I','T'],['D','O','I']])
solutii_posibile(_, _) :- false.

% rezolvare/2
% rezolvare(+Integ, -Solutie)
% Rezolvare produce în Solutie soluția integramei Integ. Soluția este
% reprezentată ca o listă de perechi de stringuri, fiecare pereche
% conținând textul unei întrebări și cuvântul (ca string) care este
% răspunsul la întrebare.
rezolvare(_, _) :- false.
