# Certificat de branch-safety des arbres + et × — 2026-11-19

Fichiers : `certify_branch_safety.py` (certificat), `test_certify_branch_safety.py` (26 essais).
Rien n'a été modifié dans l'article, le README ni les scripts existants ; rien n'a été poussé sur GitHub.

## Résultat en une phrase

L'hypothèse publiée (« aucun nœud ne traverse la coupure ») est **fausse** pour l'arbre ×, mais la propriété dont
le Théorème 4.3 a réellement besoin (les arbres + et × calculent des fonctions **continues** sur K) est maintenant
**certifiée par ordinateur** sur K = [−3, 3] + i[−3, 3] privé du petit carré |Re|, |Im| < 0,1 autour de 0, pour x et y,
moyennant un lemme de trois lignes prouvé à la main (lemme A ci-dessous).

## 1. Ce qui était faux

- Dans `E_mul`, le sous-arbre `E_ln(E_ln(x))` prend le logarithme de ln x. Pour x réel dans (0, 1), ln x < 0 :
  ce log est pris **sur la coupure**. Mesuré : E_ln(ln(0,7 ± 10⁻¹² i)) = −1,031 ± 3,1416 i (saut de 2πi).
  En x = 1, il prend ln 0.
- Le résultat final reste juste (E_mul(0,7 ± 10⁻¹² i ; 1,3) = 0,91 ± 1,3·10⁻¹² i) : le saut est effacé plus haut par
  une exponentielle.
- L'Assumption de `eml_star_branch_safety_lemma_proof.tex` (« Im de chaque argument intermédiaire dans (−π, π] ») peut
  rester vraie (ln(ln x) vaut Im = +π juste au-dessus et ≈ −π juste en dessous) : ce n'est pas la bonne condition.
  La phrase « every composite tree … evaluates without branch-cut crossings on K » est fausse.
- Le test du README (« 1000 random points, zero violations ») vérifie que la valeur finale est juste, pas l'absence
  de traversée (mpmath calcule sur la coupure sans prévenir).

## 2. La bonne règle : les deux portes

`eml(x, y) = exp(x) − ln(y)`. Un log qui traverse l'axe réel négatif fait sauter la valeur de son nœud de 2πi·k.

- Nœud en **porte x** (premier argument du parent) : le parent calcule exp(nœud) ; exp(u + 2πik) = exp(u) : effacé.
- Nœud en **porte y** (son log est repris) ou à la **racine** : le saut reste : refusé.
- Log d'une valeur qui peut valoir **0** : refusé (sauf constantes, §3).

**Pourquoi un certificat réussi prouve la continuité.** Sur chaque rectangle accepté, le programme calcule avec une
branche du log continue sur ce rectangle (principale si le rectangle évite (−∞, 0] ; tournée, ln(−y) + iπ, si le
rectangle est à cheval sur l'axe négatif et loin de 0). Aux nœuds en porte x, la valeur ainsi calculée diffère de la
valeur « principale » de 2πik, ce que l'exponentielle du parent efface ; ailleurs elle est égale. La fonction de la
racine coïncide donc, sur chaque rectangle fermé, avec une fonction continue. Un nombre fini de rectangles fermés
recouvre K : la fonction est continue sur K (lemme de recollement). La correction repose sur l'arrondi vers
l'extérieur de mpmath.iv (50 chiffres).

## 3. Les constantes

Les sous-arbres sans variable sont calculés exactement sur les **réels étendus** [−∞, +∞] (pas la sphère de Riemann :
il faut distinguer exp(−∞) = 0 de exp(+∞) = +∞) : ln 0 = −∞, exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞ ; ∞ − ∞ et ln(−∞)
refusés. E_zero = e − ln(e^e) = 0 est décidé exactement par SymPy. Chaîne de E_neg : eml(1, 0) = +∞,
eml(+∞, 1) = +∞, eml(1, +∞) = −∞, puis exp(−∞) = 0.

## 4. Le lemme A (à la main) et le point x = 1

**Lemme A.** Pour tout w ∈ ℂ \ {0}, exp(E_ln(w)) = w.
*Preuve.* E_ln(w) = e − Log(exp(e − Log w)). Pour tout u ∈ ℂ, Log(exp u) = u − 2πi·k(u) avec k(u) ∈ ℤ. Donc
E_ln(w) = Log w + 2πi·k et exp(E_ln(w)) = exp(Log w) = w. Pour w = 0, avec les conventions du §3 :
eml(1, 0) = +∞, eml(+∞, 1) = +∞, eml(1, +∞) = −∞, exp(−∞) = 0 = w. ∎

**Conséquence.** eml(E_ln(a), c) = a − Log c pour tout a. Dans E_mul, le seul passage par ln(ln x) est
eml(E_ln(E_ln(x)), c) = E_ln(x) − Log c : on obtient l'arbre `E_mul_lemme`, égal à E_mul en tout point (vérifié aussi
en des points, jusqu'à 10⁻³⁰, dont x = 1 + 10⁻⁶ + 10⁻⁷ i). C'est `E_mul_lemme` que le programme certifie : x = 1 n'a plus
besoin d'être exclu, et le saut de ln x pour x réel négatif est effacé par l'exponentielle de la racine.

## 5. Mesures (`python3 certify_branch_safety.py` et essais)

| Arbre | Domaine (pour x et pour y) | Résultat |
|---|---|---|
| E_add | [0,1 ; 3] + i[−1, 1] (domaine du README) | certifié |
| E_mul (arbre d'origine) | x ∈ [0,9 ; 1,1] + i[−0,1 ; 0,1] | échec (ln(ln x)), attendu |
| E_mul (arbre d'origine) | x ∈ [0,5 ; 1,5] + i[−0,3 ; 0,3] privé du disque \|x − 1\| < 0,1 | certifié (traversées absorbées) |
| E_mul (arbre d'origine) | [0,1 ; 3] + i[−1, 1] privé du disque \|x − 1\| < 0,05 | certifié : 75 592 rectangles, 22 784 traversées absorbées, 163 s |
| E_mul_lemme | [0,1 ; 3] + i[−1, 1] (disque compris) | certifié en un rectangle |
| E_add | [−3, 3] + i[−1, 1] privé de \|Re\|, \|Im\| < 0,1 | certifié : 16 rectangles, 8 traversées absorbées |
| E_mul_lemme | idem | certifié : 16 rectangles, 24 traversées absorbées |
| E_add | [−3, 3] + i[−3, 3] privé de \|Re\|, \|Im\| < 0,1 | certifié : 2 056 rectangles, 584 traversées absorbées, 2 s |
| E_mul_lemme | idem | certifié : 16 rectangles, 24 traversées absorbées |

Les essais vérifient aussi que le programme **refuse** ce qui est vraiment discontinu (traversée en porte y, à la
racine, log d'un rectangle contenant 0, valeur qui saute utilisée sous un log).

## 6. Énoncé proposé pour remplacer l'Assumption (à décider par Anthony, rien n'est modifié)

> **Proposition (continuité certifiée des arbres arithmétiques).** Soit K = [−3, 3] + i[−3, 3] privé de
> {|Re z| < 0,1, |Im z| < 0,1}. Les arbres E_add et E_mul d'Odrzywołek, évalués avec le logarithme principal et les
> conventions des réels étendus pour leurs sous-arbres constants, définissent des fonctions continues sur K × K, égales
> à x + y et x·y. Des nœuds intermédiaires traversent la coupure (−∞, 0] (par exemple ln(ln x) pour x ∈ (0, 1)), mais
> chaque traversée est suivie d'une exponentielle qui en efface le saut de 2πi (règle des deux portes). Preuve : lemme A
> et certificat par intervalles `certify_branch_safety.py`.

À corriger aussi : la phrase « evaluates without branch-cut crossings on K », la ligne du tableau du README
(« numerically confirmed … »), et la section « Branch-Safety Verification » du README.

## 7. Limites

- K exclut un voisinage de 0 : ln 0 y est atteint par la variable, la multiplication passe par ln x.
- |Im| ≤ 3 < π : au-delà, ln(exp(y)) ≠ y et E_neg change de valeur ; non essayé.
- Tailles mesurées des arbres de `branch_safety_final.py` (feuilles + opérateurs) : E_add 27, E_mul 41. Le README
  annonce « RPN-length 19 pour +, 17 pour × » et le script lui-même « K=13 et K=20 » : ces nombres ne désignent pas les
  mêmes arbres. À vérifier contre `eml_compiler_v4.py` (Zenodo) avant de citer des longueurs.
- Sortie brute des calculs : `CONSOLE_20261119_branch_safety_certificate.txt`.

## 8. Axe orthogonal (même jour) : ce que toutes nos analyses gardaient fixe

`axe_orthogonal.py` est le script de REPRODUCTION de ce paragraphe (remis en fichier le soir même : la trace d'origine
venait des mêmes lignes tapées sans fichier). Il affiche ses résultats sans contrôler de verdict : son code de sortie
ne signale pas un échec. La commande qui contrôle les six verdicts et sort en erreur en cas d'échec reste
`python certify_branch_safety.py`. Sorties : `CONSOLE_20261119_axe_orthogonal.txt` (origine, avec la durée) et
`CONSOLE_20261119_axe_orthogonal_rejeu.txt` (rejeu ; mêmes résultats affichés). Portée précisée après relecture : ce
que le script OBSERVE et ce qui est DÉMONTRÉ sont séparés.

- **Partie réelle et saut de 2πi.**
  *Démontré (une ligne, sans le script)* : si les deux limites latérales de la valeur d'un nœud diffèrent de 2πi·k,
  k entier, leurs parties réelles sont égales, car 2πi·k est imaginaire pur ; de même pour u + conj(u) = 2 Re u. Les
  parties imaginaires diffèrent de 2πk ; la conjugaison change l'écart en −2πi·k sans l'effacer.
  *Observé par le script, en UN seul couple de points* (x = 0,7 ± 10⁻¹² i, nœud ln(ln x)) : différence des deux valeurs
  proche de 2πi (affichée `6.2831853j`), différence des parties réelles affichée 0, de u + ū affichée 0, des parties
  imaginaires 6,2831853. Le saut exact de 2πi est celui des limites latérales ; deux points à distance 10⁻¹² n'en
  donnent qu'une valeur approchée. C'est une illustration, pas une vérification sur un domaine. Remarque : Re n'est
  pas un opérateur de la grammaire {eml, eml★, 1} ; cette « porte » concerne la base alternative {eml, Re, 1} ou une
  partie réelle prise à l'extérieur de l'arbre.
- **Addition jusqu'à zéro : seulement pour la FORME RÉCRITE par le lemme.** Notons P = [−3, 3] + i[−3, 3] et K = P privé
  du petit carré ouvert |Re|, |Im| < 0,1.
  Cas 6 des certificats : la forme récrite x − Log(exp(E_neg(y))) est certifiée sur P × P, zéro compris (256 rectangles).
  Cas 4 : l'arbre d'origine E_add est certifié sur K × K (le petit carré est retiré de chacun des deux domaines, celui
  de x et celui de y).
  Le lien entre les deux formes est le lemme d'absorption, eml(E_ln(x), c) = x − Log c, y compris en x = 0 avec
  l'évaluation étendue : il repose sur la preuve, pas sur le programme.
- **Échantillon de points** (indépendant de la logique d'intervalles, même bibliothèque mpmath) : 6 000 couples (x, y),
  trois comparaisons par couple (E_mul et sa forme récrite contre x·y, E_add contre x + y), soit 18 000 écarts absolus,
  au log principal ordinaire. x est collé au segment (0, 1) ou à l'axe négatif (jusqu'à 10⁻³⁰), ou pris dans
  [−3, 3] + i([−3, −0,1] ∪ [0,1 ; 3]) ; y vérifie |Re y| ∈ [0,1 ; 3] et Im y ∈ [−1, 1]. Maximum affiché sur cet
  échantillon : 1,09·10⁻⁴⁹. Un échantillon ne prouve rien hors des points tirés.

Pistes notées à cette date, avec ce qu'elles sont devenues :
- Multiplication en zéro (« porte de l'infini »). Le produit admet un prolongement continu de valeur zéro. Mais
  évaluer l'arbre brut en x = 0 demanderait une convention supplémentaire ou une réécriture justifiée : E_ln(0) = −∞,
  puis le sous-arbre E_ln(E_ln(0)) demande log(−∞), indéfini selon la convention retenue ; élargir les bornes des
  intervalles ne résout pas cette opération.
- Invariance de branche : démontrée pour des décalages simultanés arbitraires des huit logarithmes à argument
  variable de l'arbre ×, les logarithmes constants gardant leur convention (§12, §13). Certains logarithmes constants
  peuvent aussi être décalés isolément sans effet (le n° 6, §15) ; le n° 19 ne le peut pas.
- Vérification par un outil indépendant : faite au §9, en Julia.

## 9. Deuxième avis : le même certificat récrit en Julia (`julia_check/certify.jl`, sortie `julia_check/sortie_julia.txt`)

Autre langage (Julia 1.12), autre bibliothèque d'intervalles (IntervalArithmetic.jl 1.0.14 sur MPFR, 200 bits), exp et log
complexes refaits à la main à partir des fonctions réelles, constantes exactes par un petit calcul symbolique maison (sans
SymPy). Ce qui reste commun aux deux versions : la règle (deux portes, lemme A) et son auteur.

| Certificat | Python (mpmath.iv) | Julia (IntervalArithmetic) |
|---|---|---|
| × d'origine, [0,1 ; 3] + i[−1, 1] hors \|x − 1\| < 0,05 | 75 592 rect., 22 784 absorbées, 98 304 exclues | **identique** : 75 592, 22 784, 98 304 (43 s) |
| + sur [−3, 3] + i[−3, 3] privé du carré autour de 0 | 2 056 rect., 584 absorbées | **identique** : 2 056, 584 |
| × forme lemme, même domaine | certifié, 16 rect. | certifié, 256 rect. |
| × forme lemme, [0,1 ; 3] + i[−1, 1] | certifié, 1 rect. | certifié, 16 rect. |
| + forme lemme, carré plein, zéro compris | certifié, 256 rect. | certifié, 256 rect. |

Les deux lignes où le nombre de rectangles diffère : même verdict, mais l'exp complexe refait à la main encadre un peu plus
large que celui de mpmath, donc Julia découpe une fois de plus. Les contrôles de refus (porte y, racine, log de zéro, × près
de x = 1 sans disque) passent aussi en Julia.

## 10. Le mur et l'infini : deux formules exactes (preuve à la main + mutant de branche)

Script `mur_et_infini.py`, sortie `CONSOLE_20261119_mur_et_infini.txt`.

Trois faits seulement servent : (i) exp(log w) = w pour toute branche du log, w ≠ 0 ; (ii) log(exp u) = u + 2πi·(entier)
pour toute branche ; (iii) exp(u + 2πik) = exp(u).

- E_ln(w) = e − log(exp(e − log w)) = log w + 2πi·k (w ≠ 0).
- E_neg(b) = 0 − log(exp(b)) = −b + 2πi·m (la constante de tête vaut −∞, exp(−∞) = 0). Aucune condition sur b.
- E_sub(a, b) = exp(E_ln(a)) − log(exp(b)) = a − b + 2πi·n (lemme A ; a = 0 compris, par les réels étendus).
- Donc **E_add(a, b) = a + b + 2πi·N**, N entier, et E_mul(x, y) = exp(E_add(E_ln x, E_ln y)) = exp(log x + log y + 2πi·N').

**Théorème B (multiplication).** Pour tous x, y ∈ ℂ \ {0}, et quelle que soit la branche choisie pour chacun des
logarithmes de l'arbre, E_mul(x, y) = x·y. (En x = 1 avec le log principal, E_ln(x) = 0 et l'arbre passe par les réels
étendus : lemme A.) Conséquences : aucune condition sur les parties imaginaires ; le mur est une pure convention pour × ;
la fonction se prolonge continûment par 0 en x = 0 et en y = 0 (le « trou » autour de 0 est une singularité effaçable :
à x = 0 exactement, l'arbre d'origine demanderait ln(−∞), que les réels étendus refusent).

**Théorème C (addition).** Avec le log principal, pour tout x ∈ ℂ et tout y ∈ ℂ, E_add(x, y) = x + enroule(y), où
enroule(y) est y ramené dans la bande −π < Im ≤ π par un nombre entier de tours 2πi (bord Im y = ±π à préciser).
En particulier E_add(x, y) = x + y si et seulement si y est dans la bande |Im y| < π : la condition porte sur y seul,
jamais sur x. Avec une autre branche, la bande se déplace : l'addition, elle, dépend du mur.

Mesures (80 chiffres) :
- Mutant de branche (chaque log reçoit 2πi·k, k tiré au hasard entre −3 et 3), 4 000 points, |Im| jusqu'à 12 :
  × : pire écart relatif avec x·y = 5,6·10⁻⁷⁹ ; + : (résultat − (x + y)) / 2πi toujours entier (écart 1,9·10⁻⁷⁹), entiers de −5 à 5.
- Log principal, 4 000 points, |Im| jusqu'à 40 (3 722 hors de la bande) : × : 1,6·10⁻⁸⁰ ; + contre x + enroule(y) : 4,1·10⁻⁷⁹.
  Exemple : E_add(1, 5i) = 1 − 1,2832 i = 1 + (5 − 2π) i.
- Vers zéro : |x| = 10⁻⁶⁰ : E_mul(x, y) = x·y à 3·10⁻⁷⁹ près (relatif), dans les deux positions.

Ce que cela change pour le Théorème 4.3 : la continuité de × sur tout compact de (ℂ \ {0})² et de + sur tout compact de
ℂ × {|Im y| < π} découle des théorèmes B et C (preuve de quelques lignes, à faire relire) ; les certificats par
intervalles (Python et Julia) en sont la confirmation par machine sur les domaines des §5, §8 et §9. Pour Stone–Weierstrass,
les fonctions manipulées devront rester dans la bande |Im| < π au moment d'être additionnées en position y : c'est la
vraie hypothèse, à écrire à la place de l'ancienne.

## 11. Ta barre face au mur (formule du miroir, Théorème 3.1)

Script `miroir_et_mur.py`, sortie `CONSOLE_20261119_miroir_et_mur.txt` (60 chiffres, 3 000 points, Im jusqu'à 15).

- **Le Théorème 3.1 est exact** dans sa bande −π ≤ Im z < π (674 points, écart 8·10⁻⁶²), bords compris (Im z = −π donne
  bien 0,5 + πi). Le domaine annoncé [−π, π) est le bon, y compris le côté fermé.
- **Hors de la bande, formule exacte nouvelle** : 1 − eml★(0, eml(z, 1)) = enroule(z̄), c'est-à-dire z̄ ramené dans la
  bande par des tours entiers (exemple z = 1 + 5i → 1 + 1,2832 i = 1 + (2π − 5) i). Preuve : la formule vaut
  Log(conj(exp z)) = Log(exp(z̄)).
- **Ce qui ne dépend jamais du mur** : exp(formule) = exp(z̄) partout, et la partie réelle de la formule = Re z partout
  (écart 8·10⁻⁶²). Avec le mutant de branche, formule − z̄ est toujours un nombre entier de tours 2πi (−5 … 5).
- **Conséquences pour ce qu'on construit avec** : Re z = (z + formule)/2 a toujours la bonne partie réelle, mais reçoit hors
  de la bande une partie imaginaire parasite (jusqu'à 2π ici) ; |z|² = z · formule n'est juste que dans la bande
  (z = 1 + 5i donne −5,42 + 6,28 i au lieu de 26).

**Défaut dans le code public : `eml_toolkit/core.py`.**
- `modulus_squared(z) = eml(z, eml★(z, 1))` vaut exp(z) − z, pas |z|² (z = 3 + 4i → −16,13 − 12,92 i au lieu de 25 ;
  `modulus(3 + 4i)` = 1,51 − 4,29 i au lieu de 5). Ce fichier est dans le dépôt GitHub (commit a4e2202).
- `real_part` et `imag_part` de core.py appellent mp.re / mp.im : pas des arbres eml★.
- `eml_toolkit/test_suite.py` ne teste PAS core.py : son import échoue (core.py n'a pas `conjugate_formula`, `eml_exp`…)
  et la suite bascule en silence sur ses propres définitions de secours (justes). D'où « 5 passed » malgré le défaut.
- Correction proposée (non appliquée) : modulus_squared(z) = z * E_conj(z) (ou E_mul(z, E_conj(z))), avec la
  condition −π ≤ Im z < π écrite dans la docstring, et une suite d'essais qui importe vraiment core.py (échec bruyant).

**Correction appliquée en local le 2026-11-19 (non commitée, non poussée)** : `eml_toolkit/core.py` récrit (original gardé :
`core.py.bak_2026-11-19`). L'ancien fichier était faux aussi pour E_sub (exp(a) − b au lieu de a − b), E_add et E_mul
(E_mul(1,3 + 0,1i ; 0,7 + 0,2i) donnait −1,7·10¹⁵ + 3·10¹⁶ i au lieu de 0,89 + 0,33i). Nouveaux essais
`eml_toolkit/test_core.py` (8, import direct) ; `test_suite.py` importe maintenant core.py pour de bon (plus de définitions
de secours ; original gardé en `.bak_2026-11-19`).

## 12. Relecture externe (ChatGPT, même jour) : cinq corrections acceptées

Vérifié sur la machine avant d'être accepté (`branches_autorisees.py`, `CONSOLE_20261119_branches_autorisees.txt`).

1. **« Pour toute branche de chaque logarithme » était FAUX (§10, théorème B).** Contre-exemple confirmé : ln 1 = 2πi au
   seul dernier nœud eml(u, 1) change le résultat d'exactement −2πi. Mon mutant ne le voyait pas : il exemptait ln 1.
   Mesure, un logarithme décalé à la fois sur les 20 de l'arbre × : les 8 logarithmes à argument variable ne changent
   rien ; ceux de la feuille 1 changent le résultat (décalage, valeur fausse ou indéfinie).
   Énoncé corrigé : log principal ; l'invariance ne vaut que pour les logarithmes dont l'argument dépend de x ou y.
2. **x = 1** : la convention des réels étendus était accordée aux seuls sous-arbres constants, alors que le nœud en cause
   contient la variable. Corrigé : évaluation étendue définie nœud par nœud, en chaque point.
3. **Bord de l'addition** : E_add(x, y) = x − P(−P(y)), P = Log∘exp ; le représentant a Im ∈ [−π, π) (vérifié : écart 0).
   E_add(x, iπ) = x − iπ ; E_add(x, −iπ) = x − iπ (juste).
4. **Théorème 4.3** : énoncé pour A_K, la sous-algèbre complexe unitaire de C(K, ℂ) engendrée par les fonctions d'arbres
   continues sur K ; elle contient ℂ[z, z̄] : la question de la constante i est réglée pour cet énoncé. La parenthèse
   « couper le second terme en N parts » est retirée (il faudrait un arbre pour g/N, y compris aux zéros de g).
5. **Portée** : ce qui est retiré est la conclusion « sans traversée de coupure », pas la borne sur Im (qu'un réel négatif
   ne réfute pas) ; la règle des deux portes est une condition SUFFISANTE (une traversée en seconde position n'est pas
   absorbée par le parent immédiat : cela ne prouve pas une discontinuité à la racine).

Les énoncés à jour sont ceux de `eml_star_branch_safety_lemma_proof.tex` et du README ; le §10 ci-dessus est conservé tel
qu'écrit, avec son erreur, pour l'historique.

## 13. Décalages simultanés (question de ChatGPT) et un argument de Gemini réfuté

- **Question** : tester chaque logarithme séparément ne prouve pas qu'on peut en déplacer plusieurs ensemble.
- **Réponse** : la preuve algébrique couvre les 8 décalages à la fois (récurrence des feuilles vers la racine, chaque valeur
  connue « à 2πi près ») ; c'est maintenant écrit dans la preuve de la proposition ×, avec le cas des arguments nuls
  (a' = Log x + 2πik = 0 impose x = 1 et k = 0).
- **Argument de Gemini réfuté par mesure** : « chaque logarithme variable entre directement dans une exponentielle ».
  Faux pour le log n° 2 (le log extérieur de E_ln(x)) : son nœud est en porte y, à l'intérieur de E_ln(E_ln(x)). Son décalage
  n'est pas effacé par le parent immédiat : il est transporté tel quel (lemme d'absorption) jusqu'à l'exponentielle de la
  racine. Les 7 autres sont en porte x. Conséquence : la règle des deux portes est bien seulement suffisante.
- **Mesure** (`branches_autorisees.py`) : décalages simultanés au hasard (−3…3 tours) sur les 8 logs, 3 000 points :
  pire écart relatif 4,6·10⁻⁵⁹ ; en x = 1 exactement avec 0, +1, −2 tours partout : écart ≤ 10⁻⁵⁹.

## 14. Deuxième relecture de ChatGPT : trois retouches acceptées

1. Logarithmes décalés définis aux valeurs exceptionnelles : L_j(w) = Log w + 2πi k_j pour w fini non nul, L_j(0) = −∞,
   L_j(+∞) = +∞. La phrase « la valeur de chaque nœud à 2πi près » était inexacte (un log appliqué à une valeur décalée
   change aussi de partie réelle) : la congruence porte sur les sorties des sous-arbres E_ln, E_neg, E_sub, E_add.
   Erreur de comptage corrigée : E_ln contient UN log de la feuille 1, pas deux.
2. Remarque après le Théorème 4.3 : représenter exactement tout élément de A_K par un seul arbre fini est IMPOSSIBLE
   (les arbres sont en nombre dénombrable, A_K contient toutes les constantes complexes). La densité uniforme des arbres
   uniques est une question distincte, ouverte.
3. Ligne du README alignée : « A_K is uniformly dense » ; ce qui est retiré est l'affirmation « sans traversée ».

Avis de ChatGPT après ces retouches : favorable au commit documentaire local ; la revue de core.py et des scripts reste distincte.

## 15. Relecture du code (auto-relecture puis ChatGPT) : défauts corrigés avant le second commit

Aucun de ces défauts ne change un verdict. Les deux implémentations annoncent le succès sur les six domaines testés.
Leurs compteurs concordent dans les cas 1, 2, 4 et 6. Les cas 3 et 5 aboutissent au même verdict avec des découpages
différents (cas 3 : 1 rectangle en Python, 16 en Julia ; cas 5 : 16 rectangles et 24 traversées comptées en Python,
256 et 400 en Julia).

1. **Mutant, log n° 6 classé « indéfini » à tort** (le script renvoyait `None` sans calculer). Calcul exact : décaler le
   ln 1 le plus profond de E_zero donne e − 2πi, puis exp(e − 2πi) = exp(e), et E_zero reste 0 : aucun effet. Le n° 7
   (ln 1 du milieu de E_zero) donne un résultat faux. Classement corrigé des 20 logarithmes : 8 à argument variable sans
   effet ; constants 8, 9, 11 sans effet ; feuille 1 : n° 6 et 10 sans effet, n° 19 décale de −2πi, n° 1, 4, 7, 13, 15, 17
   faussent le résultat. Régressions : `test_branches_autorisees.py`.
2. **Bord du miroir.** Valeur exacte de la formule : Log(exp(z̄)), donc ½ + iπ pour z = ½ + iπ ET pour z = ½ − iπ.
   L'arbre brut en flottant rend ½ − iπ au bord supérieur : sur le bord, exp(z) est à un cheveu de l'axe négatif et le
   log choisit un côté au hasard de l'arrondi. `core.py` : `conjugate_formula` reste l'arbre brut (documenté comme non
   fiable SUR le bord) ; nouveau `conjugate_folded`, repliement explicite avec convention garantie aux bords ; essais
   aux deux bords et de chaque côté.
3. **Description de l'ancien défaut** : l'ancien `modulus_squared` valait exp(z) − Log(exp(z)) ; « exp(z) − z » (§11) n'est
   vrai que dans la bande (pour z = 3 + 4i les deux diffèrent, la console le montrait).
4. **Constantes des certificats** : la marge absolue de 1e−55 était avalée par l'arrondi à 50 chiffres (intervalle de
   largeur nulle). Remplacé par un recalcul nœud par nœud en arithmétique d'intervalles, avec refus si un logarithme de
   constante touche la coupure. Un compteur prouve qu'AUCUNE constante finie non nulle n'intervient dans les six
   certificats : le compteur nul signifie précisément qu'aucune constante n'a emprunté le chemin de recalcul des
   constantes finies non nulles hors feuille 1. Les résultats publiés n'étaient pas touchés.
5. **Bornes décimales** (0,1 ; 0,05) : données en texte et arrondies vers l'extérieur ; disque exclu testé en intervalles
   (majorant de la distance² < minorant du rayon²).
6. **Succès à vide** : un certificat qui n'accepte aucun rectangle est maintenant un échec (Python et Julia).
7. **Reproductibilité** : `python3 certify_branch_safety.py` et `julia certify.jl` rejouent chacun les six certificats du
   README et sortent avec un code non nul si l'un échoue.

53 essais passent (certificat, mutant, toolkit).

## 16. Un miroir global de neuf symboles (nuit du 19 au 20 novembre 2026)

Origine : en comparant eml★ à une proposition publique de Quadruplay (un opérateur ternaire
T(x, y, z) = exp(x) − Log(y) + Re(z) Im(z), discussion sur le dépôt VA00/SymbolicRegressionPackage), ChatGPT a appliqué
le lemme d'absorption au DERNIER nœud et obtenu l'arbre ci-dessous ; une recherche exhaustive locale, lancée avant,
avait trouvé un autre arbre de même longueur, juste seulement dans la bande.

**Proposition (miroir global ; distincte du Théorème 3.1).** L'arbre de 9 symboles RPN, feuilles z et 1,

    F(z) = eml( eml( 1, eml( eml★(1, z), 1 ) ), 1 )        postfixé : 1 1 z eml★ 1 eml eml 1 eml

vérifie F(z) = conj(z) pour tout z ∈ ℂ \ {0} (log principal), et F(0) = 0 avec l'évaluation étendue
(sorties des quatre nœuds : +∞, +∞, −∞, 0). Sans évaluation étendue, son domaine est ℂ \ {0}.

*Preuve.* Pour z ≠ 0, posons w = Log(conj z) et v = exp(e − w), qui n'est jamais nul. Les deux derniers nœuds donnent
exp(e − Log v) = e^e / v = exp(w) = conj(z). Seul exp(Log u) = u (u ≠ 0) est utilisé, jamais Log(exp u) = u : aucune
bande n'intervient.

*Contrôle « une seule étoile ».* La même forme avec eml à la place de l'unique eml★ donne z.

Le Théorème 3.1 reste exact pour SA formule, 1 − eml★(0, eml(z, 1)), dans sa bande −π ≤ Im z < π ; hors de la bande
elle rend conj(z) ramené modulo 2πi (§11). C'est une autre formule qui réalise la conjugaison partout.

**Statuts.**

| Résultat | Statut | Pièces |
|---|---|---|
| F(z) = conj(z) pour z fini non nul ; contrôle sans étoile = z ; en i : −i et +i ; en 0 : trace +∞, +∞, −∞, 0 | ÉTABLI : preuve courte + juge exact SymPy 1.14.0, sans tolérance (4 PASS, rejoué par Anthony) | `juge_miroir_exact.py`, `CONSOLE_20261120_juge_miroir_exact.txt` |
| Mesure : 5 013 points jusqu'à \|Im\| = 200, axe négatif et bords ±π compris, écart relatif 7·10⁻⁵¹ | OBSERVÉ (échantillon) | `miroir_global.py`, `CONSOLE_20261119_miroir_global.txt` |
| Aucun arbre {eml, eml★} à feuilles {1, z} de longueur RPN ≤ 7 ne réalise le miroir : 9 est minimal | CERTIFIÉ par machine (à relire) : les 714 arbres (2 + 8 + 64 + 640) sont encadrés au point rationnel exact z0 = (37 + 61i)/100 en arithmétique d'intervalles ; chacun a une coordonnée de T(z0) − conj(z0) séparée de 0 (713 par la partie réelle, 1 par la partie imaginaire), 0 indécis ; les constantes (zéro, infini, réel négatif) sont décidées exactement par SymPy. Contrôle positif : F reste indécis en z0, comme attendu. L'unique élimination par la partie imaginaire est la feuille z (prédit en relecture : z0 − conj(z0) = 61i/50). Contrôles du code inclus dans le script : 714 arbres produits = 714 distincts = dénombrement (formes de Catalan × feuilles × opérateurs) ; l'intervalle de z0 contient strictement 37/100 et 61/100 (comparaison exacte des bornes binaires) ; la fonction de décision réellement utilisée n'élimine jamais un écart dont les deux intervalles contiennent 0 ; le code de sortie n'est 0 que si les cinq attendus sont vérifiés (714 arbres, 0 indécis, 713 par la partie réelle, la feuille z seule par la partie imaginaire, F non éliminé avec le motif exact) ; une constante dont SymPy ne sait pas décider le signe rend l'arbre indécis au lieu d'interrompre le calcul. Convention : log principal, évaluation étendue. La borne vaut déjà sur ℂ \ {0} (z0 ≠ 0) | `certificat_minimalite.py`, `CONSOLE_20261120_certificat_minimalite.txt` (et `minimalite_miroir.py`, version flottante antérieure) |
| Corollaire : pour tout compact K ⊂ ℂ, l'algèbre A_K est uniformément dense dans C(K, ℂ) (évaluation étendue si 0 ∈ K ; sinon tout compact de ℂ \ {0}) | DÉRIVÉ, à relire : la feuille donne z, F donne conj(z), continu sur K ; ℂ[z, z̄]\|_K ⊂ A_K et l'argument de Stone–Weierstrass déjà relu s'applique. Ne dépend pas de la minimalité | — |

**Comparaison avec l'opérateur ternaire de Quadruplay** (même feuilles z et 1, même convention). Son extraction de Re et
Im par translation est exacte (vérifié sur 2 000 points). Le témoin ternaire examiné ici (construit par ChatGPT, déplié
et vérifié sur 313 points) construit la constante i et compte 811 symboles une fois déplié (270 nœuds, 541 feuilles),
contre 9 pour F. La longueur minimale dans la grammaire ternaire reste indéterminée.

**Résultats de la recherche exhaustive** (`recherche_exhaustive_miroir.py`, `recherche_exhaustive_etoile_seule.py`,
`CONSOLE_20261120_recherche_exhaustive.txt` ; flottants, comparaison en 3 points, tolérance 1e−9 : HEURISTIQUE) :
- |z|² par un seul arbre de 17 symboles, e(e(1,e(e(e(1,e(e(1,s(1,z)),1)),z),1)),1) (e = eml, s = eml★), mesuré égal à
  |z|² sur 5 013 points (écart relatif 3·10⁻⁵⁰) ;
- Re(z), Im(z) et la constante i : aucun arbre jusqu'à la longueur 17 ;
- eml★ SEUL (sans eml), feuilles {1, z} : ni conj(z), ni Re, Im, −z, 1/z, i jusqu'à la longueur 19 (4,3 millions de
  fonctions au dernier niveau). Question ouverte : {eml★, 1} est-il à lui seul une base ?

Toujours ouvert : la densité uniforme des arbres INDIVIDUELS (§ 4.3 du manuscrit).

## 17. eml★ seule suffit : un miroir strict de 101 symboles (20 novembre 2026) — non relu par un mathématicien

Ce paragraphe répond à la question laissée ouverte à la fin du § 16 (« {eml★, 1} est-il à lui seul une base ? »).

**Auteur de la découverte : ChatGPT**, le 20 novembre 2026 (texte `PROOFS.md`, « EML orthogonal-axis study », reçu
dans l'archive `EML_orthogonal_20261120.zip`). Il a d'abord donné un miroir de 297 symboles, puis, après le retour de
Claude sur cette longueur, la version de 101 symboles ci-dessous. Anthony Monnerot a rejoué le juge sur son
ThinkCentre ; Claude a rebâti les arbres indépendamment (sans exécuter le script reçu), les a jugés avec ses propres
évaluateurs, puis a écrit le vérificateur autonome `miroir_etoile_seule_101.py`, qu'Anthony a rejoué. Aucun mathématicien n'a relu. Aucune antériorité n'a été cherchée.

**Énoncé (phrase proposée par ChatGPT).** « Un arbre de 101 symboles en eml★ seule réalise la conjugaison complexe en
évaluation stricte sur ℂ \ {0}. »

*Évaluation stricte* : chaque nœud de l'arbre a une valeur finie (aucun passage par ln 0 = −∞). *eml★ seule* : toutes
les portes sont eml★, les feuilles sont 1 et z ; l'arbre est compté entièrement déplié, comme au § 16.

**Construction.** On note S = eml★, S(a, b) = exp(a) − Log(conj b), log principal (Arg dans ]−π, π]). Chaque ligne
est un arbre ; la colonne de droite est la valeur démontrée, pas une syntaxe supplémentaire.

| Arbre | Valeur | Symboles |
|---|---|---|
| X(t) = S(t, 1) | exp(t) | 3 |
| C = X(X(1)) | e^e | 5 |
| O = S(1, C) | 0 | 7 |
| Q(t) = X(S(1, t)) | C / conj(t), t ≠ 0 | 5 |
| L(t) = S(1, Q(t)) | e − Log(C/t) : un logarithme de t, exp(L(t)) = t, t ≠ 0 | 7 |
| d(t, r) = S(L(t), X(r)) | t − r, pour t ≠ 0 et r réel | — |
| n = d(S(1, X(1)), X(1)) | −1 | 17 |
| k = d(S(O, n), 1) | −iπ (car S(O, n) = 1 − Log(−1) = 1 − iπ) | 35 |
| h = S(L(k), n) | −2πi : un tour complet | 59 |
| b = L(h) | un logarithme de h : exp(b) = −2πi | 65 |
| a = L(S(1, n)) | exp(a) = e − iπ, donc exp(exp(a)) = −C | 25 |
| P(t) = X(S(a, Q(t))) | −C / conj(Q(t)) = −t, t ≠ 0 : changement de signe | 33 |
| **J(z) = X(P(S(b, z)))** | **conj(z), z ≠ 0** | **101** |

**Pourquoi c'est strict sur tout ℂ \ {0}.** L'entrée donnée à P est w = S(b, z) = −2πi − Log(conj z). Sa partie
imaginaire est dans [−3π, −π[ : w n'est jamais nul, même en z = 1 et sur l'axe réel négatif. Le changement de signe
est donc licite, et J(z) = exp(−w) = exp(2πi) · conj(z) = conj(z). Les calculs n'utilisent que exp(Log u) = u ; ils
ne remplacent jamais Log(exp u) par u. C'est le « tour complet » (décaler un logarithme intermédiaire de −2πi, que
l'exponentielle suivante efface) qui évite le zéro intermédiaire.

**Conséquence : la porte eml reconstruite.** Pour y ≠ 0, eml★(x, J(y)) = exp(x) − Log(conj(conj y)) = eml(x, y) :
103 symboles. En évaluation stricte, toute porte eml d'un arbre mixte peut donc être remplacée par un montage en
eml★ seule ; les deux grammaires calculent les mêmes fonctions strictes.

| Résultat | Statut | Pièces |
|---|---|---|
| Miroir de 297 symboles en eml★ seule (première version, avec un multiplicateur général), J(z) = conj(z) pour z ≠ 0 | ÉTABLI par juge exact : SymPy 1.14.0 sur un symbole complexe non nul et cas particuliers ; encadrements Arb (python-flint 0.9.0) aux points listés ; rejoué par Anthony | `~/EML_replay_20261120_Vw6P7N/` : `replay.py` (sha256 6507109b…), `PROOFS.md` (928217ad…), `replay-console.txt` (ec77fe7a…) |
| Miroir de 101 symboles : uniquement eml★, feuilles {1, z}, longueur dépliée 101, constantes exactes, identité J(z) = conj(z) sur un symbole complexe non nul et en 9 points exacts (dont −1, −50, e^e), refus strict en 0, eml★(x, J(y)) = eml(x, y), contrôle sans étoile | ÉTABLI par juge exact, deux fois : (1) juge de ChatGPT rejoué par Anthony (sortie « PASS … 297 -> 101 », collée dans la conversation, non archivée en fichier) ; (2) script autonome écrit par Claude sans reprendre ce juge, rejoué par Anthony sur son ThinkCentre : Python 3.12.3, SymPy 1.14.0, 24 contrôles sur 24, code de sortie 0 (trace personnelle `~/EML_101_replay_NSwy5F/`, console sha256 d48f62c4…, identique octet pour octet à celle du dépôt). Ce sont deux implémentations distinctes, toutes deux vérifiées avec SymPy | `miroir_etoile_seule_101.py`, `CONSOLE_20261120_miroir_etoile_seule_101.txt` |
| Reconstruction indépendante des deux miroirs : 297 et 101 symboles, uniquement eml★, aucun nœud infini ; égalité à conj(z) sur 89 points difficiles (écart < 4·10⁻⁴⁷ ; axe négatif, bords ±π, \|Im\| jusqu'à 100) et en calcul exact SymPy en 5 points dont deux sur l'axe négatif ; refus strict en 0 | OBSERVÉ (échantillon + points exacts) | dossier de travail local `boite_double/` (NON versé au dépôt) : `etoile_seule.py`, `test_boite.py` (37 tests) |
| Domaine strict = tout ℂ \ {0} (Im w ∈ [−3π, −π[) | DÉRIVÉ : preuve courte ci-dessus, relue par Claude ; à faire relire | ce paragraphe |
| eml(x, y) = eml★(x, J(y)) pour y ≠ 0 | DÉRIVÉ (une ligne) + juge exact sur des symboles et en 3 points dont y réel négatif (script autonome) | `miroir_etoile_seule_101.py` |
| Produit strict, partie réelle stricte sur tout le plan (2 489 symboles), porte NAND sur des bits (535 symboles), en eml★ seule | OBSERVÉ sur échantillons (415 couples, 93 points, 4 entrées) ; preuves dans `PROOFS.md`, à relire | dossier local `boite_double/` (non versé), `PROOFS.md` (non versé) |
| Densité des arbres INDIVIDUELS sur les compacts ; obstruction stricte quand 0 ∈ K et que exp n'est pas injective sur K (exemple K = {0, 2πi}) | NON VÉRIFIÉ ici : preuve écrite de ChatGPT (`PROOFS.md` § 5–6), en attente de relecture indépendante | `PROOFS.md` (non versé au dépôt) |
| 101 est-il minimal en eml★ seule ? | OUVERT. La recherche exhaustive du § 16 n'a rien trouvé jusqu'à la longueur 19 (heuristique, flottants) : le minimum est donc vraisemblablement entre 21 et 101 | — |

**Précisions de portée.**
- *En zéro.* Les arbres de ce paragraphe refusent z = 0 en évaluation stricte. Ce n'est pas une impossibilité au seul
  point 0 (la feuille z y donne la bonne valeur). L'énoncé de ChatGPT, non vérifié ici, est qu'aucun arbre strict ne
  peut réaliser le miroir à la fois en 0 et en 2πi. Le miroir mixte de 9 symboles du § 16 atteint 0 en évaluation
  ÉTENDUE.
- *Le tour complet.* Retirer le décalage −2πi DANS CETTE CONSTRUCTION fait passer l'entrée 1 par ln 0 (contrôle local,
  non versé). Cela ne prouve pas qu'aucun autre montage ne puisse s'en passer.
- *Flottants sur la coupure.* En virgule flottante, eml★(x, J(y)) peut différer de eml(x, y) de 2πi quand y est réel
  négatif (observé en y = −50) : J(y) retombe à 10⁻⁴⁹ près du mauvais côté de l'axe. Le calcul exact donne zéro
  d'écart. C'est la fragilité numérique déjà signalée sur les bords pour `conjugate_formula` (`eml_toolkit/core.py`), pas un défaut de l'identité.
- *Les deux miroirs ne répondent pas à la même question.* 9 symboles : le plus court en grammaire mixte, minimal,
  juste jusqu'en 0 (étendu). 101 symboles : un seul type de porte, strict, hors de 0, minimalité inconnue.

Toujours ouvert : la longueur minimale du miroir en eml★ seule ; la relecture de la preuve de densité ; une relecture
humaine de l'ensemble.
