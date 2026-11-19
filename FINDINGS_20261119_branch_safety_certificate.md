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

Mesuré (sortie dans `CONSOLE_20261119_axe_orthogonal.txt`) :
- **Porte « partie réelle ».** Le saut 2πi est purement imaginaire : autour de x = 0,7, ln(ln x) saute de 2πi, sa partie
  réelle de 0, et u + conj(u) de 0. Prendre Re (ou u + ū) efface donc le saut, comme l'exponentielle. La partie imaginaire,
  elle, garde le saut (2π) ; la conjugaison seule le retourne (+2πi → −2πi) sans l'effacer.
- **Le trou autour de 0 disparaît pour l'addition.** Avec le lemme A (cas w = 0 compris), E_add(x, y) = x − Log(exp(−y)) :
  certifié sur le carré PLEIN [−3, 3] + i[−3, 3] pour x et y, zéro compris (256 rectangles).
- **Témoin indépendant de la logique d'intervalles.** 18 000 valeurs calculées au log principal ordinaire, serrées contre
  le segment (0, 1) et l'axe négatif (jusqu'à 10⁻³⁰) : pire écart avec x·y et x + y = 1,1·10⁻⁴⁹.

Hypothèses ouvertes (non testées) :
- Porte « infini » : avec des intervalles à bornes infinies, la multiplication serait certifiable en x = 0 et y = 0
  (ln 0 = −∞ … exp(−∞) = 0 = x·y), ce qui fermerait le dernier trou.
- Invariance de branche : un arbre qui respecte la règle des portes calcule la même fonction quel que soit l'endroit où
  l'on place la coupure des logs en porte x.
- Vérification par un outil vraiment indépendant (Julia IntervalArithmetic ou PARI, installés) : les deux exécutions du
  19 novembre utilisent le même programme sur la même machine, ce n'est pas deux avis.

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
