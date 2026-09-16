# Analyse numérique du transport de sel en aquifère côtier

Stage de fin d'études — Master 2 Analyse Appliquée et Modélisation
(UPJV / laboratoires IMATH, CPT, MIO — Université de Toulon).

Modélisation de l'intrusion saline en zone de décharge d'eaux souterraines :
construction d'un schéma volumes finis pour le système couplé écoulement /
transport, démonstration de ses propriétés, et vérification numérique.

📄 **[Rapport complet](Rapport_de_stage complet.pdf)**

---

## Le modèle

Système couplé sous l'approximation de Boussinesq, en charge équivalente eau
douce $h_f$ et concentration en sel $C$ :

$$\nabla\cdot\vec q = 0, \qquad \vec q = -K_f\left(\nabla h_f + b(C)\,\vec e_z\right)$$

$$n_e\,\partial_t C + \nabla\cdot(\vec q\,C) - \nabla\cdot\left(n_e\,\mathbf{D}\,\nabla C\right) = 0$$

avec $b(C) = \dfrac{\rho(C)-\rho_f}{\rho_f}$ la flottabilité et la loi d'état
linéaire $\partial\rho/\partial C = 0{,}7143$.

---

## Choix de discrétisation

| Choix | Justification |
|---|---|
| **Volumes finis, flux à deux points (TPFA)** | Localement conservatif ; consistant d'ordre 2 sur maillage cartésien |
| **Moyenne harmonique pour $K$ aux faces** | $K$ est discontinue (plage / aquifère / mer) ; la moyenne arithmétique y est *inconsistante*, avec un facteur $\kappa \approx 2{,}5\times 10^{10}$ sur la face plage/aquifère |
| **Moyenne arithmétique pour $\rho$ aux faces** | $\rho$ est une variable d'état continue, pas un coefficient discontinu |
| **Pénalisation volumique de la zone marine** | La pente ne coïncide pas avec les faces du maillage cartésien ; évite un maillage conforme qui ferait perdre la K-orthogonalité |
| **Cible $h^{\text{ref}}(z) = -b_s z$** | Charge d'une colonne d'eau de mer au repos ; une cible nulle créerait une vitesse verticale fictive dans toute la mer |
| **Décentrement amont (upwind)** | Assure le signe des coefficients hors-diagonaux, donc le principe du maximum. Ordre 1 : c'est le prix de la monotonie (théorème de Godunov) |
| **Euler implicite** | Stabilité inconditionnelle du transport |
| **Couplage séquentiel explicite** | Densité gelée au pas précédent ; ordre 1 en temps, pas de temps limité par une CFL de *précision* et non de stabilité |

### Le résultat central

Le schéma en charge et la reconstruction du champ de vitesse utilisent
**exactement les mêmes coefficients de face** $K_\sigma$ et $b_\sigma$ :

$$\sum_{\sigma} \tau_\sigma (h_A - h_P) = -\sum_{\sigma} K_\sigma b_\sigma (\vec e_z\cdot\vec n_{P\sigma})\,|\sigma|
\qquad\text{et}\qquad
q_\sigma = -K_\sigma\left(\frac{h_A - h_P}{d_{PA}} + b_\sigma (\vec e_z\cdot\vec n_{P\sigma})\right)$$

Il en découle une chaîne de propriétés :

$$\text{mêmes coefficients} \;\Longrightarrow\; \sum_\sigma F_\sigma = 0 \;\Longrightarrow\; \text{somme des lignes} = n_e m_P \;\Longrightarrow\; \text{M-matrice} \;\Longrightarrow\; C \in [C_f, C_s]$$

La divergence discrète nulle est une propriété **algébrique**, valable sur tout
maillage — et c'est elle qui, transmise par récurrence sur les pas de temps,
donne la stabilité du schéma couplé, sans condition sur $\Delta t$.

---

## Les codes

Dépendances :

```bash
pip install fipy numpy sympy matplotlib scipy
```

### `mms_stationnaire.py`

Vérification par **solutions manufacturées** en régime permanent. Les champs
$h^{\text{ex}}$ et $C^{\text{ex}}$ sont choisis analytiquement, les résidus
calculés par `sympy` et injectés comme termes sources. Résout le système couplé
(itérations de Picard) et les deux équations découplées, sur cinq maillages
emboîtés.

Ordres mesurés à $N = 256$ :

| test | ordre |
|---|---|
| charge seule (flottabilité exacte) | **2,00** |
| transport seul (champ de Darcy exact) | **1,12** $\to 1$ |
| système couplé | **1,18** $\to 1$ |

L'ordre 2 de la charge confirme la consistance du flux à deux points ; l'ordre 1
du transport est celui du décentrement amont ; le couplé tend vers l'ordre du
plus faible, la contamination s'exerçant du transport vers la charge via la
flottabilité. Le résidu de continuité $|\mathrm{div}_h q + S_h|$ reste entre
$10^{-22}$ et $10^{-19}$ sans suivre la décroissance en $\Delta x$ : c'est le
résidu du solveur linéaire, pas une erreur de schéma.

### `mms_transitoire.py`

Même démarche avec le **terme d'évolution** $n_e\partial_t C$ et sa
discrétisation d'Euler implicite. Champs manufacturés modulés en temps,
conditions de Dirichlet dépendantes du temps, raffinement simultané
$\Delta t \propto \Delta x$. Ordre global mesuré : **1,14** à $N = 256$,
conforme à la valeur asymptotique attendue.

### `simulation_robinson.py`

Simulation d'intrusion complète : domaine $200 \times 32$ m, pente inclinée,
zone marine pénalisée, flux d'eau douce imposé — reproduisant le cas sans marée
de Robinson, Li & Barry (2007).

Diagnostics en cours de calcul : divergence discrète dans l'aquifère, bornes de
la salinité, flux échangés à l'interface mer/aquifère, position du pied du
biseau.

**Résultat :** masse de sel à l'équilibre $M = 14\,358$ kg contre
$14\,710$ kg dans l'article, soit un écart de **2,4 %** — réduit à **1,0 %**
après raffinement du maillage, le résidu restant de l'ordre du paramètre de
Boussinesq $\Delta\rho/\rho_f \approx 2{,}5\,\%$.

---

## Référence

C. Robinson, L. Li, D. A. Barry,
*Effect of tidal forcing on a subterranean estuary*,
**Advances in Water Resources** 30(4):851–865, 2007.
[doi:10.1016/j.advwatres.2006.07.006](https://doi.org/10.1016/j.advwatres.2006.07.006)
