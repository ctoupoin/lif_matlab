Protocole pour une expérience de LIF
====================================

Calibration spatiale et distorsions
-----------------------------------

On prépare la calibration spatiale avant de fermer la cellule, avec la cellule remplie d'eau (pas de
colorant).

* On utilise une mire (empruntée à Facundo) que l'on déplace dans la cuve remplie d'eau.
  On fait 20 plans, donc au moins 1 plan dans le plan de la nappe laser.

* La routine `camera_calibrationv11.m`, adaptée du code de calibration de
  `Jean-Yves Bouguet <http://www.vision.caltech.edu/bouguetj>`_ permet de déduire les
  paramètres intrinsèques et extrinsèques de la caméra, et de corriger les distorsions.
  Après avoir corrigé les distorsions, on travaille dans le plan de la nappe laser et on calcule
  la matrice de transformation (dans le même script).

Préparation de la rhodamine
---------------------------

* On prépare deux flacons de :math:`200~\mathrm{mL}`, l'un pour la Rhodamine 110, et l'autre pour la
  Rhodamine B. Attention,
  la Rhodamine B tend à précipiter donc il faut remélanger avant d'utiliser le flacon. On utilise une
  concentration de :math:`4.40\times 10^{-4}~\mathrm{mol/L}`, obtenue à partir des Rhodamine en poudre.

* Pour préparer la cellule expérimentale, on introduit quelques mL dans la cuve. Pour cela, on utilise
  une seringue avec un petit tuyau flexible qui permet d'introduire la solution de rhodamine directement
  dans le cœur. On fait cette opération avec la convection en route pour aider au mélange. Il faut
  une concentration de Rhodamine B plus élevée que de Rhodamine 110 à cause du cross-talk entre les
  spectres d'absorption/émission. On utilise :math:`9~\mathrm{mL}` du flacon de Rhodamine B, et
  :math:`1.5~\mathrm{mL}` du flacon de Rhodamine 110, ce qui correspond à des concentrations finales
  dans la cellule de :math:`4\times 10^{-6}~\mathrm{mol/L}` et :math:`6.6\times 10^{-7}~\mathrm{mol/L}`
  pour la Rhodamine B et la Rhodamine 110.

* On laisse mélanger pendant typiquement une journée.

Calibration en température
--------------------------

* On attend que le bain thermique du laser soit stable à la consigne (20°C), on allume le laser et on 
  le laisse
  tirer avec shutter fermé pendant au moins 2 minutes pour qu'il se stabilise en température.

* On place la cellule dans une configure de température homogène (possible aux températures supérieures
  à la température ambiante), et on place la sonde PT-100 au centre de la cellule.
  On prend 100 images (10 secondes). Mieux: prendre sur une durée d'un run
  classique.

* On sélectionne notamment une température pour établir la fonction de transfert (généralement après avoir
  bien relaxé pendant un week-end).

* On utilise directement le logiciel CamWare pour l'acquisition, avec déclenchement par le GBF en mode
  pulse.


Acquisition des images
----------------------

* Bien mettre un tissu noir derrière les caméras pour éviter les réflexions sur le tableau blanc, et
  obtenir un fond totalement noir. Bien positionner les cartons noirs autours des caméras et sur les
  profilés alu.

* On fait une image dans le noir (sans laser) pour vérifier qu'il n'y a pas de lumière parasite et
  déterminer l'offet des caméras. L'image d'offset est obtenue en moyennant 10 images. Autre possibilité:
  utiliser plutôt la médiane.

Traitement des images
---------------------

Pour obtenir une image calibrée à partir des images brutes, on utilise le script `processing.m`.

Pour chaque acquisition, les caractéristiques sont renseignées dans un tableau (fichier tsv): date, numéro,
nombre d'images, puissance du laser (4W en général), concentration de Rhodamine, temps d'exposition (95 ms
en général), fréquence d'acquisition (10 Hz en général), numéro de la database du suivi des sondes,
temps en heures depuis le début du suivi des sondes, switch (1 ou 0) pour définir si cette acquisition doit
être incluse pour le calcul.

Il y a deux fichiers tsv: un pour les runs avec convection, l'autre pour les runs dans convection.
Le script `processing.m` va chercher automatiquement les bons fichiers.


