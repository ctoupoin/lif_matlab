Protocole pour une expérience de LIF
====================================

Calibration spatiale et distorsions
-----------------------------------

On prépare la calibration spatiale avant de fermer la cellule, avec la cellule remplie d'eau (pas de
colorant).

* On utilise une mire (empruntée à Facundo, faite par HexisGroup) que l'on déplace dans la cuve remplie d'eau.
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


Protocole point à point
---------------------
* Allumer le bain thermique du laser, le laser et les deux caméras
* Appuyer sur "start" sur le bain thermique du laser pour lancer le refroidissement
* Sur le PC de manip (plume), lancer CamWare, OPSL (logiciel de contrôle du laser), Notepad++ (ou autre éditeur de texte), 
  et Anaconda Prompt
* Naviguer dans anaconda Prompt vers D:\Toupoint\thermal_plumes\data
* Ouvrir avec Notepadd++ D:\Toupoint\thermal_plumes\data\trigger_burst.py, y régler le nombre d'images et la fraquence d'acuisition
* Lorsque le bain thermique du laser a atteint la température de consigne, noter la température de la pièce indiquée
  par la sonde PT-100 externe
* Plonger la sonde de température externe dans la cellule, via le deuxièmme tube en partant de la gauche
* Laisser la sponde de température acquérir quelques points, puis noter la température obtenue dans la spreadhseet correspondante
* Noter le nom du fichier .db contenant les données des sondes de température, ainsi que le temps 
* Remplir le reste des entrées de la spreadsheet
* Retirer le cache en papier des caméras, accrocher la blouse "fond noir" derrière les caméras
* Eteindre la lumière dans la salle, allumer le voyant signalant l'utilisation du laser dans la pièce, et fermer le rideau laser
* Mettre les caméras en attente d'un trigger en cliquant sur Record
* Eteindre l'écran du PC de'acquisition de la température (Joe)
* Mettre les lunettes laser, régler dans OPSL la puissance laser voulue (4W généralement), Tourner la clé du boîter de
  contrôle laser sur ON, tourner la clé, mais garder le shutter du laser fermé
* Attendre 2min alors que le laser tire
* Ouvrir le shutter du laser
* Via Anaconda Prompt, exécuter python trigger_burst.py
* Eteindre l'écran du PC de manip (Plume)
* On peut rallumer les écrans une fois l'acquisition terminée
* Une fois l'acquisition terminée, vérifier sur la partie Recorder de Camware (en bas à gauche) que les images sont synchronisées
  sur les deux caméras
* Etindre le laser via OPSL, puis tourner la clé sur OFF sur le boîter de contrôle du laser, et enfin fermer le shutter
* Si c'est OK, File -> Save Raw Recorder Sequence -> D:\Toupoint\thermal_plumes\conv_[on/off]\[exp_date]\run[run_nulmber]\
* Vérifier que le format des fichiers est 16 bits TIFF file
* Pendant l'enregistrement des images, on peut remettre les caches sur les caméras, rallumer la lumière dans la pièce et
  éteidre la lumière de signalisation du laser )à l'extérieur de la salle
* Copier depuis Joe le fichier .db contenant les informations en température du run effectué
* Exporter les données de "log" en .tsv via DB Browser for SQ Lite dans D:\Toupoint\thermal_plumes\data\temperature_probes\
