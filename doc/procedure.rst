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

Protocole d'acuiqisition des images point par point
---------------------
* Allumer le laser, les caméras et le bain thermique du laser
* Attendre que le bain thermique du laser atteigne la température de consigne (20°C)
* Pendant que le bain thermique du laser refroidit, sur le PC de manip (Plume), lancer:
  Camware, OPSL (logiciel de contrôle laser), Anaconda Promptn et Notepadd++ (ou autre édtieur de texte)
* Vérifier les paramètres des caméras via Camware: Delay 0, Exposure 95ms, Trigger mode Ext. Exp. Start, Pixel Clock 100MHz
* Sur OPSL, cliquer sur la clé (qui est ON par défaut), puis sur Control -> Local
* Dans Anaconda Prompt, naviguer vers D:\Toupoint\thermal_plumes\data
* Dans Notepadd++, ouvir D:\Toupoint\thermal_plumes\data\trigger_burst.py et sélectionner le nombre d'images voulul
* Lorsque le bain thermique du laser a atteint la consigne, relever la température de la pièce, donnée par la sonde PT100 externe
* Insérer la sonde PT100 externe dans la cellule parle deuxième tube débouchant en partant de la gauche
  attendre quelques points pour obtenir la mesure de température
* Retirer la sonde PT100 externe, reboucher le tube débouchant
* Reporter sur la spreadheet correspondante (conv on ou off) les valeurs de la température de la pièce,
  à l'intérieur de la cellule, et dans les plaques, et les autres informations requises
* Eteindre la lumière dans la salle, allumer la lumière "laser" à l'extérieur de la salle
* Fermer le rideau laser, retirer les caches en papier des caméras et accrocher la blouse "fond noir" derrière les caméras
* Eteindre l'écran du PC contrôlant les sondes de température et la convection (Joe)
* Sur Camware, mettre les caméras en attente d'un trigger en cliquant sur Record (l'interface deviendra rouge)
* Si premier run de la journée, lancer trigger_burst.py avec 10 images sans illumination laser pour acquérir l'offset
* Eteindre l'écran du PC de manip durant l'acquisition pour diminuer la luminosité dans la salle
* Sinon, tourner la clé sur ON sur le boîtier de contrôle laser, définir une puissance de 4W sur OPSL et cliquer sur la clé
* Lancer un chrono de 2min
* A l'issue de ces deux minutes, ouvrir le shutter, et exécuter trigger_burst.py avec le nombre d'images voulu
* Eteindre l'écran du PC de manip durant l'acquisition
* Une fois l'acquisition terminée, couper le laser via OPSL, puis en tournant la clé sur le boîtier, et enfin fermer le shutter
* Dans Camware, vérifier dans le Recorder (en bas à gauche) les timestamps des images, qui doivent être synchronisées entre les deux caméras
(en ignorant les 5 premières images, qui ne le sont jamais).
* Si elles le sont, File -> Save Raw Recorder Sequence (ou Ctrl + S), et enregistrer les images dans
  D:\Toupoint\thermal_plumes\data\conv_[on/off]\[expdate]\run[run_number]\[exp_date]_run[run_number]
* Valider le prompt de changement de nom, faire Esc à la demande de commentaire
* Remettre les caches en papier sur les caméras, rallumer la lumière, éteinre la lumière "laser"
* Copier le fichier .db contenant les données de température du run effecuté de Joe vers Plume
* Ouvrir ce fichier (via DB Brower for SQlite par exmeple), et exporter le conetnu de "log" vers
  D:\Toupoint\thermal_plumes\data\temperature_probes\[exp_date].tsv
