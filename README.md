[![Documentation Status](https://readthedocs.org/projects/plumex-lif/badge/?version=latest)](https://plumex-lif.readthedocs.io/en/latest/?badge=latest)

# Scripts de traitement de LIF

Dans ce dépôt se trouvent les scripts écrits par Clément Toupoint pendant son
[séjour post-doctoral](https://plumexlyon.github.io/arrivee-de-clement-toupoint.html)
dans le projet [PlumeX](https://plumexlyon.github.io/pages/a-propos.html).

La documentation est disponible [en ligne](https://plumex-lif.readthedocs.io/en/latest/).

Le code s'appuie sur celui de [Jean-Yves Bouguet](http://www.vision.caltech.edu/bouguetj/)
pour la calibration, modifié pour notre usage (dans `processing/calibration`).
L'adaptation pour des grilles de point utilise les scripts de
[George Vogiatzis et Carlos Hernández](http://george-vogiatzis.org/calib/).
La réduction de bruit est une version modifiée de filtre Kalman disponible sur
[Matlab Central](https://fr.mathworks.com/matlabcentral/fileexchange/26334-kalman-filter-for-noisy-movies).

Tous les codes sont en accès libre. Ceux de Jean-Yves Bouguet et ceux de George Vogiatzis et Carlos
Hernández ne précisent pas la licence. Le filtre Kalman est soumis à la licence de Mathworks.
Les codes écrits par Clément Toupoint dans le cadre de ce projet sont sous
licence [CeCILL-B](https://cecill.info/licences/Licence_CeCILL-B_V1-en.html).
