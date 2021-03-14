# Raytracing

Le but de ce rapport est de détailler les méthodes mises en place pour réaliser le développement d’un algorithme de raytracing. L’objectif est de créer et de simuler le comportement de la lumière afin de réaliser un rendu le plus réaliste (prise en compte de reflets, d’ombres, de réfraction etc …). Ce TP se décompose en plusieurs phase :

- 1 : La mise en place de méthodes de calcul d'intersection rayons-sphère et génération de lumière
- 2 : La mise en place d’ombres portées, de correction gamma, de surfaces spéculaires et transparentes.
- 3 :La mise en place de l’équation du rendu, d’intégration de Monte-Carlo et d’éclairage indirect.
- 4 : La mise en place d’anti-Aliasing, d’ombres douces, de différents modèles de caméra, depth-of-field.
- 5 : La mise en place d’intersection rayon-plan, rayon-triangle, rayon-boite englobante, de gestion des maillages. Le but de toutes ces implémentations est de pouvoir travailler avec des formes plus complexes (maillages) et de diminuer les temps de calcul.

L’ensemble du code à été réalisé en C++. Ce code peut être trouvé en pièce jointe de ce rapport.

## 1. Mise en place du projet

La première étape consiste à créer les classes élémentaires de notre projet :

- L’objet Vector qui permet de stocker triplet et qui est surchargé d’un ensemble de méthodes utiles par la suite.
- L’objet Sphère qui permet de modéliser une sphère dans une scène. Une sphère est définie par un centre, un rayon et une couleur.
- L’objet Ray qui permet de définir un rayon. Cet objet est une demi droite définie par une origine et une direction.
- L’objet Scene qui stockera la collection d’objets qui seront stockés dans notre scène.

Une fois ces objets créés (cf.code), il est possible de créer notre première scène qui consistera en 5 sphères (4 pour les murs et 1 pour l’objet à observer) conformément au schéma proposé en cours :

![Structure de la cene](images/scene_structure.png)

La lumière de la scène est ponctuelle.

Il est important de noter que pour pouvoir calculer les intéractions de la lumière, il est nécessaire de calculer l'ensemble des intersections entre les sphères et les rayons lumineux.

On a alors l'équation suivante à résoudre pour chaque rayon et chaque sphère :

![Equation intesection shère - rayon](images/intersect_sph_ray.png)

Cette equation d'ordre 2 peut être résolue et ainsi détecter les points de rebond des rayons.

On implémente alors des fonctions intersect dans les classes Sphere et Scene afin de calculer toutes ces intersections

Une fois l'ensemble de ces objets initialisés dans la scene, on obtinent le résultat suivant :

![Première scène](images/first_scene.png)

## 2. Mise en place d'ombres, de surfaces mirroirs et transparentes
