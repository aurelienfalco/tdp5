###############
# Compilation #
###############
make

#############
# Execution #
#############
make exec n=a bs=b
avec n le nombre de processus
avec bs la taille des blocs de découpage de la matrice.

##############
# Validation #
##############
make test N=n e=eps
N correspond aux tailles de matrices utilisées pour les tests.
e est la précision au delà de laquelle deux valeurs sont considérées comme fausses.

###########
# Courbes #
###########
make plot
Cette commande trace les courbes du fichier perfs.dat. Vous pouvez choisir de n'imprimer qu'une partie du code, en affectant "fox", "scatter" ou "gather" à la variable stat.
Par exemple:
make plot stat=scatter

Pour imprimer la courbe de speedup, à partir du fichier spup.dat :
make plot-sp

