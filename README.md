###############
# Compilation #
###############
make

#############
# Execution #
#############
make exec n=a m=size seq=bool p=bool
Paramètres (optionnels):
- n : le nombre de processus ;
- m : la taille des blocs de découpage de la matrice ;
- seq : un entier. Si égal à 0 : parallèle, sinon séquentiel ;
- p : un entier. Si différent de 0, imprime les matrices L, U, et A finale.


##############
# Validation #
##############
make test N=n e=eps p=bool
Paramètres (optionnels):
- N : longueur d'un côté des matrices utilisées pour les tests ;
- e : précision au delà de laquelle deux valeurs sont considérées comme fausses ;
- p : un entier. Si différent de 0, imprime les matrices calculées.


#########
# Stats #
#########
make stat
ou 
make qsub
Cette dernière option seulement si vous êtes sur plafrim, mais attention au dossier contenant les sources. Dans ce cas, il faudra peut-être modifier le batch.
La série de speedup ainsi affichée par ligne correspond au nombre de processeur suivi du speedup calculé. 

###########
# Courbes #
###########
make plot
et
make plot-sp
Par défaut, la première commande trace les courbes du fichier sp-proc.data.
La deuxième, celles du fichier sp-size.data.
Vous pouvez choisir le fichier contenant les données par la variable stat.
Par exemple:
make plot stat=scatter
