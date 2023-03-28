L'algorithme se trouve dans le fichier projetNW_MEDINASolene.py et utilise argparse qui est une bibliothèque standard de python.

## Description

Il consiste à renvoyer l'alignement optimal de deux séquences données en :
    - lisant les fichiers fasta donnés pour obtenir les séquences
    - créant une matrice de substitutions qui donne les scores de matches et mismatches,
    - créant une matrice de score correspondant aux séquences et à la matrice de substitution
    - créant une matrice de traceback à partir de cette matrice de score
    - déduction de l'alignement optimal grâce à la matrice traceback
    - renvoi de cet alignement avec ses statistiques


Pour réaliser tout cela nous prenons en arguments (comment les renseigner est expliqué dans "Consignes d'utilisation" ): 
    - 2 chemins de fichiers existants dans lequel se trouve une sequence fasta 
    - le type des sequences, nucleotidique par defaut ou sinon pour des sequences proteiques
    - le type d'alignement, global comme la methode Needleman-Wunsch ou local comme Smith et Waterman
    - les valeurs associee au score d'une ouverture de gap (gap) et d'extension de gap(extension) pour la matrice de score. Par defaut le score de gap=-10 et extension=-1
    - les valeurs associee au score d'un match, mismatch entre purines ou entre pyrimidines et mismatch(mismatchpupy) entre purine et pyrimidines. Par defaut match=2, mismatchpupy=1, mismatch=-1.
    Affiche : les 2 sequences, les matrices de score et de traceback, l'alignement optimal trouve ainsi que le score total et les nombres de gaps, mismatches et matches


## Consignes d'utilisation

Se positionner dans le bon dossier puis :

python3.9  projetNW_MEDINASolene.py + arguments cités ci-dessous 

Voici les arguments utilisés dans cette pipeline : (présentés avec le raccourci, le nom complet, la valeur par defaut si besoin, le type de la variable et une description pour aider)

'-F' ou  '--fasta1', type=str, help='chemin du premier fichier fasta' OBLIGATOIRE

'-f' ou  '--fasta2', type=str, help='chemin du deuxieme fichier fasta'  OBLIGATOIRE

'-s' ou  '--typeseq', default="nuc", type=str, help="le type des sequences, --typeseq 'nuc' par defaut et --typeseq 'aa' pour des sequences proteiques"

'-a' ou  '--typeali', default="global", type=str, help="le type d'alignement, --typeali 'global' pour alignement global et 'local' pour alignement local"

'-g' ou  '--gap',default=-10, type=int, help="valeur associee au score d'une ouverture de gap"

'-e' ou  '--extension', default=-1,  type=int, help="valeur associee au score d'extension de gap"

'-M' ou  '--match', default=2, type=int, help="valeur associee au score d'un match"

'-p' ou  '--mismatchpupy', default=1, type=int, help="valeur associee au score d'un match entre purines ou entre pyrimidines"

'-m' ou  '--mismatch', default=-1, type=int, help="valeur associee au score d'un mismatch entre 1 purine et 1 pyrimidine"


## Exemples et résultats associés

 # Exemple 1 : séquences de nucléotides en global
 
 projetNW_MEDINASolene.py -F "C:/Users/Nouveau/Documents/.Licence/s6bioinfo/fastant1.fa.txt" -f "C:/Users/Nouveau/Documents/.Licence/s6bioinfo/fastant2.fa.txt"

 Résultat  :

 SEQUENCE 1 =  AGTCGATC 
 SEQUENCE 2 =  AGTACCG
 
 La matrice score est :
 
 0,      -10,    -11,    -12,    -13,    -14,    -15,    -16,    -17
 -10,    2,      -8,     -9,     -10,    -11,    -12,    -13,    -14
 -11,    -8,     4,      3,      2,      1,      0,      -1,     -2
 -12,    -12,    -6,     6,      5,      4,      3,      2,      1
 -13,    -10,    -11,    -4,     5,      6,      6,      5,      4
 -14,    -14,    -11,    -10,    -2,     4,      5,      7,      7
 -15,    -15,    -15,    -10,    -8,     -3,     3,      6,      9
 -16,    -14,    -13,    -16,    -11,    -6,     -2,     2,      5
 
 La matrice traceback est :
 
 done,   left,   left,   left,   left,   left,   left,   left,   left
 up,     diag,   left,   left,   left,   left,   diag,   left,   left
 up,       up,   diag,   left,   left,   left,   left,   left,   left
 up,     diag,     up,   diag,   left,   left,   left,   diag,   left
 up,     diag,   diag,     up,   diag,   diag,   diag,   left,   left
 up,     diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag
 up,     diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag
 up,     diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag
 AFFICHAGE DE L'ALIGNEMENT GLOBAL OPTIMAL :
  AGTCGATC
  ||| ::::
  AGT-ACCG
 Le score final de l'alignement est :  5
 L'alignement trouve :  1  gaps ( ),  3  matchs (|) et  4  mismatchs (:).



 # Exemple 2 : séquence d'acides aminés en local

 python3.9 projetNW_MEDINASolene.py -F "C:\Users\Nouveau\Documents\.Licence\s6bioinfo\fastaaa1.txt" -f "C:\Users\Nouveau\Documents\.Licence\s6bioinfo\fastaaa2.txt"  -a 'local' -s 'aa' 

 #Résultat :

 SEQUENCE 1 =  ARLGMRPLAA 
 SEQUENCE 2 =  AMRLRPLA
 
 La matrice score est :
 
 0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0
 0,      4,      0,      0,      0,      0,      0,      0,      0,      4,      4
 0,      0,      3,      2,      0,      5,      0,      0,      2,      0,      3
 0,      0,      5,      1,      0,      0,      10,     0,      1,      1,      2
 0,      0,      0,      9,      0,      2,      9,      8,      7,      6,      5
 0,      0,      5,      0,      7,      6,      8,      7,      6,      6,      5
 0,      0,      0,      2,      0,      5,      4,      15,     5,      5,      5
 0,      0,      0,      4,      0,      2,      3,      5,      19,     18,     17
 0,      4,      0,      0,      4,      0,      1,      2,      9,      23,     22
 
 La matrice traceback est :
 
 done,   left,   left,   left,   left,   left,   left,   left,   left,   left,   left
 up,     diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag
 up,     diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag,   diag
 up,     diag,   diag,   diag,   diag,   diag,   diag,   left,     up,   diag,     up
 up,     diag,   diag,   diag,   left,   diag,     up,   left,   left,   left,   left
 up,     diag,   diag,     up,   diag,   left,     up,   diag,   diag,   diag,   diag
 up,     diag,   diag,   diag,   diag,   diag,   diag,   diag,   left,   diag,   diag
 up,     diag,   diag,   diag,   diag,   diag,   diag,     up,   diag,   left,   left
 up,     diag,   diag,   diag,   diag,   diag,   diag,   diag,     up,   diag,   diag
 
    AFFICHAGE DE L'ALIGNEMENT LOCAL OPTIMAL :
  MR--PLA
  ||  |||
  MRLRPLA
 Le score final de l'alignement est :  23
 L'alignement trouve :  2  gaps ( ),  5  matchs (|) et  0  mismatchs (:).



 
 

