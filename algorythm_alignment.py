import argparse

###ARGUMENTS

parser = argparse.ArgumentParser()
parser.add_argument('-F', '--fasta1', type=str, help=' chemin du premier fichier fasta')
parser.add_argument('-f', '--fasta2', type=str, help='chemin du deuxieme fichier fasta')
parser.add_argument('-s', '--typeseq', default="nuc", type=str, help="le type des sequences, --typeseq 'nuc' par defaut et --typeseq 'aa' pour des sequences proteiques")
parser.add_argument('-a', '--typeali', default="global", type=str, help="le type d'alignement, --typeali 'global' pour alignement global et 'local' pour alignement local")
parser.add_argument('-g', '--gap',default=-10, type=int, help="valeur associee au score d'une ouverture de gap")
parser.add_argument('-e', '--extension', default=-1,type=int, help="valeur associee au score d'extension de gap")
parser.add_argument('-M', '--match',default=2, type=int, help="valeur associee au score d'un match")
parser.add_argument('-p', '--mismatchpupy', default=1,type=int, help="valeur associee au score d'un match entre purines ou entre pyrimidines")
parser.add_argument('-m', '--mismatch', default=-1,type=int, help="valeur associee au score d'un mismatch entre 1 purine et 1 pyrimidine")
args = parser.parse_args()

if not args.fasta1 or not args.fasta2 :
    raise Exception("Fichier(s) fastas manquant(s), merci de remplir ces 2 valeurs.")

### LECTURE FASTA

def lit_fasta(fasta):
    """ Prend en entree le chemin d'un fichier existant dans lequel se trouve une sequence fasta.
    Retourne la sequence presente dans le fichier sous forme de chaine de caractere en majuscules."""
    dico={}
    nom="" 
    with open(fasta,'r') as f:
        long=len(f.readlines())
        f.seek(0)
        for l in range(0,long):
            line = f.readline()
            if line[0] == '>':
                nom=line[1:-1]
                dico[nom]=""
            else:
                if line[-1]=='\n':
                    dico[nom]=dico[nom]+line[:-1].upper()
                else: dico[nom]=dico[nom]+line[:].upper()
    return dico[nom]

### AFFICHAGE

def affichematrice(matrice):
    """Prend en entree une liste de listes et l'affiche sous forme de matrice."""
    for i in range(0,len(matrice)):
        ligne=""
        for j in range(0,len(matrice[0])):
            ligne=(ligne+', \t'+str(matrice[i][j]))
        print(ligne[3:])


### MATRICES DE SUBSTITUTION

def matsubstitutionAa():
    """Retourne la matrice de subsitution blosum62 sous forme de dictionnaire de dictionnaires."""
    blosum62 = {
    'C':{'C':9, 'S':-1, 'T':-1, 'P':-3, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-4, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':-2, 'Y':-2, 'W':-2},
    'S':{'C':-1,'S':4,  'T':1,  'P':-1, 'A':1,  'G':0,  'N':1,  'D':0,  'E':0,  'Q':0,  'H':-1, 'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
    'T':{'C':-1,'S':1,  'T':4,  'P':1,  'A':-1, 'G':1,  'N':0,  'D':1,  'E':0,  'Q':0,  'H':0,  'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
    'P':{'C':-3,'S':-1, 'T':1,  'P':7,  'A':-1, 'G':-2, 'N':-1, 'D':-1, 'E':-1, 'Q':-1, 'H':-2, 'R':-2, 'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-4, 'Y':-3, 'W':-4},
    'A':{'C':0, 'S':1,  'T':-1, 'P':-1, 'A':4,  'G':0,  'N':-1, 'D':-2, 'E':-1, 'Q':-1, 'H':-2, 'R':-1, 'K':-1, 'M':-1, 'I':-1, 'L':-1, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
    'G':{'C':-3,'S':0,  'T':1,  'P':-2, 'A':0,  'G':6,  'N':-2, 'D':-1, 'E':-2, 'Q':-2, 'H':-2, 'R':-2, 'K':-2, 'M':-3, 'I':-4, 'L':-4, 'V':0,  'F':-3, 'Y':-3, 'W':-2},
    'N':{'C':-3,'S':1,  'T':0,  'P':-2, 'A':-2, 'G':0,  'N':6,  'D':1,  'E':0,  'Q':0,  'H':-1, 'R':0,  'K':0,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-4},
    'D':{'C':-3,'S':0,  'T':1,  'P':-1, 'A':-2, 'G':-1, 'N':1,  'D':6,  'E':2,  'Q':0,  'H':-1, 'R':-2, 'K':-1, 'M':-3, 'I':-3, 'L':-4, 'V':-3, 'F':-3, 'Y':-3, 'W':-4},
    'E':{'C':-4,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':2,  'E':5,  'Q':2,  'H':0,  'R':0,  'K':1,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
    'Q':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':0,  'E':2,  'Q':5,  'H':0,  'R':1,  'K':1,  'M':0,  'I':-3, 'L':-2, 'V':-2, 'F':-3, 'Y':-1, 'W':-2},
    'H':{'C':-3,'S':-1, 'T':0,  'P':-2, 'A':-2, 'G':-2, 'N':1,  'D':1,  'E':0,  'Q':0,  'H':8,  'R':0,  'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-1, 'Y':2,  'W':-2},
    'R':{'C':-3,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-2, 'N':0,  'D':-2, 'E':0,  'Q':1,  'H':0,  'R':5,  'K':2,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
    'K':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':-1, 'E':1,  'Q':1,  'H':-1, 'R':2,  'K':5,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
    'M':{'C':-1,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':0,  'H':-2, 'R':-1, 'K':-1, 'M':5,  'I':1,  'L':2,  'V':-2, 'F':0,  'Y':-1, 'W':-1},
    'I':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':1,  'I':4,  'L':2,  'V':1,  'F':0,  'Y':-1, 'W':-3},
    'L':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-4, 'E':-3, 'Q':-2, 'H':-3, 'R':-2, 'K':-2, 'M':2,  'I':2,  'L':4,  'V':3,  'F':0,  'Y':-1, 'W':-2},
    'V':{'C':-1,'S':-2, 'T':-2, 'P':-2, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-2, 'Q':-2, 'H':-3, 'R':-3, 'K':-2, 'M':1,  'I':3,  'L':1,  'V':4,  'F':-1, 'Y':-1, 'W':-3},
    'F':{'C':-2,'S':-2, 'T':-2, 'P':-4, 'A':-2, 'G':-3, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-1, 'R':-3, 'K':-3, 'M':0,  'I':0,  'L':0,  'V':-1, 'F':6,  'Y':3,  'W':1},
    'Y':{'C':-2,'S':-2, 'T':-2, 'P':-3, 'A':-2, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':-1, 'H':2,  'R':-2, 'K':-2, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':3,  'Y':7,  'W':2},
    'W':{'C':-2,'S':-3, 'T':-3, 'P':-4, 'A':-3, 'G':-2, 'N':-4, 'D':-4, 'E':-3, 'Q':-2, 'H':-2, 'R':-3, 'K':-3, 'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':1,  'Y':2,  'W':11}
    }
    return blosum62

def matsubstitutionNt(match=2,mismatchpupy=1,mismatch=-1):
    """Prend en entree facultatives le score affecte au match, mismatch entre purines ou entre pyrimidines et mismatch(mismatchpupy) entre purine et pyrimidines. 
    Par defaut match=2, mismatchpupy=1, mismatch=-1.
    Retourne la matrice de substituion des nucleotides sous forme de dictionnaire de dictionnaires."""
    res={}
    liste=['C','T','A','G']
    purine=['A','G']
    for i in liste: #pour chaque nucleotide
        dico={}
        for j in liste: 
            if i==j :  #si nucleotide identique : match
                dico[j]=match
            elif i in purine and j in purine or i not in purine and j not in purine: #si nucleotide tous les 2 purines ou pyrimides : mismatch purine/pyrimidine
                dico[j]=mismatchpupy
            else: dico[j]=mismatch #sinon mismatch complet
        res[i]=dico
    return res


### MATRICES DE SCORE ET TRACEBACK GLOBAL

def initGlobal(seq1,seq2,gap=-10,extension=-1): #Initialisation des matrices
    """ Prend en entree les deux sequences, la valeur associee au score d'une ouverture de gap (gap) et d'extension de gap(extension).
    Par defaut le score de gap=-10 et extension=-1.
    Retourne les matrices de score et de traceback initialisees : avec la premiere ligne et colonne remplies de la façon Needleman et Wunsch."""
    score=[]
    traceback=[]
    ncol=len(seq1)
    nlignes=len(seq2)
    ##première ligne score
    lignescore=[0]
    lignescore.append(gap)
    for c in range(1,ncol): 
        lignescore.append(lignescore[c]+extension)
    score.append(lignescore)
    ##première ligne traceback
    lignetraceback=['done']
    for c in range(0,ncol): 
        lignetraceback.append('left')
    traceback.append(lignetraceback)
    ##première colonne des matrices et remplissage du reste avec '-'
    for l in range (0,nlignes):
        lignescore=[]
        lignetraceback=[]
        for c in range (0,ncol):
            if l==0 and c==0:
                lignescore.append(gap)
                lignetraceback.append('up')
            elif c==0 :
                lignescore.append(score[l][0]+extension)
                lignetraceback.append('up')
            lignescore.append('-')
            lignetraceback.append('-')
        score.append(lignescore)
        traceback.append(lignetraceback)
    return score,traceback


def matscoreGlobal(seq1,seq2, matricesub,gap=-10,extension=-1):
    """ Prend en entree les deux sequences, la matrice de substitution correspondante et la valeur associee au score d'une ouverture de gap (gap) et d'extension de gap(extension).
    Par defaut le score de gap=-10 et extension=-1.
    Retourne les matrices de score et de traceback completees avec la methode Needleman-Wunsch."""
    res,traceback=initGlobal(seq1,seq2,gap,extension) #matrice finale de score et matrice finale de traceback
    ncol=len(seq1)
    nlignes=len(seq2)
    for l in range (1,nlignes+1): #lignes
        gapOuvert=False
        for c in range (1,ncol+1): #colonnes
            qdiag=res[l-1][c-1]+matricesub[seq2[l-1]][seq1[c-1]]

            ###affectation des variables qup, qdiag et qleft
            if gapOuvert==True: #verification si il y a un gap avant
                qup=res[l-1][c]+extension
                qleft=res[l][c-1]+extension
            else :
                qup=res[l-1][c]+gap
                qleft=res[l][c-1]+gap
            
            ###remplissage des matrices avec la valeur maximale
            if qdiag >= qup and qdiag >= qleft:
                res[l][c]=qdiag
                traceback[l][c]='diag'
            elif qdiag <= qup and qup >= qleft:
                res[l][c]=qup
                traceback[l][c]='  up'
                gapOuvert = True
            else : 
                res[l][c]=qleft
                traceback[l][c]='left'
                gapOuvert=True
    return res,traceback


### MATRICES DE SCORE ET TRACEBACK LOCAL

def initLocal(seq1,seq2,gap=-10,extension=-1): #Initialisation des matrices
    """ Prend en entree les deux sequences, la valeur associee au score d'une ouverture de gap (gap) et d'extension de gap(extension).
    Par defaut le score de gap=-10 et extension=-1.
    Retourne les matrices de score et de traceback initialisees : avec la premiere ligne et colonne remplies de la façon Smith et Waterman."""
    score=[]
    traceback=[]
    ncol=len(seq1)
    nlignes=len(seq2)
    #première ligne score
    lignescore=[0]
    lignescore.append(gap)
    for c in range(1,ncol): 
        lignescore.append(lignescore[c]+extension)
    score.append(lignescore)
    #première ligne traceback
    lignetraceback=['done']
    for c in range(0,ncol): 
        lignetraceback.append('left')
    traceback.append(lignetraceback)
    #première colonne des matrices et remplissage du reste avec '-'
    for l in range (0,nlignes):
        lignescore=[]
        lignetraceback=[]
        for c in range (0,ncol):
            if l==0 and c==0:
                lignescore.append(gap)
                lignetraceback.append('up')
            elif c==0 :
                lignescore.append(score[l][0]+extension)
                lignetraceback.append('up')
            lignescore.append('-')
            lignetraceback.append('-')
        score.append(lignescore)
        traceback.append(lignetraceback)
    return score,traceback


def matscoreLocal(seq1,seq2, matricesub,gap=-10,extension=-1):
    """ Prend en entree les deux sequences, la matrice de substitution correspondante et la valeur associee au score d'une ouverture de gap (gap) et d'extension de gap(extension).
    Par defaut le score de gap=-10 et extension=-1.
    Retourne les matrices de score et de traceback completees avec la methode Smith-Waterman."""
    e=0
    res,traceback=initLocal(seq1,seq2,gap,extension) #matrice finale de score et matrice finale de traceback
    while e < len (res[0]):
        for i in res[0]:
            if i<0:
                res[0][e]=0
            e+=1
    ncol=len(seq1)
    nlignes=len(seq2)
    for l in range (1,nlignes+1): #lignes
        gapOuvert=False
        if res[l][0]<0:
            res[l][0]=0
        for c in range (1,ncol+1): #colonnes
            qdiag=res[l-1][c-1]+matricesub[seq2[l-1]][seq1[c-1]]
            if gapOuvert==True:
                qup=res[l-1][c]+extension
                qleft=res[l][c-1]+extension
            else :
                qup=res[l-1][c]+gap
                qleft=res[l][c-1]+gap
            
            ##remplissage de la matrice socre
            maximum=max(qup, qdiag, qleft)
            if maximum < 0 :
                maximum = 0
            res[l][c]=maximum 

            ##remplissage de la matrice traceback
            if qdiag >= qup and qdiag >= qleft:
                traceback[l][c]='diag'
            elif qdiag <= qup and qup >= qleft:
                traceback[l][c]='  up'
                gapOuvert = True
            else : 
                traceback[l][c]='left'
                gapOuvert=True
    return res,traceback


### ALIGNEMENTS

def alignementGlobal(matricescore,traceback,seq1,seq2):
    """Prend en entree les matrice de score, de traceback et les deux sequences etudiees.
    Affiche l'alignement optimal de ces sequences, le score total et les nombres de gaps, mismatches et matches."""
    ngap,nmismatch,nmatch=0,0,0
    sequence1=''
    alignement=''
    sequence2=''
    j=len(matricescore[0])-1 #colonnes
    i=len(matricescore)-1 #lignes
    score=matricescore[i][j] #score final

    ## lecture de la matrice traceback pour faire le chemin d'alignement
    while traceback[i][j] != 'done': #tant que l'alignement n'est pas fini
        if traceback[i][j] =='diag':
            sequence1=seq1[j-1]+sequence1
            sequence2=seq2[i-1]+sequence2
            i=i-1
            j=j-1
        elif traceback[i][j] =='left':
            ngap+=1
            sequence1=seq1[j-1]+sequence1
            sequence2='-'+sequence2
            j=j-1
        else:
            ngap+=1
            sequence1='-'+sequence1
            sequence2=seq2[i-1]+sequence2
            i=i-1
        if sequence1[0]=='-' or sequence2[0]=='-': #si en face d'un gap
            alignement=' '+alignement
        elif sequence1[0]==sequence2[0]: #si en face de la meme chose
            alignement='|'+alignement
            nmatch+=1
        else : 
            nmismatch+=1
            alignement=':'+alignement #si en face d'un autre
    print("AFFICHAGE DE L'ALIGNEMENT GLOBAL OPTIMAL : \n",sequence1, "\n",alignement,"\n",sequence2)
    print("Le score final de l'alignement est : ",score)
    print("L'alignement trouve : ",ngap, " gaps ( ), ",nmatch," matchs (|) et ",nmismatch," mismatchs (:).")


def alignementLocal(matricescore,traceback,seq1,seq2):
    ngap,nmismatch,nmatch=0,0,0
    sequence1=''
    alignement=''
    sequence2=''
    dicomax={'max':0}
    for l in range (len(matricescore)): #récuperer le score maximum
            m=max(matricescore[l])
            if int(m) > dicomax['max']:
                dicomax={'max':0, 'ligne':l,'colonne':matricescore[l].index(m)}
    i=dicomax['ligne']
    j=dicomax['colonne']
    score=matricescore[i][j] #score final

    ## lecture de la matrice traceback pour faire le chemin d'alignement
    while matricescore[i][j] != 0: #tant que l'alignement n'est pas fini
        if traceback[i][j] =='diag':
            sequence1=seq1[j-1]+sequence1
            sequence2=seq2[i-1]+sequence2
            i=i-1
            j=j-1
        elif traceback[i][j] =='left':
            ngap+=1
            sequence1=seq1[j-1]+sequence1
            sequence2='-'+sequence2
            j=j-1
        else:
            ngap+=1
            sequence1='-'+sequence1
            sequence2=seq2[i-1]+sequence2
            i=i-1
        if sequence1[0]=='-' or sequence2[0]=='-': #si en face d'un gap
            alignement=' '+alignement
        elif sequence1[0]==sequence2[0]: #si en face de la meme chose
            alignement='|'+alignement
            nmatch+=1
        else : 
            nmismatch+=1
            alignement=':'+alignement #si en face d'un autre
    print("\n   AFFICHAGE DE L'ALIGNEMENT LOCAL OPTIMAL : \n",sequence1, "\n",alignement,"\n",sequence2)
    print("Le score final de l'alignement est : ",score)
    print("L'alignement trouve : ",ngap, " gaps ( ), ",nmatch," matchs (|) et ",nmismatch," mismatchs (:).")



### FONCTION PRINCIPALE PIPELINE

def execAlignement(fasta1=args.fasta1,fasta2=args.fasta2,typeseq=args.typeseq, typeali=args.typeali,gap=args.gap,extension=args.extension,match=args.match,mismatchpupy=args.mismatchpupy,mismatch=args.mismatch):
    """Prend en arguments : 
    - 2 chemins de fichiers existants dans lequel se trouve une sequence fasta 
    - le type des sequences, typeseq='nuc' par defaut et typeseq='aa' pour des sequences proteiques
    - les valeurs associees au score d'une ouverture de gap (gap) et d'extension de gap(extension) pour la matrice de score. Par defaut le score de gap=-10 et extension=-1
    - les valeurs associees au score d'un match, mismatch entre purines ou entre pyrimidines et mismatch(mismatchpupy) entre purine et pyrimidine. Par defaut match=2, mismatchpupy=1, mismatch=-1.
    Affiche : les 2 sequences, les matrices de score et de traceback, l'alignement optimal trouve ainsi que le score total et les nombres de gaps, mismatches et matches"""
    seq1=lit_fasta(fasta1)
    seq2=lit_fasta(fasta2)
    print('SEQUENCE 1 = ',seq1, '\nSEQUENCE 2 = ', seq2)
    ###Alignement nucleotidique
    if typeseq=='nuc':
        ## Alignement global
        if typeali=='global':
            matricesub=matsubstitutionNt(match,mismatchpupy,mismatch)
            matricescore,traceback=matscoreGlobal(seq1,seq2, matricesub,gap,extension)
            print('\nLa matrice score est : \n')
            affichematrice(matricescore)
            print('\nLa matrice traceback est :\n')
            affichematrice(traceback)
            alignementGlobal(matricescore,traceback,seq1,seq2)
        ## Alignement local
        elif typeali=='local': 
            matricesub=matsubstitutionNt(match,mismatchpupy,mismatch)
            matricescore,traceback=matscoreLocal(seq1,seq2, matricesub,gap,extension)
            print('\nLa matrice score est : \n')
            affichematrice(matricescore)
            print('\nLa matrice traceback est :\n')
            affichematrice(traceback)
            alignementLocal(matricescore,traceback,seq1,seq2)
        else:
            raise Exception("Argument typeali invalide. Renseignez --typeseq 'global' pour un alignement global ou 'local' si local.")
    ###Alignement proteique
    elif typeseq=='aa':
        ## Alignement global
        if typeali=='global':
            matricesub=matsubstitutionAa()
            matricescore,traceback=matscoreGlobal(seq1,seq2,matricesub,gap,extension)
            print('\nLa matrice score est :\n')
            affichematrice(matricescore)
            print('\nLa matrice traceback est :\n')
            affichematrice(traceback)
            alignementGlobal(matricescore,traceback,seq1,seq2)
        ## Alignement local
        elif typeali=='local': 
            matricesub=matsubstitutionAa()
            matricescore,traceback=matscoreLocal(seq1,seq2,matricesub,gap,extension)
            print('\nLa matrice score est :\n')
            affichematrice(matricescore)
            print('\nLa matrice traceback est :\n')
            affichematrice(traceback)
            alignementLocal(matricescore,traceback,seq1,seq2)
    else:
        raise Exception("Argument typeseq invalide. Renseignez --typeseq 'nuc' pour une séquence nucléotidique ou 'aa' si protéique.")


###EXECUTION

if __name__ == '__main__': 
    execAlignement()
